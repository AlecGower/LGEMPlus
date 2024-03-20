# %%
from os import cpu_count
from typing import TYPE_CHECKING, Optional
import sys

if TYPE_CHECKING:
    from cobra.core import Model, Solution, Gene

from cobra.io import read_sbml_model
from cobra.core import Metabolite, Reaction
from cobra.flux_analysis import (
    single_gene_deletion,
    single_reaction_deletion,
    double_gene_deletion,
    double_reaction_deletion,
)

from cobra.flux_analysis.moma import moma, add_moma

import numpy as np
import pandas as pd

from tqdm.auto import tqdm
from copy import deepcopy
import time

import multiprocessing as mp
from functools import partial

# def log_message(string):
#     sys.stderr.write(string + "\n")

# %%
num_workers = int(np.floor(mp.cpu_count() / 2))
run_moma = True


# %%
def track_job(pool, job, update_interval=10):
    while job._number_left > 0:
        print(
            "Tasks remaining = {}".format(job._number_left * job._chunksize), sys.stderr
        )
        time.sleep(update_interval)


num_genes = len(fba_model.genes)


# %%
def load_list(fp):
    mets = []
    with open(fp) as fi:
        for ln in fi:
            try:
                mets.append(ln.rstrip().split("\t")[3])
            except IndexError:
                pass
    return ynb


# %%
def fba_deletion(model, gene):
    print("{} - calculating FBA deletion...".format(gene.id), sys.stderr)
    with model:
        model.genes.get_by_id(gene.id).knock_out()
        solution = model.optimize()
        print(
            "{} - done! ({}, {})".format(
                gene.id, solution.objective_value, solution.status
            ),
            sys.stderr,
        )
        return gene.id, solution.objective_value, solution.status


def moma_deletion(
    model: "Model", wt_solution: Optional["Solution"], gene: "Gene", linear: bool = True
) -> "Growth":
    print("{} - calculating MOMA deletion...".format(gene.id), sys.stderr)
    with model:
        model.genes.get_by_id(gene.id).knock_out()
        print("{} - adding MOMA constraints...".format(gene.id), sys.stderr)
        add_moma(model=model, solution=wt_solution, linear=linear)
        print("{} - MOMA constraints added...".format(gene.id), sys.stderr)
        # test tm_lim
        model.solver.configuration._smcp.tm_lim = 120000

        print("{} - optimising...".format(gene.id), sys.stderr)
        solution = model.optimize()
        print("{} - optimised.".format(gene.id), sys.stderr)

        try:
            if "moma_old_objective" in model.solver.variables:
                model.slim_optimize()
                growth = model.solver.variables.moma_old_objective.primal
            else:
                print("{} - MOMA growth not found.".format(gene.id), sys.stderr)
                growth = model.slim_optimize()
        except SolverError:
            growth = float("nan")

        print(
            "{} - done! ({}, {})".format(gene.id, growth, model.solver.status),
            sys.stderr,
        )
        return gene.id, growth, model.solver.status


def compile_results(res):
    print("All tasks completed!", sys.stderr)
    deletion_results = res
    deletion_results = pd.DataFrame(
        deletion_results, columns=["ids", "objective_value", "status"]
    )

    # deletion_results.ids = deletion_results.ids.apply(lambda s: list(s)[0])
    deletion_results = deletion_results.set_index("ids")
    print("Deletion results compiled.", sys.stderr)


# %%

ynb = load_list("model-files/ynb-compounds-yeastGEM.tsv")
ynb = [
    m
    for m in fba_model.metabolites
    if m.name[: m.name.find("[")].strip() in ynb and m.compartment == "e"
]
print("Yeast Nitrogen Base metabolites loaded.", sys.stderr)

ubiquitous = load_list("model-files/ubiquitous-compounds-yeastGEM.tsv")
ubiquitous = [
    m for m in fba_model.metabolites if m.name[: m.name.find("[")].strip() in ubiquitous
]
print("Ubiquitous metabolites loaded.", sys.stderr)

with fba_model:
    # Set medium
    medium = []
    for r in fba_model.exchanges:
        if (
            r.reactants[0] in ynb
            or r.reactants[0].id == "s_0565[e]"
            or r.reactants[0] in ubiquitous
        ):
            medium.append(r.id)
    fba_model.medium = {rid: 1000.0 for rid in medium}

    if run_moma:
        wt_soln = fba_model.optimize()
        print(wt_soln, sys.stderr)

    # # Not parallel
    # deletion_results = single_gene_deletion(fba_model, gene_list=fba_model.genes[:50])

    # Parallel
    print("Defining parallel functions.", sys.stderr)
    if run_moma:
        fba_model_deletion = partial(moma_deletion, fba_model, wt_soln)
    else:
        fba_model_deletion = partial(fba_deletion, fba_model)

    if __name__ == "__main__":
        print("Starting pool...", sys.stderr)
        with mp.Pool(processes=num_workers, maxtasksperchild=16) as pool:
            res = pool.map_async(
                fba_model_deletion,
                fba_model.genes[:num_genes],
                callback=compile_results,
            )
            pool.close()
            print("Starting simulations...", sys.stderr)
            track_job(pool, res)
            pool.join()
