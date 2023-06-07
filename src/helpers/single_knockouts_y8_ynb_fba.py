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
        print("Tasks remaining = {}".format(job._number_left * job._chunksize))
        time.sleep(update_interval)


# %%
path_to_model = "ModelFiles/yeastGEM.xml"
y8 = read_sbml_model(path_to_model)
# y8.solver.configuration.verbosity = 0
print("Yeast8 model loaded.")

num_genes = len(y8.genes)

# %%
ynb = []
with open("ModelFiles/ynb-compounds-yeastGEM.tsv") as fi:
    for ln in fi:
        try:
            ynb.append(ln.rstrip().split('\t')[3])
        except IndexError:
            pass

ynb = [m for m in y8.metabolites if m.name[:m.name.find('[')].strip() in ynb and m.compartment =='e']
print("Yeast Nitrogen Base metabolites loaded.")

# %%
ubiquitous = []
with open("ModelFiles/ubiquitous-compounds-yeastGEM.tsv") as fi:
    for ln in fi:
        try:
            ubiquitous.append(ln.rstrip().split('\t')[3])
        except IndexError:
            pass

ubiquitous = [m for m in y8.metabolites if m.name[:m.name.find('[')].strip() in ubiquitous]
print("Ubiquitous metabolites loaded.")

# %%
def fba_deletion(model, gene):
    print("{} - calculating FBA deletion...".format(gene.id))
    with model:
        model.genes.get_by_id(gene.id).knock_out()
        solution = model.optimize()
        print(
            "{} - done! ({}, {})".format(
                gene.id, solution.objective_value, solution.status
            )
        )
        return gene.id, solution.objective_value, solution.status


def moma_deletion(
    model: "Model", wt_solution: Optional["Solution"], gene: "Gene", linear: bool = True
) -> "Growth":
    print("{} - calculating MOMA deletion...".format(gene.id))
    with model:
        model.genes.get_by_id(gene.id).knock_out()
        print("{} - adding MOMA constraints...".format(gene.id))
        add_moma(model=model, solution=wt_solution, linear=linear)
        print("{} - MOMA constraints added...".format(gene.id))
        # test tm_lim
        model.solver.configuration._smcp.tm_lim = 120000
        
        print("{} - optimising...".format(gene.id))
        solution = model.optimize()
        print("{} - optimised.".format(gene.id))

        try:
            if "moma_old_objective" in model.solver.variables:
                model.slim_optimize()
                growth = model.solver.variables.moma_old_objective.primal
            else:
                print("{} - MOMA growth not found.".format(gene.id))
                growth = model.slim_optimize()
        except SolverError:
            growth = float("nan")

        print("{} - done! ({}, {})".format(gene.id, growth, model.solver.status))
        return gene.id, growth, model.solver.status


def compile_results(res):
    print("All tasks completed!")
    ynb_deletion_results = res
    ynb_deletion_results = pd.DataFrame(
        ynb_deletion_results, columns=["ids", "objective_value", "status"]
    )

    # ynb_deletion_results.ids = ynb_deletion_results.ids.apply(lambda s: list(s)[0])
    ynb_deletion_results = ynb_deletion_results.set_index("ids")
    print("YNB deletion results compiled.")

    # deletion_results = pd.read_csv("../results/deletion_results.csv")
    # dr = deletion_results
    # dr.ids = dr.ids.apply(lambda s: s[2:-2])
    # dr = dr.set_index("ids")
    # # print(dr.head())

    results_path = "experiments/results/fba_results"
    if run_moma:
        results_path += "_moma"
    results_path += ".csv"

    results = ynb_deletion_results #.join(dr, lsuffix=".ynb")
    results.to_csv(results_path)
    print("YNB deletion results written to '{}'.\nProgram terminated.".format(results_path))


# %%
with y8:
    # cut off exchange reactions
    for r in y8.exchanges:

        # except ynb exchanges
        if r.reactants[0] in ynb:
            # print(r.name)
            r.upper_bound = 1000.0

        # and except glucose
        if r.reactants[0].id == "s_0565[e]":
            # print(r.name)
            r.upper_bound = 1000.0

        # and except ubiquitous metabolites
        elif r.reactants[0] in ubiquitous:
            # print(r.name)
            r.upper_bound = 1000.0

        else:
            pass
            # print(r.lower_bound, r.name)
            r.upper_bound = 1.0

    if run_moma:
        wt_soln = y8.optimize()
        print(wt_soln)

    # # Not parallel
    # ynb_deletion_results = single_gene_deletion(y8, gene_list=y8.genes[:50])

    # Parallel
    print("Defining parallel functions.")
    if run_moma:
        y8_deletion = partial(moma_deletion, y8, wt_soln)
    else:
        y8_deletion = partial(fba_deletion, y8)

    if __name__ == "__main__":
        print("Starting pool...")
        with mp.Pool(processes=num_workers, maxtasksperchild=16) as pool:
            res = pool.map_async(y8_deletion, y8.genes[:num_genes], callback=compile_results)
            pool.close()
            print("Starting simulations...")
            track_job(pool, res)
            pool.join()

