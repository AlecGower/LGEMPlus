import sys
from tqdm.auto import tqdm
import os
from getopt import getopt
import re
from cobra.core.reaction import Reaction
from cobra.io import read_sbml_model
from cobra.core import Model, Solution

# from genesis.io import logical_model_from_sbml
import numpy as np

from typing import Optional, Tuple, Union, Iterable

# from multiprocessing import Pool, cpu_count
from joblib import Parallel, delayed, cpu_count
from functools import partial

from sympy import N

import pandas as pd

debug = False

cpus = cpu_count()
if debug:
    print("CPUs - {}".format(cpus), file=sys.stderr)

## Define helpful functions
def trim_reaction(
    model,
    rid,
    thresh: float = 10e-5,
    boundary_ids: Iterable[str] = [],
    ignore: Iterable[Union[str, Reaction]] = None,
    verbose: bool = False,
):

    rxn = model.reactions.get_by_id(rid)
    if (
        # any(
        #     [
        #         # rxn in model.exchanges, # possibly these lines that are taking lots of time
        #         # rxn in model.sinks,
        #         rxn in model.boundary, # this is the line that is taking lots of time
        #     ]
        # )
        rid in boundary_ids
        or rid in ignore
        or rxn in ignore
    ):
        if rid in boundary_ids and verbose:
            print("{:<6}: BOUNDARY REACTION {}.".format(rid, rxn.name), file=sys.stderr)
        if verbose:
            print("{:<6}: Ignoring reaction {}.".format(rid, rxn.name), file=sys.stderr)
        pass
    else:
        if verbose:
            print(
                "{:<6}: Setting reaction bounds to [{:>10},{:>10}] for {}.".format(
                    rid, (-1 * thresh), thresh, rxn.name
                ),
                file=sys.stderr,
            )
        rxn.lower_bound, rxn.upper_bound = (-1 * thresh), thresh


def promote_reaction(
    model,
    rid,
    reversed,
    thresh: float = 1e-9,
    boundary_ids: Iterable[str] = [],
    ignore: Iterable[Union[str, Reaction]] = None,
    verbose: bool = False,
):

    rxn = model.reactions.get_by_id(rid)
    if (
        # any(
        #     [
        #         # rxn in model.exchanges,
        #         rxn in model.sinks,
        #         rxn in model.boundary,
        #     ]
        # )
        rid in boundary_ids
        or "diffusion" in rxn.name.lower()
        or rid in ignore
        or rxn in ignore
    ):
        if verbose:
            print(
                "{:<6}: Ignoring reaction '{}'.".format(rid, rxn.name), file=sys.stderr
            )
        pass
    else:
        if reversed:
            if verbose:
                print(
                    "{:<6}: Setting reaction bounds to [{:>10},{:>10}] for '{}'.".format(
                        rid, -1000.0, (-1 * thresh), rxn.name
                    ),
                    file=sys.stderr,
                )
            rxn.lower_bound = -1000.0
            rxn.upper_bound = -1 * thresh
        else:
            if verbose:
                print(
                    "{:<6}: Setting reaction bounds to [{:>10},{:>10}] for '{}'.".format(
                        rid, thresh, 1000.0, rxn.name
                    ),
                    file=sys.stderr,
                )
            rxn.upper_bound = 1000.0
            rxn.lower_bound = thresh


def reactions_list_to_solution(
    model: Model,
    reaction_ids: Iterable[Tuple[str, bool]],
    mode: str,
    thresh: float = 1e-9,
    boundary_ids: Iterable[str] = [],
    verbose: bool = False,
    **kwargs
) -> Solution:

    with model:

        if mode == "trim":
            keep = [t[0] for t in reaction_ids]
            print("Trimming Reactions({}) ...".format(reactions_fp), file=sys.stderr)
            for rxn in tqdm(model.reactions):
                if rxn.id not in keep:
                    trim_reaction(
                        model,
                        rxn.id,
                        thresh=thresh,
                        boundary_ids=boundary_ids,
                        verbose=verbose,
                        **kwargs
                    )

        elif mode == "promote":
            print("Promoting Reactions({}) ...".format(reactions_fp), file=sys.stderr)
            for rid, reversed in tqdm(reaction_ids):
                # # future TO DO # rids processed (track of reverse and forward) [[[dumb solution for now]]]
                # print(
                #     "Promoting reaction: {}{}".format(
                #         rid, (" (reversed)" if reversed else "")
                #     )
                # )
                promote_reaction(
                    model,
                    rid,
                    boundary_ids=boundary_ids,
                    reversed=reversed,
                    thresh=thresh,
                    verbose=verbose,
                    **kwargs
                )
                # ## Reaction debugging
                # with model:
                #     soln = model.optimize()
                #     if soln.status == "infeasible":
                #         print(rid)
                #         breakpoint()

        else:
            pass
        print("Solving FBA({}) ...".format(reactions_fp), file=sys.stderr)
        return model.optimize()


# def hyp_to_soln(
#     hyp: str,
#     model: Model,
#     mode: str,
#     reaction_dicts: dict,
#     fp: str,
#     reaction_id_pattern: re.Pattern,
#     verbose: Optional[bool] = True,
#     return_soln: Optional[bool] = False,
#     inplace: Optional[bool] = False,
# ) -> Union[None, Solution]:
#     reaction_dicts[fp][hyp] = {
#         "ids": [
#             (reaction_id_pattern.match(rxn).group(1), "reverse" in rxn.lower())
#             for rxn in reaction_dicts[fp][hyp]
#         ],
#         "raw": reaction_dicts[fp][hyp],
#     }

#     # Solve and store solution
#     print("Solving FBA for {} ({})".format(hyp, fp), file=sys.stderr)
#     soln = reactions_list_to_solution(
#         model=model,
#         mode=mode,
#         boundary_ids=[r.id for r in model.boundary],
#         reaction_ids=reaction_dicts[fp][hyp]["ids"],
#         ignore=ignore_rxns,
#     )

#     # Print output
#     if verbose:
#         print('"{}";"{}";"{}"'.format(fp, hyp, soln.objective_value))

#     if inplace:
#         reaction_dicts[fp][hyp]["solution"] = soln

#     if return_soln:
#         return fp, hyp, soln


## Process arguments
arguments, values = getopt(
    sys.argv[1:],
    None,
    [
        "reactions=",
        "model=",
        "mode=",
        "threshold=",
        "verbose=",
        "ignore_rxn=",
        "ko_gene=",
    ],
)
reactions_fp = [a[1] for a in arguments if a[0].strip() == "--reactions"][0]
model_fp = [a[1] for a in arguments if a[0].strip() == "--model"][0]
mode = [a[1] for a in arguments if a[0].strip() == "--mode"][0]
try:
    verbose_output = [a[1] for a in arguments if a[0].strip() == "--verbose"][0]
    verbose_output = eval(verbose_output.capitalize())
except IndexError:
    verbose_output = False
try:
    thresh = float([a[1] for a in arguments if a[0].strip() == "--threshold"][0])
except IndexError:
    thresh = 1e-9
ignore_rxns = [a[1] for a in arguments if a[0].strip() == "--ignore_rxn"]
genes_ko = [a[1] for a in arguments if a[0].strip() == "--ko_gene"]
assert mode in ["trim", "promote"]
# try:
#     verbose_output = bool(verbose_output)
# except TypeError:
#     verbose_output = False
# assert verbose_output in [True, False]

# print(reactions_fp)

script_mode_msg = (
    "Running script in {} mode with absolute flux threshold {}. Ignoring the following reactions:".format(
        mode, thresh
    )
    + "\n"
    + "\n".join(["\t{}".format(r) for r in ignore_rxns])
)
print(script_mode_msg, file=sys.stderr)

## FBA model
# Load FBA model
fba_model = read_sbml_model(model_fp)
# Get reference FBA solution
print("Solving reference solution({}) ...".format(reactions_fp), file=sys.stderr)
with fba_model:
    for gn in genes_ko:
        try:
            fba_model.genes.get_by_id(gn).knock_out()
        except KeyError:
            [g for g in fba_model.genes if g.name.upper() == gn.upper()][0].knock_out()
    reference_solution = fba_model.optimize()
    # print(reference_solution.status, reference_solution.objective_value)
print(
    "REFERENCE",
    reactions_fp,
    reference_solution.status,
    reference_solution.objective_value,
    genes_ko,
)


# ## Load reactions for each hypothesis in each file
# # Get filepaths in directory
# if reaction_dir.endswith("/"):
#     reaction_dir = reaction_dir[:-1]
# reaction_dicts_fps = [
#     reaction_dir + "/" + fp for fp in os.listdir(reaction_dir) if fp.endswith(".json")
# ]

# # Create empty dictionary for each separate file
# reaction_dicts = dict()
# # Create regex for reaction ID
reaction_id_pattern = re.compile(r"^(r_[0-9]+)")
# # Load dictionaries from file and extract IDs
# for fp in reaction_dicts_fps:
# print("Evaluating:", fp)
# with open(fp, "r") as fi:
#     reaction_dicts[fp] = eval(fi.read())

# for hyp in reaction_dicts[fp]:
#     hyp_to_soln(hyp=hyp, reaction_dicts=reaction_dicts, fp=fp)

## Multiprocessing approach
# if __name__ == "__main__":
#     with Pool(processes=(cpus - 4)) as pool:
#         hypotheses = reaction_dicts[fp].keys()
#         # print(hypotheses)
#         f_ = partial(
#             hyp_to_soln,
#             model=fba_model,
#             mode=mode,
#             reaction_dicts=reaction_dicts,
#             fp=fp,
#             reaction_id_pattern=reaction_id_pattern,
#             verbose=False,
#             return_soln=True,
#             inplace=False,
#         )
#         outputs_async = pool.map_async(f_, hypotheses)
#         solutions = outputs_async.get()

#     # Print solutions
#     # print("Solutions:")
#     # print(solutions)
#     print(
#         *map(
#             lambda t: '"{}";"{}";"{}";"{}"'.format(*t),
#             zip(
#                 [fp] * len(hypotheses),
#                 hypotheses,
#                 [soln.status for soln in solutions],
#                 [soln.objective_value for soln in solutions]
#             ),
#         ),
#         sep="\n"
#     )

with open(reactions_fp) as fi:
    reactions = [
        # (reaction_id_pattern.match(rxn).group(1), "reverse" in rxn.lower())
        ## Different reaction input
        # (rxn.split(";")[0], rxn.split(";")[1] == "-1")
        (rxn.rstrip(), False)
        for rxn in fi.readlines()
    ]
# breakpoint()
# print("Solving FBA({}) ...".format(reactions_fp))
with fba_model:
    for gn in genes_ko:
        try:
            fba_model.genes.get_by_id(gn).knock_out()
        except KeyError:
            [g for g in fba_model.genes if g.name.upper() == gn.upper()][0].knock_out()
    soln = reactions_list_to_solution(
        model=fba_model,
        mode=mode,
        thresh=thresh,
        reaction_ids=reactions,
        boundary_ids=[r.id for r in fba_model.boundary],
        ignore=ignore_rxns,
        verbose=verbose_output,
    )

print("RESULT_WITH_KO", reactions_fp, soln.status, soln.objective_value, genes_ko)

genes_ko = []

with fba_model:
    for gn in genes_ko:
        try:
            fba_model.genes.get_by_id(gn).knock_out()
        except KeyError:
            [g for g in fba_model.genes if g.name.upper() == gn.upper()][0].knock_out()
    soln = reactions_list_to_solution(
        model=fba_model,
        mode=mode,
        thresh=thresh,
        reaction_ids=reactions,
        boundary_ids=[r.id for r in fba_model.boundary],
        ignore=ignore_rxns,
        verbose=verbose_output,
    )

print("RESULT_WITHOUT_KO", reactions_fp, soln.status, soln.objective_value, genes_ko)
# print(np.mean(soln.fluxes.abs() == thresh))
# print(soln.fluxes.abs().head())

# ## JobLib approach
# f_ = partial(
#     hyp_to_soln,
#     model=fba_model,
#     mode=mode,
#     reaction_dicts=reaction_dicts,
#     reaction_id_pattern=reaction_id_pattern,
#     verbose=False,
#     return_soln=True,
#     inplace=False,
# )


# def joblib_solve(fp_hyp: Tuple[str, str]):

#     return f_(hyp=fp_hyp[1], fp=fp_hyp[0])


# task = (
#     delayed(f_)(fp=fp, hyp=hyp)
#     for fp in reaction_dicts.keys()
#     for hyp in reaction_dicts[fp].keys()
# )

# with Parallel(n_jobs=cpus - 8) as parallel:
#     solutions = parallel(task)

# print(
#     *map(
#         lambda t: '"{}";"{}";"{}";"{}"'.format(
#             t[0], t[1], t[2].status, t[2].objective_value
#         ),
#         solutions,
#     ),
#     sep="\n"
# )

fluxes = pd.concat(
    [reference_solution.fluxes.rename("reference"), soln.fluxes.rename("result")],
    axis=1,
).join(pd.DataFrame(reactions, columns=["id", "reversed"]).set_index("id"), how="left")
fluxes["absDiff"] = (fluxes.reference - fluxes.result).abs()
print("Sum of absolute differences between fluxes: {}".format(fluxes.absDiff.sum()))
