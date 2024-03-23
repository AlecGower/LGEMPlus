# Import packages
import os, subprocess
from argparse import ArgumentParser
from pathlib import Path
import multiprocessing as mp
from functools import partial
from tqdm.auto import tqdm
from typing import Optional, Iterable
import re
from pprint import pprint
import pickle
from datetime import datetime, timezone

start_time = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H%M%SZ")

NUM_PROCESSES = mp.cpu_count() - 2
MY_ENV = os.environ.copy()

ap = ArgumentParser(
    description="Conducts multiple-gene deletions using iProver given a directory containing theory files, and outputs activated genes."
)
ap.add_argument("theory_root")
ap.add_argument(
    "-p", "--additional_problem_file", action="append", dest="additional_problem_files"
)
ap.add_argument("-g", "--gene_ko", action="append", nargs="*")

arguments = ap.parse_args()
theory_root = Path(arguments.theory_root)
ko_genes = arguments.gene_ko
try:
    ko_genes = [g[0] for g in ko_genes]  # fix this later
except TypeError:
    ko_genes = []

with open(theory_root / "info.txt") as fi:
    for ln in fi:
        if ln.startswith("Base deletants:"):
            base_deletants = ln[len("Base deletants:") :].strip().split(" ")
        elif ln.startswith("Base GEM ID:"):
            model_id = ln[len("Base GEM ID:") :].strip()

orf_name_extract = re.compile(r"g_([A-Z0-9_]+)_([A-Z0-9_]+)__in_genome")
genes = []


def calculate_pathway(
    ko_genes: Iterable[tuple],
    theory_root: os.PathLike,
    additional_problem_files: Optional[Iterable] = None,
):
    if additional_problem_files is None:
        additional_problem_files = []

    # print(ko_genes)
    # print(
    #     ";".join(
    #         [
    #             f"/{g[0]}/ s/gn(g_/~gn(g_/;/{g[0]}/ s/_in_genome/_deletion/"
    #             for g in ko_genes
    #         ]
    #     )
    # )
    # print(
    #     *[
    #         "sed",
    #         ";".join(
    #             [
    #                 f"/{g[0]}/ s/gn(g_/~gn(g_/;/{g[0]}/ s/_in_genome/_deletion/"
    #                 for g in ko_genes
    #             ]
    #         ),
    #         str(theory_root / "genes.p"),
    #     ],
    #     sep=" ",
    # )
    genome = subprocess.run(
        [
            "sed",
            ";".join(
                [
                    f"/{g}/ s/gn(g_/~gn(g_/;/{g}/ s/_in_genome/_deletion/"
                    for g in ko_genes
                ]
            ),
            str(theory_root / "genes.p"),
        ],
        check=True,
        capture_output=True,
    )
    # print("Genome constructed...")
    theory = subprocess.run(
        [
            "cat",
            "-",
            *[
                str(theory_root / fname)
                for fname in [
                    "reactions.p",
                    "compound_synonyms.p",
                    "ubiquitous_compounds.p",
                    "media_compounds.p",
                    "abduced_extra_compounds.p",
                    # "query.p",
                ]
            ],
            *additional_problem_files,
        ],
        input=genome.stdout,
        capture_output=True,
    )
    # print("Theory constructed...")
    proof = subprocess.run(
        [
            # "/usr/bin/time",
            # "gtime",
            # # "-f%U\t%S\t%P\t%e",
            # "bash",
            # "./helpers/gene_knockout_simple.sh",
            MY_ENV["IPROVER_HOME"] + "/iproveropt",
            # "/Users/alexander/workspace/iprover/iproveropt",
            "--stdin",
            "true",
            "--proof_out",
            "false",
            "--sat_out_model",
            "none",
            "--sat_out_clauses",
            "false",
            "--sat_out_model",
            "none",
        ],
        input=theory.stdout,
        capture_output=True,
    )
    # print("Simple proof found...")
    # print(proof.stdout)
    deletion_result = subprocess.run(
        ["grep", "-c", "% SZS status Satisfiable"],
        input=proof.stdout,
        capture_output=True,
    )
    # print("Result calculated...")
    if deletion_result.stdout.rstrip().decode("utf-8") == "0":
        proof = subprocess.run(
            [
                # "/usr/bin/time",
                # "gtime",
                # "-f%U\t%S\t%P\t%e",
                # "bash",
                # "./helpers/gene_knockout_simple.sh",
                MY_ENV["IPROVER_HOME"] + "/iproveropt",
                # "/Users/alexander/workspace/iprover/iproveropt",
                "--stdin",
                "true",
            ],
            input=theory.stdout,
            capture_output=True,
        )
        # print("Proof found...")
        activations = subprocess.run(
            ["grep", "__in_genome"],
            input=proof.stdout,
            capture_output=True,
        )
        # print("Pathway extracted...")

        return (
            ko_genes,
            [
                orf_name_extract.findall(g)[0]
                for g in activations.stdout.rstrip().decode("utf-8").split("\n")
            ],
            deletion_result.stdout.rstrip().decode("utf-8") == "0",
        )
    else:
        return (ko_genes, [], deletion_result.stdout.rstrip().decode("utf-8") == "0")


# results = calculate_pathway(
#     ko_genes=ko_genes,
#     theory_root=theory_root,
#     additional_problem_files=arguments.additional_problem_files,
# )


## Credit for this function to leimao@github
def run_imap_multiprocessing(func, argument_list, num_processes, **tqdm_kwargs):
    with mp.Pool(processes=num_processes) as pool:
        result_list_tqdm = []
        for result in tqdm(
            pool.imap(func=func, iterable=argument_list),
            total=len(argument_list),
            **tqdm_kwargs,
        ):
            result_list_tqdm.append(result)

        return result_list_tqdm


# Final function
def strain_to_pathway(g, ko_genes, theory_root, additional_problem_files):
    strain = ko_genes + [g[0]]
    new = calculate_pathway(
        ko_genes=strain,
        theory_root=theory_root,
        additional_problem_files=additional_problem_files,
    )
    return new


## Multiprocessing bit
pathway_map_base = partial(
    calculate_pathway,
    theory_root=theory_root,
    additional_problem_files=arguments.additional_problem_files,
)


def print_with_lock(s: str, l: mp.Lock):
    l.acquire()
    try:
        print(s)
    finally:
        l.release()


# Loop


def worker(ko_genes, qnext, results_dict, results_lock, printl):
    # ko_genes = q.get()
    results_key = frozenset(ko_genes)

    with results_lock:
        new_strain = results_key not in results_dict

    if new_strain:
        ko_orfs = [g[0] for g in ko_genes]
        # printl(f"Calculating pathway for strain: {ko_genes}")
        result = pathway_map_base(ko_orfs)
        if result[2]:
            # printl(f"VIABLE\tStrain:{result[0]}\tPathway Length: {len(result[1]):>3}")
            for g in result[1]:
                # qnext.put(ko_genes + [g])
                qnext.append(ko_genes + [g])
            # q.task_done()
        else:
            pass
            # printl(f"LETHAL\tStrain:{result[0]}")
        with results_lock:
            results_dict[results_key] = {
                "result": result,
                "pathway_hash": hash(frozenset(result[1])),
            }


if __name__ == "__main__":
    print(f"Starting multiprocessing with {NUM_PROCESSES} processes...")
    mp.set_start_method("spawn")
    with mp.Manager() as manager:
        # q = manager.Queue()
        q = manager.list()
        # qnext = manager.Queue()
        qnext = manager.list()
        results_dict = manager.dict()

        printlock = manager.Lock()
        results_lock = manager.Lock()
        printl = partial(print_with_lock, l=printlock)

        # lethal_combos = []
        # viable_pathways_untested = []
        # viable_pathways = []
        # initial = calculate_pathway(
        #     ko_genes=ko_genes,
        #     theory_root=theory_root,
        #     additional_problem_files=arguments.additional_problem_files,
        # )

        # qnext.put(ko_genes)
        qnext.append(ko_genes)
        # printl(f"q:{q}\tqnext{qnext}")
        # while not qnext.empty():
        ko_count = 0
        while len(qnext) > 0:
            q = qnext
            # qnext = manager.Queue()
            qnext = manager.list()
            # printl(f"q: {q}\tqnext: {qnext}")
            workern = partial(
                worker,
                qnext=qnext,
                results_dict=results_dict,
                results_lock=results_lock,
                printl=printl,
            )
            run_imap_multiprocessing(
                workern,
                q,
                num_processes=NUM_PROCESSES,
                desc=f"Strain KO Count - {ko_count}",
            )

            with open(f"experiments/results/glycolysis_pathways_{start_time}.pickle", "wb") as f:
                # Pickle the 'data' dictionary using the highest protocol available.
                pickle.dump(results_dict, f, pickle.HIGHEST_PROTOCOL)
            # print(f"Strains up to {ko_count} knockouts calculated. Next level")
            ko_count += 1

            # pprint(dict(results_dict))