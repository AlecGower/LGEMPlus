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


def worker(ko_genes, qnext, printl):
    # ko_genes = q.get()
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
    return result


if __name__ == "__main__":
    print(f"Starting multiprocessing with {NUM_PROCESSES} processes...")
    mp.set_start_method("spawn")
    with mp.Manager() as manager:
        # q = manager.Queue()
        q = manager.list()
        # qnext = manager.Queue()
        qnext = manager.list()
        results = manager.list()

        l = manager.Lock()
        printl = partial(print_with_lock, l=l)

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
            workern = partial(worker, qnext=qnext, printl=printl)
            results.extend(
                run_imap_multiprocessing(
                    workern,
                    q,
                    num_processes=NUM_PROCESSES,
                    desc=f"Strain KO Count - {ko_count}",
                )
            )
            # print(f"Strains up to {ko_count} knockouts calculated. Next level")
            ko_count += 1

        print(len(results))
        pprint(results[:5])


# while len(viable_pathways_untested) > 0:
#     print(
#         "{} pathways tested, {} pathways in queue. {} lethal mutants found.".format(
#             len(viable_pathways),
#             len(viable_pathways_untested),
#             len(lethal_combos),
#         )
#     )
#     ko_genes, pathway = viable_pathways_untested.pop(0)
#     pathway_map = partial(pathway_map_base, ko_genes=ko_genes)
#     print(pathway_map(pathway[0]))
#     if __name__ == "__main__":
#         pathway_results = run_imap_multiprocessing(
#             pathway_map, pathway, mp.cpu_count(), desc=f"Strain{ko_genes}", leave=False
#         )
#     # for g in pathway:
#     #     strain = ko_genes + [g[0]]
#     #     new = calculate_pathway(
#     #         ko_genes=strain,
#     #         theory_root=theory_root,
#     #         additional_problem_files=arguments.additional_problem_files,
#     #     )
#     # print(pathway_results)

#     # Add results to lists
#     lethal_combos.extend([new[0] for new in pathway_results if not new[2]])
#     viable_pathways_untested.extend([new[:2] for new in pathway_results if not new[2]])
#     viable_pathways.append((ko_genes, pathway))
#     print("Lethal:")
#     pprint(lethal_combos)
#     print("Untested:")
#     pprint(viable_pathways_untested)
#     print("Tested:")
#     pprint(viable_pathways)


# pprint(lethal_combos)
# pprint(viable_pathways)


# pathway_results = # list of tuples '(orf, name, essential, time)'

# ## Output
# print("orf", "name", "activatedorf", "activatedname", sep="\t")
# for orf, name, l in pathway_results:
#     for actorf, actname in l:
#         print(orf, name, actorf, actname, sep="\t")
