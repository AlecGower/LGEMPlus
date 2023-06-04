# Import packages
import os, subprocess
from argparse import ArgumentParser
from pathlib import Path
import multiprocessing as mp
from functools import partial
from tqdm.auto import tqdm

ap = ArgumentParser(
    description="Conducts single-gene deletions using iProver given a directory containing theory files."
)
ap.add_argument("theory_root")
ap.add_argument(
    "-g",
    "--gene_list",
    dest="gene_list",
    default="gene_list.txt",
    help="Absolute path or path relative to theory directory that contains list"
    + " of genes to be evaluated. Default = 'gene_list.txt'.",
)

arguments = ap.parse_args()
theory_root = Path(arguments.theory_root)
try:
    gene_list = Path(arguments.gene_list)
    with open(gene_list) as fi:
        pass
except FileNotFoundError:
    gene_list = theory_root / arguments.gene_list

# Create directory for results
results_root = list(theory_root.parts)
results_root[results_root.index("theories")] = "results"
results_root = Path("/".join(results_root))
results_root.mkdir(parents=True, exist_ok=True)

with open(theory_root / "info.txt") as fi:
    for ln in fi:
        if ln.startswith("Base deletants:"):
            base_deletants = ln[len("Base deletants:") :].strip().split(" ")
        elif ln.startswith("Base GEM ID:"):
            model_id = ln[len("Base GEM ID:") :].strip()

genes = []
with open(gene_list) as fi:
    for ln in fi:
        orf, name = ln.rstrip().split("\t")
        if orf not in base_deletants and name not in base_deletants:
            genes.append((orf, name))


def calculate_gene_essentiality(gene: tuple, theory_root: os.PathLike):
    # Calculate deletion
    orf, name = gene
    proof = subprocess.run(
        [
            "bash",
            "./helpers/gene_knockout_simple.sh",
            orf,
            str(theory_root / "genes.p"),
            *[
                str(theory_root / fname)
                for fname in [
                    "reactions.p",
                    "ubiquitous_compounds.p",
                    "media_compounds.p",
                    "abduced_extra_compounds.p",
                    "query.p",
                ]
            ],
        ],
        check=True,
        capture_output=True,
    )
    deletion_result = subprocess.run(
        ["grep", "-c", "% SZS status Satisfiable"],
        input=proof.stdout,
        capture_output=True,
    )
    return orf, name, deletion_result.stdout.rstrip().decode("utf-8")


## Multiprocessing bit
pool = mp.Pool(mp.cpu_count())
essentiality_results = []
essentiality_map = partial(calculate_gene_essentiality, theory_root=theory_root)

## Credit for this function to leimao@github
def run_imap_multiprocessing(func, argument_list, num_processes, **tqdm_kwargs):

    with mp.Pool(processes=num_processes) as pool:

        result_list_tqdm = []
        for result in tqdm(
            pool.imap(func=func, iterable=argument_list),
            total=len(argument_list),
            **tqdm_kwargs
        ):
            result_list_tqdm.append(result)

        return result_list_tqdm


if __name__ == "__main__":
    essentiality_results = run_imap_multiprocessing(
        essentiality_map, genes, mp.cpu_count(), desc=model_id
    )

# essentiality_results = # list of tuples '(orf, name, essential)'

## Output
with open(results_root / "single_gene_deletions.txt", "w") as fo:
    # print("orf", "name", "essential", sep="\t", file=fo)
    print("orf", "name", "growth", sep="\t", file=fo)
    for orf, name, essential in essentiality_results:
        # print(orf, name, essential, sep="\t", file=fo)
        print(orf, name, 1 - int(essential), sep="\t", file=fo)
