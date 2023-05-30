# Import packages
import os, sys, subprocess
import multiprocessing as mp
from functools import partial

sys.path.append(os.getcwd())

from tqdm.auto import tqdm

from genesis.io import logical_model_from_sbml
from genesis.logicalmet import LogicalClause, LogicalLiteral
from pathlib import Path
from datetime import datetime

theory_root = Path(sys.argv[-1])
results_root = list(theory_root.parts)
results_root[results_root.index("theories")] = "results"
results_root = Path("/".join(results_root))
results_root.mkdir(parents=True, exist_ok=True)

stream = os.popen('grep "Base deletants:" {}'.format(theory_root / "info.txt"))
base_deletants = stream.read()[len("Base deletants:") :].strip().split(" ")
genes = []
with open(theory_root / "gene_list.txt") as fi:
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
def run_imap_multiprocessing(func, argument_list, num_processes):

    with mp.Pool(processes=num_processes) as pool:

        result_list_tqdm = []
        for result in tqdm(
            pool.imap(func=func, iterable=argument_list), total=len(argument_list)
        ):
            result_list_tqdm.append(result)

        return result_list_tqdm


if __name__ == "__main__":
    essentiality_results = run_imap_multiprocessing(essentiality_map, genes, mp.cpu_count())
        
# essentiality_results = # list of tuples '(orf, name, essential)'

## Output
with open(results_root / "single_gene_deletions.txt", "w") as fo:
    print("orf", "name", "essential", sep="\t", file=fo)
    for orf, name, essential in essentiality_results:
        print(orf, name, essential, sep="\t", file=fo)
