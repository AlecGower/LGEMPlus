# Import packages
import os, sys, subprocess

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
with open(results_root / "single_gene_deletions.txt", "w") as fo:
    print("orf", "name", "essential", sep="\t", file=fo)
    n_genes = int(
        os.popen("wc {}".format(theory_root / "gene_list.txt")).read().split()[0]
    )
    with open(theory_root / "gene_list.txt") as fi:
        for ln in tqdm(fi, total=n_genes):
            orf, name = ln.rstrip().split("\t")
            if orf not in base_deletants and name not in base_deletants:
                # Calculate deletion
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
                print(
                    orf,
                    name,
                    deletion_result.stdout.rstrip().decode("utf-8"),
                    sep="\t",
                    file=fo,
                )

