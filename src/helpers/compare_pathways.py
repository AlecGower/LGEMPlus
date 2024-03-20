# Import packages
from argparse import ArgumentParser
from pathlib import Path
import pandas as pd

ap = ArgumentParser(
    description="Compares pathways given a directory containing results files, and outputs ordered list of divergent pathways."
)
ap.add_argument("results_root")

arguments = ap.parse_args()
results_root = Path(arguments.results_root)

with open(results_root / "wt_pathway.txt") as fi:
    wt = list(map(str.strip, fi.readlines()))
wt = set([g.split()[0] for g in wt])

pwys = pd.read_csv(results_root / "single_gene_deletions_pathways.txt", delimiter="\t")
pwdict = {}
for row in pwys.iterrows():
    if not row[1]["orf"] in pwdict:
        pwdict[row[1]["orf"]] = {"name": row[1]["name"], "pathway": []}
    pwdict[row[1]["orf"]]["pathway"].append(row[1]["activatedorf"])
pwdict = {
    k: {"name": v["name"], "pathway": set(v["pathway"])} for k, v in pwdict.items()
}

for k, d in pwdict.items():
    d["wt_overlap"] = d["pathway"] & wt
    d["pathway - wt"] = d["pathway"] - wt
    d["wt - pathway"] = wt - d["pathway"]

print(
    "orf",
    "name",
    "pathway - wt",
    "wt - pathway",
    "pathway",
    sep="\t",
)
for orf in filter(
    lambda k: len(pwdict[k]["pathway - wt"]) + len(pwdict[k]["wt - pathway"]) > 0,
    sorted(
        pwdict.keys(),
        key=lambda k: len(pwdict[k]["pathway - wt"]) + len(pwdict[k]["wt - pathway"]),
        reverse=True,
    ),
):
    print(
        orf,
        pwdict[orf]["name"],
        len(pwdict[orf]["pathway - wt"]),
        len(pwdict[orf]["wt - pathway"]),
        len(pwdict[orf]["pathway"]),
        sep="\t",
    )
