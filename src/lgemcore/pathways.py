# Import packages
import os, subprocess
from argparse import ArgumentParser
from pathlib import Path
import multiprocessing as mp
from functools import partial
from tqdm.auto import tqdm
from typing import Optional, Iterable, Union
import re
import yaml

from random import sample, choice

NUM_PROCESSES = mp.cpu_count() - 2
MY_ENV = os.environ.copy()
ORF_NAME_EXTRACT = re.compile(r"g_([A-Z0-9_]+)_([A-Z0-9_]+)__in_genome")
REACTION_NAME_EXTRACT = re.compile(r"(r_[0-9]+)_([A-z0-9_]+(?:_reverse)?)\)")
METABOLITE_EXTRACT = re.compile(r"[^~] met\(([A-z0-9_]+),([a-z_]+)\)")

def theory_to_genome(theory_root, gene_list: Union[str, os.PathLike] = "gene_list.txt"):
    theory_root = Path(theory_root)
    try:
        genelst = Path(gene_list)
        with open(gene_list) as fi:
            pass
    except FileNotFoundError:
        gene_list = theory_root / gene_list

    with open(theory_root / "info.txt") as fi:
       info = yaml.safe_load(fi)
    #    for ln in fi:
    #        if ln.startswith("Base deletants:"):
    #            base_deletants = ln[len("Base deletants:") :].strip().split(" ")
    #        elif ln.startswith("Base GEM ID:"):
    #            model_id = ln[len("Base GEM ID:") :].strip()
    base_deletants = info["Base deletants"]
    model_id = info["Base GEM ID"]

    genes = []
    with open(gene_list) as fi:
        for ln in fi:
            orf, name = ln.rstrip().split()
            if orf not in base_deletants and name not in base_deletants:
                genes.append((orf, name))

    return genes

def calculate_pathway(
    ko_genes: Iterable[tuple],
    ko_reactions: Iterable[str],
    theory_root: Union[str, os.PathLike],
    additional_problem_files: Optional[Iterable] = None,
    pathway_format: Optional[str] = "genes",
):
    if not isinstance(theory_root, os.PathLike):
        theory_root = Path(theory_root)

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

    reactions = subprocess.run(
        [
            "sed",
            ";".join(
                [
                    f"/{rxn}/ s/axiom, ( reaction_toggle(r_/axiom, ( ~reaction_toggle(r_/; "
                    for rxn in ko_reactions ### Turn Off Reaction
                ]
            ),
            str(theory_root / "reactions.p"),
        ],
        check=True,
        capture_output=True,
    )
    # print(reactions.stdout.rstrip().decode("utf-8")[:2500])
    theory = subprocess.run(
        [
            "cat",
            "-",
            *[
                str(theory_root / fname)
                for fname in [
                    # "reactions.p",
                    "compound_synonyms.p",
                    "ubiquitous_compounds.p",
                    "media_compounds.p",
                    "abduced_extra_compounds.p",
                    "query.p",
                ]
            ],
            *additional_problem_files,
        ],
        input=genome.stdout + reactions.stdout,
        capture_output=True,
    )
    # print("Theory constructed...")
    # print(theory.stdout.rstrip().decode("utf-8"))
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
    # print(proof.stdout.rstrip().decode("utf-8"))
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
        # print(proof.stdout.rstrip().decode("utf-8"))

        results = []

        if "genes" in pathway_format:
            activations = subprocess.run(
                ["grep", "__in_genome"],
                input=proof.stdout,
                capture_output=True,
            )
            results.extend([
                ("gene", ORF_NAME_EXTRACT.findall(g)[0])
                for g in activations.stdout.rstrip().decode("utf-8").split("\n")
            ])
        if "reactions" in pathway_format:
            reactions = subprocess.run(
                ["grep", "_toggled"],
                input=proof.stdout,
                capture_output=True,
            )
            results.extend([
                ("reaction", REACTION_NAME_EXTRACT.findall(rxn)[0])
                for rxn in reactions.stdout.rstrip().decode("utf-8").split("\n")
            ])
        if "metabolites" in pathway_format:
            metabolites = subprocess.run(
                ["grep", "met("],
                input=proof.stdout,
                capture_output=True
            )
            results.extend([
                ("metabolite", t) 
                for m in metabolites.stdout.rstrip().decode("utf-8").split("\n")
                for t in METABOLITE_EXTRACT.findall(m)
            ])
        # print("Pathway extracted...")

        return (
            ko_genes,
            ko_reactions,
            results,
            deletion_result.stdout.rstrip().decode("utf-8") == "0",
        )
    else:
        return (ko_genes, ko_reactions, [], deletion_result.stdout.rstrip().decode("utf-8") == "0")


def disruption_reaction(
    ko_genes,
    potential_ko_reactions: Iterable[str],
    min_reactions: int,
    max_reactions: int,
    theory_root: Union[str, os.PathLike],
    additional_problem_files: Optional[Iterable] = None,
    pathway_format: Optional[str] = "genes",
):
    return calculate_pathway(
        ko_genes=ko_genes,
        ko_reactions=sample(
            potential_ko_reactions,
            choice(range(min_reactions, max_reactions + 1))
        ),
        theory_root=theory_root,
        additional_problem_files=additional_problem_files,
        pathway_format=pathway_format,
    )


def results_to_diff_df(results_list, wt_filter, allMets):
    rFilteredBinary = [[m in r[2] for m in allMets] for r in results_list]
    rdf = pd.DataFrame(data=rFilteredBinary, columns=list(allMets))
    rDiff = []
    for i, row in enumerate(rdf[rdf.any(axis=1)].iterrows()):
        rDiff.append(np.array(wt_filter) ^ row[1])
    return pd.DataFrame(data=rDiff, columns=list(allMets))
