# Import packages
import os, sys
from typing import Iterable

sys.path.append(os.getcwd())
from lgemcore.io import logical_model_from_sbml
from lgemcore.logicalmet import LogicalClause, LogicalLiteral
from helpers.abduce_hypotheses import abduce_hypotheses
from pathlib import Path
from datetime import datetime

gen_time = datetime.utcnow()

model = sys.argv[-2]
media_name = sys.argv[-1]

model_xml = "model-files/{}.xml".format(model)
query_compounds = "model-files/essential-compounds-{}.tsv".format(model)
ubiquitous_compounds = "model-files/ubiquitous-compounds-{}.tsv".format(model)
# media_compounds = "model-files/ynb-compounds-{}.tsv".format(model)
# media_name = "ynb"
# media_name = "scgal"
media_compounds = "model-files/{}-compounds-{}.tsv".format(media_name, model)

# # base_knocked_out = "HIS3 LEU2 LYS2 MET17 URA3".split(" ")
# base_knocked_out = "YOR202W YCL018W YBR115C YLR303W YEL021W".split(" ") # BY4743
# base_knocked_out = "HIS3 LEU2 MET17 URA3".split(" ")
# base_knocked_out = "YOR202W YCL018W YLR303W YEL021W".split(" ")  # BY4741
base_knocked_out = "YOR202W YCL018W YBR115C YEL021W".split(" ")  # BY4742

GENE_ACTIVATIONS = True
REVERSE_REACTIONS = True
COMPARTMENTLESS = False

# Load model
lm = logical_model_from_sbml(model_file=model_xml, include_reverse=REVERSE_REACTIONS)
base_knocked_out = [
    g
    for g in lm.genes
    if any(
        [
            g.name.upper() in base_knocked_out,
            g.orf.upper() in base_knocked_out,
            any([g.orf.endswith(nm) for nm in base_knocked_out]),
        ]
    )
]


# Load lists of query and ubiquitous compounds
def read_compounds(file: os.PathLike) -> list:
    compound_list = []
    with open(file, "r") as fi:
        for lnno, row in enumerate(fi):
            if lnno > 0:
                try:
                    compound = row.split("\t")[3].rstrip()
                except IndexError:
                    continue
                if compound != "":
                    compound_list.append(compound)
    return compound_list


compounds = {
    ky: read_compounds(eval(ky + "_compounds"))
    for ky in ["query", "ubiquitous", "media"]
}

# Load compound synonyms
compound_synonyms = {}
try:
    with open("model-files/compound-synonyms-{}.tsv".format(model)) as fi:
        for lnno, row in enumerate(fi):
            if lnno > 0:
                try:
                    compound1, compound2 = row.rstrip().split("\t")
                    if compound1 != "":
                        compound_synonyms[compound1] = compound2
                except ValueError:
                    pass
except FileNotFoundError:
    pass
# Create theories directory with following structure
# theories
# └── <GEM>
#     └── <date_time>
#         ├── info.txt (contains information about when and how the theories were created)
#         ├── reactions.p
#         ├── genes.p
#         ├── <media_name>.p
#         ├── ubiquitous.p
#         ├── abduced_extra_metabolites.p
#         └── queries
#             └── ess.p

# create theories directory if doesn't exist
# create GEM directory if doesn't exist
# create new <date_time> directory
theory_root = Path(
    "./experiments/theories/{}/{}".format(
        lm.model_id, gen_time.strftime("%Y-%m-%dT%H%M%SZ")
    )
)
theory_root.mkdir(parents=True)

# create info.txt file
template = """
Theory files information
========================

Base GEM ID:         {}
Base GEM filepath:   {}
Generation time:     {}

Base deletants:      {}
Growth medium:       {}

Reverse reactions included:         {}
Compartmentless metabolites:        {}
""".strip()
with open(theory_root / "info.txt", "w") as fo:
    fo.write(
        template.format(
            lm.model_id,
            model_xml,
            gen_time.strftime("%Y-%m-%dT%H:%M:%SZ"),
            " ".join(sorted([g.orf for g in base_knocked_out])),
            media_name,
            REVERSE_REACTIONS,
            COMPARTMENTLESS,
        )
    )

# Write out reactions
with open(theory_root / "reactions.p", "w") as fo:
    lm.print_reactions(
        compartmentless=COMPARTMENTLESS, skip_reaction_ids=["R_r_0767"], out_fo=fo
    )

# Write out compound synonyms
with open(theory_root / "compound_synonyms.p", "w") as fo:
    for c1, c2 in compound_synonyms.items():
        c1 = [m for m in lm.metabolites if m.species == c1][0]
        c2 = [m for m in lm.metabolites if m.species == c2][0]
        for c in lm.compartments.values():
            fo.write(
                "{}\n".format(
                    LogicalClause(
                        f"synonym_{c1.word}_{c2.word}_{c}",
                        literals=[
                            LogicalLiteral("met", [c1.word, c.word], negated=False),
                            LogicalLiteral("met", [c2.word, c.word], negated=True),
                        ],
                        type="axiom",
                    )
                )
            )
            fo.write(
                "{}\n".format(
                    LogicalClause(
                        f"synonym_{c2.word}_{c1.word}_{c}",
                        literals=[
                            LogicalLiteral("met", [c2.word, c.word], negated=False),
                            LogicalLiteral("met", [c1.word, c.word], negated=True),
                        ],
                        type="axiom",
                    )
                )
            )

# Write out gene activations
with open(theory_root / "genes.p", "w") as fo:
    lm.print_genes(knocked_out=base_knocked_out, out_fo=fo)

# Write out media
with open(theory_root / "media_compounds.p", "w") as fo:
    for m in lm.metabolites:
        if m.species in compounds["media"]:
            if not COMPARTMENTLESS:
                fo.write(
                    m.cnf_presence(
                        compartment=lm.compartments.get("e"),
                        presence_type="MEDIA",
                    )
                )
            else:
                fo.write(
                    m.cnf_presence(
                        compartment=False,
                        presence_type="MEDIA",
                    )
                )

# Write out ubiquitous metabolites
with open(theory_root / "ubiquitous_compounds.p", "w") as fo:
    for m in lm.metabolites:
        if m.species in compounds["ubiquitous"]:
            if not COMPARTMENTLESS:
                for c in lm.compartments.values():
                    fo.write(
                        m.cnf_presence(
                            compartment=c,
                            presence_type="UBIQUITOUS",
                        )
                    )
            else:
                fo.write(m.cnf_presence(compartment=False, presence_type="UBIQUITOUS"))

# Write query
query = LogicalClause(
    "query_metabolites",
    literals=[],
    type="negated_conjecture",
)
cytosol = lm.compartments.get("c")
for m in lm.metabolites:
    if m.species in compounds["query"]:
        if not COMPARTMENTLESS:
            if cytosol in m.compartments:
                compound_compartment = cytosol
            else:
                compound_compartment = list(m.compartments)[0]
            query.literals.append(
                LogicalLiteral("met", [m, compound_compartment.word], negated=True)
            )
        else:
            query.literals.append(
                LogicalLiteral("met", [m, "no_compartment"], negated=True)
            )


with open(theory_root / "query.p", "w") as fo:
    fo.write("{}".format(query))

# Abduce extra
hypotheses = abduce_hypotheses(
    theory_root,
    excluded_predicates=["rxn", "pro", "gn", "enz"],
    # excluded_compounds=["glucose", "galactose"],
    # excluded_compounds=["_glucose", "_ADP"],
    excluded_compounds=["_glucose"],
)

# Choose the hypothesis with the lowest number of additional compounds
selected = sorted(hypotheses, key=lambda h: len(h))[0]

# Write abduced extra to file
template = "cnf(abduced_extra_compound_{:02d}, axiom, {} ).\n"
with open(theory_root / "abduced_extra_compounds.p", "w") as fo:
    for i, s in enumerate(selected):
        fo.write(template.format(i, s[2:-1]))

# Write gene list to file for easier deletant
with open(theory_root / "gene_list.txt", "w") as fo:
    for gene in lm.genes:
        fo.write("{}\t{}\n".format(gene.orf, gene.name))

# Write KEGG synonyms to file
with open(theory_root / "kegg_synonyms.p", "w") as fo:
    for species in lm.metabolites:
        for compartment in lm.compartments.values():
            fo.write(species.cnf_kegg_synonym(compartment=compartment))

# Write KEGG synonyms to file
with open(theory_root / "kegg_synonyms.txt", "w") as fo:
    for species in lm.metabolites:
        fo.write("{}\t{}\n".format(species.word, species.kegg_word))
