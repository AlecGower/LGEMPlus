# Import packages
import os, sys

sys.path.append(os.getcwd())
from genesis.io import logical_model_from_sbml
from genesis.logicalmet import LogicalClause, LogicalLiteral
from pathlib import Path
from datetime import datetime

gen_time = datetime.utcnow()

# Set parameters
model_xml = "ModelFiles/yeastGEM.xml"
query_compounds = "ModelFiles/essentialAberyeastGEM.txt"
ubiquitous_compounds = "ModelFiles/ubiquitousAberyeastGEM.txt"
media_compounds = "ModelFiles/ynbAberyeastGEM.txt"
media_name = "ynb"

model_xml = "ModelFiles/iMM904.xml"
query_compounds = "ModelFiles/essentialAberiMM904.txt"
ubiquitous_compounds = "ModelFiles/ubiquitousAberiMM904.txt"
media_compounds = "ModelFiles/ynbAberiMM904.txt"
media_name = "ynb"

model_xml = "ModelFiles/iND750.xml"
query_compounds = "ModelFiles/essentialAberiND750.txt"
ubiquitous_compounds = "ModelFiles/ubiquitousAberiND750.txt"
media_compounds = "ModelFiles/ynbAberiND750.txt"
media_name = "ynb"

base_knocked_out = "HIS3 LEU2 LYS2 MET17 URA3".split(" ")

GENE_ACTIVATIONS = True
REVERSE_REACTIONS = True
COMPARTMENTLESS = False

# Load model
lm = logical_model_from_sbml(model_file=model_xml, include_reverse=REVERSE_REACTIONS)
base_knocked_out = [
    g
    for g in lm.genes
    if g.name.upper() in base_knocked_out or g.orf.upper() in base_knocked_out
]

# Load lists of query and ubiquitous compounds
compounds = {"query": [], "ubiquitous": [], "media": []}
with open(query_compounds, "r") as fi:
    for (lnno, row) in enumerate(fi):
        if lnno > 0:
            try:
                compound = row.split("\t")[3].rstrip()
            except IndexError:
                continue
            if compound != "":
                compounds["query"].append(compound)

with open(ubiquitous_compounds, "r") as fi:
    for (lnno, row) in enumerate(fi):
        if lnno > 0:
            try:
                compound = row.split("\t")[3].rstrip()
            except IndexError:
                continue
            if compound != "":
                compounds["ubiquitous"].append(compound)

with open(media_compounds, "r") as fi:
    for (lnno, row) in enumerate(fi):
        if lnno > 0:
            try:
                compound = row.split("\t")[3].rstrip()
            except IndexError:
                continue
            if compound != "":
                compounds["media"].append(compound)

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
            " ".join(sorted([g.name for g in base_knocked_out])),
            media_name,
            REVERSE_REACTIONS,
            COMPARTMENTLESS,
        )
    )

# Write out reactions
with open(theory_root / "reactions.p", "w") as fo:
    lm.print_reactions(compartmentless=COMPARTMENTLESS, out_fo=fo)

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
                        compartment=lm.compartments.get("e"), presence_type="MEDIA",
                    )
                )
            else:
                fo.write(m.cnf_presence(compartment=False, presence_type="MEDIA",))

# Write out ubiquitous metabolites
with open(theory_root / "ubiquitous_compounds.p", "w") as fo:
    for m in lm.metabolites:
        if m.species in compounds["ubiquitous"]:
            if not COMPARTMENTLESS:
                for c in lm.compartments.values():
                    fo.write(m.cnf_presence(compartment=c, presence_type="UBIQUITOUS",))
            else:
                fo.write(m.cnf_presence(compartment=False, presence_type="UBIQUITOUS"))

# Write query
query = LogicalClause("query_metabolites", literals=[], type="negated_conjecture",)
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
stream = os.popen(
    "|".join(
        [
            '$IPROVER_HOME"/iproveropt" $(ls {}/*.p)'.format(theory_root),
            'grep "{~("',
            'grep -v "rxn("',
            'grep -v "enz("',
            'grep -v "pro("',
            'grep -v "gn("',
        ]
    )
)
hypotheses = [h[h.find("{") + 1 : h.find("}")].split(";") for h in stream.readlines()]
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
