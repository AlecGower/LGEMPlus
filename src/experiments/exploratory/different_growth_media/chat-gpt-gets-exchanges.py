import libsbml

# Load the SBML file
sbml_file = 'model-files/yeastGEM.xml'
reader = libsbml.SBMLReader()
document = reader.readSBML(sbml_file)
model = document.getModel()

# Dictionary of compound names and their corresponding identifiers
yeast8_identifiers = {
    "Oxygen": "https://identifiers.org/bigg.reaction/EX_o2",
    "NH3": "https://identifiers.org/bigg.reaction/EX_nh4",
    "SLF": "https://identifiers.org/bigg.reaction/EX_so4",
    "PI": "https://identifiers.org/bigg.reaction/EX_pi",
    "Glucose": "https://identifiers.org/bigg.reaction/EX_glc",
    "Glycerol": "https://identifiers.org/bigg.reaction/EX_glyc",
    "Lactate": "https://identifiers.org/bigg.reaction/EX_lac_L",
    "Galactose": "https://identifiers.org/bigg.reaction/EX_gal",
    "Raffinose": "https://identifiers.org/bigg.reaction/EX_raffin",
    "Ethanol": "https://identifiers.org/bigg.reaction/EXtoh",
    "Acetate": "https://identifiers.org/bigg.reaction/EX_ac",
    "Ergosterol": "https://identifiers.org/bigg.reaction/EXrgst",
    "Zymosterol": "https://identifiers.org/bigg.reaction/EX_zymst",
    "Palmitoleate": "https://identifiers.org/bigg.reaction/EX_palme",
    "Stearate": "https://identifiers.org/bigg.reaction/EX_stear",
    "Oleate": "https://identifiers.org/bigg.reaction/EX_ole",
    "Linoleate": "https://identifiers.org/bigg.reaction/EX_lnlc",
    "ALA": "https://identifiers.org/bigg.reaction/EX_ala_L",
    "ARG": "https://identifiers.org/bigg.reaction/EX_arg_L",
    "ASN": "https://identifiers.org/bigg.reaction/EX_asn_L",
    "ASP": "https://identifiers.org/bigg.reaction/EX_asp_L",
    "CYS": "https://identifiers.org/bigg.reaction/EX_cys_L",
    "GLU": "https://identifiers.org/bigg.reaction/EX_glu_L",
    "GLN": "https://identifiers.org/bigg.reaction/EX_gln_L",
    "GLY": "https://identifiers.org/bigg.reaction/EX_gly",
    "HIS": "https://identifiers.org/bigg.reaction/EX_his_L",
    "ILE": "https://identifiers.org/bigg.reaction/EX_ile_L",
    "LEU": "https://identifiers.org/bigg.reaction/EX_leu_L",
    "LYS": "https://identifiers.org/bigg.reaction/EX_lys_L",
    "MET": "https://identifiers.org/bigg.reaction/EX_met_L",
    "PHE": "https://identifiers.org/bigg.reaction/EX_phe_L",
    "PRO": "https://identifiers.org/bigg.reaction/EX_pro_L",
    "SER": "https://identifiers.org/bigg.reaction/EX_ser_L",
    "THR": "https://identifiers.org/bigg.reaction/EX_thr_L",
    "TRP": "https://identifiers.org/bigg.reaction/EX_trp_L",
    "TYR": "https://identifiers.org/bigg.reaction/EX_tyr_L",
    "VAL": "https://identifiers.org/bigg.reaction/EX_val_L",
    "Potassium": "https://identifiers.org/bigg.reaction/EX_k",
    "Sodium": "https://identifiers.org/bigg.reaction/EX_na1",
    "Biotin": "https://identifiers.org/bigg.reaction/EX_btn",
    "Choline": "https://identifiers.org/bigg.reaction/EX_chol",
    "Riboflavin": "https://identifiers.org/bigg.reaction/EX_ribflv",
    "Thiamine": "https://identifiers.org/bigg.reaction/EX_thm",
    "Inositol": "https://identifiers.org/bigg.reaction/EX_inost",
    "Thymidine": "https://identifiers.org/bigg.reaction/EX_thymd",
    "Nicotinate": "https://identifiers.org/bigg.reaction/EX_nac",
    "4-Aminobenzoate": "https://identifiers.org/bigg.reaction/EX_4abz",
    "(R)-Pantothenate": "https://identifiers.org/bigg.reaction/EX_pnto_R",
    "Pyridoxine": "https://identifiers.org/bigg.reaction/EX_pydxn",
    "Uracile": "https://identifiers.org/bigg.reaction/EX_ura",
    "Adenine": "https://identifiers.org/bigg.reaction/EX_adn",
    "Antimycin A": "https://identifiers.org/bigg.reaction/EX_antim"
}

# Function to look up reactions for each compound
def get_reactions_for_compounds(compound_identifiers, model):
    compound_reactions = {}
    for compound, identifier in compound_identifiers.items():
        for i in range(model.getNumReactions()):
            reaction = model.getReaction(i)
            annotation = reaction.getAnnotation()
            if annotation and identifier in annotation.toXMLString():
                compound_reactions[compound] = reaction.getId()
                break
            if compound in compound_reactions:
                break
        else:
            compound_reactions[compound] = "Not found"
    return compound_reactions

# Perform the lookup
compound_reactions = get_reactions_for_compounds(yeast8_identifiers, model)

# Print the results
for compound, reaction in compound_reactions.items():
    print(f'"{compound}";"{reaction[2:]}"')