# %%
import sys

sys.path.append("/usr/src/logmet")
from genesis.logicalmet import *
import re
from bs4 import BeautifulSoup


# helpers (from cobrapy)
pattern_from_sbml = re.compile(r"__(\d+)__")
species_name_pattern_y8 = re.compile(r"^(.+)\s\[(.+)\]$")


def _escape_non_alphanum(nonASCII: str) -> str:
    """converts a non alphanumeric character to a string representation of
    its ascii number"""
    return "__" + str(ord(nonASCII.group())) + "__"


def _number_to_chr(numberStr: str) -> str:
    """converts an ascii number to a character """
    return chr(int(numberStr.group(1)))


def logical_model_from_sbml(
    model_file: str, include_reverse: bool = False
) -> "LogicalMetabolicModel":

    # load model file
    with open(model_file, "r") as fi:
        b = BeautifulSoup(fi.read(), "xml")

    # Create model
    M = LogicalMetabolicModel()

    # Load compartments
    for c in b.find_all("compartment"):
        cid = c.get("id")
        M.compartments[cid] = LogicalCompartment(model=M, cid=cid, name=c.get("name"))

    # Load species as metabolites
    for s in b.find_all("species"):
        if b.find("model").get("id").startswith("M_yeastGEM_v8"):
            species_name, species_comp = species_name_pattern_y8.match(
                s.get("name")
            ).groups()
        else:
            species_name = s.get("name")

        # comp = list(filter(lambda c: c.id == s.get("compartment"), M.compartments))[0]
        comp = M.compartments.get(s.get("compartment"))
        s_chem_form = s.get("fbc:chemicalFormula")
        try:
            s_ids = list(
                map(
                    lambda li: li.get("rdf:resource"),
                    s.find("bqbiol:is").findAll("rdf:li"),
                )
            )
        except AttributeError:
            s_ids = []

        print

        try:
            existing_mets = list(
                filter(lambda m: m.species == species_name, M.metabolites)
            )
            m = existing_mets[0]
            assert len(existing_mets) == 1
            m.compartments.add(comp)
            m.sbml_ids[s.get("id")] = comp

            if m.chemical_formula is None:
                m.chemical_formula = s_chem_form

            try:
                m.identifiers.extend(s_ids)
            except AttributeError:
                m.identifiers = s_ids if s_ids else []

        except IndexError:
            M.metabolites.append(
                LogicalMetabolite(
                    model=M,
                    species=species_name,
                    chemical_formula=s_chem_form,
                    compartments=[comp],
                    sbml_ids={s.get("id"): comp},
                    identifiers=(s_ids if s_ids else []),
                )
            )

    # Load genes
    for g in b.find_all("fbc:geneProduct"):
        M.genes.append(
            LogicalGene(
                model=M,
                orf=g.get("fbc:label").replace("G_", ""),
                name=g.get("fbc:name"),
                sbml_id=g.get("fbc:id"),
            )
        )

    # Load reactions
    for r in b.find_all("reaction"):
        # Process GPR
        r_complexes = []
        gpr = r.find("fbc:geneProductAssociation")
        if gpr is None:
            pass
        else:
            # print(gpr)
            if gpr.find("fbc:or") is not None:
                gpr = gpr.find("fbc:or")

            for comp in filter(lambda child: child != "\n", gpr.children):
                # Make complex
                if comp.name == "and":
                    comp_genes = comp.find_all("fbc:geneProductRef")
                else:
                    comp_genes = [comp]

                complex = LogicalEnzymeComplex(
                    model=M,
                    genes=[
                        lg
                        for g in comp_genes
                        for lg in M.genes
                        if g.get("fbc:geneProduct") == lg.sbml_id
                    ],
                )
                # Check if exists: if yes - link it, if no - add it
                try:
                    r_complexes.append(
                        M.enzyme_complexes[
                            [c.word for c in M.enzyme_complexes].index(complex.word)
                        ]
                    )
                except ValueError:
                    M.enzyme_complexes.append(complex)
                    r_complexes.append(M.enzyme_complexes[-1])

        try:
            reactants = [
                (m.sbml_ids[sid], m)
                for sid in [
                    s.get("species")
                    for s in r.find("listOfReactants").find_all("speciesReference")
                ]
                for m in M.metabolites
                if sid in m.sbml_ids
            ]
        except AttributeError:
            reactants = []

        try:
            products = [
                (m.sbml_ids[sid], m)
                for sid in [
                    s.get("species")
                    for s in r.find("listOfProducts").find_all("speciesReference")
                ]
                for m in M.metabolites
                if sid in m.sbml_ids
            ]
        except AttributeError:
            products = []

        M.reactions.append(
            LogicalReaction(
                model=M,
                name=r.get("name"),
                reactants=reactants,
                products=products,
                enzymes=r_complexes,
                sbml_id=r.get("id"),
            )
        )
        if include_reverse and r.get("reversible") == "true":
            M.reactions.append(
                LogicalReaction(
                    model=M,
                    name=r.get("name"),
                    reactants=products,
                    products=reactants,
                    enzymes=r_complexes,
                    sbml_id=r.get("id"),
                    reverse=True,
                )
            )
    return M


def theory_files_from_logical_model(
    model: "LogicalMetabolicModel",
    root_directory: str,
    COMPARTMENTLESS: bool = False,
    GENE_ACTIVATIONS: bool = True,
    SPLIT_ESSENTIAL: bool = False,
) -> None:
    return
