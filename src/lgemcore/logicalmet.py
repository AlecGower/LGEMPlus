# %%
from ast import arguments
from genericpath import exists
from itertools import permutations
from lib2to3.pgen2 import literals
from os import name, PathLike
from pyexpat import model

# from ctypes import Union
import re
from typing import Dict, Iterable, Union, Optional, Tuple, TextIO
from copy import copy, deepcopy

from pyparsing import originalTextFor
from cobra.core import Metabolite, Gene, Reaction


# %%
class Constant:
    def __init__(self, word: str):
        self.word = word

    def __str__(self):
        return self.word


# %%
class LogicalMetabolicModel:
    def __init__(
        self,
        model_id: str = "unnamed_model",
        compartments: dict = {},
        metabolites: Iterable = [],
        genes: Iterable = [],
        enzyme_complexes: Optional[Iterable[Iterable["LogicalGene"]]] = None,
        reactions: Iterable = [],
    ) -> None:
        # self.MET_SPECIES_PATTERN = re.compile(r"^(?P<SPECIES>.+)\s\[[^\[]+\]$")
        # self.MET_NAME_PATTERN = re.compile(
        #     r"^(?P<SPECIES>.+)\s\[(?P<COMPARTMENT>[^\[]+)\]$"
        # )
        # self.MET_ID_PATTERN = re.compile(r"^(?P<ID>s_\d{4})\[\w+\]$")
        self.model_id = model_id
        self.clauses = set()

        self.compartments = {
            cid: LogicalCompartment(model=self, cid=cid, name=cname)
            for cid, cname in compartments.items()
        }

        self.metabolites = [
            LogicalMetabolite(model=self, cobra_metabolite=m) for m in metabolites
        ]
        self.genes = [LogicalGene(model=self, cobra_gene=g) for g in genes]

        if enzyme_complexes is not None:
            self.enzyme_complexes = [
                LogicalEnzymeComplex(model=self, genes=glist)
                for glist in enzyme_complexes
            ]
        else:
            self.enzyme_complexes = []

        self.reactions = []
        for r in reactions:
            self.reactions.append(LogicalReaction(model=self, cobra_reaction=r))
            if r.reversibility:
                self.reactions.append(
                    LogicalReaction(model=self, cobra_reaction=r, reverse=True)
                )

    def add_metabolite(self, metabolite: "Metabolite") -> None:
        self.metabolites.append(
            LogicalMetabolite(model=self, cobra_metabolite=metabolite)
        )
        pass

    def add_gene(self, gene: "Gene") -> None:
        pass

    def add_reaction(self, reaction: "Reaction") -> None:
        pass

    def get_gene(self, gene_orf: str):
        return [g for g in self.genes if g.orf == gene_orf][0]

    def print_reactions(
        self, compartmentless: bool = False, skip_reaction_ids: Iterable[str] = [], out_fo: Optional[TextIO] = None
    ):
        for r in self.reactions:
            if not r.sbml_id in skip_reaction_ids:
                print(r.cnf_lines(compartments=(not compartmentless)), file=out_fo)

    def print_genes(self, knocked_out: Iterable = [], out_fo: Optional[TextIO] = None):
        for g in self.genes:
            if g in knocked_out:
                print(
                    g.cnf_presence(present=False, presence_type="_deletion"),
                    file=out_fo,
                )
            else:
                print(
                    g.cnf_presence(present=True, presence_type="_in_genome"),
                    file=out_fo,
                )


# %%
class LogicalLiteral:
    def __init__(
        self, predicate_symbol: str, arguments: Iterable, negated: bool = False
    ) -> None:
        self.SKELETON = "{}{}({})"
        self.NEGATOR = "~"
        self.predicate_symbol = predicate_symbol
        self.arguments = arguments
        self.negated = negated

    def __str__(self) -> str:
        return self.SKELETON.format(
            self.NEGATOR if self.negated else "",
            self.predicate_symbol,
            ", ".join(map(str, self.arguments)),
        )

    def negate(self, inplace: bool = False) -> Union[None, "LogicalLiteral"]:
        if inplace:
            self.negated = not self.negated
        else:
            _ = copy(self)
            _.negated = not _.negated
            return _


class LogicalClause:
    def __init__(
        self, name: str, literals: Iterable["LogicalLiteral"] = [], type: str = "axiom"
    ) -> None:
        self.SKELETON = "cnf({}, {}, ( {} ))."
        self.name = name
        self.literals = list(literals)
        self.type = type

    def __str__(self) -> str:
        return self.SKELETON.format(
            self.name, self.type, " | ".join([l.__str__() for l in self.literals])
        )


class LogicalCompartment:
    def __init__(
        self,
        model: "LogicalMetabolicModel",
        cid: str,
        name: str,
        word: Optional[str] = None,
    ) -> None:
        self.SKELETON = "{}"
        self.NW = re.compile(r"\W")

        self.model = model
        self.id = cid
        self.name = name
        if word is not None:
            self.word = word
        else:
            self.word = "c_" + self.NW.sub("_", self.id)

    def __str__(self) -> str:
        return self.SKELETON.format(self.word)

    def __repr__(self) -> str:
        return "LogicalCompartment(cid='{}', name='{}', word='{}')".format(
            self.id, self.name, self.word
        )


# %%
class LogicalMetabolite:
    def __init__(
        self,
        model: "LogicalMetabolicModel",
        word: Optional[str] = None,
        species: Optional[str] = None,
        chemical_formula: Optional[str] = None,
        compartments: Optional[Iterable[LogicalCompartment]] = None,
        sbml_ids: Optional[dict] = None,
        identifiers: Iterable = [],
        cobra_metabolite: Optional["Metabolite"] = None,
        kegg_id: Optional[str] = None,
    ) -> None:

        self.model = model
        # TODO check if exists

        assert any(
            [
                (cobra_metabolite is not None),
                all([species is not None, compartments is not None]),
            ]
        )

        self.NW = re.compile(r"\W")

        self.sbml_ids = sbml_ids
        self._cobra_metabolite = cobra_metabolite

        if self._cobra_metabolite is not None:
            self.MET_SPECIES_PATTERN = re.compile(r"^(?P<SPECIES>[^\[]+)\s\[[^\[]+\]$")
            self._from_cobra_metabolite(self._cobra_metabolite)
        else:
            self.species = species
            self.chemical_formula = chemical_formula
            if hasattr(compartments, "__iter__"):
                self.compartments = set(compartments)
            else:
                self.compartments = set([compartments])
            self.identifiers = identifiers
            self.kegg_id = kegg_id

        if not word:
            self.word = "m" + "_" + self.NW.sub("_", self.species)
        else:
            self.word = word

        if kegg_id is not None:
            self.kegg_word = "kegg_{}".format(self.kegg_id)
        else:
            self.kegg_word = None

    def __str__(self) -> str:
        return self.word

    def __repr__(self) -> str:
        return "LogicalMetabolite(word='{}', species='{}', chemical_formula='{}', compartments={})".format(
            self.word, self.species, self.chemical_formula, self.compartments,
        )

    def _from_cobra_metabolite(self, cobra_metabolite: "Metabolite") -> None:
        self.MET_SPECIES_PATTERN = re.compile(r"^(?P<SPECIES>.+)\s\[[^\[]+\]$")
        self.species = self.MET_SPECIES_PATTERN.match(cobra_metabolite.name).group(1)
        self.compartments = set(
            [self.model.compartments.get(cid) for cid in [cobra_metabolite.compartment]]
        )
        self.identifiers = {
            "sbml.id": cobra_metabolite.id,
            "sbml.name": cobra_metabolite.name,
        }
        self.chemical_formula = cobra_metabolite.formula
        self.kegg_id = cobra_metabolite.annotation.get("kegg.compound")

    def as_literal(
        self,
        compartment: Union["LogicalCompartment", bool] = False,
        negated: bool = False,
    ) -> Union[None, "LogicalLiteral"]:
        if compartment in self.compartments:
            return LogicalLiteral(
                predicate_symbol="met",
                arguments=[self.word, compartment.word],
                negated=negated,
            )
        elif not compartment:
            return LogicalLiteral(
                predicate_symbol="met",
                arguments=[self.word, "no_compartment"],
                negated=negated,
            )
        else:
            return None

    def cnf_presence(
        self,
        compartment: Union["LogicalCompartment", bool] = False,
        present: bool = True,
        presence_type: str = "medium",
    ) -> str:
        if (compartment in self.compartments) or (not compartment):
            if not compartment:
                title_line = "%----- METABOLITE PRESENCE ({}) {}".format(
                    presence_type.upper(), self.species
                )
                clause_name = self.word + "_" + presence_type
            else:
                title_line = "%----- METABOLITE PRESENCE ({}) {} ({})".format(
                    presence_type.upper(), self.species, compartment.name
                )
                clause_name = self.word + "_" + presence_type + "_" + compartment.word

            return "\n".join(
                [
                    title_line,
                    str(
                        LogicalClause(
                            name=clause_name,
                            literals=[
                                self.as_literal(
                                    compartment=compartment, negated=(not present)
                                )
                            ],
                            type="axiom",
                        )
                    ),
                    "",
                ]
            )
        else:
            return ""

    def cnf_kegg_synonym(self, compartment: LogicalCompartment) -> str:
        if self.kegg_word is not None:
            return "\n".join(
                [
                    str(
                        LogicalClause(
                            name="synonym_{}_{}_{}_{}".format(
                                self.kegg_word, self.word, compartment.word, i
                            ),
                            literals=[
                                LogicalLiteral(
                                    predicate_symbol="met",
                                    arguments=[w1, compartment.word],
                                    negated=False,
                                ),
                                LogicalLiteral(
                                    predicate_symbol="met",
                                    arguments=[w2,  compartment.word],
                                    negated=True,
                                ),
                            ],
                        )
                    )
                    for i, (w1, w2) in enumerate(
                        permutations([self.word, self.kegg_word])
                    )
                ]
                + [""]
            )
        else:
            return ""


class LogicalGene:
    def __init__(
        self,
        model: "LogicalMetabolicModel",
        orf: Optional[str] = None,
        name: Optional[str] = None,
        sbml_id: Optional[str] = None,
        cobra_gene: Optional["Gene"] = None,
    ) -> None:

        self.model = model
        # TODO check if exists

        assert any(
            [(cobra_gene is not None), all([orf is not None, name is not None]),]
        )

        self.NW = re.compile(r"\W")
        if cobra_gene is not None:
            self._cobra_gene = cobra_gene
            self.orf = self._cobra_gene.id
            self.name = self._cobra_gene.name
        else:
            self.orf = orf
            self.name = name
            self.sbml_id = sbml_id

        self.word = (
            "g" + "_" + self.NW.sub("_", self.orf) + "_" + self.NW.sub("_", self.name)
        )

    def __str__(self) -> str:
        return self.word

    def as_literal(self, negated: bool = False) -> "LogicalLiteral":
        return LogicalLiteral(
            predicate_symbol="gn", arguments=[self.word], negated=negated
        )

    def cnf_presence(
        self, present: bool = True, presence_type: str = "in_genome"
    ) -> str:
        return "\n".join(
            [
                "%----- GENE PRESENCE ({})  {} ({})".format(
                    presence_type, self.orf, self.name
                ),
                str(
                    LogicalClause(
                        name=self.word + "_" + presence_type,
                        literals=[self.as_literal(negated=(not present))],
                        type="axiom",
                    )
                ),
                "",
            ]
        )


class LogicalReaction:
    def __init__(
        self,
        model: "LogicalMetabolicModel",
        name: Optional[str] = None,
        reactants: Optional[Iterable[Tuple[str, "Metabolite"]]] = None,
        products: Optional[Iterable[Tuple[str, "Metabolite"]]] = None,
        enzymes: Optional[Iterable["LogicalEnzymeComplex"]] = None,
        sbml_id: Optional[str] = None,
        cobra_reaction: Optional["Reaction"] = None,
        reverse: bool = False,
    ) -> None:
        self.model = model
        # TODO check if exists
        self.NW = re.compile(r"\W")
        self.REVERSE = reverse
        if cobra_reaction is not None:
            self._from_cobra_reaction(cobra_reaction)
            self._compile_inputs()
            self._compile_outputs()
        else:
            self.name = name
            self.inputs = reactants
            self.outputs = products
            self.sbml_id = sbml_id
            self.word = (
                "r_"
                + self.sbml_id
                + "_"
                + self.NW.sub("_", self.name)
                + ("_reverse" if self.REVERSE else "")
            )
            self._enzyme_base_word = "e_" + self.word
            if enzymes is not None:
                self.enzymes = list(enzymes)
            else:
                self.enzymes = []

        self.metabolites = set(self.inputs + self.outputs)
        self.compartments = set(
            self.model.compartments.get(m[0].id) for m in self.metabolites
        )

    def _from_cobra_reaction(self, cobra_reaction: "Reaction") -> None:
        self._cobra_reaction = cobra_reaction
        self.name = self._cobra_reaction.name
        if not self.REVERSE:
            self._inputs = self._cobra_reaction.reactants
            self._outputs = self._cobra_reaction.products
        else:
            self._inputs = self._cobra_reaction.products
            self._outputs = self._cobra_reaction.reactants

        self.word = (
            cobra_reaction.id
            + "_"
            + self.NW.sub("_", cobra_reaction.name)
            + ("_reverse" if self.REVERSE else "")
        )
        self._enzyme_base_word = "e_" + self.word
        self.enzymes = self._process_cobra_gene_reaction_rule(
            self._cobra_reaction.gene_reaction_rule
        )

    def __str__(self) -> str:
        return self.word

    def _compile_inputs(self) -> None:
        self.inputs = [
            (
                self.model.compartments.get(m.compartment),
                LogicalMetabolite(model=self.model, cobra_metabolite=m),
            )
            for m in self._inputs
        ]

    def _compile_outputs(self) -> None:
        self.outputs = [
            (
                self.model.compartments.get(m.compartment),
                LogicalMetabolite(model=self.model, cobra_metabolite=m),
            )
            for m in self._outputs
        ]

    def as_literal(self, negated: bool = False) -> "LogicalLiteral":
        return LogicalLiteral(
            predicate_symbol="rxn", arguments=[self.word], negated=negated
        )

    def _build_input_clause(self, compartments: bool = False) -> None:
        self._input_clause = LogicalClause(
            name=self.word + "_in",
            literals=map(
                lambda input: input[1].as_literal(
                    compartment=(input[0] if compartments else False), negated=True
                ),
                self.inputs,
            ),
            type="axiom",
        )
        self._input_clause.literals.insert(0, self.as_literal())

    def _build_output_clauses(self, compartments: bool = False) -> None:
        self._output_clauses = [
            LogicalClause(
                name=self.word + "_out_{}".format(i),
                literals=[
                    output[1].as_literal(
                        compartment=(output[0] if compartments else False)
                    ),
                    self.as_literal(negated=True),
                ],
                type="axiom",
            )
            for (i, output) in enumerate(self.outputs)
        ]

    def _build_activation_clauses(self) -> None:
        self.activation_word = "p_" + self.word
        # self._enzyme_base_word = "e_" + self.word
        # self.isoenzymes = self._process_cobra_gene_reaction_rule(
        #     self._cobra_reaction.gene_reaction_rule
        # )
        self._activation_clauses = [
            LogicalClause(
                name=self.activation_word + "_{}".format(i),
                literals=[
                    LogicalLiteral(
                        predicate_symbol="pro", arguments=[self.activation_word]
                    ),
                    e.as_literal(negated=True),
                ],
                type="axiom",
            )
            for (i, e) in enumerate(self.enzymes)
        ]
        ## Enzyme complex activation causes now can be elsewhere
        self._activation_clauses.extend(
            [
                LogicalClause(
                    name=e.word,
                    literals=[e.as_literal()]
                    + [g.as_literal(negated=True) for g in e.genes],
                    type="axiom",
                )
                for e in self.enzymes
            ]
        )

        if len(self._activation_clauses) > 0:
            self._input_clause.literals.append(
                LogicalLiteral(
                    predicate_symbol="pro",
                    arguments=[self.activation_word],
                    negated=True,
                )
            )

    def _process_cobra_gene_reaction_rule(
        self, rule_string: str
    ) -> Iterable[Iterable["LogicalGene"]]:
        if " or " in rule_string:
            complexes = rule_string.split(" or ")
            complexes = [
                LogicalEnzymeComplex(
                    model=self.model,
                    word=self._enzyme_base_word + "_{}".format(i),
                    genes=[
                        self.model.get_gene(g.replace("(", "").replace(")", "").strip())
                        for g in c.split(" and ")
                    ],
                )
                for (i, c) in enumerate(complexes)
            ]
            return complexes
        elif " and " in rule_string:
            return [
                LogicalEnzymeComplex(
                    model=self.model,
                    word=self._enzyme_base_word + "_{}".format(0),
                    genes=[
                        self.model.get_gene(g.replace("(", "").replace(")", "").strip())
                        for g in rule_string.split(" and ")
                    ],
                )
            ]
        else:
            try:
                return [
                    LogicalEnzymeComplex(
                        model=self.model,
                        word=self._enzyme_base_word + "_{}".format(0),
                        genes=[
                            self.model.get_gene(
                                g.replace("(", "").replace(")", "").strip()
                            )
                            for g in rule_string.split(" and ")
                        ],
                    )
                ]
            except IndexError:
                return []

    def cnf_lines(self, compartments: bool = False) -> str:
        self._build_input_clause(compartments=compartments)
        self._build_output_clauses(compartments=compartments)
        self._build_activation_clauses()
        cnf_lines = [
            "%----- REACTION: {} ({})".format(
                self.name + (" - reverse" if self.REVERSE else ""),
                ",".join(sorted(map(lambda c: c.name, self.compartments))),
            ),
            "%----- input",
            str(self._input_clause),
        ]

        if len(self._activation_clauses) > 0:
            cnf_lines.extend(
                ["%----- activation", *map(str, self._activation_clauses),]
            )

        cnf_lines.extend(
            ["%----- outputs", *map(str, self._output_clauses),]
        )

        return "\n".join(cnf_lines) + "\n" * 2

    # def _print_enzyme_cnf(self) -> str:
    #     es = {
    #         "e_{}_{i}".format(self.word, i=i): [
    #             "cnf({e}, axiom, ( e({e}) | {rhs}) ).".format(
    #                 e="e_{}_{i}".format(self.word, i=i),
    #                 rhs=(
    #                     " | ".join(["~g({})".format(process_name(g, names)) for g in c])
    #                 ),
    #             )
    #         ]
    #         for (i, c) in enumerate([c for c in cs if len(c) > 1])
    #     }

    #     # if es != {}:
    #     #     print("")
    #     #     print(type(es), es)

    #     # generate "enzyme_presence" cnfs, one for each complex or gene

    #     enzyme_presence = [
    #         "cnf({p}_{i}, axiom, ( p({p}) | ~e({e}) )).".format(i=i, p=eid, e=p)
    #         for (i, p) in enumerate(es)
    #     ] + [
    #         "cnf({p}_{i}, axiom, ( p({p}) | ~g({g}) )).".format(
    #             i=i + len(es), p=eid, g=process_name(c[0], names)
    #         )
    #         for (i, c) in enumerate(cs)
    #         if len(c) == 1
    #     ]


# %%
class LogicalEnzymeComplex:
    def __init__(
        self,
        model: "LogicalMetabolicModel",
        genes: Iterable["LogicalGene"],
        word: Optional[str] = None,
    ):

        self.model = model
        self.genes = set(genes)
        self.word = word

        if self.word is None:
            self.word = "ec_" + "_".join(
                map(str, sorted(self.genes, key=lambda g: g.orf))
            )

    def __str__(self):
        try:
            return " and ".join([g._cobra_gene.id for g in self.genes])
        except AttributeError:
            return " and ".join([g.word for g in self.genes])

    def __len__(self) -> int:
        return len(self.genes)

    def _build_clause(self):
        self._clause = LogicalClause(
            name=self.word,
            literals=[self.as_literal()]
            + [g.as_literal(negated=True) for g in self.genes],
            type="axiom",
        )

    def cnf_lines(self):
        self._build_clause()
        cnf_lines = [
            "%----- ENZYME COMPLEX: {{{}}} ({{{}}})".format(
                " and ".join(sorted(map(lambda g: g.name, self.genes))),
                " and ".join(sorted(map(lambda g: g.orf, self.genes))),
            ),
            str(self._clause),
        ]
        return "\n".join(cnf_lines) + "\n" * 2

    def as_literal(self, negated: bool = False) -> "LogicalLiteral":
        return LogicalLiteral(
            predicate_symbol="enz", arguments=[self.word], negated=negated
        )
