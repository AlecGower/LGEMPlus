import os
from pathlib import Path
from typing import Iterable, Optional


def abduce_hypotheses(
    theory_directory: os.PathLike,
    excluded_predicates: Iterable = ["rxn", "pro", "gn"],
    knockouts: Optional[Iterable[str]] = None,
):
    theory_directory = Path(theory_directory)

    if knockouts is not None:
        deletion_command = [
            "sed",
            '"/{0}/ s/gn(g_/~gn(g_/;/{0}/ s/_in_genome/_deletion/"'.format(
                "|".join(knockouts)
            ),
            str(theory_directory / "genes.p"),
        ]
    theory_command = [
        "cat",
        "-",
        *[
            str(theory_directory / fp)
            for fp in [
                "reactions.p",
                "abduced_extra_compounds.p",
                "media_compounds.p",
                "ubiquitous_compounds.p",
                "query.p",
            ]
        ],
    ]
    proof_command = ['$IPROVER_HOME"/iproveropt"', "--stdin", "true"]
    hypotheses_command = [
        "grep",
        '"{~("',
    ]
    filtered_command = [
        "grep",
        "-vE",
        '"({})\("'.format("|".join(excluded_predicates)),
    ]

    if knockouts is not None:
        commands = [deletion_command]
    else:
        commands = []
    commands.extend(
        [theory_command, proof_command, hypotheses_command, filtered_command,]
    )

    filtered = os.popen("|".join([" ".join(command) for command in commands]))
    hypotheses = [
        h[h.find("{") + 1 : h.find("}")].split(";") for h in filtered.readlines()
    ]

    return hypotheses
