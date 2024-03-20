import sys

sys.path.append("experiments/exploratory/different_growth_media")
from load_media import load_media_dataframe
import re
from cobra.io import read_sbml_model

ygem = read_sbml_model("model-files/yeastGEM.xml")

media_comp = load_media_dataframe(
    "experiments/exploratory/different_growth_media/media-compositions-Snitkin-et-al-2008.tsv",
    "experiments/exploratory/different_growth_media/exchange-fluxes-dictionary.csv",
)

NW = re.compile(r"\W")
MET_SPECIES_PATTERN = re.compile(r"^(?P<SPECIES>.+)\s\[[^\[]+\]$")

media_comp.columns = list(map(lambda md: NW.sub("_", md), media_comp.columns))
media_comp["Compound Name"] = list(
    map(
        lambda rxnid: MET_SPECIES_PATTERN.match(
            ygem.reactions.get_by_id(rxnid).reactants[0].name
        ).group(1),
        media_comp.index,
    )
)
media_comp["Compound"] = ""
media_comp["yeastGEM"] = media_comp["Compound Name"]

for medium in media_comp.columns[:-3]:
    exchanges = media_comp.loc[:, ["Compound", "Compound Name", "yeastGEM", medium]]
    exchanges[exchanges[medium] > 0.0].reset_index(drop=True)[
        ["Compound", "Compound Name", "yeastGEM"]
    ].to_csv(f"model-files/{medium}-compounds-yeastGEM.tsv", sep="\t")
