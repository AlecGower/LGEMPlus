# %%
# This code aims to set up and run FBA and MOMA experiments for yeast-GEM, based
# on the methodology in Snitkin, Evan S et al. (2008)

# %% Imports
from cobra.io import read_sbml_model
import pandas as pd
from pprint import pprint
import sys

sys.path.append("experiments/exploratory/different_growth_media")
from load_media import load_media_dataframe

# %% Load yeast-GEM model
ygem = read_sbml_model("model-files/yeastGEM.xml")

# %% Define growth media
media_comp = load_media_dataframe(
    "experiments/exploratory/different_growth_media/media-compositions-Snitkin-et-al-2008.tsv",
    "experiments/exploratory/different_growth_media/exchange-fluxes-dictionary.csv",
)


# %% Run base simulation
def implement_growth_media_fba(model, medium, media_df):
    assert medium in media_df.columns
    medium_dict = model.medium
    for rxnid in media_df.index:
        # rxn = model.reactions.get_by_id(rxnid)
        bound = media_df.loc[rxnid, medium]
        # rxn.lower_bound, rxn.upper_bound = -1 * bound, bound
        medium_dict[rxnid] = bound
    model.medium = medium_dict
    # pprint({medium : {r : {'name' : ygem.reactions.get_by_id(r).name, 'bound' : v} for r,v in medium_dict.items()}})
    # pprint({medium : {r : {'name' : ygem.reactions.get_by_id(r).name, 'bound' : v} for r,v in model.medium.items()}})
    return model


def get_solution(model, medium, media_df):
    with model:
        solution = implement_growth_media_fba(model, medium, media_df).optimize()

    # return solution.fluxes.to_frame().join(media_df['exchange'], how='inner').to_dict()
    return solution.objective_value, solution.status, solution.fluxes["r_1714"]


media_solutions = {
    med: get_solution(ygem, med, media_comp)
    for med in list(media_comp.columns)
    if med != "exchange"
}
pprint(media_solutions)

# %% Run mutant simulations
