import pandas as pd


def load_media_dataframe(filepath, exchange_flux_dictionary_filepath):
    media_comp = (
        pd.read_csv(
            filepath,
            sep="\t",
        )
        .fillna(0.0)
        .rename(columns={"Unnamed: 0": "exchange"})
    )
    media_comp["exchange"] = media_comp["exchange"].str.strip()
    with open(exchange_flux_dictionary_filepath) as fi:
        exchange_flux_dictionary = dict(map(lambda l: l.rstrip().split(",")[:2], fi))
    exchange_flux_dictionary = {
        k: v for k, v in exchange_flux_dictionary.items() if v != ""
    }
    media_comp["ygem_reaction"] = media_comp["exchange"].apply(
        exchange_flux_dictionary.get
    )
    media_comp = media_comp.dropna(subset=["ygem_reaction"])
    media_comp = media_comp.set_index("ygem_reaction")

    return media_comp.drop(columns="exchange").astype(float)
