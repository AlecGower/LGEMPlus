from argparse import ArgumentParser
import pandas as pd
from copy import copy
from sklearn.metrics import (
    classification_report,
    confusion_matrix,
    f1_score,
    jaccard_score,
)

ap = ArgumentParser(description="Takes a result list and compares with empirical data.")

ap.add_argument("results_file")
ap.add_argument(
    "--fba",
    default=False,
    type=bool,
    help="Boolean value specifying whether results are from an FBA simulation. Default = True",
)
ap.add_argument(
    "--fba_threshold",
    default=1e-6,
    type=float,
    help="growth rate threshold below which strain is judged inviable. Default = 1e-6",
)
arguments = ap.parse_args()

ess = pd.read_csv("data/sgd_ess.csv", index_col=0)
res = (
    pd.read_csv(arguments.results_file, delimiter="\t", names=["orf", "name", "growth"])
    .set_index("orf")
    .drop(columns="name")
)

if arguments.fba:
    print("Results:", len(res))
    preds = res.join(ess).loc[:, ["inviable", "growth"]].fillna(0)
    # # Changed because of train test split
    # preds = res.join(ess).loc[:, ["inviable", "growth"]].dropna()
    preds.growth = preds.growth < arguments.fba_threshold
    print("Preds:", len(preds))
else:
    print("Results:", len(res))
    preds = res.join(ess).loc[:, ["inviable", "growth"]].fillna(0)
    # preds = res.join(ess).loc[:, ["inviable", "growth"]].dropna()
    print("Preds:", len(preds))

majority = copy(preds)
majority["growth"] = 0

print("Majority class ('viable') classifier")

print(classification_report(majority.inviable, majority.growth, zero_division=0))

cm = confusion_matrix(majority.inviable, majority.growth)
# print("Confusion matrix:\n\r", cm)
print(
    pd.DataFrame(
        cm,
        columns=pd.MultiIndex.from_product([["Predicted"], ["G", "NG"]]),
        index=pd.MultiIndex.from_product([["Actual"], ["G", "NG"]]),
    ).to_markdown()
)

f1 = f1_score(majority.inviable, majority.growth, zero_division=0)
print("\nF1 score:", f1)

jaccard = jaccard_score(majority.inviable, majority.growth,)
print("Jaccard score: {:.3f}".format(jaccard))

print("\nOur classifier")

print(classification_report(preds.inviable, preds.growth,))

cm = confusion_matrix(preds.inviable, preds.growth,)
# print("Confusion matrix:\n\r", cm)
print(
    pd.DataFrame(
        cm,
        columns=pd.MultiIndex.from_product([["Predicted"], ["G", "NG"]]),
        index=pd.MultiIndex.from_product([["Actual"], ["G", "NG"]]),
    ).to_markdown()
)

f1 = f1_score(preds.inviable, preds.growth,)
print("\nF1 score: {:.3f}".format(f1))

jaccard = jaccard_score(preds.inviable, preds.growth,)
print("Jaccard score: {:.3f}".format(jaccard))

ngg = preds[(preds.inviable == 0) & (preds.growth == 1)].join(
    ess.loc[:, ["name", "full_name"]], how="left", rsuffix=".ess"
)
ngng = preds[(preds.inviable == 1) & (preds.growth == 1)].join(
    ess.loc[:, ["name", "full_name"]], how="left", rsuffix=".ess"
)
gng = preds[(preds.inviable == 1) & (preds.growth == 0)].join(
    ess.loc[:, ["name", "full_name"]], how="left", rsuffix=".ess"
)
gg = preds[(preds.inviable == 0) & (preds.growth == 0)].join(
    ess.loc[:, ["name", "full_name"]], how="left", rsuffix=".ess"
)


# fproot = os.path.split(results_file)[0]
# fproot += "/{}".format(datetime.now().strftime("%Y%m%dT%H%M%S"))
# if not os.path.exists(fproot):
#     os.mkdir(fproot)
# fproot += "/fba" if fba else "/log"

# for suffix in ["ngg", "gng", "ngng", "gg"]:
#     df = eval(suffix)
#     with open(fproot + "_{}.txt".format(suffix), "w") as fo:
#         for g, row in df.sort_values("name").iterrows():
#             print(g, row["name"], file=fo)