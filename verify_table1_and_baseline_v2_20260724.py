"""
Reproduce Table 1 and the exploratory isolation-date-only baseline in the Letter.

Required public inputs:
  1. Bandaru et al., J Clin Invest. 2026;136(13):e196284, Supplemental Table 1.
     Repository path:
     Supplementary_Tables/ST1/RAW_HNSCC_METADATA_NEW_v10.csv
  2. epifluidlab/headneck, commit
     dcc7a36a46f3d36e34992b08db8a57f0d24aef85
     Utils/Lists/cv_ids.txt
     Utils/Lists/holdout_ids.txt

Run from the checked-out repository root:
    python verify_table1_and_baseline_v2_20260724.py

The encoder and logistic regression are fitted only on the training patients in
each round. The date-only analysis is exploratory and unpaired with the rMDS
results because the exact published resampling indices were not available.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import sklearn
from sklearn.compose import ColumnTransformer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder

COMMIT = "dcc7a36a46f3d36e34992b08db8a57f0d24aef85"
META = Path("Supplementary_Tables/ST1/RAW_HNSCC_METADATA_NEW_v10.csv")
CV_IDS = Path("Utils/Lists/cv_ids.txt")
HO_IDS = Path("Utils/Lists/holdout_ids.txt")


def read_ids(path):
    return [value.strip() for value in path.read_text().splitlines() if value.strip()]


def to_patients(samples):
    """Collapse longitudinal samples after confirming patient-level consistency."""
    inconsistent_response = samples.groupby("OGID")["resp"].nunique()
    inconsistent_date = samples.groupby("OGID")["iso"].nunique()
    assert int((inconsistent_response > 1).sum()) == 0
    assert int((inconsistent_date > 1).sum()) == 0
    return (
        samples.groupby("OGID", as_index=False)
        .agg(resp=("resp", "first"), iso=("iso", "first"))
        .sort_values("OGID")
        .reset_index(drop=True)
    )


def make_model():
    preprocessing = ColumnTransformer(
        [("isolation_date", OneHotEncoder(handle_unknown="ignore"), ["iso"])],
        remainder="drop",
    )
    return Pipeline(
        [
            ("preprocessing", preprocessing),
            (
                "classifier",
                LogisticRegression(
                    C=1.0,
                    solver="lbfgs",
                    max_iter=5000,
                ),
            ),
        ]
    )


def preholdout_auc(patients, rounds=100, test_size=10, seed=42):
    """Repeat random 10-patient holdout evaluation with fold-contained encoding."""
    rng = np.random.RandomState(seed)
    y = patients["resp"].to_numpy()
    x = patients[["iso"]]
    all_indices = np.arange(len(patients))
    aucs = []

    for _ in range(rounds):
        test_indices = rng.choice(all_indices, test_size, replace=False)
        train_indices = np.setdiff1d(all_indices, test_indices)
        if np.unique(y[test_indices]).size < 2 or np.unique(y[train_indices]).size < 2:
            continue

        model = make_model()
        model.fit(x.iloc[train_indices], y[train_indices])
        probabilities = model.predict_proba(x.iloc[test_indices])[:, 1]
        aucs.append(roc_auc_score(y[test_indices], probabilities))

    return float(np.mean(aucs)), float(np.std(aucs)), len(aucs)


def main():
    for path in (META, CV_IDS, HO_IDS):
        if not path.exists():
            raise FileNotFoundError(f"Run from the repository root; missing: {path}")

    metadata = pd.read_csv(META)
    cv_ids = read_ids(CV_IDS)
    holdout_ids = read_ids(HO_IDS)
    assert len(cv_ids) == 151 and len(holdout_ids) == 25

    metadata["resp"] = (
        metadata["Treatment Response"].astype(str).str.strip() == "Responder"
    ).astype(int)
    metadata["iso"] = metadata["cfDNA Isolation Date"].astype(str)

    analysis = to_patients(metadata[metadata["ID"].isin(cv_ids)].copy())
    pooled = to_patients(metadata[metadata["ID"].isin(cv_ids + holdout_ids)].copy())
    assert len(analysis) == 58
    assert len(pooled) == 68

    table = pd.crosstab(analysis["iso"], analysis["resp"]).rename(
        columns={0: "NR", 1: "R"}
    )
    for column in ("NR", "R"):
        if column not in table:
            table[column] = 0
    table = table[["NR", "R"]]
    table["Total"] = table.sum(axis=1)

    single_class = table[(table[["NR", "R"]] > 0).sum(axis=1) == 1]
    multi_patient_single_class = single_class[single_class["Total"] >= 2]

    print(f"Repository commit: {COMMIT}")
    print(
        "Versions: "
        f"numpy={np.__version__}; pandas={pd.__version__}; "
        f"scikit-learn={sklearn.__version__}"
    )
    print("\nTable 1 — isolation batch × response, 58-patient analysis set")
    print(table.to_string())
    print(
        f"\nSingle-class batches: {len(single_class)} of {len(table)}; "
        f"{len(multi_patient_single_class)} multi-patient batches cover "
        f"{int(multi_patient_single_class['Total'].sum())} patients."
    )

    pooled_mean, pooled_sd, pooled_rounds = preholdout_auc(pooled)
    analysis_mean, analysis_sd, analysis_rounds = preholdout_auc(analysis)
    print(
        "\nDate-only baseline, 68 patients: "
        f"AUC {pooled_mean:.2f} ± {pooled_sd:.2f} "
        f"[{pooled_rounds} evaluable rounds]"
    )
    print(
        "Date-only baseline, 58-patient analysis set: "
        f"AUC {analysis_mean:.2f} ± {analysis_sd:.2f} "
        f"[{analysis_rounds} evaluable rounds]"
    )


if __name__ == "__main__":
    main()
