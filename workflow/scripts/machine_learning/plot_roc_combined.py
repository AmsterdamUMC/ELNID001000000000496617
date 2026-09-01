#!/usr/bin/env python3
"""
plot_roc_combined.py
====================
Plot ROC curves for all 4 modalities on a SINGLE figure, each in a different
colour, showing the mean curve ± 95% CI band across seeds at the best (model × k).

Output: results/plots/roc_curves_combined.png
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_curve, auc

BASE_DIR = Path(__file__).resolve().parent.parent
RAW_DIR  = BASE_DIR / "results" / "raw"
OUT_DIR  = BASE_DIR / "results"
PLOT_DIR = OUT_DIR / "plots"
PLOT_DIR.mkdir(parents=True, exist_ok=True)

MODALITIES = ["PBL_PROT", "Tissue_PROT", "PBL_RNA", "Tissue_RNA"]

# Colour palette — one per modality
COLORS = {
    "PBL_PROT":    "#2196F3",   # blue
    "Tissue_PROT": "#E91E63",   # pink/red
    "PBL_RNA":     "#4CAF50",   # green
    "Tissue_RNA":  "#FF9800",   # orange
}


def load_data():
    dfs = [pd.read_csv(f) for f in RAW_DIR.glob("*__*.csv") if "_FAILED" not in f.name]
    all_df = pd.concat(dfs, ignore_index=True)
    best_df = pd.read_csv(OUT_DIR / "summary_best_k.csv")
    return all_df, best_df


def compute_mean_roc_with_ci(sub_df, fpr_grid):
    """
    Per-seed ROC → interpolate to common FPR grid → mean ± 95% CI.
    Returns mean_tpr, ci_lower, ci_upper, mean_auc, std_auc, n_seeds.
    """
    tprs, seed_aucs = [], []

    for seed in sub_df["Seed"].unique():
        sdf    = sub_df[sub_df["Seed"] == seed]
        y_true = sdf["TrueLabel"].values
        y_score= sdf["PredScore"].values

        if len(np.unique(y_true)) < 2:
            continue

        fpr_s, tpr_s, _ = roc_curve(y_true, y_score)
        seed_aucs.append(auc(fpr_s, tpr_s))

        tpr_i    = np.interp(fpr_grid, fpr_s, tpr_s)
        tpr_i[0] = 0.0
        tprs.append(tpr_i)

    tprs_arr = np.array(tprs)
    mean_tpr = tprs_arr.mean(axis=0)
    mean_tpr[-1] = 1.0

    n      = tprs_arr.shape[0]
    se     = tprs_arr.std(axis=0) / np.sqrt(n)
    t_crit = stats.t.ppf(0.975, df=max(n - 1, 1))
    ci_lo  = np.clip(mean_tpr - t_crit * se, 0, 1)
    ci_hi  = np.clip(mean_tpr + t_crit * se, 0, 1)

    return mean_tpr, ci_lo, ci_hi, float(np.mean(seed_aucs)), float(np.std(seed_aucs)), n


def main():
    all_df, best_df = load_data()

    # Best combo per modality (highest mean AUC across all models × k)
    best_per_mod = (
        best_df
        .sort_values("AUC_mean", ascending=False)
        .groupby("Modality")
        .first()
        .reset_index()
    )

    fpr_grid = np.linspace(0, 1, 300)

    fig, ax = plt.subplots(figsize=(8, 7))

    for mod in MODALITIES:
        row = best_per_mod[best_per_mod["Modality"] == mod]
        if row.empty:
            continue

        best_model = row["ModelName"].values[0]
        best_k     = int(row["k_requested"].values[0])
        color      = COLORS[mod]

        sub = all_df[
            (all_df["Modality"]    == mod) &
            (all_df["ModelName"]   == best_model) &
            (all_df["k_requested"] == best_k)
        ]

        mean_tpr, ci_lo, ci_hi, mean_auc, std_auc, n = compute_mean_roc_with_ci(sub, fpr_grid)

        label = (
            f"{mod}  [{best_model}, k={best_k}]\n"
            f"  AUC = {mean_auc:.3f} ± {std_auc:.3f}  (n={n} seeds)"
        )

        ax.plot(fpr_grid, mean_tpr, color=color, linewidth=2.2, label=label)
        ax.fill_between(fpr_grid, ci_lo, ci_hi, color=color, alpha=0.12)

    # Chance line
    ax.plot([0, 1], [0, 1], "k--", linewidth=1, alpha=0.5, label="Chance (AUC = 0.50)")

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.02)
    ax.set_xlabel("False Positive Rate", fontsize=13)
    ax.set_ylabel("True Positive Rate", fontsize=13)
    ax.set_title(
        "ROC Curves — Best model × k per modality\n"
        "Mean ± 95% CI across 100 seeds (LOOCV per-fold predictions)",
        fontsize=12
    )
    ax.legend(loc="lower right", fontsize=9.5, framealpha=0.9)
    ax.grid(True, alpha=0.25, linestyle="--")

    plt.tight_layout()
    out = PLOT_DIR / "roc_curves_combined.png"
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved → {out}")


if __name__ == "__main__":
    main()
