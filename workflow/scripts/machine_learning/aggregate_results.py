#!/usr/bin/env python3
"""
aggregate_results.py
====================
Run after all SLURM jobs complete.
Reads every per-job CSV from results/raw/, builds:

  results/
  ├── aggregated_fold_level.csv
  ├── summary_by_k.csv
  ├── summary_best_k.csv
  ├── failed_jobs.csv
  └── plots/
      ├── auc_vs_k_<modality>.png
      ├── heatmap_best_k.png
      ├── seed_stability_violin.png
      └── roc_curves_best_combo.png   ← NEW: one panel per modality, 95% CI
"""

import glob
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from sklearn.metrics import roc_curve, auc

BASE_DIR   = Path(__file__).resolve().parent.parent
RAW_DIR    = BASE_DIR / "results" / "raw"
OUT_DIR    = BASE_DIR / "results"
PLOT_DIR   = OUT_DIR / "plots"
PLOT_DIR.mkdir(parents=True, exist_ok=True)

K_VALUES    = [10, 25, 50, 100, 200, 300, 400, 500]
MODALITIES  = ["PBL_PROT", "Tissue_PROT", "PBL_RNA", "Tissue_RNA"]
MODEL_NAMES = ["LogisticRegression", "GaussianNB", "XGBoost", "KNeighbors"]


# ── 1. Load all raw CSVs ──────────────────────────────────────
def load_all_results():
    ok_files     = sorted(RAW_DIR.glob("*__*.csv"))
    failed_files = sorted(RAW_DIR.glob("*_FAILED.csv"))

    fail_rows = []
    for f in failed_files:
        try:
            fail_rows.append(pd.read_csv(f))
        except Exception:
            pass
    failed_df = pd.concat(fail_rows, ignore_index=True) if fail_rows else pd.DataFrame()

    dfs = []
    for f in ok_files:
        if "_FAILED" in f.name:
            continue
        try:
            dfs.append(pd.read_csv(f))
        except Exception as e:
            print(f"Warning: could not read {f}: {e}")

    if not dfs:
        raise RuntimeError(f"No result CSVs found in {RAW_DIR}")

    all_df = pd.concat(dfs, ignore_index=True)
    print(f"Loaded {len(dfs):,} job files → {len(all_df):,} fold rows")
    if len(failed_df):
        print(f"  ⚠  {len(failed_df)} failed jobs — see results/failed_jobs.csv")
    return all_df, failed_df


# ── 2. Summary tables ─────────────────────────────────────────
def build_summary(all_df):
    job_df = (
        all_df
        .groupby(["Modality", "ModelName", "k_requested", "Seed"], as_index=False)
        ["LOOCV_Score"].first()
    )
    summary = (
        job_df
        .groupby(["Modality", "ModelName", "k_requested"])
        .agg(
            AUC_mean   = ("LOOCV_Score", "mean"),
            AUC_std    = ("LOOCV_Score", "std"),
            AUC_median = ("LOOCV_Score", "median"),
            AUC_min    = ("LOOCV_Score", "min"),
            AUC_max    = ("LOOCV_Score", "max"),
            N_seeds    = ("LOOCV_Score", "count"),
        )
        .reset_index()
        .sort_values(["Modality", "ModelName", "k_requested"])
    )
    return summary, job_df


def best_k_per_combo(summary):
    idx  = summary.groupby(["Modality", "ModelName"])["AUC_mean"].idxmax()
    best = summary.loc[idx].copy().sort_values("AUC_mean", ascending=False)
    return best


# ── 3. Standard plots ─────────────────────────────────────────
def plot_auc_vs_k(summary):
    for mod in MODALITIES:
        sub = summary[summary["Modality"] == mod]
        if sub.empty:
            continue
        fig, ax = plt.subplots(figsize=(8, 5))
        for model in MODEL_NAMES:
            msub = sub[sub["ModelName"] == model].sort_values("k_requested")
            if msub.empty:
                continue
            ax.errorbar(
                msub["k_requested"], msub["AUC_mean"],
                yerr=msub["AUC_std"],
                marker="o", linewidth=1.8, capsize=4, label=model,
            )
        ax.set_xscale("log")
        ax.set_xticks(K_VALUES)
        ax.set_xticklabels(K_VALUES)
        ax.set_xlabel("k (SelectKBest features)")
        ax.set_ylabel("Mean LOOCV ROC-AUC (± std)")
        ax.set_title(f"{mod} — AUC vs k")
        ax.axhline(0.5, color="gray", linestyle="--", alpha=0.6, label="Chance")
        ax.legend(fontsize=9)
        ax.set_ylim(0, 1.05)
        plt.tight_layout()
        plt.savefig(PLOT_DIR / f"auc_vs_k_{mod}.png", dpi=150)
        plt.close()
        print(f"Saved → plots/auc_vs_k_{mod}.png")


def plot_best_k_heatmap(best_df):
    pivot_auc = best_df.pivot(index="Modality", columns="ModelName", values="AUC_mean")
    pivot_k   = best_df.pivot(index="Modality", columns="ModelName", values="k_requested")
    pivot_auc = pivot_auc.reindex(index=MODALITIES, columns=MODEL_NAMES)
    pivot_k   = pivot_k.reindex(index=MODALITIES, columns=MODEL_NAMES)
    annot = pivot_auc.round(3).astype(str) + "\n(k=" + pivot_k.astype(str) + ")"

    fig, ax = plt.subplots(figsize=(10, 5))
    sns.heatmap(pivot_auc, annot=annot, fmt="", cmap="YlGnBu",
                vmin=0.4, vmax=1.0, ax=ax, linewidths=0.5,
                linecolor="white", annot_kws={"size": 9})
    ax.set_title("Best mean LOOCV AUC per modality × model (annotation: AUC | best k)")
    plt.tight_layout()
    plt.savefig(PLOT_DIR / "heatmap_best_k.png", dpi=150)
    plt.close()
    print("Saved → plots/heatmap_best_k.png")


def plot_seed_distribution(job_df, best_df):
    key_cols = ["Modality", "ModelName", "k_requested"]
    merged   = job_df.merge(best_df[key_cols], on=key_cols, how="inner")
    merged["combo"] = merged["Modality"] + "\n" + merged["ModelName"]
    combos = merged["combo"].unique()

    fig, ax = plt.subplots(figsize=(max(10, len(combos) * 1.2), 6))
    data_by_combo = [merged[merged["combo"] == c]["LOOCV_Score"].values for c in combos]
    parts = ax.violinplot(data_by_combo, positions=range(len(combos)),
                          showmedians=True, showextrema=True)
    for pc in parts["bodies"]:
        pc.set_alpha(0.7)
    ax.axhline(0.5, color="red", linestyle="--", alpha=0.5, label="Chance")
    ax.set_xticks(range(len(combos)))
    ax.set_xticklabels(combos, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("LOOCV AUC across seeds (at best k)")
    ax.set_title("Seed stability at best k per combo")
    ax.legend()
    plt.tight_layout()
    plt.savefig(PLOT_DIR / "seed_stability_violin.png", dpi=150)
    plt.close()
    print("Saved → plots/seed_stability_violin.png")


# ── 4. ROC curves — all modalities on one figure, coloured ──
def plot_roc_curves(all_df, best_df):
    """
    All 4 modalities on a SINGLE figure, each in a different colour.
    Mean ROC curve + 95% CI band across seeds at the best (model × k).
    Legend states the model and k used for each modality.
    """
    COLORS = {
        "PBL_PROT":    "#2196F3",   # blue
        "Tissue_PROT": "#E91E63",   # pink/red
        "PBL_RNA":     "#4CAF50",   # green
        "Tissue_RNA":  "#FF9800",   # orange
    }

    best_per_mod = (
        best_df
        .sort_values("AUC_mean", ascending=False)
        .groupby("Modality")
        .first()
        .reset_index()
    )

    fpr_grid = np.linspace(0, 1, 300)
    fig, ax  = plt.subplots(figsize=(8, 7))

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

        tprs, seed_aucs = [], []
        for seed in sub["Seed"].unique():
            sdf    = sub[sub["Seed"] == seed]
            y_true = sdf["TrueLabel"].values
            y_scr  = sdf["PredScore"].values
            if len(np.unique(y_true)) < 2:
                continue
            fpr_s, tpr_s, _ = roc_curve(y_true, y_scr)
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

        mean_auc = float(np.mean(seed_aucs))
        std_auc  = float(np.std(seed_aucs))

        label = (
            f"{mod}  [{best_model}, k={best_k}]\n"
            f"  AUC = {mean_auc:.3f} ± {std_auc:.3f}  (n={n} seeds)"
        )
        ax.plot(fpr_grid, mean_tpr, color=color, linewidth=2.2, label=label)
        ax.fill_between(fpr_grid, ci_lo, ci_hi, color=color, alpha=0.12)
        print(f"  {mod}: {best_model}, k={best_k}, AUC={mean_auc:.3f}±{std_auc:.3f}, n={n}")

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
    out_path = PLOT_DIR / "roc_curves_combined.png"
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print("Saved → plots/roc_curves_combined.png")


# ── Main ─────────────────────────────────────────────────────
def main():
    print("Loading results...")
    all_df, failed_df = load_all_results()

    fold_path = OUT_DIR / "aggregated_fold_level.csv"
    all_df.to_csv(fold_path, index=False)
    print(f"Saved fold-level table → {fold_path}  ({len(all_df):,} rows)")

    if not failed_df.empty:
        failed_df.to_csv(OUT_DIR / "failed_jobs.csv", index=False)

    summary, job_df = build_summary(all_df)
    summary.to_csv(OUT_DIR / "summary_by_k.csv", index=False)
    print(f"Saved → summary_by_k.csv")

    best_df = best_k_per_combo(summary)
    best_df.to_csv(OUT_DIR / "summary_best_k.csv", index=False)
    print(f"\nBest k per combo:")
    print(best_df[["Modality","ModelName","k_requested","AUC_mean","AUC_std","N_seeds"]].to_string(index=False))

    print("\nGenerating plots...")
    plot_auc_vs_k(summary)
    plot_best_k_heatmap(best_df)
    plot_seed_distribution(job_df, best_df)

    print("\nGenerating ROC curves with 95% CI...")
    plot_roc_curves(all_df, best_df)

    print("\nDone. Results written to:", OUT_DIR)


if __name__ == "__main__":
    main()
