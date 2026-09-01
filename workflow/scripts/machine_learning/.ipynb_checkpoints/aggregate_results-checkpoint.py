#!/usr/bin/env python3
"""
aggregate_results.py
====================
Run after all SLURM jobs complete.
Reads every per-job CSV from results/raw/, builds:

  results/
  ├── raw/                         ← per-job CSVs (written by worker.py)
  ├── aggregated_fold_level.csv    ← every fold row, all jobs merged
  ├── summary_by_k.csv             ← mean/std LOOCV AUC per modality×model×k
  ├── summary_best_k.csv           ← best k per modality×model (across seeds)
  ├── failed_jobs.csv              ← jobs that errored
  └── plots/
      ├── auc_vs_k_<modality>.png  ← AUC vs k curve per modality
      └── heatmap_best_k.png       ← best k heatmap
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

BASE_DIR   = Path(__file__).resolve().parent.parent
RAW_DIR    = BASE_DIR / "results" / "raw"
OUT_DIR    = BASE_DIR / "results"
PLOT_DIR   = OUT_DIR / "plots"
PLOT_DIR.mkdir(parents=True, exist_ok=True)

K_VALUES     = [10, 25, 50, 100, 200, 300, 400, 500]
MODALITIES   = ["PBL_PROT", "Tissue_PROT", "PBL_RNA", "Tissue_RNA"]
MODEL_NAMES  = ["LogisticRegression", "GaussianNB", "XGBoost", "KNeighbors"]


# ── 1. Load all raw CSVs ──────────────────────────────────────
def load_all_results():
    ok_files     = sorted(RAW_DIR.glob("*__*.csv"))
    failed_files = sorted(RAW_DIR.glob("*_FAILED.csv"))

    # Failed jobs
    fail_rows = []
    for f in failed_files:
        try:
            fail_rows.append(pd.read_csv(f))
        except Exception:
            pass
    failed_df = pd.concat(fail_rows, ignore_index=True) if fail_rows else pd.DataFrame()

    # Successful fold-level rows
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


# ── 2. Build LOOCV summary table ─────────────────────────────
def build_summary(all_df: pd.DataFrame) -> pd.DataFrame:
    """
    One row per (Modality, ModelName, k_requested, Seed) = LOOCV_Score.
    Then aggregate across seeds.
    """
    # LOOCV_Score is repeated for every fold row — take first (all identical)
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


def best_k_per_combo(summary: pd.DataFrame) -> pd.DataFrame:
    """For each (Modality, ModelName), find the k with the highest mean AUC."""
    idx = summary.groupby(["Modality", "ModelName"])["AUC_mean"].idxmax()
    best = summary.loc[idx].copy()
    best = best.sort_values("AUC_mean", ascending=False)
    return best


# ── 3. Plots ─────────────────────────────────────────────────
def plot_auc_vs_k(summary: pd.DataFrame):
    """One figure per modality: AUC ± std vs k, one line per model."""
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
        ax.set_ylabel("Mean LOOCV ROC-AUC (± std across seeds)")
        ax.set_title(f"{mod} — AUC vs k")
        ax.axhline(0.5, color="gray", linestyle="--", alpha=0.6, label="Chance")
        ax.legend(fontsize=9)
        ax.set_ylim(0, 1.05)

        path = PLOT_DIR / f"auc_vs_k_{mod}.png"
        plt.tight_layout()
        plt.savefig(path, dpi=150)
        plt.close()
        print(f"Saved → {path}")


def plot_best_k_heatmap(best_df: pd.DataFrame):
    """Heatmap: rows = modality, cols = model, cell = best AUC."""
    pivot_auc = best_df.pivot(index="Modality", columns="ModelName", values="AUC_mean")
    pivot_k   = best_df.pivot(index="Modality", columns="ModelName", values="k_requested")

    # Reorder to consistent order
    pivot_auc = pivot_auc.reindex(index=MODALITIES, columns=MODEL_NAMES)
    pivot_k   = pivot_k.reindex(index=MODALITIES, columns=MODEL_NAMES)

    annot = pivot_auc.round(3).astype(str) + "\n(k=" + pivot_k.astype(str) + ")"

    fig, ax = plt.subplots(figsize=(10, 5))
    sns.heatmap(
        pivot_auc, annot=annot, fmt="", cmap="YlGnBu",
        vmin=0.4, vmax=1.0, ax=ax,
        linewidths=0.5, linecolor="white",
        annot_kws={"size": 9},
    )
    ax.set_title("Best mean LOOCV AUC per modality × model\n(annotation: AUC | best k)")
    plt.tight_layout()
    path = PLOT_DIR / "heatmap_best_k.png"
    plt.savefig(path, dpi=150)
    plt.close()
    print(f"Saved → {path}")


def plot_seed_distribution(job_df: pd.DataFrame, best_df: pd.DataFrame):
    """
    Violin plot of per-seed LOOCV AUC at the best k for each combo.
    Shows stability across TPOT's internal randomness.
    """
    # Filter job_df to best k only
    key_cols = ["Modality", "ModelName", "k_requested"]
    best_keys = best_df[key_cols].copy()
    merged = job_df.merge(best_keys, on=key_cols, how="inner")

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
    path = PLOT_DIR / "seed_stability_violin.png"
    plt.savefig(path, dpi=150)
    plt.close()
    print(f"Saved → {path}")

def plot_best_roc_curves(all_df: pd.DataFrame, best_df: pd.DataFrame):
    """
    Plot final ROC curves for the best model × k for each Modality.
    Assumes fold-level table contains:
      Modality, ModelName, k_requested, y_true, y_score
    """

    required = ["Modality", "ModelName", "k_requested", "y_true", "y_score"]
    missing = [c for c in required if c not in all_df.columns]

    if missing:
        raise ValueError(
            "Cannot plot ROC curves. Missing columns in aggregated_fold_level.csv: "
            + ", ".join(missing)
        )

    roc_rows = []

    for _, row in best_df.iterrows():

        mod   = row["Modality"]
        model = row["ModelName"]
        k     = row["k_requested"]

        sub = all_df[
            (all_df["Modality"] == mod) &
            (all_df["ModelName"] == model) &
            (all_df["k_requested"] == k)
        ].copy()

        sub = sub.dropna(subset=["y_true", "y_score"])

        if sub.empty:
            print(f"Skipping ROC for {mod} / {model} / k={k}: no prediction rows")
            continue

        if sub["y_true"].nunique() < 2:
            print(f"Skipping ROC for {mod} / {model} / k={k}: only one class")
            continue

        fpr, tpr, thresholds = roc_curve(sub["y_true"], sub["y_score"])
        roc_auc = auc(fpr, tpr)

        for x, y, thr in zip(fpr, tpr, thresholds):
            roc_rows.append({
                "Modality": mod,
                "ModelName": model,
                "k_requested": k,
                "AUC": roc_auc,
                "FPR": x,
                "TPR": y,
                "Threshold": thr,
            })

    roc_df = pd.DataFrame(roc_rows)

    if roc_df.empty:
        print("No ROC curves generated.")
        return

    roc_path = OUT_DIR / "best_roc_curves_data.csv"
    roc_df.to_csv(roc_path, index=False)
    print(f"Saved ROC curve data → {roc_path}")

    n_panels = roc_df["Modality"].nunique()
    n_cols = 2
    n_rows = int(np.ceil(n_panels / n_cols))

    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=(6 * n_cols, 5 * n_rows),
        squeeze=False
    )

    axes = axes.flatten()

    for ax, mod in zip(axes, roc_df["Modality"].unique()):

        sub = roc_df[roc_df["Modality"] == mod]
        first = sub.iloc[0]

        label = (
            f"{first['ModelName']} | k={int(first['k_requested'])} | "
            f"AUC={first['AUC']:.3f}"
        )

        ax.plot(sub["FPR"], sub["TPR"], linewidth=2, label=label)
        ax.plot([0, 1], [0, 1], linestyle="--", color="gray", linewidth=1)

        ax.set_title(mod)
        ax.set_xlabel("False positive rate")
        ax.set_ylabel("True positive rate")
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1.05)
        ax.legend(loc="lower right", fontsize=8)

    for ax in axes[n_panels:]:
        ax.axis("off")

    plt.suptitle("Best ROC curve per modality", fontsize=14, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.96])

    path = PLOT_DIR / "best_roc_curves_by_modality.png"
    plt.savefig(path, dpi=200)
    plt.close()

    print(f"Saved → {path}")

# ── Main ─────────────────────────────────────────────────────
def main():
    print("Loading results...")
    all_df, failed_df = load_all_results()

    # Save fold-level master table
    fold_path = OUT_DIR / "aggregated_fold_level.csv"
    all_df.to_csv(fold_path, index=False)
    print(f"Saved fold-level table → {fold_path}  ({len(all_df):,} rows)")

    if not failed_df.empty:
        failed_df.to_csv(OUT_DIR / "failed_jobs.csv", index=False)

    # Summary tables
    summary, job_df = build_summary(all_df)
    summary.to_csv(OUT_DIR / "summary_by_k.csv", index=False)
    print(f"Saved summary_by_k → {OUT_DIR / 'summary_by_k.csv'}")

    best_df = best_k_per_combo(summary)
    best_df.to_csv(OUT_DIR / "summary_best_k.csv", index=False)
    print(f"\nBest k per combo:")
    print(best_df[["Modality","ModelName","k_requested","AUC_mean","AUC_std","N_seeds"]].to_string(index=False))

    # Plots
    plot_auc_vs_k(summary)
    plot_best_k_heatmap(best_df)
    plot_seed_distribution(job_df, best_df)
    plot_best_roc_curves(all_df, best_df)
    print("\nDone. Results written to:", OUT_DIR)


if __name__ == "__main__":
    main()
