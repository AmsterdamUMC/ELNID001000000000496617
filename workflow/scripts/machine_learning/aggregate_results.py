#!/usr/bin/env python3
"""Aggregate repeated-CV ML results and generate summary plots under Snakemake."""

import sys
import traceback
from pathlib import Path

# Start logging before third-party imports.
log_path = Path(snakemake.log[0])
log_path.parent.mkdir(parents=True, exist_ok=True)

log_fh = log_path.open("w", buffering=1)
sys.stdout = log_fh
sys.stderr = log_fh

print("[checkpoint] aggregation script started", flush=True)

try:
    print("[checkpoint] importing packages", flush=True)

    from itertools import combinations

    import matplotlib
    matplotlib.use("Agg")

    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import seaborn as sns
    from scipy import stats
    from sklearn.metrics import auc, roc_auc_score, roc_curve

    print("[checkpoint] imports completed", flush=True)

except Exception:
    print("[error] imports failed", flush=True)
    traceback.print_exc()
    log_fh.flush()
    raise


REQUIRED_RESULT_COLUMNS = {
    "Modality", "ModelName", "k_requested", "Seed", "CV_AUC",
    "TestSample", "TrueLabel", "PredScore",
}


def read_many(files, kind, required=None):
    frames = []
    for filename in map(Path, files):
        try:
            frame = pd.read_csv(filename)
        except Exception as exc:
            raise RuntimeError(f"Cannot read {kind} file {filename}: {exc}") from exc
        if required:
            missing = required.difference(frame.columns)
            if missing:
                raise ValueError(
                    f"{kind} file {filename} is missing columns: {sorted(missing)}"
                )
        if not frame.empty:
            frames.append(frame)
    if not frames:
        raise RuntimeError(f"No non-empty {kind} files were supplied")
    return pd.concat(frames, ignore_index=True)


def per_repeat_auc(results):
    keys = ["Modality", "ModelName", "k_requested", "Seed"]
    # Every sample row within a repeat should carry the same pooled CV_AUC.
    inconsistent = results.groupby(keys)["CV_AUC"].nunique(dropna=False)
    if (inconsistent > 1).any():
        bad = inconsistent[inconsistent > 1].index.tolist()[:5]
        raise ValueError(f"Inconsistent CV_AUC within repeat(s), for example: {bad}")
    return results.groupby(keys, as_index=False)["CV_AUC"].first()


def build_summaries(results):
    repeat = per_repeat_auc(results)
    summary = (
        repeat.groupby(["Modality", "ModelName", "k_requested"], as_index=False)
        .agg(
            AUC_mean=("CV_AUC", "mean"),
            AUC_std=("CV_AUC", "std"),
            AUC_median=("CV_AUC", "median"),
            AUC_min=("CV_AUC", "min"),
            AUC_max=("CV_AUC", "max"),
            N_repeats=("CV_AUC", "count"),
        )
        .sort_values(["Modality", "ModelName", "k_requested"])
    )
    idx = summary.groupby(["Modality", "ModelName"])["AUC_mean"].idxmax()
    best = summary.loc[idx].sort_values("AUC_mean", ascending=False).reset_index(drop=True)
    return repeat, summary, best


def paired_comparison(results, best, n_boot, seed):
    best_modality = (
        best.sort_values("AUC_mean", ascending=False)
        .groupby("Modality", as_index=False)
        .first()
    )
    predictions = {}
    for row in best_modality.itertuples(index=False):
        subset = results[
            (results["Modality"] == row.Modality)
            & (results["ModelName"] == row.ModelName)
            & (results["k_requested"] == row.k_requested)
        ]
        predictions[row.Modality] = subset.groupby("TestSample").agg(
            TrueLabel=("TrueLabel", "first"), Score=("PredScore", "mean")
        )

    modalities = list(predictions)
    if len(modalities) < 2:
        return pd.DataFrame(), pd.DataFrame()
    shared = sorted(set.intersection(*(set(predictions[m].index) for m in modalities)))
    if len(shared) < 4:
        print(f"Warning: only {len(shared)} samples are shared; skipping paired bootstrap")
        return pd.DataFrame(), pd.DataFrame()

    labels = {m: predictions[m].loc[shared, "TrueLabel"].to_numpy() for m in modalities}
    scores = {m: predictions[m].loc[shared, "Score"].to_numpy() for m in modalities}
    reference = labels[modalities[0]]
    for modality in modalities[1:]:
        if not np.array_equal(reference, labels[modality]):
            raise ValueError(f"True labels disagree across modalities for {modality}")

    rng = np.random.default_rng(seed)
    boot = {m: [] for m in modalities}
    for _ in range(n_boot):
        idx = rng.integers(0, len(shared), size=len(shared))
        y = reference[idx]
        if np.unique(y).size < 2:
            continue
        for modality in modalities:
            boot[modality].append(roc_auc_score(y, scores[modality][idx]))

    metadata = best_modality.set_index("Modality")
    rows = []
    for modality in modalities:
        arr = np.asarray(boot[modality])
        lo, hi = np.percentile(arr, [2.5, 97.5]) if arr.size else (np.nan, np.nan)
        rows.append({
            "Modality": modality,
            "Best_Model": metadata.loc[modality, "ModelName"],
            "Best_k": int(metadata.loc[modality, "k_requested"]),
            "AUC_shared": roc_auc_score(reference, scores[modality]),
            "CI95_low": lo,
            "CI95_high": hi,
            "N_shared": len(shared),
            "N_boot_valid": arr.size,
        })
    paired = pd.DataFrame(rows).sort_values("AUC_shared", ascending=False)

    pair_rows = []
    for first, second in combinations(modalities, 2):
        a, b = np.asarray(boot[first]), np.asarray(boot[second])
        pair_rows.append({
            "Modality_A": first,
            "Modality_B": second,
            "P_AUC_A_gt_AUC_B": np.mean(a > b) if a.size else np.nan,
            "Mean_AUC_difference_A_minus_B": np.mean(a - b) if a.size else np.nan,
            "N_boot_valid": a.size,
        })
    return paired, pd.DataFrame(pair_rows)


def make_ranking(best, paired):
    ranking = (
        best.sort_values("AUC_mean", ascending=False)
        .groupby("Modality", as_index=False)
        .first()[["Modality", "ModelName", "k_requested", "AUC_mean", "AUC_std"]]
        .rename(columns={"ModelName": "Best_Model", "k_requested": "Best_k"})
    )
    if not paired.empty:
        ranking = ranking.merge(
            paired[["Modality", "AUC_shared", "CI95_low", "CI95_high", "N_shared"]],
            on="Modality", how="left",
        )
    order = "AUC_shared" if "AUC_shared" in ranking else "AUC_mean"
    return ranking.sort_values(order, ascending=False).reset_index(drop=True)


def palette(modalities):
    colors = sns.color_palette("tab10", n_colors=max(len(modalities), 1))
    return dict(zip(modalities, colors))


def plot_auc_by_k(summary, modalities, models, k_values, plot_dir):
    for modality in modalities:
        sub = summary[summary["Modality"] == modality]
        fig, ax = plt.subplots(figsize=(8, 5))
        for model in models:
            data = sub[sub["ModelName"] == model].sort_values("k_requested")
            if not data.empty:
                ax.errorbar(data["k_requested"], data["AUC_mean"], yerr=data["AUC_std"],
                            marker="o", capsize=4, label=model)
        if len(k_values) > 1 and min(k_values) > 0:
            ax.set_xscale("log")
        ax.set_xticks(k_values); ax.set_xticklabels(k_values)
        ax.axhline(0.5, color="gray", linestyle="--", alpha=0.7)
        ax.set(xlabel="k (SelectKBest features)", ylabel="Mean repeated-CV ROC-AUC",
               title=f"{modality}: AUC versus k", ylim=(0, 1.05))
        ax.legend(fontsize=8); fig.tight_layout()
        fig.savefig(plot_dir / f"auc_vs_k_{modality}.pdf", dpi=180)
        plt.close(fig)


def plot_heatmap(best, modalities, models, plot_dir):
    auc_matrix = best.pivot(index="Modality", columns="ModelName", values="AUC_mean").reindex(
        index=modalities, columns=models
    )
    k_matrix = best.pivot(index="Modality", columns="ModelName", values="k_requested").reindex(
        index=modalities, columns=models
    )
    annotations = auc_matrix.applymap(lambda x: f"{x:.3f}" if pd.notna(x) else "")
    for row in modalities:
        for col in models:
            if pd.notna(k_matrix.loc[row, col]):
                annotations.loc[row, col] += f"\n(k={int(k_matrix.loc[row, col])})"
    fig, ax = plt.subplots(figsize=(max(8, 1.5 * len(models)), max(4, 0.8 * len(modalities))))
    sns.heatmap(auc_matrix, annot=annotations, fmt="", cmap="YlGnBu", vmin=0.4, vmax=1,
                linewidths=0.5, ax=ax)
    ax.set_title("Best mean CV AUC (and best k)")
    fig.tight_layout(); fig.savefig(plot_dir / "heatmap_best_k.pdf", dpi=180); plt.close(fig)


def plot_repeat_stability(repeat, best, plot_dir):
    keys = ["Modality", "ModelName", "k_requested"]
    data = repeat.merge(best[keys], on=keys, how="inner")
    data["Configuration"] = data["Modality"] + "\n" + data["ModelName"]
    fig, ax = plt.subplots(figsize=(max(10, 1.2 * data["Configuration"].nunique()), 6))
    sns.violinplot(data=data, x="Configuration", y="CV_AUC", inner="quartile", ax=ax)
    ax.axhline(0.5, color="red", linestyle="--", alpha=0.6)
    ax.tick_params(axis="x", rotation=45)
    ax.set_title("Repeated-CV stability at the best k per model")
    fig.tight_layout(); fig.savefig(plot_dir / "repeat_stability_violin.pdf", dpi=180); plt.close(fig)

def plot_roc_curve(results, best, modalities, plot_dir):
    best_modality = best.sort_values("AUC_mean", ascending=False).groupby("Modality").first()
    colors = palette(modalities)
    grid = np.linspace(0, 1, 300)
    fig, ax = plt.subplots(figsize=(7, 7))
    for modality in modalities:
        row = best_modality.loc[modality]
        subset = results[
            (results["Modality"] == modality)
            & (results["ModelName"] == row["ModelName"])
            & (results["k_requested"] == row["k_requested"])
        ]
        curves, aucs = [], []
        for _, seed_data in subset.groupby("Seed"):
            if seed_data["TrueLabel"].nunique() < 2:
                continue
            fpr, tpr, _ = roc_curve(seed_data["TrueLabel"], seed_data["PredScore"])
            curves.append(np.interp(grid, fpr, tpr)); aucs.append(auc(fpr, tpr))
        curves = np.asarray(curves)
        if not curves.size:
            continue
        mean = curves.mean(axis=0); mean[0] = 0; mean[-1] = 1
        sem = curves.std(axis=0, ddof=1) / np.sqrt(len(curves)) if len(curves) > 1 else np.zeros_like(mean)
        critical = stats.t.ppf(0.975, max(len(curves) - 1, 1))
        ax.plot(grid, mean, color=colors[modality], linewidth=2,
                label=f"{modality}: {row['ModelName']}, k={int(row['k_requested'])}, AUC={np.mean(aucs):.3f}")
        ax.fill_between(grid, np.clip(mean-critical*sem, 0, 1),
                        np.clip(mean+critical*sem, 0, 1), color=colors[modality], alpha=0.13)
    ax.plot([0, 1], [0, 1], "k--", alpha=0.5)
    ax.set(xlabel="FPR", ylabel="TPR",
           title="ROC curves for the best configuration per modality", xlim=(0, 1), ylim=(0, 1))
    ax.legend(fontsize=8, loc="lower right"); fig.tight_layout()
    fig.savefig(plot_dir / "roc_curves_combined.pdf", dpi=180); plt.close(fig)


def plot_top_features(features, best, top_n, plot_dir):
    if features.empty:
        return
    for row in best.itertuples(index=False):
        data = features[
            (features["Modality"] == row.Modality)
            & (features["ModelName"] == row.ModelName)
            & (features["k_requested"] == row.k_requested)
        ].copy()
        frequency = "ModelUsedFrequency" if (
            "ModelUsedFrequency" in data and data["ModelUsedFrequency"].fillna(0).sum() > 0
        ) else "SelectionFrequency"
        data = data.nlargest(top_n, frequency).sort_values(frequency)
        if data.empty:
            continue
        fig, ax = plt.subplots(figsize=(8, max(4, 0.3 * len(data))))
        ax.barh(data["Feature"], data[frequency], color="#3F51B5")
        ax.set(xlabel=frequency, title=f"{row.Modality}: {row.ModelName}, k={int(row.k_requested)}",
               xlim=(0, 1))
        fig.tight_layout()
        fig.savefig(plot_dir / f"top_features_{row.Modality}_{row.ModelName}.pdf", dpi=180)
        plt.close(fig)


def run():
    outputs = snakemake.output
    plot_dir = Path(outputs.plots)
    plot_dir.mkdir(parents=True, exist_ok=True)

    results = read_many(snakemake.input.folds, "fold-result", REQUIRED_RESULT_COLUMNS)
    features = read_many(snakemake.input.features, "feature-usage")
    modalities = list(snakemake.params.modalities)
    models = list(snakemake.params.models)
    k_values = sorted(map(int, snakemake.params.k_values))

    unknown_modalities = sorted(set(results["Modality"]) - set(modalities))
    unknown_models = sorted(set(results["ModelName"]) - set(models))
    if unknown_modalities or unknown_models:
        raise ValueError(f"Results/config mismatch: modalities={unknown_modalities}, models={unknown_models}")

    repeat, summary, best = build_summaries(results)
    paired, pairwise = paired_comparison(
        results, best, int(snakemake.params.n_boot), int(snakemake.params.seed)
    )
    ranking = make_ranking(best, paired)

    results.to_csv(outputs.aggregated, index=False)
    repeat.to_csv(outputs.repeat_auc, index=False)
    summary.to_csv(outputs.summary_by_k, index=False)
    best.to_csv(outputs.summary_best_k, index=False)
    paired.to_csv(outputs.paired, index=False)
    pairwise.to_csv(outputs.pairwise, index=False)
    ranking.to_csv(outputs.ranking, index=False)
    features.to_csv(outputs.feature_summary, index=False)

    plot_auc_by_k(summary, modalities, models, k_values, plot_dir)
    plot_heatmap(best, modalities, models, plot_dir)
    plot_repeat_stability(repeat, best, plot_dir)
    plot_roc_curve(results, best, modalities, plot_dir)
    plot_top_features(features, best, int(snakemake.params.top_n_features), plot_dir)
    print(f"Aggregated {len(snakemake.input.folds)} training files and {len(results)} prediction rows")
    print(ranking.to_string(index=False))

try:
    run()

except Exception:
    print("[error] aggregation failed", flush=True)
    traceback.print_exc()
    log_fh.flush()
    raise

finally:
    log_fh.close()
