#!/usr/bin/env python3
"""
worker.py
=========
One SLURM job = one (modality × model × k), running all 100 seeds internally.
Results saved to results/raw/{modality}__{model}__k{k}.csv

Usage (called by sbatch):
    python worker.py --modality PBL_PROT --model_idx 0 --k 100
"""

import argparse
import os
import sys
import traceback

import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.metrics import roc_auc_score, accuracy_score
from sklearn.model_selection import LeaveOneOut
from sklearn.preprocessing import LabelEncoder
from tpot import TPOTClassifier

# ── Paths ─────────────────────────────────────────────────────
BASE_DIR   = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR   = os.path.join(BASE_DIR, "..", "data", "model")
RESULT_DIR = os.path.join(BASE_DIR, "results", "raw")
os.makedirs(RESULT_DIR, exist_ok=True)

# ── Seeds: 100 fixed seeds, generated once ────────────────────
SEEDS = np.random.default_rng(2025).choice(2000, size=100, replace=False).tolist()

# ── TPOT configs ──────────────────────────────────────────────
TPOT_CONFIG_LIST = [
    {
        'sklearn.linear_model.LogisticRegression': {
            'C':       [0.001, 0.01, 0.1],
            'penalty': ['l1', 'l2'],
            'solver':  ['liblinear'],
        }
    },
    {
        'sklearn.naive_bayes.GaussianNB': {}
    },
    {
        'xgboost.XGBClassifier': {
            'max_depth':        [2, 3],
            'learning_rate':    [0.05, 0.1],
            'n_estimators':     [50],
            'subsample':        [0.8],
            'colsample_bytree': [0.8],
            'reg_lambda':       [1.0, 5.0],
            'use_label_encoder':[False],
            'eval_metric':      ['logloss'],
        }
    },
    {
        'sklearn.neighbors.KNeighborsClassifier': {
            'n_neighbors': [3, 4, 5],
            'weights':     ['uniform', 'distance'],
            'metric':      ['euclidean', 'manhattan'],
        }
    },
]

MODEL_NAMES = ["LogisticRegression", "GaussianNB", "XGBoost", "KNeighbors"]
K_VALUES    = [10, 25, 50, 100, 200, 300, 400, 500]


# ── Helpers ───────────────────────────────────────────────────
class SafeSelectKBest(BaseEstimator, TransformerMixin):
    """SelectKBest that silently caps k to the available number of features."""
    def __init__(self, k=100, score_func=f_classif):
        self.k = k
        self.score_func = score_func

    def fit(self, X, y):
        k_safe = min(self.k, X.shape[1])
        self.selector_ = SelectKBest(score_func=self.score_func, k=k_safe)
        self.selector_.fit(X, y)
        return self

    def transform(self, X):
        return self.selector_.transform(X)


def load_data(modality: str):
    """Load X, y for a modality using its own available samples."""
    file_map = {
        "PBL_PROT":    "prot_pbl_data.csv",
        "Tissue_PROT": "prot_tissue_data.csv",
        "PBL_RNA":     "rna_pbl_data.csv",
        "Tissue_RNA":  "rna_tissue_data.csv",
    }
    labels_df = pd.read_csv(
        os.path.join(DATA_DIR, "labels.csv"),
        index_col="ClientGroupName", sep='\t'
    )
    data = pd.read_csv(
        os.path.join(DATA_DIR, file_map[modality]),
        index_col=0, sep='\t'
    )
    le = LabelEncoder()
    le.fit(labels_df['Response'])

    common = data.index.intersection(labels_df.index)
    x = data.loc[common]
    y = pd.Series(le.transform(labels_df.loc[common, 'Response']), index=common)
    return x, y


def run_loocv_one_seed(x, y, model_idx, seed, k, metric='roc_auc'):
    """Run full LOOCV for a single seed. Returns per-fold DataFrame + AUC."""
    config     = TPOT_CONFIG_LIST[model_idx]
    model_name = MODEL_NAMES[model_idx]
    selector   = SafeSelectKBest(k=k, score_func=f_classif)

    loo   = LeaveOneOut()
    X_arr = x.values
    y_arr = y.values
    ids   = x.index.tolist()
    rows  = []

    for fold, (train_idx, test_idx) in enumerate(loo.split(X_arr)):
        X_train, X_test = X_arr[train_idx], X_arr[test_idx]
        y_train, y_test = y_arr[train_idx], y_arr[test_idx]

        # Feature selection fitted on train only — no leakage
        X_train_fs = selector.fit_transform(X_train, y_train)
        X_test_fs  = selector.transform(X_test)

        tpot = TPOTClassifier(
            generations     = 5,
            population_size = 10,
            cv              = 5,
            scoring         = metric,
            config_dict     = config,
            verbosity       = 0,
            random_state    = seed,
        )

        tpot.fit(X_train_fs, y_train)

        if metric == 'roc_auc':
            score = tpot.predict_proba(X_test_fs)[0, 1]
        else:
            score = float(tpot.predict(X_test_fs)[0])

        rows.append({
            "Fold":             fold,
            "TestSample":       ids[test_idx[0]],
            "TrueLabel":        int(y_test[0]),
            "PredScore":        score,
            "BestPipeline":     str(tpot.fitted_pipeline_),
            "InternalCV_Score": tpot._optimized_pipeline_score,
            "k_requested":      k,
            "k_actual":         X_train_fs.shape[1],
            "ModelName":        model_name,
            "Seed":             seed,
        })

    fold_df = pd.DataFrame(rows)
    y_true  = fold_df["TrueLabel"].values
    y_score = fold_df["PredScore"].values

    if metric == 'roc_auc':
        loocv_score = roc_auc_score(y_true, y_score)
    else:
        loocv_score = accuracy_score(y_true, (y_score > 0.5).astype(int))

    fold_df["LOOCV_Score"] = loocv_score
    fold_df["Metric"]      = metric
    return fold_df, loocv_score


# ── Main ──────────────────────────────────────────────────────
def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--modality",  required=True,
                   choices=["PBL_PROT", "Tissue_PROT", "PBL_RNA", "Tissue_RNA"])
    p.add_argument("--model_idx", required=True, type=int, choices=[0, 1, 2, 3])
    p.add_argument("--k",         required=True, type=int, choices=K_VALUES)
    p.add_argument("--metric",    default="roc_auc")
    return p.parse_args()


def main():
    args = parse_args()
    model_name = MODEL_NAMES[args.model_idx]

    out_path = os.path.join(
        RESULT_DIR,
        f"{args.modality}__{model_name}__k{args.k}.csv"
    )
    if os.path.exists(out_path):
        print(f"Already done, skipping: {out_path}")
        sys.exit(0)

    x, y = load_data(args.modality)
    print(
        f"Loaded {args.modality}: n={len(y)} | "
        f"model={model_name}, k={args.k} | running {len(SEEDS)} seeds"
    )

    all_folds = []
    failed_seeds = []

    for seed_i, seed in enumerate(SEEDS):
        print(f"  seed {seed_i+1}/{len(SEEDS)} (seed={seed})", flush=True)
        try:
            fold_df, loocv_score = run_loocv_one_seed(
                x, y, args.model_idx, seed, args.k, args.metric
            )
            fold_df["Modality"] = args.modality
            fold_df["Status"]   = "OK"
            all_folds.append(fold_df)
            print(f"    LOOCV AUC = {loocv_score:.4f}")
        except Exception:
            failed_seeds.append({"seed": seed, "error": traceback.format_exc()})
            print(f"    FAILED — {traceback.format_exc()}", flush=True)

    if not all_folds:
        print("ERROR: all seeds failed, no output written.")
        pd.DataFrame(failed_seeds).to_csv(
            out_path.replace(".csv", "_ALL_FAILED.csv"), index=False
        )
        sys.exit(1)

    # Save all seeds in one CSV
    result_df = pd.concat(all_folds, ignore_index=True)
    result_df.to_csv(out_path, index=False)
    n_ok = len(all_folds)
    print(
        f"\nDone: {n_ok}/{len(SEEDS)} seeds OK | "
        f"{len(result_df)} total fold rows → {out_path}"
    )

    # Log any partial failures alongside the main output
    if failed_seeds:
        fail_path = out_path.replace(".csv", "_partial_failures.csv")
        pd.DataFrame(failed_seeds).to_csv(fail_path, index=False)
        print(f"  {len(failed_seeds)} failed seeds logged → {fail_path}")


if __name__ == "__main__":
    main()
