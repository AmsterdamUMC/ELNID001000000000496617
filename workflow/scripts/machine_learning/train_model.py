#!/usr/bin/env python3
"""
train_model.py  (Snakemake script for rule `ml_train`)
========================================================
One Snakemake job = one (modality x model x k).
Outer resampling = Stratified k-fold CV, REPEATED `n_seeds` times.
The `n_seeds` seeds drive the OUTER split (not TPOT's internal search).

Results -> {output.folds_csv}
Feature usage -> {output.feature_usage_csv}

Assumes `input.data` and `input.labels` are ALREADY restricted to the donors
shared across all modalities (done upstream by rule `ml_prepare_data`), so no
complete-case merging is done here — this script simply intersects the one
modality's sample index with the labels' donor index as a sanity check.

Design
------
* Outer loop: StratifiedKFold(n_splits=cv, shuffle=True, random_state=split_seed),
  run once per split_seed in SPLIT_SEEDS (`n_seeds` repeats). Repeat r therefore
  has a genuinely different train/test partition -> the spread of CV_AUC across
  the repeats reflects partition variability, i.e. a real (if in-sample)
  stability estimate, not just TPOT search noise.

* TPOT random_state is FIXED (TPOT_SEED) so it does not contribute variance
  across repeats.

* Per repeat: each sample is held out in exactly one fold, predicted by a model
  trained on the other folds. Pool those out-of-fold predictions and compute
  ONE AUC per repeat (CV_AUC).

* Leakage: zero-variance drop + SelectKBest are fit on the TRAIN split only.
"""

import sys
import traceback

# ── Log FIRST, before any heavy imports or touching snakemake.input/output/params ──
# A crash in `import tpot`/`import xgboost` (missing or mismatched conda env), a bad
# config key, or a missing wildcard must all land in the log file instead of
# vanishing into Snakemake's own stderr capture, which is easy to lose track of.
log_fh = open(snakemake.log[0], "w")
sys.stdout = log_fh
sys.stderr = log_fh

import faulthandler
faulthandler.enable(file=log_fh)  # dump a trace to the log even on a hard crash (e.g. segfault)


def _excepthook(exc_type, exc_value, exc_tb):
    traceback.print_exception(exc_type, exc_value, exc_tb, file=log_fh)
    log_fh.flush()


sys.excepthook = _excepthook

try:
    import warnings
    warnings.filterwarnings(
        "ignore",
        message=r"passing a class to None is deprecated.*",
        category=FutureWarning,
        module=r"sklearn\.base",
    )
    from collections import Counter
    import numpy as np
    import pandas as pd
    from sklearn.base import BaseEstimator, TransformerMixin
    from sklearn.feature_selection import SelectKBest, f_classif
    from sklearn.metrics import roc_auc_score, accuracy_score
    from sklearn.model_selection import StratifiedKFold
    from tpot import TPOTClassifier

except Exception:
    print("FAILED while importing dependencies (check the conda env):", flush=True)
    traceback.print_exc(file=log_fh)
    log_fh.flush()
    raise

try:
    # ── Snakemake I/O ─────────────────────────────────────────────
    data_csv   = snakemake.input.data
    labels_csv = snakemake.input.labels
    folds_csv          = snakemake.output.folds_csv
    feature_usage_csv  = snakemake.output.feature_usage_csv
    _ml  = snakemake.params._ml
    modality   = snakemake.params.modality
    model_name = snakemake.params.model
    k          = int(snakemake.params.k)
    metric     = snakemake.params.metric
    n_seeds    = int(snakemake.params.n_seeds)
    seed       = int(snakemake.params.seed)
    tpot_cfg   = snakemake.params.tpot
    tpot_generations  = int(tpot_cfg.get("generations", 5))
    tpot_population = int(tpot_cfg.get("population_size", 10))
    donor_col      = _ml["donor_column"]
    label_col      = _ml["label_column"]
    classes        = _ml["classes"]
    positive_class = _ml["positive_class"]
except Exception:
    print("FAILED while reading snakemake input/output/params:", flush=True)
    traceback.print_exc(file=log_fh)
    log_fh.flush()
    raise

# ── Resampling config ─────────────────────────────────────────
N_SPLITS  = int(tpot_cfg.get("cv", 5))
TPOT_SEED = 42

# n_seeds fixed OUTER-split seeds, generated once (unique, reproducible from `seed`).
SPLIT_SEEDS = np.random.default_rng(seed).choice(2000, size=n_seeds, replace=False).tolist()

# ── TPOT configs, keyed by model name (matches config["machine_learning"]["models"]) ──
TPOT_CONFIG_BY_MODEL = {
    "Lasso": {
        # Lasso = L1-penalised logistic regression (sparse feature signature).
        'sklearn.linear_model.LogisticRegression': {
            'penalty':  ['l1'],
            'solver':   ['liblinear'],
            'C':        [0.001, 0.01, 0.1, 1.0],
            'max_iter': [5000],
        }
    },
    "Ridge": {
        # Ridge = L2-penalised logistic regression (shrinks, keeps all features).
        'sklearn.linear_model.LogisticRegression': {
            'penalty':  ['l2'],
            'solver':   ['liblinear'],
            'C':        [0.001, 0.01, 0.1, 1.0],
            'max_iter': [5000],
        }
    },
    "ElasticNet": {
        # Elastic Net = L1+L2 penalised logistic regression (saga supports it).
        'sklearn.linear_model.LogisticRegression': {
            'penalty':  ['elasticnet'],
            'solver':   ['saga'],
            'C':        [0.001, 0.01, 0.1, 1.0],
            'l1_ratio': [0.15, 0.5, 0.85],
            'max_iter': [5000],
        }
    },
    "XGBoost": {
        'xgboost.XGBClassifier': {
            'max_depth':        [2, 3],
            'learning_rate':    [0.05, 0.1],
            'n_estimators':     [50],
            'subsample':        [0.8],
            'colsample_bytree': [0.8],
            'reg_lambda':       [1.0, 5.0],
            'eval_metric':      ['logloss'],
        }
    },
    "GaussianNB": {
        'sklearn.naive_bayes.GaussianNB': {}
    },
}

if model_name not in TPOT_CONFIG_BY_MODEL:
    raise ValueError(
        f"Unknown model '{model_name}'. Must be one of {list(TPOT_CONFIG_BY_MODEL)}."
    )
config = TPOT_CONFIG_BY_MODEL[model_name]

# ── Helpers ──────────────────────────
class SafeSelectKBest(BaseEstimator, TransformerMixin):
    """SelectKBest that caps k to the available number of features and records
    WHICH features (by name) it picked, so callers can track feature usage."""
    def __init__(self, k=100, score_func=f_classif, feature_names=None):
        self.k = k
        self.score_func = score_func
        self.feature_names = feature_names

    def fit(self, X, y):
        k_safe = min(self.k, X.shape[1])
        self.selector_ = SelectKBest(score_func=self.score_func, k=k_safe)
        self.selector_.fit(X, y)
        support = self.selector_.get_support()
        if self.feature_names is not None:
            names = np.asarray(self.feature_names)
            self.selected_features_ = names[support].tolist()
        else:
            self.selected_features_ = np.where(support)[0].tolist()
        return self

    def transform(self, X):
        return self.selector_.transform(X)


def drop_zero_variance(X_train, X_test, feature_names):
    """Drop columns constant within the TRAIN split (fit on train only)."""
    keep = X_train.std(axis=0) > 0
    if keep.all():
        return X_train, X_test, feature_names
    names = np.asarray(feature_names)[keep].tolist()
    return X_train[:, keep], X_test[:, keep], names


def positive_class_scores(fitted, X):
    """Class-1 scores for a batch of samples: predict_proba -> decision_function
    -> predict. roc_auc_score accepts any of these."""
    if hasattr(fitted, "predict_proba"):
        return fitted.predict_proba(X)[:, 1]
    if hasattr(fitted, "decision_function"):
        return np.ravel(fitted.decision_function(X))
    return fitted.predict(X).astype(float)


def model_used_features(fitted_pipeline, selected_names):
    """Best-effort: of the SelectKBest-selected features, which ones the fitted model
    actually relies on — non-zero coefficient (Lasso/ElasticNet) or non-zero
    importance (XGBoost). Models without per-feature weights (e.g. GaussianNB) use all
    selected features. Fully guarded: any failure falls back to 'all selected', so it
    can never break a fold."""
    try:
        est = fitted_pipeline.steps[-1][1] if hasattr(fitted_pipeline, "steps") else fitted_pipeline
        if hasattr(est, "coef_"):
            w = np.abs(np.asarray(est.coef_)).ravel()
        elif hasattr(est, "feature_importances_"):
            w = np.asarray(est.feature_importances_).ravel()
        else:
            return list(selected_names)                 # no sparsity (e.g. GaussianNB)
        if len(w) != len(selected_names):
            return list(selected_names)                 # structure mismatch -> can't map
        return [n for n, wi in zip(selected_names, w) if wi != 0]
    except Exception:
        return list(selected_names)


# ── Data loading ─────────────────────────────────────────────────
def load_data():
    """Load X, y for this modality. `data_csv` and `labels_csv` are already
    restricted to the shared donor cohort by ml_prepare_data; here we just
    intersect indices as a sanity check (handles any residual mismatch, e.g. a
    donor dropped by the low-count filter on this modality's assay)."""
    labels_df = pd.read_csv(labels_csv, index_col=0)
    labels_df = labels_df[labels_df[label_col].isin(classes)]
    data = pd.read_csv(data_csv, index_col=0)
    common = sorted(set(data.index) & set(labels_df.index))
    if len(common) < N_SPLITS * 2:
        raise ValueError(
            f"Only {len(common)} donors shared between data ({data_csv}) and "
            f"labels ({labels_csv}) — too few for {N_SPLITS}-fold CV."
        )
    x = data.loc[common]
    y = pd.Series(
        (labels_df.loc[common, label_col] == positive_class).astype(int).values,
        index=common,
    )
    return x, y


def run_one_repeat(x, y, split_seed):
    """One repeat = one stratified k-fold split (seeded by split_seed).
    Returns per-sample DataFrame (n rows), the pooled repeat AUC, a Counter of
    feature selections, a Counter of model-used features, and the number of
    SelectKBest fits (= folds run)."""

    feature_names = x.columns.tolist()

    X_arr = x.values
    y_arr = y.values.copy()
    ids   = x.index.tolist()
    # Guard: stratified k-fold needs each class to have >= n_splits members.
    counts = np.bincount(y_arr)
    if counts.min() < N_SPLITS:
        raise ValueError(
            f"{model_name}: smallest class has {counts.min()} samples "
            f"< n_splits={N_SPLITS}; reduce cv for this modality."
        )
    skf = StratifiedKFold(n_splits=N_SPLITS, shuffle=True, random_state=split_seed)
    rows = []
    feat_counter = Counter()        # SelectKBest picks (filter stability)
    model_counter = Counter()       # features the fitted model actually uses
    n_fits = 0

    for fold, (train_idx, test_idx) in enumerate(skf.split(X_arr, y_arr)):

        X_train, X_test = X_arr[train_idx], X_arr[test_idx]
        y_train, y_test = y_arr[train_idx], y_arr[test_idx]

        X_train, X_test, fold_feats = drop_zero_variance(X_train, X_test, feature_names)
        selector   = SafeSelectKBest(k=k, score_func=f_classif, feature_names=fold_feats)
        X_train_fs = selector.fit_transform(X_train, y_train)
        X_test_fs  = selector.transform(X_test)

        selected = selector.selected_features_
        feat_counter.update(selected)
        n_fits += 1

        tpot = TPOTClassifier(
            generations     = tpot_generations,
            population_size = tpot_population, 
            cv              = N_SPLITS,
            scoring         = metric,
            config_dict     = config,
            verbosity       = 0,
            random_state    = TPOT_SEED,
        )
        tpot.fit(X_train_fs, y_train)

        # Features the model actually relies on (non-zero weight), among those selected.
        model_counter.update(model_used_features(tpot.fitted_pipeline_, selected))

        if metric == 'roc_auc':
            scores = positive_class_scores(tpot, X_test_fs)
        else:
            scores = tpot.predict(X_test_fs).astype(float)

        internal = getattr(tpot, "_optimized_pipeline_score", np.nan)
        pipe_str = str(tpot.fitted_pipeline_)

        for j, s_idx in enumerate(test_idx):
            rows.append({
                "Fold":             fold,
                "TestSample":       ids[s_idx],
                "TrueLabel":        int(y_test[j]),
                "PredScore":        float(scores[j]),
                "BestPipeline":     pipe_str,
                "InternalCV_Score": internal,
                "k_requested":      k,
                "k_actual":         X_train_fs.shape[1],
                "SelectedFeatures": ";".join(map(str, selected)),
                "ModelName":        model_name,
                "Seed":             split_seed,   # outer-split seed for this repeat
            })

    repeat_df = pd.DataFrame(rows)
    y_true  = repeat_df["TrueLabel"].values
    y_score = repeat_df["PredScore"].values

    if metric == 'roc_auc':
        cv_auc = roc_auc_score(y_true, y_score)         # pooled over the folds
    else:
        cv_auc = accuracy_score(y_true, (y_score > 0.5).astype(int))

    repeat_df["CV_AUC"] = cv_auc
    repeat_df["Metric"] = metric
    return repeat_df, cv_auc, feat_counter, model_counter, n_fits


# ── Main ──────────────────────────────────────────────────────────
def main():
    x, y = load_data()

    all_repeats = []
    failed = []
    feat_counter = Counter()        # SelectKBest selection tallies
    model_counter = Counter()       # model-used (non-zero weight) tallies
    n_fits_ok = 0

    for r, split_seed in enumerate(SPLIT_SEEDS):
        try:
            rep_df, cv_auc, rep_feats, rep_model_feats, n_fits = run_one_repeat(x, y, split_seed)
            rep_df["Repeat"]   = r
            rep_df["Modality"] = modality
            rep_df["Status"]   = "OK"
            all_repeats.append(rep_df)
            feat_counter.update(rep_feats)
            model_counter.update(rep_model_feats)
            n_fits_ok += n_fits
            print(f"Modality: {modality}, Model {model_name}, Repeat: {r}, CV AUC = {cv_auc:.4f}")
        except Exception:
            failed.append({"split_seed": split_seed, "error": traceback.format_exc()})
            print(f"    FAILED — {traceback.format_exc()}", flush=True)

    print(f"[checkpoint] main: repeat loop finished, "
          f"{len(all_repeats)} OK / {len(failed)} failed", flush=True)

    if not all_repeats:
        print("ERROR: all repeats failed, no output written.", flush=True)
        pd.DataFrame(failed).to_csv(folds_csv.replace(".csv", "_ALL_FAILED.csv"), index=False)
        # Still create the declared outputs (empty) so Snakemake doesn't error on a
        # missing-output check; the failure is visible in the log and the *_ALL_FAILED file.
        pd.DataFrame().to_csv(folds_csv, index=False)
        pd.DataFrame().to_csv(feature_usage_csv, index=False)
        sys.exit(1)

    result_df = pd.concat(all_repeats, ignore_index=True)
    result_df.to_csv(folds_csv, index=False)
    print(
        f"\nDone: {len(all_repeats)}/{len(SPLIT_SEEDS)} repeats OK | "
        f"{len(result_df)} rows -> {folds_csv}", flush=True
    )

    # Feature usage across the repeats x folds (real run only). Two frequencies:
    #   SelectionFrequency  = fraction of fold-fits where SelectKBest picked it.
    #   ModelUsedFrequency  = fraction of fold-fits where the fitted model gave it a
    #                         non-zero weight (the sparse signature for Lasso/ElasticNet,
    #                         non-zero importance for XGBoost; = SelectionFrequency for NB).

    if feat_counter and n_fits_ok:
        feature_df = pd.DataFrame(feat_counter.most_common(),
                                   columns=["Feature", "TimesSelected"])
        feature_df["TimesModelUsed"]     = feature_df["Feature"].map(
            lambda f: model_counter.get(f, 0))
        feature_df["TotalRuns"]          = n_fits_ok          # = repeats x folds
        feature_df["SelectionFrequency"] = feature_df["TimesSelected"] / n_fits_ok
        feature_df["ModelUsedFrequency"] = feature_df["TimesModelUsed"] / n_fits_ok
        feature_df["Modality"]           = modality
        feature_df["ModelName"]          = model_name
        feature_df["k_requested"]        = k
        feature_df = feature_df.sort_values(
            ["ModelUsedFrequency", "TimesSelected"], ascending=False)
        feature_df.to_csv(feature_usage_csv, index=False)
    else:
        # Output is declared in the rule, so it must exist even if nothing was tallied.
        pd.DataFrame(columns=[
            "Feature", "TimesSelected", "TimesModelUsed", "TotalRuns",
            "SelectionFrequency", "ModelUsedFrequency", "Modality", "ModelName", "k_requested",
        ]).to_csv(feature_usage_csv, index=False)

    if failed:
        fail_path = folds_csv.replace(".csv", "_partial_failures.csv")
        pd.DataFrame(failed).to_csv(fail_path, index=False)
        print(f"  {len(failed)} failed repeats logged -> {fail_path}")

try:
    main()
except Exception:
    traceback.print_exc(file=log_fh)
    log_fh.flush()
    raise
finally:
    log_fh.close()