#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import itertools
import math
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import ndtr
from sklearn.exceptions import ConvergenceWarning
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import (
    ConstantKernel,
    Matern,
    WhiteKernel,
)
from sklearn.preprocessing import StandardScaler

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

except ImportError:
    plt = None

# Parameter and objective definitions

OPTIMIZE_COLS = [
    "hspthresh",
    "gappedthresh",
    "ydrop",
    "inner",
    "step",
    "masking",
    "notransition",
    "chain",
]

TARGET_COLS = [
    "length_ratio",
    "num_ratio",
]

CATEGORICAL_COLS = {
    "masking",
    "notransition",
    "chain",
}

NUMERIC_COLS = [
    col for col in OPTIMIZE_COLS
    if col not in CATEGORICAL_COLS
]

PARAM_RULES = {
    "hspthresh": {
        "band": 1000,
        "step": 50,
        "min": 2200,
        "max": 4500,
    },
    "gappedthresh": {
        "band": 1000,
        "step": 200,
        "min": 3000,
        "max": 10000,
    },
    "ydrop": {
        "band": 1000,
        "step": 200,
        "min": 3400,
        "max": 15000,
    },
    "inner": {
        "band": 1000,
        "step": 200,
        "min": 2000,
        "max": 3000,
    },
    "step": {
        "band": 5,
        "step": 1,
        "min": 1,
    },
    "masking": {
        "categorical": True,
        "min": 0,
    },
    "notransition": {
        "categorical": True,
        "min": 0,
    },
    "chain": {
        "categorical": True,
        "min": 0,
    },
}
DEFAULT_LEN_TOL = 0.003
DEFAULT_NUM_TOL = 0.005

DEFAULT_LEN_WEIGHT = 0.60
DEFAULT_NUM_WEIGHT = 0.40

DEFAULT_MAX_CANDIDATES = 200000

# Data loading and validation

def read_data(csv_path: Path) -> pd.DataFrame:
    """
    Read historical experimental data.

    Required columns:
        8 optimization parameters
        length_ratio
        num_ratio

    time_min is no longer required or used.
    """

    df = pd.read_csv(csv_path)

    required_cols = OPTIMIZE_COLS + TARGET_COLS

    missing = [
        col for col in required_cols
        if col not in df.columns
    ]

    if missing:
        raise ValueError(
            f"CSV missing required columns: {missing}"
        )

    data = df[required_cols].copy()

    for col in required_cols:
        data[col] = pd.to_numeric(
            data[col],
            errors="coerce",
        )

    if data.empty:
        raise ValueError(
            "History CSV contains no rows"
        )

    if data.isna().any().any():
        bad_cols = data.columns[
            data.isna().any()
        ].tolist()

        raise ValueError(
            "History CSV contains NaN or "
            f"non-numeric values in: {bad_cols}"
        )

    if not np.isfinite(
        data.to_numpy(dtype=float)
    ).all():
        raise ValueError(
            "History CSV contains inf or -inf"
        )

    if (
        data["length_ratio"] < 0
    ).any() or (
        data["length_ratio"] > 1
    ).any():
        raise ValueError(
            "length_ratio must be between 0 and 1"
        )

    if (
        data["num_ratio"] < 0
    ).any() or (
        data["num_ratio"] > 1
    ).any():
        raise ValueError(
            "num_ratio must be between 0 and 1"
        )

    # Categorical parameters must be integers.
    for col in CATEGORICAL_COLS:
        values = data[col].to_numpy(
            dtype=float
        )

        rounded = np.rint(values)

        if not np.allclose(
            values,
            rounded,
        ):
            raise ValueError(
                f"Categorical parameter {col} "
                "contains non-integer values"
            )

        data[col] = rounded.astype(int)

    # Numeric parameters are converted to integers
    # because the candidate generation is grid-based.
    for col in NUMERIC_COLS:
        data[col] = np.rint(
            data[col].to_numpy(dtype=float)
        ).astype(int)

    return data


# Historical ranking and local anchor selection

def rank_observed(
    df: pd.DataFrame,
) -> pd.DataFrame:
    """
    Rank historical observations using lexicographic ranking:

        1. maximize length_ratio
        2. maximize num_ratio

    Runtime is no longer considered.
    """

    return df.sort_values(
        by=[
            "length_ratio",
            "num_ratio",
        ],
        ascending=[
            False,
            False,
        ],
        kind="mergesort",
    ).reset_index(drop=True)

def build_seen_set(
    df: pd.DataFrame,
):
    """
    Build a set containing all historically tested
    parameter combinations.
    """

    seen = set()

    for row in df[
        OPTIMIZE_COLS
    ].itertuples(
        index=False,
        name=None,
    ):
        key = tuple(
            int(round(float(value)))
            for value in row
        )

        seen.add(key)

    return seen

def topk_center(
    df: pd.DataFrame,
    topk: int,
):
    """
    Select Top-K historical observations and
    construct the local search center.

    Numerical parameters:
        arithmetic mean of Top-K values

    Categorical parameters:
        mode of Top-K values
    """

    if topk <= 0:
        raise ValueError(
            "anchor-topk must be positive"
        )

    top = rank_observed(
        df
    ).head(topk)

    if top.empty:
        raise ValueError(
            "No anchor rows available"
        )

    center = {}

    for col in OPTIMIZE_COLS:
        values = top[col].astype(
            float
        ).to_numpy()

        if col in CATEGORICAL_COLS:
            mode = pd.Series(
                values
            ).mode().iloc[0]

            center[col] = float(mode)

        else:
            center[col] = float(
                np.mean(values)
            )

    return center, top


# Local candidate generation

def make_numeric_values(
    center: float,
    band: int,
    step: int,
    min_value: int = 0,
    max_value: int = None,
):
    """
    Generate grid values around a local center,
    restricted by min_value and max_value.
    """

    if step <= 0:
        step = 1

    raw_low = int(math.floor(center - band))
    raw_high = int(math.ceil(center + band))

    low = max(
        min_value,
        raw_low,
    )

    high = raw_high

    if max_value is not None:
        high = min(
            max_value,
            high,
        )

    if low > high:
        clipped_center = int(round(center))
        clipped_center = max(
            min_value,
            clipped_center,
        )

        if max_value is not None:
            clipped_center = min(
                max_value,
                clipped_center,
            )

        return [clipped_center]

    start = (low // step) * step

    if start < low:
        start += step

    values = list(
        range(
            start,
            high + 1,
            step,
        )
    )

    if not values:
        clipped_center = int(round(center))
        clipped_center = max(
            min_value,
            clipped_center,
        )

        if max_value is not None:
            clipped_center = min(
                max_value,
                clipped_center,
            )

        values = [clipped_center]

    rounded_center = int(round(center))

    rounded_center = max(
        min_value,
        rounded_center,
    )

    if max_value is not None:
        rounded_center = min(
            max_value,
            rounded_center,
        )

    if rounded_center not in values:
        values.append(rounded_center)

    return sorted(
        set(
            int(value)
            for value in values
            if value >= min_value
            and (
                max_value is None
                or value <= max_value
            )
        )
    )

def make_candidate_values(
    df: pd.DataFrame,
    center: dict,
    top: pd.DataFrame,
):
    """
    Generate a local candidate value set for
    each optimization parameter.
    """

    values = {}

    for col in OPTIMIZE_COLS:
        rule = PARAM_RULES[col]

        # Categorical parameters are restricted to
        # categories appearing in the Top-K anchors.
        if col in CATEGORICAL_COLS:
            top_values = sorted(
                set(
                    int(round(float(value)))
                    for value in top[col].tolist()
                )
            )

            values[col] = top_values
            continue

        # The parameter "step" is fixed at 1 when
        # the local center is already at or below 1.
        if col == "step":
            if int(round(center[col])) <= 1:
                values[col] = [1]
                continue

        values[col] = make_numeric_values(
            center=center[col],
            band=rule["band"],
            step=rule["step"],
            min_value=rule["min"],
            max_value=rule.get("max"),
        )

    return values

def generate_candidates(
    value_map: dict,
):
    """
    Generate the Cartesian product of all
    parameter-specific candidate value sets.
    """

    grids = [
        value_map[col]
        for col in OPTIMIZE_COLS
    ]

    if any(
        len(grid) == 0
        for grid in grids
    ):
        return pd.DataFrame(
            columns=OPTIMIZE_COLS
        )

    rows = []

    for combo in itertools.product(*grids):
        row = {
            col: int(value)
            for col, value in zip(
                OPTIMIZE_COLS,
                combo,
            )
        }

        rows.append(row)

    if not rows:
        return pd.DataFrame(
            columns=OPTIMIZE_COLS
        )

    candidates = pd.DataFrame(
        rows,
        columns=OPTIMIZE_COLS,
    )

    return candidates.drop_duplicates(
        subset=OPTIMIZE_COLS,
    ).reset_index(drop=True)

def limit_candidates(
    candidates: pd.DataFrame,
    max_candidates: int,
    seed: int,
):
    """
    Randomly limit the number of candidates
    when the Cartesian product is too large.
    """

    if max_candidates <= 0:
        raise ValueError(
            "max-candidates must be positive"
        )

    if len(candidates) <= max_candidates:
        return candidates

    sampled = candidates.sample(
        n=max_candidates,
        random_state=seed,
    )

    return sampled.reset_index(drop=True)


# Feature encoding

def get_category_values(
    df: pd.DataFrame,
):
    """
    Obtain the historical category values used
    for one-hot encoding.
    """

    categories = {}

    for col in sorted(
        CATEGORICAL_COLS
    ):
        categories[col] = sorted(
            set(
                int(value)
                for value in df[col].astype(
                    int
                ).tolist()
            )
        )

    return categories

def build_feature_matrix(
    df: pd.DataFrame,
    categories: dict,
):
    """
    Build the GP feature matrix.

    Numerical parameters are retained as numerical features.
    Categorical parameters are one-hot encoded.
    """

    parts = []

    numeric_part = df[
        NUMERIC_COLS
    ].astype(
        float
    ).to_numpy()

    if len(NUMERIC_COLS) > 0:
        parts.append(numeric_part)

    for col in sorted(
        CATEGORICAL_COLS
    ):
        values = df[col].astype(
            int
        ).to_numpy()

        for category in categories[col]:
            encoded = (
                values == int(category)
            ).astype(
                float
            ).reshape(
                -1,
                1,
            )

            parts.append(encoded)

    if not parts:
        raise ValueError(
            "No feature columns available"
        )

    return np.hstack(parts)


# Gaussian-process surrogate model

def fit_gp(
    X: np.ndarray,
    y: np.ndarray,
    seed: int,
):
    """
    Fit one Gaussian-process regression model.

    Kernel:
        ConstantKernel × Matern(nu=2.5)
        + WhiteKernel
    """

    kernel = (
        ConstantKernel(
            1.0,
            (
                1e-4,
                1e4,
            ),
        )
        * Matern(
            length_scale=np.ones(
                X.shape[1]
            ),
            length_scale_bounds=(
                1e-3,
                1e6,
            ),
            nu=2.5,
        )
        + WhiteKernel(
            noise_level=1e-5,
            noise_level_bounds=(
                1e-10,
                1e0,
            ),
        )
    )

    gp = GaussianProcessRegressor(
        kernel=kernel,
        normalize_y=True,
        n_restarts_optimizer=5,
        random_state=seed,
        alpha=1e-8,
    )

    with warnings.catch_warnings():
        warnings.simplefilter(
            "always",
            ConvergenceWarning,
        )

        gp.fit(X, y)

    return gp


# Local distance

def nearest_anchor_distance(
    candidates: pd.DataFrame,
    anchor_rows: pd.DataFrame,
):
    """
    Calculate the nearest-anchor distance for
    each candidate configuration.
    """

    scales = {}

    for col in OPTIMIZE_COLS:
        values = anchor_rows[
            col
        ].astype(
            float
        ).to_numpy()

        span = float(
            np.max(values)
            - np.min(values)
        )

        scales[col] = max(
            1.0,
            span,
        )

    distances = []

    for _, candidate in candidates.iterrows():
        best_distance = None

        for _, anchor in anchor_rows.iterrows():
            distance = 0.0

            for col in OPTIMIZE_COLS:
                distance += (
                    abs(
                        float(candidate[col])
                        - float(anchor[col])
                    )
                    / scales[col]
                )

            if (
                best_distance is None
                or distance < best_distance
            ):
                best_distance = distance

        distances.append(
            0.0
            if best_distance is None
            else best_distance
        )

    return np.asarray(
        distances,
        dtype=float,
    )


# Expected improvement

def expected_improvement(
    mu: np.ndarray,
    sigma: np.ndarray,
    best: float,
    xi: float = 0.01,
):
    """
    Expected improvement for a maximization objective.
    """

    mu = np.asarray(
        mu,
        dtype=float,
    )

    sigma = np.asarray(
        sigma,
        dtype=float,
    )

    sigma = np.maximum(
        sigma,
        1e-12,
    )

    improvement = (
        mu
        - best
        - xi
    )

    z = improvement / sigma

    normal_pdf = (
        np.exp(-0.5 * z * z)
        / math.sqrt(2.0 * math.pi)
    )

    normal_cdf = ndtr(z)

    ei = (
        improvement * normal_cdf
        + sigma * normal_pdf
    )

    return np.maximum(
        ei,
        0.0,
    )


# Robust scaling and score construction

def robust_scale(
    values: np.ndarray,
):
    """
    Robust scale based primarily on the interquartile range.
    """

    values = np.asarray(
        values,
        dtype=float,
    )

    q25, q75 = np.percentile(
        values,
        [
            25,
            75,
        ],
    )

    scale = float(
        q75 - q25
    )

    if scale <= 1e-12:
        scale = float(
            np.std(values)
        )

    return max(
        scale,
        1e-4,
    )

def normalize_weights(
    length_weight: float,
    num_weight: float,
):
    """
    Normalize the two objective weights.
    """

    weights = np.asarray(
        [
            length_weight,
            num_weight,
        ],
        dtype=float,
    )

    if np.any(
        weights < 0
    ):
        raise ValueError(
            "Objective weights cannot be negative"
        )

    total = float(
        np.sum(weights)
    )

    if total <= 0:
        raise ValueError(
            "At least one objective weight "
            "must be positive"
        )

    return weights / total

# Pareto analysis

def pareto_mask(
    length_values: np.ndarray,
    num_values: np.ndarray,
):
    """
    Return the non-dominated points for:

        maximize length_ratio
        maximize num_ratio
    """

    length_values = np.asarray(
        length_values,
        dtype=float,
    )

    num_values = np.asarray(
        num_values,
        dtype=float,
    )

    n = len(length_values)

    is_pareto = np.ones(
        n,
        dtype=bool,
    )

    for i in range(n):
        dominated = (
            (length_values >= length_values[i])
            & (num_values >= num_values[i])
            & (
                (length_values > length_values[i])
                | (
                    num_values > num_values[i]
                )
            )
        )

        if np.any(dominated):
            is_pareto[i] = False

    return is_pareto


# Diagnostic plotting

def plot_optimization_summary(
    history: pd.DataFrame,
    candidates: pd.DataFrame,
    anchors: pd.DataFrame,
    recommended: pd.Series,
    output_path: Path,
):
    """
    Generate a three-panel optimization diagnostic figure.

    Panels:
        1. Historical objective space
        2. Predicted candidate objective space
        3. Balanced score versus local-anchor distance
    """

    if plt is None:
        print(
            "\nWarning: matplotlib is not installed. "
            "Skipping plot generation."
        )

        return

    if candidates.empty:
        print(
            "\nWarning: no candidates available. "
            "Skipping plot generation."
        )

        return

    output_path = Path(
        output_path
    )

    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    history_length = history[
        "length_ratio"
    ].to_numpy(
        dtype=float
    )

    history_num = history[
        "num_ratio"
    ].to_numpy(
        dtype=float
    )

    observed_pareto = pareto_mask(
        history_length,
        history_num,
    )

    candidate_length = candidates[
        "pred_length_ratio"
    ].to_numpy(
        dtype=float
    )

    candidate_num = candidates[
        "pred_num_ratio"
    ].to_numpy(
        dtype=float
    )

    candidate_pareto = pareto_mask(
        candidate_length,
        candidate_num,
    )

    fig, axes = plt.subplots(
        1,
        3,
        figsize=(17, 5.5),
        constrained_layout=True,
    )

    ax_observed = axes[0]
    ax_candidates = axes[1]
    ax_score = axes[2]


    # Panel 1: observed objective space

    observed_scatter = ax_observed.scatter(
        history_length,
        history_num,
        color="#9e9e9e",
        s=62,
        alpha=0.85,
        edgecolors="white",
        linewidths=0.6,
        label="Historical observations",
    )

    # Highlight Top-K local anchors.
    ax_observed.scatter(
        anchors["length_ratio"],
        anchors["num_ratio"],
        color="#d62728",
        s=85,
        marker="o",
        edgecolors="black",
        linewidths=0.7,
        label="Top-K local anchors",
        zorder=5,
    )

    if np.any(observed_pareto):
        pareto_order = np.argsort(
            history_length[
                observed_pareto
            ]
        )

        ax_observed.plot(
            history_length[
                observed_pareto
            ][pareto_order],
            history_num[
                observed_pareto
            ][pareto_order],
            color="#8c2d04",
            linewidth=1.8,
            marker="o",
            markersize=4,
            label="Observed Pareto front",
        )

    best_observed_index = int(
        rank_observed(
            history
        ).index[0]
    )

    best_observed_row = rank_observed(
        history
    ).iloc[0]

    ax_observed.scatter(
        [
            best_observed_row[
                "length_ratio"
            ]
        ],
        [
            best_observed_row[
                "num_ratio"
            ]
        ],
        marker="*",
        s=220,
        color="#e41a1c",
        edgecolors="black",
        linewidths=0.8,
        label="Best observed",
        zorder=6,
    )

    ax_observed.set_title(
        "(a) Historical observations"
    )

    ax_observed.set_xlabel(
        "length_ratio"
    )

    ax_observed.set_ylabel(
        "num_ratio"
    )

    ax_observed.grid(
        alpha=0.25
    )

    ax_observed.legend(
        fontsize=8,
        loc="best",
    )


    # Panel 2: predicted candidate objective space

    candidate_scatter = ax_candidates.scatter(
        candidate_length,
        candidate_num,
        color="#2ca25f",
        s=45,
        alpha=0.68,
        edgecolors="none",
        label="Candidate predictions",
    )

    if np.any(candidate_pareto):
        candidate_pareto_order = np.argsort(
            candidate_length[
                candidate_pareto
            ]
        )

        ax_candidates.plot(
            candidate_length[
                candidate_pareto
            ][candidate_pareto_order],
            candidate_num[
                candidate_pareto
            ][candidate_pareto_order],
            color="#1f77b4",
            linewidth=1.5,
            marker="o",
            markersize=3,
            label="Predicted Pareto front",
        )

    ax_candidates.scatter(
        [
            recommended[
                "pred_length_ratio"
            ]
        ],
        [
            recommended[
                "pred_num_ratio"
            ]
        ],
        marker="*",
        s=260,
        color="#ff7f00",
        edgecolors="black",
        linewidths=0.9,
        label="Recommended candidate",
        zorder=6,
    )

    ax_candidates.set_title(
        "(b) Predicted candidate objectives"
    )

    ax_candidates.set_xlabel(
        "predicted length_ratio"
    )

    ax_candidates.set_ylabel(
        "predicted num_ratio"
    )

    ax_candidates.grid(
        alpha=0.25
    )

    ax_candidates.legend(
        fontsize=8,
        loc="best",
    )


    # Panel 3: score versus local distance

    score_values = candidates[
        "balanced_score"
    ].to_numpy(
        dtype=float
    )

    distance_values = candidates[
        "nearest_anchor_dist"
    ].to_numpy(
        dtype=float
    )

    ax_score.scatter(
        distance_values,
        score_values,
        color="#756bb1",
        s=48,
        alpha=0.72,
        edgecolors="none",
        label="Candidate configurations",
    )

    recommended_distance = float(
        recommended[
            "nearest_anchor_dist"
        ]
    )

    recommended_score = float(
        recommended[
            "balanced_score"
        ]
    )

    ax_score.scatter(
        [
            recommended_distance
        ],
        [
            recommended_score
        ],
        marker="*",
        s=260,
        color="#e41a1c",
        edgecolors="black",
        linewidths=0.9,
        label="Recommended candidate",
        zorder=6,
    )

    ax_score.set_title(
        "(c) Acquisition score and local distance"
    )

    ax_score.set_xlabel(
        "nearest-anchor distance"
    )

    ax_score.set_ylabel(
        "balanced acquisition score"
    )

    ax_score.grid(
        alpha=0.25
    )

    ax_score.legend(
        fontsize=8,
        loc="best",
    )

    fig.suptitle(
        "Local Bayesian Optimization Based on Two Quality Objectives",
        fontsize=15,
        fontweight="bold",
    )

    fig.savefig(
        output_path,
        dpi=180,
        bbox_inches="tight",
    )

    plt.close(fig)

    print(
        f"Saved optimization plot to: {output_path}"
    )

# Main procedure

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Local Bayesian optimization for LASTZ parameters "
            "using length_ratio and num_ratio objectives"
        )
    )

    parser.add_argument(
        "--history",
        required=True,
        help=(
            "History CSV containing the eight optimization "
            "parameters, length_ratio, and num_ratio"
        ),
    )

    parser.add_argument(
        "--output",
        default="next_params.csv",
        help=(
            "Output CSV for ranked candidates"
        ),
    )

    parser.add_argument(
        "--topk",
        type=int,
        default=5,
        help=(
            "Number of top candidate configurations to print"
        ),
    )

    parser.add_argument(
        "--anchor-topk",
        type=int,
        default=5,
        help=(
            "Number of best historical observations "
            "used as local anchors"
        ),
    )

    parser.add_argument(
        "--len-tol",
        type=float,
        default=DEFAULT_LEN_TOL,
        help=(
            "Allowed predicted drop from the best historical "
            "length_ratio before progressive relaxation"
        ),
    )

    parser.add_argument(
        "--num-tol",
        type=float,
        default=DEFAULT_NUM_TOL,
        help=(
            "Allowed predicted drop from the local best "
            "num_ratio before progressive relaxation"
        ),
    )

    parser.add_argument(
        "--len-weight",
        type=float,
        default=DEFAULT_LEN_WEIGHT,
        help=(
            "Weight of length_ratio in the balanced score"
        ),
    )

    parser.add_argument(
        "--num-weight",
        type=float,
        default=DEFAULT_NUM_WEIGHT,
        help=(
            "Weight of num_ratio in the balanced score"
        ),
    )

    parser.add_argument(
        "--max-candidates",
        type=int,
        default=DEFAULT_MAX_CANDIDATES,
        help=(
            "Maximum number of candidates after local "
            "Cartesian-product generation"
        ),
    )

    parser.add_argument(
        "--plot",
        default=None,
        help=(
            "Output path for the optimization plot. "
            "Defaults to the output CSV path with .png suffix"
        ),
    )

    parser.add_argument(
        "--no-plot",
        action="store_true",
        help=(
            "Disable plot generation"
        ),
    )

    parser.add_argument(
        "--seed",
        type=int,
        default=0,
        help=(
            "Random seed used when limiting candidates"
        ),
    )

    args = parser.parse_args()

    if args.len_tol < 0:
        raise ValueError(
            "len-tol cannot be negative"
        )

    if args.num_tol < 0:
        raise ValueError(
            "num-tol cannot be negative"
        )

    history_path = Path(
        args.history
    )

    if not history_path.exists():
        raise FileNotFoundError(
            f"History CSV not found: {history_path}"
        )


    # 1. Read historical data

    df = read_data(
        history_path
    )

    if len(df) < 5:
        raise ValueError(
            "At least 5 historical observations "
            "are recommended"
        )


    # 2. Identify previously tested configurations

    seen_set = build_seen_set(
        df
    )


    # 3. Select Top-K anchors and local search center

    center, top = topk_center(
        df,
        args.anchor_topk,
    )


    # 4. Generate local candidate values

    value_map = make_candidate_values(
        df,
        center,
        top,
    )

    candidates = generate_candidates(
        value_map
    )

    if len(candidates) == 0:
        raise RuntimeError(
            "No candidates generated. "
            "Check the local parameter rules."
        )

    generated_candidate_count = len(
        candidates
    )

    candidates = limit_candidates(
        candidates,
        args.max_candidates,
        args.seed,
    )


    # 5. Remove previously tested configurations

    keep_mask = []

    for row in candidates[
        OPTIMIZE_COLS
    ].itertuples(
        index=False,
        name=None,
    ):
        key = tuple(
            int(round(float(value)))
            for value in row
        )

        keep_mask.append(
            key not in seen_set
        )

    candidates = candidates.loc[
        keep_mask
    ].drop_duplicates(
        subset=OPTIMIZE_COLS,
    ).reset_index(
        drop=True
    )

    if len(candidates) == 0:
        raise RuntimeError(
            "All generated candidate configurations "
            "have already been tested"
        )


    # 6. Feature encoding and standardization

    categories = get_category_values(
        df
    )

    X_raw = build_feature_matrix(
        df[
            OPTIMIZE_COLS
        ],
        categories,
    )

    C_raw = build_feature_matrix(
        candidates[
            OPTIMIZE_COLS
        ],
        categories,
    )

    scaler = StandardScaler()

    X = scaler.fit_transform(
        X_raw
    )

    C = scaler.transform(
        C_raw
    )


    # 7. Prepare the two objective variables

    y_len = df[
        "length_ratio"
    ].astype(
        float
    ).to_numpy()

    y_num = df[
        "num_ratio"
    ].astype(
        float
    ).to_numpy()


    # 8. Fit two independent GP models

    gp_len = fit_gp(
        X,
        y_len,
        seed=args.seed + 1,
    )

    gp_num = fit_gp(
        X,
        y_num,
        seed=args.seed + 2,
    )


    # 9. Predict candidate objective values

    mu_len, std_len = gp_len.predict(
        C,
        return_std=True,
    )

    mu_num, std_num = gp_num.predict(
        C,
        return_std=True,
    )

    out = candidates.copy()

    out[
        "pred_length_ratio"
    ] = mu_len

    out[
        "pred_num_ratio"
    ] = mu_num

    out[
        "std_length_ratio"
    ] = std_len

    out[
        "std_num_ratio"
    ] = std_num

    out[
        "nearest_anchor_dist"
    ] = nearest_anchor_distance(
        out,
        top,
    )


    # 10. Calculate expected improvement

    best_len = float(
        np.max(y_len)
    )

    # For num_ratio, use the best observed number ratio
    # among historical observations whose length_ratio is
    # close to the best observed length_ratio.
    near_best_length = (
        y_len >= best_len - args.len_tol
    )

    if np.any(
        near_best_length
    ):
        best_num_at_len = float(
            np.max(
                y_num[
                    near_best_length
                ]
            )
        )

    else:
        best_num_at_len = float(
            np.max(y_num)
        )

    out[
        "ei_length_ratio"
    ] = expected_improvement(
        mu=mu_len,
        sigma=std_len,
        best=best_len,
        xi=0.01,
    )

    out[
        "ei_num_ratio"
    ] = expected_improvement(
        mu=mu_num,
        sigma=std_num,
        best=best_num_at_len,
        xi=0.005,
    )


    # 11. Construct the two-objective acquisition score

    len_scale = robust_scale(
        y_len
    )

    num_scale = robust_scale(
        y_num
    )

    length_gain = (
        out[
            "ei_length_ratio"
        ].to_numpy()
        / len_scale
    )

    num_gain = (
        out[
            "ei_num_ratio"
        ].to_numpy()
        / num_scale
    )

    weights = normalize_weights(
        args.len_weight,
        args.num_weight,
    )

    length_weight = float(
        weights[0]
    )

    num_weight = float(
        weights[1]
    )

    out[
        "length_gain"
    ] = length_gain

    out[
        "num_gain"
    ] = num_gain

    out[
        "balanced_score"
    ] = (
        length_weight * length_gain
        + num_weight * num_gain
    )


    # 12. Progressive quality constraints

    filtered = pd.DataFrame()

    for multiplier in [
        1.0,
        2.0,
        4.0,
    ]:
        length_floor = (
            best_len
            - args.len_tol * multiplier
        )

        num_floor = (
            best_num_at_len
            - args.num_tol * multiplier
        )

        trial = out[
            (
                out[
                    "pred_length_ratio"
                ] >= length_floor
            )
            & (
                out[
                    "pred_num_ratio"
                ] >= num_floor
            )
        ].copy()

        if len(trial) > 0:
            filtered = trial
            break

    # If no candidate satisfies the progressively relaxed
    # quality constraints, retain all candidates.
    if len(filtered) == 0:
        filtered = out.copy()


    # 13. Rank candidate configurations

    filtered = filtered.sort_values(
        by=[
            "balanced_score",
            "pred_length_ratio",
            "pred_num_ratio",
            "nearest_anchor_dist",
        ],
        ascending=[
            False,
            False,
            False,
            True,
        ],
        kind="mergesort",
    ).reset_index(
        drop=True
    )


    # 14. Save ranked candidates

    filtered.to_csv(
        args.output,
        index=False,
    )


    # 15. Print optimization summary

    print(
        "Local optimization summary:"
    )

    print(
        f"History rows: {len(df)}"
    )

    print(
        "Generated candidates before "
        f"limiting: {generated_candidate_count}"
    )

    print(
        f"Candidates after filtering tested configurations: "
        f"{len(candidates)}"
    )

    print(
        f"Filtered candidates: {len(filtered)}"
    )

    print(
        "Objectives: maximize length_ratio "
        "and num_ratio"
    )

    print(
        "Weights: "
        f"length={length_weight:.3f}, "
        f"num={num_weight:.3f}"
    )

    print(
        "\nGP kernels:"
    )

    print(
        "length_ratio: "
        f"{gp_len.kernel_}"
    )

    print(
        "num_ratio:    "
        f"{gp_num.kernel_}"
    )


    # 16. Save diagnostic plot

    best = filtered.iloc[0]

    if not args.no_plot:
        plot_path = (
            Path(args.plot)
            if args.plot
            else Path(args.output).with_suffix(
                ".png"
            )
        )

        plot_optimization_summary(
            history=df,
            candidates=filtered,
            anchors=top,
            recommended=best,
            output_path=plot_path,
        )


    # 17. Print recommended configuration

    print(
        "\nNext recommended parameter set:"
    )

    print(
        best[
            OPTIMIZE_COLS
            + [
                "pred_length_ratio",
                "pred_num_ratio",
                "std_length_ratio",
                "std_num_ratio",
                "ei_length_ratio",
                "ei_num_ratio",
                "balanced_score",
                "nearest_anchor_dist",
            ]
        ].to_string()
    )

    print(
        "\nTop candidates:"
    )

    print(
        filtered.head(
            args.topk
        )[
            OPTIMIZE_COLS
            + [
                "pred_length_ratio",
                "pred_num_ratio",
                "std_length_ratio",
                "std_num_ratio",
                "ei_length_ratio",
                "ei_num_ratio",
                "balanced_score",
                "nearest_anchor_dist",
            ]
        ].to_string(
            index=False
        )
    )

    print(
        f"\nSaved ranked candidates to: "
        f"{args.output}"
    )

if __name__ == "__main__":
    main()