"""Local-randomization utilities for the RD application."""

from itertools import combinations
from math import comb

import numpy as np
import pandas as pd
from scipy.stats import t as student_t

from .discrete import (
    compute_discrete_gammas,
    gamma_from_block_pvalues,
    pooled_upper_tail_pvalues,
)

def compute_statistic(values, left_n, right_n, assignment_indices, method="diff_means"):
    """Compute the observed statistic and its assignment distribution."""
    values = np.asarray(values, dtype=float).reshape(-1)
    if left_n <= 0 or right_n <= 0:
        raise ValueError("Both sides of the cutoff must contain observations.")

    left_sums = np.sum(values[assignment_indices], axis=1)
    right_sums = np.sum(values) - left_sums
    differences = right_sums / right_n - left_sums / left_n
    observed_difference = values[left_n:].mean() - values[:left_n].mean()

    if method == "diff_means":
        return observed_difference, differences

    if left_n <= 1 or right_n <= 1:
        raise ValueError(f"{method} requires at least two observations per side.")

    squared = values**2
    left_squares = np.sum(squared[assignment_indices], axis=1)
    right_squares = np.sum(squared) - left_squares
    left_variances = np.maximum(
        (left_squares - left_sums**2 / left_n) / (left_n - 1),
        0,
    )
    right_variances = np.maximum(
        (right_squares - right_sums**2 / right_n) / (right_n - 1),
        0,
    )
    left_terms = left_variances / left_n
    right_terms = right_variances / right_n
    standard_errors = np.sqrt(left_terms + right_terms)
    permutation_t = np.divide(
        differences,
        standard_errors,
        out=np.full_like(differences, np.nan),
        where=standard_errors > 0,
    )

    observed_left_variance = np.var(values[:left_n], ddof=1)
    observed_right_variance = np.var(values[left_n:], ddof=1)
    observed_left_term = observed_left_variance / left_n
    observed_right_term = observed_right_variance / right_n
    observed_se = np.sqrt(observed_left_term + observed_right_term)
    observed_t = observed_difference / observed_se

    if method == "t_stat":
        return observed_t, permutation_t
    if method != "welch_p":
        raise ValueError(f"Unknown method: {method}")

    permutation_df_denominator = (
        left_terms**2 / (left_n - 1)
        + right_terms**2 / (right_n - 1)
    )
    permutation_df = np.divide(
        (left_terms + right_terms) ** 2,
        permutation_df_denominator,
        out=np.full_like(permutation_df_denominator, np.nan),
        where=permutation_df_denominator > 0,
    )
    permutation_pvalues = student_t.sf(permutation_t, df=permutation_df)

    observed_df = (
        (observed_left_term + observed_right_term) ** 2
        / (
            observed_left_term**2 / (left_n - 1)
            + observed_right_term**2 / (right_n - 1)
        )
    )
    observed_pvalue = student_t.sf(observed_t, df=observed_df)

    # Negating one-sided p-values makes a smaller p-value a larger statistic.
    return -observed_pvalue, -permutation_pvalues


def assignment_indices(
    left_n,
    right_n,
    *,
    draws,
    rng,
    exact=False,
    exact_limit=200_000,
):
    """Return assignment indices and whether they enumerate the full design."""
    n_total = left_n + right_n
    number_assignments = comb(n_total, left_n)
    if exact and number_assignments <= exact_limit:
        indices = np.asarray(
            list(combinations(range(n_total), left_n)),
            dtype=np.intp,
        )
        return indices, True

    random_scores = rng.random((draws, n_total))
    # Keep the original notebook's full-sort algorithm.  With RandomState,
    # this reproduces both the sampled assignments and their index order.
    indices = np.argsort(random_scores, axis=1)[:, :left_n]
    return indices, False



def compute_pvalues_window(
    data,
    outcome_var,
    *,
    level="town",
    alpha=0.01,
    test_stat="diff_means",
    dcp_draws=70,
    rng,
    distance_column="dzag",
    bw_l=0.1,
    bw_r=0.1,
    gamma_hat=False,
    return_structures=False,
):
    """Compute exact local and pooled Monte Carlo permutation p-values."""
    control_statistics = {}
    observed_statistics = {}
    local_pvalues = {}
    exact_statistics = {}

    for level_id, group in data.groupby(level):
        working = group.copy()
        left = working[
            (working[distance_column] < 0)
            & (working[distance_column] >= -bw_l)
        ]
        right = working[
            (working[distance_column] > 0)
            & (working[distance_column] <= bw_r)
        ]
        left_n, right_n = len(left), len(right)
        values = pd.concat([left, right])[outcome_var].to_numpy()

        dcp_indices, _ = assignment_indices(
            left_n,
            right_n,
            draws=dcp_draws,
            rng=rng,
            exact=False,
        )
        observed, controls = compute_statistic(
            values,
            left_n,
            right_n,
            dcp_indices,
            method=test_stat,
        )
        observed_statistics[level_id] = observed
        control_statistics[level_id] = controls

        local_indices, is_exact = assignment_indices(
            left_n,
            right_n,
            draws=20_000,
            rng=rng,
            exact=True,
        )
        observed_local, local_distribution = compute_statistic(
            values,
            left_n,
            right_n,
            local_indices,
            method=test_stat,
        )
        if is_exact:
            local_pvalue = np.mean(local_distribution >= observed_local)
        else:
            local_pvalue = (
                np.sum(local_distribution >= observed_local) + 1
            ) / (len(local_distribution) + 1)
        local_pvalues[level_id] = local_pvalue
        exact_statistics[level_id] = local_distribution

    if not observed_statistics:
        raise ValueError("No groups remain after applying the window restrictions.")

    level_ids = list(observed_statistics)
    block_matrix = np.column_stack(
        [
            np.append(
                control_statistics[level_id],
                observed_statistics[level_id],
            )
            for level_id in level_ids
        ]
    )
    sorted_pool = np.sort(block_matrix.ravel())
    observed_array = np.array(
        [observed_statistics[level_id] for level_id in level_ids]
    )
    compound_pvalues = pooled_upper_tail_pvalues(
        sorted_pool,
        observed_array,
    )
    block_pvalues = pooled_upper_tail_pvalues(
        sorted_pool,
        block_matrix,
    )

    results = pd.DataFrame(
        {
            level: level_ids,
            "local_p": [local_pvalues[level_id] for level_id in level_ids],
            "compound_p": compound_pvalues,
        }
    )
    if gamma_hat:
        gammas = compute_discrete_gammas(
            compound_pvalues,
            block_pvalues,
            alpha,
        )
    else:
        gammas = (gamma_from_block_pvalues(block_pvalues, alpha),)

    if return_structures:
        return (
            results,
            gammas,
            sorted_pool,
            block_matrix,
            block_pvalues,
            exact_statistics,
            observed_statistics,
        )
    return results, gammas

def filter_local_designs(
    data,
    *,
    level="town_z",
    distance_column="dzag",
    bandwidth=0.051,
    minimum_per_side=4,
    maximum_total=20,
):
    """Apply the balance and sample-size restrictions in the paper."""
    frame = data.copy()
    grouped = frame.groupby(level)[distance_column]
    frame["_minimum_distance"] = grouped.transform("min")
    frame["_maximum_distance"] = grouped.transform("max")
    frame = frame[
        (frame["_minimum_distance"] <= -bandwidth)
        & (frame["_maximum_distance"] >= bandwidth)
        & (frame[distance_column].abs() <= bandwidth)
    ].copy()

    grouped = frame.groupby(level)[distance_column]
    frame["_left_size"] = grouped.transform(
        lambda values: (values < 0).sum()
    )
    frame["_right_size"] = grouped.transform(
        lambda values: (values > 0).sum()
    )
    frame["_total_size"] = (
        frame["_left_size"] + frame["_right_size"]
    )
    return frame[
        (frame["_left_size"] >= minimum_per_side)
        & (frame["_right_size"] >= minimum_per_side)
        & (frame["_total_size"] <= maximum_total)
    ].copy()


def construct_placebo_design(
    data,
    cutoff,
    *,
    distance_column="dzag",
    output_column="placebo_distance",
):
    """Recenter the running variable at a placebo cutoff."""
    frame = data.copy()
    frame[output_column] = frame[distance_column] - cutoff
    return frame


def filter_placebo_designs(
    data,
    *,
    level="town_z",
    distance_column="placebo_distance",
    bandwidth=0.051,
    minimum_per_side=4,
    maximum_total=20,
):
    """Apply the placebo window and sample-size restrictions."""
    frame = data[
        data[distance_column].abs() <= bandwidth
    ].copy()
    grouped = frame.groupby(level)[distance_column]
    frame["_left_size"] = grouped.transform(
        lambda values: (values < 0).sum()
    )
    frame["_right_size"] = grouped.transform(
        lambda values: (values > 0).sum()
    )
    frame["_total_size"] = (
        frame["_left_size"] + frame["_right_size"]
    )
    return frame[
        (frame["_left_size"] >= minimum_per_side)
        & (frame["_right_size"] >= minimum_per_side)
        & (frame["_total_size"] <= maximum_total)
    ].copy()
