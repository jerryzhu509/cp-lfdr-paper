"""Pooled permutation p-values and worst-block CDF calculations."""

from itertools import combinations

import numpy as np

from .multiple_testing import compute_alpha_hat_fast

def pooled_t_statistic(group_a, group_b):
    """Equal-variance two-sample t-statistic, oriented as group_b - group_a."""
    group_a = np.asarray(group_a, dtype=float)
    group_b = np.asarray(group_b, dtype=float)
    n_a, n_b = group_a.shape[0], group_b.shape[0]
    mean_a = group_a.mean(axis=0)
    mean_b = group_b.mean(axis=0)
    ss_a = np.sum((group_a - mean_a) ** 2, axis=0)
    ss_b = np.sum((group_b - mean_b) ** 2, axis=0)
    pooled_variance = (ss_a + ss_b) / (n_a + n_b - 2)
    standard_error = np.sqrt(pooled_variance * (1 / n_a + 1 / n_b))
    return (mean_b - mean_a) / standard_error


def pooled_upper_tail_pvalues(sorted_pool, statistics):
    """Pooled p-values using the inclusive upper tail P(T >= t)."""
    statistics = np.asarray(statistics)
    tail_counts = len(sorted_pool) - np.searchsorted(
        sorted_pool,
        statistics,
        side="left",
    )
    return tail_counts / len(sorted_pool)


def gamma_from_block_pvalues(block_pvalues, level):
    """Worst blockwise conditional CDF at a specified p-value level."""
    if level <= 0:
        return 0.0
    return float(np.max(np.mean(block_pvalues <= level, axis=0)))


def compute_gene_permutation_metrics(
    observed_statistics,
    data_matrix,
    n_group_a,
    alpha=0.2,
):
    """Enumerate all assignments without double-counting the observed assignment."""
    n_total = data_matrix.shape[0]
    all_indices = np.arange(n_total)
    assignments = list(combinations(range(n_total), n_group_a))
    permutation_statistics = np.empty(
        (len(assignments), data_matrix.shape[1]),
        dtype=float,
    )

    for row, assignment in enumerate(assignments):
        group_a_indices = np.asarray(assignment)
        group_b_indices = np.setdiff1d(
            all_indices,
            group_a_indices,
            assume_unique=True,
        )
        permutation_statistics[row] = pooled_t_statistic(
            data_matrix[group_a_indices],
            data_matrix[group_b_indices],
        )

    absolute_permutations = np.abs(permutation_statistics)
    absolute_observed = np.abs(observed_statistics)
    ordinary_pvalues = np.mean(
        absolute_permutations >= absolute_observed,
        axis=0,
    )

    # The observed assignment is already one of the exhaustive assignments.
    # Therefore the pooled grid has m * choose(12, 3), not one extra row.
    sorted_pool = np.sort(absolute_permutations.ravel())
    compound_pvalues = pooled_upper_tail_pvalues(
        sorted_pool,
        absolute_observed,
    )
    block_pvalues = pooled_upper_tail_pvalues(
        sorted_pool,
        absolute_permutations,
    )

    alpha_minus, alpha_bh = compute_alpha_hat_fast(compound_pvalues, alpha)
    gamma = gamma_from_block_pvalues(block_pvalues, alpha)
    gamma_minus = gamma_from_block_pvalues(block_pvalues, np.max(alpha_minus))
    gamma_bh = gamma_from_block_pvalues(block_pvalues, alpha_bh)

    return {
        "ordinary_pvalues": ordinary_pvalues,
        "compound_pvalues": compound_pvalues,
        "permutation_statistics": permutation_statistics,
        "block_pvalues": block_pvalues,
        "gamma": gamma,
        "gamma_minus": gamma_minus,
        "gamma_bh": gamma_bh,
    }

def compute_discrete_gammas(compound_pvalues, block_pvalues, alpha):
    """Return fixed, leave-one-out, and BH worst-block CDF values."""
    alpha_minus, alpha_bh = compute_alpha_hat_fast(
        compound_pvalues,
        alpha,
    )
    gamma = gamma_from_block_pvalues(block_pvalues, alpha)
    gamma_minus = gamma_from_block_pvalues(
        block_pvalues,
        np.max(alpha_minus),
    )
    gamma_bh = gamma_from_block_pvalues(
        block_pvalues,
        alpha_bh,
    )
    return gamma, gamma_minus, gamma_bh
