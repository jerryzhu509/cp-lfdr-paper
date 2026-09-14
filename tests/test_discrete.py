from math import comb

import numpy as np

from cp_lfdr.discrete import (
    compute_gene_permutation_metrics,
    gamma_from_block_pvalues,
    pooled_t_statistic,
    pooled_upper_tail_pvalues,
)


def test_pooled_pvalues_use_inclusive_upper_tail():
    pool = np.array([1.0, 1.0, 2.0, 3.0])
    pvalues = pooled_upper_tail_pvalues(pool, np.array([1.0, 2.0, 3.0]))
    np.testing.assert_allclose(pvalues, [1.0, 0.5, 0.25])


def test_gamma_uses_block_pvalues_in_significant_tail():
    blocks = np.array(
        [
            [0.01, 0.30],
            [0.02, 0.40],
            [0.90, 0.50],
        ]
    )
    assert gamma_from_block_pvalues(blocks, 0.05) == 2 / 3


def test_gene_pool_does_not_append_observation_twice():
    data = np.array(
        [
            [0.0, 0.3, 1.0],
            [0.1, 0.2, 0.9],
            [1.0, 0.4, 0.2],
            [1.1, 0.5, 0.1],
        ]
    )
    observed = pooled_t_statistic(data[:2], data[2:])
    result = compute_gene_permutation_metrics(
        observed,
        data,
        n_group_a=2,
        alpha=0.2,
    )
    expected_assignments = comb(4, 2)
    assert result["permutation_statistics"].shape == (
        expected_assignments,
        data.shape[1],
    )
    assert result["block_pvalues"].shape[0] == expected_assignments
