import numpy as np
import pandas as pd

from cp_lfdr.rd import (
    assignment_indices,
    compute_pvalues_window,
    compute_statistic,
)


def test_difference_in_means_matches_manual_calculation():
    values = np.array([1.0, 2.0, 4.0, 6.0])
    assignments = np.array([[0, 1], [0, 2], [1, 3]])
    observed, distribution = compute_statistic(
        values,
        left_n=2,
        right_n=2,
        assignment_indices=assignments,
    )
    assert observed == 3.5
    np.testing.assert_allclose(distribution[0], 3.5)


def test_seed_zero_matches_original_legacy_assignment_sampler():
    actual, is_exact = assignment_indices(
        2,
        4,
        draws=5,
        rng=np.random.RandomState(0),
    )
    legacy_rng = np.random.RandomState(0)
    expected = np.argsort(
        legacy_rng.rand(5, 6),
        axis=1,
    )[:, :2]
    assert not is_exact
    np.testing.assert_array_equal(actual, expected)


def test_rd_pool_contains_controls_plus_one_observation_per_block():
    frame = pd.DataFrame(
        {
            "group": np.repeat(["a", "b", "c"], 4),
            "distance": np.tile([-0.2, -0.1, 0.1, 0.2], 3),
            "outcome": [
                1.0, 1.2, 1.8, 2.0,
                0.9, 1.1, 1.4, 1.7,
                1.3, 1.4, 1.6, 1.9,
            ],
        }
    )
    draws = 5
    (
        results,
        _,
        pool,
        block_matrix,
        block_pvalues,
        _,
        _,
    ) = compute_pvalues_window(
        frame,
        "outcome",
        level="group",
        alpha=0.2,
        test_stat="diff_means",
        dcp_draws=draws,
        rng=np.random.default_rng(4),
        distance_column="distance",
        bw_l=0.3,
        bw_r=0.3,
        gamma_hat=True,
        return_structures=True,
    )
    assert len(results) == 3
    assert block_matrix.shape == (draws + 1, 3)
    assert block_pvalues.shape == block_matrix.shape
    assert len(pool) == 3 * (draws + 1)
