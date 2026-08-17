import numpy as np
import pandas as pd

from cp_lfdr.rd import (
    assignment_indices,
    compute_pvalues_window,
    compute_statistic,
    construct_year_specific_midpoint_base,
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


def test_exact_local_pvalue_includes_the_observed_assignment():
    frame = pd.DataFrame(
        {
            "group": ["a"] * 11,
            "distance": [
                -0.05,
                -0.05,
                -0.02,
                -0.01,
                0.01,
                0.02,
                0.02,
                0.03,
                0.03,
                0.04,
                0.04,
            ],
            "outcome": [
                6.70,
                6.56,
                7.04,
                7.00,
                8.16,
                7.10,
                7.06,
                7.07,
                7.12,
                7.37,
                7.16,
            ],
        }
    )
    results, _ = compute_pvalues_window(
        frame,
        "outcome",
        level="group",
        alpha=0.2,
        test_stat="welch_p",
        dcp_draws=5,
        rng=np.random.RandomState(0),
        distance_column="distance",
        bw_l=0.05,
        bw_r=0.05,
    )
    np.testing.assert_allclose(results.loc[0, "local_p"], 1 / 330)


def test_midpoints_are_translation_invariant_and_pool_the_same_pair():
    rows = []
    cutoff_positions = {
        1: {10: 0.0, 3: 0.4, 7: 1.0},
        2: {8: 99.5, 10: 100.0, 3: 100.4, 7: 101.0},
    }
    student_scores = {
        1: {"a": 0.15, "b": 0.25},
        2: {"c": 100.15, "d": 100.25},
    }
    for year, scores in student_scores.items():
        for student, score in scores.items():
            for cutoff_id, cutoff in cutoff_positions[year].items():
                rows.append(
                    {
                        "town": 1,
                        "Y": year,
                        "z": cutoff_id,
                        "sid2": student,
                        "dzag": score - cutoff,
                        "bcg": float(ord(student)),
                    }
                )
    frame = pd.DataFrame(rows)

    _, year_gaps, paired = construct_year_specific_midpoint_base(frame)
    same_pair = year_gaps[year_gaps["placebo_id"] == "1_placebo3_10"]
    assert set(same_pair["gap_rank"]) == {1, 2}
    assert set(same_pair["z_low"]) == {3}
    assert set(same_pair["z_high"]) == {10}
    np.testing.assert_allclose(
        sorted(
            paired.loc[
                paired["placebo_id"] == "1_placebo3_10",
                "placebo_distance",
            ]
        ),
        [-0.05, -0.05, 0.05, 0.05],
    )
