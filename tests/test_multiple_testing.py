import numpy as np
import pytest

from cp_lfdr.multiple_testing import (
    bh_procedure,
    compute_alpha_hat_minus_fast,
    compute_alpha_hat_thresholds,
    sl_procedure,
)


def test_bh_returns_indices_and_count():
    rejected, count = bh_procedure([0.001, 0.02, 0.8, 0.9], 0.05)
    assert count == 2
    assert set(rejected) == {0, 1}


def test_support_line_has_consistent_api():
    rejected, count = sl_procedure([0.01, 0.02, 0.9], 0.1)
    assert count == len(rejected)


def test_invalid_pvalue_is_rejected():
    with pytest.raises(ValueError):
        sl_procedure([0.1, 1.2], 0.1)


def modified_support_line_threshold(values, alpha, truncation):
    ordered = np.sort(values)
    candidates = np.r_[0.0, ordered]
    ranks = np.arange(len(candidates))
    valid = candidates <= truncation
    gaps = candidates - alpha * ranks / len(values)
    minimum = np.min(gaps[valid])
    return candidates[np.flatnonzero(valid & (gaps == minimum))[-1]]


def test_fast_alpha_hat_minus_matches_dense_small_problem():
    values = np.array([0.01, 0.08, 0.2, 0.7])
    alpha = 0.1
    actual = compute_alpha_hat_minus_fast(values, alpha)
    augmented_alpha = alpha * (len(values) + 1) / len(values)
    grid = np.unique(np.r_[np.linspace(0, alpha, 101), values, alpha])
    expected = []
    for index in range(len(values)):
        remaining = np.delete(values, index)
        largest_threshold = max(
            modified_support_line_threshold(
                np.r_[remaining, first, second],
                augmented_alpha,
                alpha,
            )
            for first in grid
            for second in grid
        )
        expected.append(largest_threshold)
    np.testing.assert_allclose(actual, expected, atol=alpha / 100)


def test_bh_posthoc_threshold_is_an_upper_bound():
    values = np.array([0.01, 0.08, 0.2, 0.7])
    alpha_hat_minus, alpha_hat_bh = compute_alpha_hat_thresholds(
        values,
        alpha=0.1,
    )
    assert np.max(alpha_hat_minus) <= alpha_hat_bh


def test_bh_posthoc_threshold_is_not_capped_at_alpha():
    _, alpha_hat_bh = compute_alpha_hat_thresholds(
        np.zeros(4),
        alpha=0.1,
    )
    assert alpha_hat_bh == pytest.approx(0.125)
