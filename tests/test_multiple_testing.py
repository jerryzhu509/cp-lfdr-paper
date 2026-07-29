import numpy as np
import pytest

from cp_lfdr.multiple_testing import (
    bh_procedure,
    compute_alpha_hat_fast,
    sl_procedure,
)


def brute_support_line_count(values, alpha):
    ordered = np.sort(values)
    gaps = ordered - alpha * np.arange(1, len(values) + 1) / len(values)
    if np.all(gaps > 0):
        return 0
    return int(np.flatnonzero(gaps == gaps.min())[-1] + 1)


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


def test_fast_alpha_minus_matches_dense_small_problem():
    values = np.array([0.01, 0.08, 0.2, 0.7])
    alpha = 0.1
    actual, _ = compute_alpha_hat_fast(values, alpha)
    grid = np.unique(np.r_[np.linspace(0, 1, 5001), values])
    expected = []
    for index in range(len(values)):
        remaining = np.delete(values, index)
        largest_count = max(
            brute_support_line_count(np.r_[remaining, value], alpha)
            for value in grid
        )
        expected.append(alpha * largest_count / len(values))
    np.testing.assert_allclose(actual, expected)
