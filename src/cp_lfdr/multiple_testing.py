"""Support-line, BH, and data-dependent threshold utilities."""

import numpy as np
from numba import njit

def _validate_pvalues(p_values):
    values = np.asarray(p_values, dtype=float)
    if values.ndim != 1:
        raise ValueError("p-values must be a one-dimensional array.")
    if not np.all(np.isfinite(values)):
        raise ValueError("p-values must be finite.")
    if np.any((values < 0) | (values > 1)):
        raise ValueError("p-values must lie in [0, 1].")
    return values


def sl_procedure(p_values, alpha, lmbdas=None):
    """Return the rejected indices and count for the support-line procedure."""
    values = _validate_pvalues(p_values)
    if not 0 <= alpha <= 1:
        raise ValueError("alpha must lie in [0, 1].")

    m = len(values)
    if m == 0:
        return np.array([], dtype=int), 0

    order = np.argsort(values, kind="mergesort")
    ordered = values[order]
    thresholds = alpha * np.arange(1, m + 1) / m
    gaps = ordered - thresholds

    if np.all(gaps > 0):
        rejection_count = 0
    else:
        if lmbdas is None:
            valid = np.arange(m)
        else:
            lambda_values = np.asarray(lmbdas, dtype=float)
            if lambda_values.shape != values.shape:
                raise ValueError("lmbdas must have the same shape as p_values.")
            valid = np.flatnonzero(ordered <= lambda_values[order])

        if valid.size == 0:
            rejection_count = 0
        else:
            minimum = np.min(gaps[valid])
            minimizers = valid[gaps[valid] == minimum]
            rejection_count = int(minimizers[-1] + 1)

    rejected = order[:rejection_count]
    if lmbdas is not None and rejection_count:
        rejected = rejected[values[rejected] <= np.asarray(lmbdas)[rejected]]
        rejection_count = len(rejected)
    return rejected, rejection_count


def bh_procedure(p_values, alpha):
    """Return the rejected indices and count for Benjamini--Hochberg."""
    values = _validate_pvalues(p_values)
    if not 0 <= alpha <= 1:
        raise ValueError("alpha must lie in [0, 1].")

    m = len(values)
    if m == 0:
        return np.array([], dtype=int), 0

    order = np.argsort(values, kind="mergesort")
    ordered = values[order]
    passes = ordered <= alpha * np.arange(1, m + 1) / m
    rejection_count = int(np.flatnonzero(passes)[-1] + 1) if np.any(passes) else 0
    return order[:rejection_count], rejection_count


def BH(p_values, alpha):
    """Compatibility wrapper returning only the BH rejection count."""
    return bh_procedure(p_values, alpha)[1]


def compute_bFDR(rejected_indices, null_indices):
    """Boundary false-discovery indicator for one realization."""
    rejected = np.asarray(rejected_indices, dtype=int)
    if rejected.size == 0:
        return 0.0
    return float(rejected[-1] in set(np.asarray(null_indices, dtype=int)))


def compute_FDR(rejected_indices, null_indices):
    rejected = np.asarray(rejected_indices, dtype=int)
    if rejected.size == 0:
        return 0.0
    false_positives = np.intersect1d(rejected, null_indices)
    return len(false_positives) / len(rejected)


def compute_power(rejected_indices, null_indices, m):
    rejected = np.asarray(rejected_indices, dtype=int)
    number_nonnull = m - len(null_indices)
    if number_nonnull == 0:
        return np.nan
    true_positives = np.setdiff1d(rejected, null_indices)
    return len(true_positives) / number_nonnull

@njit
def _compute_alpha_hat_minus_sorted_numba(ps, order, alpha):
    """Compute the universal leave-one-out thresholds on sorted p-values.

    For each original index ``i``, the returned value is

        sup_{x,t} tau_{alpha (m+1) / m}^{alpha}(p_{-i}, x, t).

    The two free coordinates may be taken equal at the supremum.  The
    implementation checks each possible insertion interval using prefix and
    suffix minima of the support-line objective, giving O(m^2) total work.
    """
    m = len(ps)
    slope = alpha / m

    out = np.empty(m)
    prefix_min = np.empty(m)
    suffix_min = np.empty(m)

    INF = 1e300

    for removed in range(m):
        # Let r be ps with ps[removed] deleted, so len(r) = m - 1.
        # prefix_min[a] is the minimum objective among the first a entries
        # of r at their unchanged ranks.
        prefix_min[0] = INF
        for a in range(1, m):
            idx = a - 1
            r_value = ps[idx] if idx < removed else ps[idx + 1]
            prefix_min[a] = min(
                prefix_min[a - 1],
                r_value - slope * a,
            )

        # suffix_min[a] is the minimum objective among entries after insertion
        # interval a.  Two inserted values shift these entries by two ranks.
        suffix_min[m - 1] = INF
        for a in range(m - 2, -1, -1):
            idx = a
            r_value = ps[idx] if idx < removed else ps[idx + 1]
            rank_after_two_insertions = a + 3
            value = INF
            if r_value <= alpha:
                value = r_value - slope * rank_after_two_insertions
            suffix_min[a] = min(suffix_min[a + 1], value)

        best_threshold = 0.0
        for a in range(m):
            lower = 0.0
            if a > 0:
                idx = a - 1
                lower = ps[idx] if idx < removed else ps[idx + 1]

            upper = alpha
            if a < m - 1:
                idx = a
                next_value = ps[idx] if idx < removed else ps[idx + 1]
                upper = min(upper, next_value)

            objective_bound = min(
                0.0,
                prefix_min[a],
                suffix_min[a],
            )
            candidate = min(
                upper,
                slope * (a + 2) + objective_bound,
            )

            # Endpoints are included because alpha_hat_{-i} is a supremum.
            if candidate >= lower and candidate > best_threshold:
                best_threshold = candidate

        out[order[removed]] = max(0.0, best_threshold)

    return out


def compute_alpha_hat_minus_fast(p_values, alpha):
    """Return the exact universal threshold ``alpha_hat_{-i}``."""
    values = _validate_pvalues(p_values)
    if not 0 < alpha < 1:
        raise ValueError("alpha must lie in (0, 1).")

    order = np.argsort(values, kind="mergesort")
    sorted_values = values[order]

    if njit is None:
        raise ImportError("Install numba for this fast version: pip install numba")

    return _compute_alpha_hat_minus_sorted_numba(
        sorted_values,
        order,
        alpha,
    )


def compute_alpha_hat_thresholds(p_values, alpha):
    """Return ``alpha_hat_{-i}`` and its BH upper bound.

    The BH comparison is run on the vector obtained by deleting the largest
    p-value and inserting two zeros.  Its critical slope is ``alpha / m``,
    matching the augmented support-line construction.
    """
    values = _validate_pvalues(p_values)
    m = len(values)
    if m == 0:
        return np.array([], dtype=float), 0.0

    alpha_hat_minus = compute_alpha_hat_minus_fast(values, alpha)

    largest_idx = int(np.argmax(values))
    augmented = np.r_[np.delete(values, largest_idx), 0.0, 0.0]
    ordered = np.sort(augmented, kind="mergesort")
    passes = ordered <= alpha * np.arange(1, m + 2) / m
    rejection_count = (
        int(np.flatnonzero(passes)[-1] + 1) if np.any(passes) else 0
    )
    alpha_hat_bh = alpha * rejection_count / m

    return alpha_hat_minus, alpha_hat_bh
