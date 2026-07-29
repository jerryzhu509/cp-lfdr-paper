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
    """
    ps: sorted p-values
    order: original indices of sorted p-values

    Computes alpha_hat_minus[i] = alpha * max_x R_SL(p_{-i}, x) / m.
    """
    m = len(ps)
    a = alpha / m

    out = np.empty(m)
    right_val = np.empty(m + 2)
    right_arg = np.empty(m + 2, dtype=np.int64)

    INF = 1e300

    for t in range(m):
        # r is ps with ps[t] removed.
        # r[j] in 0-based indexing is:
        #     ps[j]     if j < t
        #     ps[j + 1] if j >= t

        # Compute suffix minima for the right side.
        #
        # For insertion position s = 1,...,m, the right-side gaps are
        #
        #     r_j - alpha * (j + 1) / m,  j >= s,
        #
        # where j is 1-based in the math.
        best_val = INF
        best_arg = -1

        for s in range(m, 0, -1):
            if s <= m - 1:
                # r_s in 1-based indexing
                idx = s - 1
                if idx < t:
                    r_s = ps[idx]
                else:
                    r_s = ps[idx + 1]

                gap = r_s - a * (s + 1)
                arg = s + 1

                # Scanning from right to left, keep largest minimizer.
                # Since new arg is smaller, update only on strict improvement.
                if gap < best_val:
                    best_val = gap
                    best_arg = arg

            right_val[s] = best_val
            right_arg[s] = best_arg

        # Scan insertion positions s from left to right.
        left_val = INF
        left_arg = -1

        R_max = 0

        for s in range(1, m + 1):
            # Non-x minimum among all positions except insertion position s.
            C = left_val
            C_arg = left_arg

            rv = right_val[s]
            ra = right_arg[s]

            if rv < C:
                C = rv
                C_arg = ra
            elif rv == C and ra > C_arg:
                C_arg = ra

            # Candidate where x is not the minimizing position.
            if C <= 0.0 and C_arg > R_max:
                R_max = C_arg

            # Candidate where x itself is the minimizing position.
            #
            # For insertion position s, the smallest feasible x is:
            #
            #     0          if s = 1
            #     r_{s - 1}  otherwise.
            #
            # x can make position s the minimizer iff
            #
            #     lower_x - alpha s / m <= min(C, 0).
            #
            if s == 1:
                lower_x = 0.0
            else:
                idx = s - 2
                if idx < t:
                    lower_x = ps[idx]
                else:
                    lower_x = ps[idx + 1]

            x_gap_lowest = lower_x - a * s

            threshold = C
            if threshold > 0.0:
                threshold = 0.0

            if x_gap_lowest <= threshold and s > R_max:
                R_max = s

            # Update left prefix for the next insertion position.
            if s <= m - 1:
                idx = s - 1
                if idx < t:
                    r_s = ps[idx]
                else:
                    r_s = ps[idx + 1]

                gap = r_s - a * s
                arg = s

                # Scanning left to right, update on ties to keep largest minimizer.
                if gap <= left_val:
                    left_val = gap
                    left_arg = arg

        original_index = order[t]
        out[original_index] = alpha * R_max / m

    return out


def compute_alpha_hat_minus_fast(cp_values, alpha):
    cp_values = np.asarray(cp_values, dtype=float)

    order = np.argsort(cp_values, kind="mergesort")
    ps = cp_values[order]

    if njit is None:
        raise ImportError("Install numba for this fast version: pip install numba")

    return _compute_alpha_hat_minus_sorted_numba(ps, order, alpha)

def compute_alpha_hat_fast(cp_values, alpha):
    cp_values = np.asarray(cp_values, dtype=float)
    m = len(cp_values)

    alpha_hat_minus = compute_alpha_hat_minus_fast(cp_values, alpha)

    # BH alpha-hat: remove largest p-value and insert 0.
    largest_idx = np.argmax(cp_values)
    p_bh = cp_values.copy()
    p_bh[largest_idx] = 0.0

    R_BH = BH(p_bh, alpha)
    alpha_hat_BH = alpha * R_BH / m

    return alpha_hat_minus, alpha_hat_BH
