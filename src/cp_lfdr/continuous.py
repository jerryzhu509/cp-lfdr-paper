"""Continuous Gaussian-means simulations from the manuscript."""

import numpy as np
from scipy import stats
from scipy.optimize import brentq

from .multiple_testing import (
    compute_FDR,
    compute_alpha_hat_fast,
    compute_bFDR,
    compute_power,
    sl_procedure,
)

def _solve_tail_cutoff(tail_probability, alpha, initial_upper=10.0):
    """Solve tail_probability(x) = alpha with an adaptively expanded bracket."""
    if alpha <= 0:
        return np.inf
    if alpha >= 1:
        return 0.0

    upper = float(initial_upper)
    objective = lambda x: float(tail_probability(x) - alpha)
    while objective(upper) > 0:
        upper *= 2
        if upper > 1e6:
            raise RuntimeError("Could not bracket the requested tail probability.")
    return brentq(objective, 0.0, upper)


def compute_gamma_t_test(sigma, u2_sum, sample_size, alpha):
    """Compute oracle and data-driven conditional null CDF values."""
    sigma = np.asarray(sigma, dtype=float)
    u2_sum = np.asarray(u2_sum, dtype=float)
    if alpha <= 0:
        zeros = np.zeros_like(sigma)
        return zeros, zeros.copy()

    oracle_tail = lambda x: np.mean(
        2 * stats.norm.cdf(-np.abs(x) / sigma)
    )
    oracle_cutoff = _solve_tail_cutoff(oracle_tail, alpha)
    gamma_oracle = 2 * stats.norm.cdf(-oracle_cutoff / sigma)

    mean_squares = u2_sum / sample_size
    data_tail = lambda x: np.mean(
        1
        - stats.beta.cdf(
            x**2 / mean_squares,
            a=0.5,
            b=(sample_size - 1) / 2,
        )
    )
    data_cutoff = _solve_tail_cutoff(data_tail, alpha)
    gamma_data = 1 - stats.beta.cdf(
        data_cutoff**2 / mean_squares,
        a=0.5,
        b=(sample_size - 1) / 2,
    )
    return gamma_oracle, gamma_data


def simulate_continuous_data(
    m,
    sample_size,
    pi0,
    effect,
    tau,
    alpha,
    rng,
    sigma_center=1.0,
):
    """Simulate the Gaussian-means setting stated in the manuscript.

    Standard deviations are sampled from U(sigma_center - tau,
    sigma_center + tau), for both null and non-null hypotheses.
    """
    if tau < 0 or sigma_center - tau <= 0:
        raise ValueError("tau must be nonnegative and sigma_center - tau positive.")

    m0 = int(pi0 * m)
    m1 = m - m0
    sigma = rng.uniform(sigma_center - tau, sigma_center + tau, size=m)
    means = np.concatenate([np.zeros(m0), np.full(m1, effect)])
    observations = rng.normal(
        loc=means[:, None],
        scale=sigma[:, None],
        size=(m, sample_size),
    )

    sample_means = observations.mean(axis=1)
    sample_variances = observations.var(axis=1, ddof=1)
    sum_squares = np.sum(observations**2, axis=1)
    t_statistics = sample_means / np.sqrt(sample_variances / sample_size)
    p_t = 2 * stats.t.sf(np.abs(t_statistics), df=sample_size - 1)

    p_oracle = np.mean(
        2
        * stats.norm.cdf(
            -np.abs(sample_means)[None, :]
            / (sigma[:, None] / np.sqrt(sample_size))
        ),
        axis=0,
    )

    p_data = np.mean(
        1
        - stats.beta.cdf(
            sample_means[None, :] ** 2
            / (sum_squares[:, None] / sample_size),
            a=0.5,
            b=(sample_size - 1) / 2,
        ),
        axis=0,
    )

    westfall_weights = m * sum_squares / np.sum(sum_squares)
    p_weighted = np.minimum(p_t / westfall_weights, 1)

    oracle_eweights = m * (sum_squares / sample_size) / np.sum(sigma**2)
    p_eweighted = np.minimum(p_t / oracle_eweights, 1)
    null_indices = np.arange(m0)

    gamma_oracle, gamma_data = compute_gamma_t_test(
        sigma,
        sum_squares,
        sample_size,
        alpha,
    )
    return {
        "p_values": {
            "t-test": p_t,
            "Oracle": p_oracle,
            "Data-driven": p_data,
            "Weighted": p_weighted,
            "E-weighted": p_eweighted,
        },
        "gamma": {
            "Oracle": gamma_oracle,
            "Data-driven": gamma_data,
        },
        "sigma": sigma,
        "sum_squares": sum_squares,
        "null_indices": null_indices,
    }

def run_continuous_setting(
    *,
    m,
    sample_size,
    pi0,
    effect,
    tau,
    alpha,
    draws,
    rng,
    alpha_grid=None,
):
    """Average rejection metrics and theoretical quantities over draws."""
    if alpha_grid is None:
        alpha_grid = np.array([alpha])

    methods = ("t-test", "Oracle", "Data-driven", "Weighted", "E-weighted")
    compound_methods = ("Oracle", "Data-driven")
    metrics = {
        name: {"bFDR": 0.0, "FDR": 0.0, "Power": 0.0}
        for name in methods
    }
    gamma_grid = {
        name: {float(level): 0.0 for level in alpha_grid}
        for name in compound_methods
    }
    gamma_minus = {name: 0.0 for name in compound_methods}
    gamma_bh = {name: 0.0 for name in compound_methods}

    for _ in range(draws):
        simulated = simulate_continuous_data(
            m=m,
            sample_size=sample_size,
            pi0=pi0,
            effect=effect,
            tau=tau,
            alpha=alpha,
            rng=rng,
        )
        null_indices = simulated["null_indices"]

        for name, values in simulated["p_values"].items():
            rejected, _ = sl_procedure(values, alpha)
            metrics[name]["bFDR"] += compute_bFDR(rejected, null_indices)
            metrics[name]["FDR"] += compute_FDR(rejected, null_indices)
            metrics[name]["Power"] += compute_power(rejected, null_indices, m)

        for level in alpha_grid:
            gamma_oracle, gamma_data = compute_gamma_t_test(
                simulated["sigma"],
                simulated["sum_squares"],
                sample_size,
                float(level),
            )
            gamma_grid["Oracle"][float(level)] += np.max(gamma_oracle)
            gamma_grid["Data-driven"][float(level)] += np.max(gamma_data)

        for name in compound_methods:
            values = simulated["p_values"][name]
            alpha_minus, alpha_bh = compute_alpha_hat_fast(values, alpha)

            gamma_oracle_bh, gamma_data_bh = compute_gamma_t_test(
                simulated["sigma"],
                simulated["sum_squares"],
                sample_size,
                alpha_bh,
            )
            gamma_bh[name] += np.max(
                gamma_oracle_bh if name == "Oracle" else gamma_data_bh
            )

            gamma_minus_values = np.empty(m)
            for level in np.unique(alpha_minus):
                gamma_oracle_level, gamma_data_level = compute_gamma_t_test(
                    simulated["sigma"],
                    simulated["sum_squares"],
                    sample_size,
                    float(level),
                )
                selected = alpha_minus == level
                gamma_minus_values[selected] = (
                    gamma_oracle_level[selected]
                    if name == "Oracle"
                    else gamma_data_level[selected]
                )
            gamma_minus[name] += np.max(gamma_minus_values)

    for name in methods:
        for metric in metrics[name]:
            metrics[name][metric] /= draws
    for name in compound_methods:
        gamma_minus[name] /= draws
        gamma_bh[name] /= draws
        for level in gamma_grid[name]:
            gamma_grid[name][level] /= draws

    return {
        "metrics": metrics,
        "gamma_grid": gamma_grid,
        "gamma_minus": gamma_minus,
        "gamma_bh": gamma_bh,
    }
