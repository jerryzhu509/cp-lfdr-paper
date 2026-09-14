"""Methods used in the density-compound p-value analyses."""

from .multiple_testing import (
    BH,
    bh_procedure,
    compute_FDR,
    compute_bFDR,
    compute_alpha_hat_minus_fast,
    compute_alpha_hat_thresholds,
    compute_power,
    sl_procedure,
)

__all__ = [
    "BH",
    "bh_procedure",
    "compute_FDR",
    "compute_bFDR",
    "compute_alpha_hat_minus_fast",
    "compute_alpha_hat_thresholds",
    "compute_power",
    "sl_procedure",
]
