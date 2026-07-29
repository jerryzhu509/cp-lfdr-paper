"""Methods used in the density-compound p-value analyses."""

from .multiple_testing import (
    BH,
    bh_procedure,
    compute_FDR,
    compute_alpha_hat_fast,
    compute_bFDR,
    compute_power,
    sl_procedure,
)

__all__ = [
    "BH",
    "bh_procedure",
    "compute_FDR",
    "compute_alpha_hat_fast",
    "compute_bFDR",
    "compute_power",
    "sl_procedure",
]
