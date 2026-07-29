import numpy as np

from cp_lfdr.continuous import (
    compute_gamma_t_test,
    simulate_continuous_data,
)


def test_simulation_uses_symmetric_sigma_interval():
    simulated = simulate_continuous_data(
        m=500,
        sample_size=5,
        pi0=0.8,
        effect=0.5,
        tau=0.3,
        alpha=0.1,
        rng=np.random.default_rng(1),
    )
    assert simulated["sigma"].min() >= 0.7
    assert simulated["sigma"].max() <= 1.3
    assert simulated["sigma"].min() < 0.8
    assert simulated["sigma"].max() > 1.2


def test_gamma_is_zero_at_zero_level():
    oracle, data_driven = compute_gamma_t_test(
        np.array([0.8, 1.2]),
        np.array([4.0, 7.0]),
        sample_size=5,
        alpha=0,
    )
    np.testing.assert_array_equal(oracle, 0)
    np.testing.assert_array_equal(data_driven, 0)
