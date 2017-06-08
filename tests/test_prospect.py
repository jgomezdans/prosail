import numpy as np

from prosail import run_prospect

def test_reflectance_prospect5():
    # runs prospect and compares to online prospect run
    w, true_refl, true_trans = np.loadtxt(
            "tests/prospect5_spectrum.txt", unpack=True)

    w, refl, trans = run_prospect(2.1, 40, 10., 0.1, 
                0.015, 0.009, prospect_version="5")
    assert np.allclose( true_refl, refl, atol=1e-4)

def test_transmittance_prospect5():
    # runs prospect and compares to online prospect run
    w, true_refl, true_trans = np.loadtxt(
            "tests/prospect5_spectrum.txt", unpack=True)

    w, refl, trans = run_prospect(2.1, 40, 10., 0.1, 
                0.015, 0.009, prospect_version="5")
    assert np.allclose( true_trans, trans, atol=1e-4)


