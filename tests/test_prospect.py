import numpy as np
import os
from prosail import run_prospect
from pytest import fixture
from distutils import dir_util


@fixture
def datadir(tmpdir, request):
    '''
    Fixture responsible for locating the test data directory and copying it
    into a temporary directory.
    Taken from  http://www.camillescott.org/2016/07/15/travis-pytest-scipyconf/
    '''
    filename = request.module.__file__
    test_dir = os.path.dirname(filename)
    data_dir = os.path.join(test_dir, 'data') 
    dir_util.copy_tree(data_dir, str(tmpdir))

    def getter(filename, as_str=True):
        filepath = tmpdir.join(filename)
        if as_str:
            return str(filepath)
        return filepath

    return getter

def test_reflectance_prospect5(datadir):
    # runs prospect and compares to online prospect run
    fname = datadir("prospect5_spectrum.txt")
    w, true_refl, true_trans = np.loadtxt(fname,
             unpack=True)

    w, refl, trans = run_prospect(2.1, 40, 10., 0.1, 
                0.015, 0.009, prospect_version="5")
    assert np.allclose( true_refl, refl, atol=1e-4)

def test_transmittance_prospect5(datadir):
    # runs prospect and compares to online prospect run
    fname = datadir("prospect5_spectrum.txt")
    w, true_refl, true_trans = np.loadtxt(fname,
            unpack=True)

    w, refl, trans = run_prospect(2.1, 40, 10., 0.1, 
                0.015, 0.009, prospect_version="5")
    assert np.allclose( true_trans, trans, atol=1e-4)


