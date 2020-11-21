import numpy as np
from scipy.io import loadmat
import os
import prosail
from pytest import fixture
from distutils import dir_util

from prosail.prospect_d import calctav


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

    w, refl, trans = prosail.run_prospect(2.1, 40, 10., 0.1, 
                0.015, 0.009, prospect_version="5")
    assert np.allclose( true_refl, refl, atol=1e-4)

def test_transmittance_prospect5(datadir):
    # runs prospect and compares to online prospect run
    fname = datadir("prospect5_spectrum.txt")
    w, true_refl, true_trans = np.loadtxt(fname,
            unpack=True)

    w, refl, trans = prosail.run_prospect(2.1, 40, 10., 0.1, 
                0.015, 0.009, prospect_version="5")
    assert np.allclose( true_trans, trans, atol=1e-4)
    
def test_calctav_prospectd(datadir):
    fname = datadir("tav_alpha40.mat")
    tav_mtlab = loadmat(fname)['tav'].squeeze()
    tav_py = calctav(40, prosail.spectral_lib.prospectd.nr)
    assert np.allclose(tav_mtlab, tav_py, atol=1.e-4)
    
def test_reflectance_prospectd(datadir):
    fname = datadir("prospect_d_test.mat")
    refl_mtlab = loadmat(fname)['LRT'][:,1]
    w, refl, trans = prosail.run_prospect(1.2, 30, 10., 0.0, 
                0.015, 0.009, ant=1., prospect_version="D")
    assert np.allclose(refl_mtlab, refl, atol=1.e-4)
    
def test_transmittance_prospectd(datadir):
    fname = datadir("prospect_d_test.mat")
    trans_mtlab = loadmat(fname)['LRT'][:,2]
    w, refl, trans = prosail.run_prospect(1.2, 30, 10., 0.0, 
                0.015, 0.009, ant=1., prospect_version="D")
    assert np.allclose(trans_mtlab, trans, atol=1.e-4)


def test_reflectance_prospectpro(datadir):
    fname = datadir("prospect_pro_test.txt")
    w, refl_mtlab, trans_mtlab = np.loadtxt(fname, unpack=True)
    w, refl, trans = prosail.run_prospect(n=1.2, cab=30, car=10., cbrown=0.0,
                                          cw=0.015, cm=0.009, ant=1.0, prot=0.001, cbc=0.009, prospect_version="PRO")
    assert np.allclose(refl_mtlab, refl, atol=1.e-4)


def test_transmittance_prospectpro(datadir):
    fname = datadir("prospect_pro_test.txt")
    w, refl_mtlab, trans_mtlab = np.loadtxt(fname, unpack=True)
    w, refl, trans = prosail.run_prospect(1.2, 30, 10., 0.0,
                                          0.015, 0.009, ant=1., prot=0.001, cbc=0.009,  prospect_version="PRO")
    assert np.allclose(trans_mtlab, trans, atol=1.e-4)