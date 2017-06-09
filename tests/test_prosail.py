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


def test_rsot_prosail5(datadir):
    fname = datadir("REFL_CAN.txt")
    w, resv, rdot, rsot, rddt, rsdt = np.loadtxt(fname,
             unpack=True)
    rr = prosail.run_prosail(1.5, 40., 8., 0.0, 0.01, 0.009, 3., -0.35, 0.01,
                        30., 10., 0., typelidf=1, lidfb=-0.15, 
                        rsoil = 1., psoil=1., factor="SDR")
    assert np.allclose(rsot, rr, atol=0.01)


def test_rdot_prosail5(datadir):
    fname = datadir("REFL_CAN.txt")
    w, resv, rdot, rsot, rddt, rsdt = np.loadtxt(fname,
             unpack=True)
    rr = prosail.run_prosail(1.5, 40., 8., 0.0, 0.01, 0.009, 3., -0.35, 0.01,
                        30., 10., 0., typelidf=1, lidfb=-0.15, 
                        rsoil = 1., psoil=1., factor="HDR")
    assert np.allclose(rdot, rr, atol=0.01)

def test_rddt_prosail5(datadir):
    fname = datadir("REFL_CAN.txt")
    w, resv, rdot, rsot, rddt, rsdt = np.loadtxt(fname,
             unpack=True)
    rr = prosail.run_prosail(1.5, 40., 8., 0.0, 0.01, 0.009, 3., -0.35, 0.01,
                        30., 10., 0., typelidf=1, lidfb=-0.15, 
                        rsoil = 1., psoil=1., factor="BHR")
    assert np.allclose(rddt, rr, atol=0.01)
    
def test_rsdt_prosail5(datadir):
    fname = datadir("REFL_CAN.txt")
    w, resv, rdot, rsot, rddt, rsdt = np.loadtxt(fname,
             unpack=True)
    rr = prosail.run_prosail(1.5, 40., 8., 0.0, 0.01, 0.009, 3., -0.35, 0.01,
                        30., 10., 0., typelidf=1, lidfb=-0.15, 
                        rsoil = 1., psoil=1., factor="DHR")
    assert np.allclose(rsdt, rr, atol=0.01)
