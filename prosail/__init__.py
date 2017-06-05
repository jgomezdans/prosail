#!/usr/bin/env python
from collections import namedtuple
import pkgutil
from StringIO import StringIO

import numpy as np
from prospect_d import run_prospect

Spectra = namedtuple('Spectra', 'prospect5 prospectd soil light')
Prospect5Spectra = namedtuple('Prospect5Spectra', 
                                'nr kab kcar kbrown kw km')
ProspectDSpectra = namedtuple('ProspectDSpectra', 
                                'nr kab kcar kbrown kw km kant')
SoilSpectra = namedtuple("SoilSpectra", "rsoil1 rsoil2")
LightSpectra = namedtuple("LightSpectra", "es ed")

spectral_lib = get_spectra()

def get_spectra():
    """Reads the spectral information and stores is for future use."""

    # PROSPECT-D
    prospect_d_spectraf = pkgutil.get_data('prosail', 'prospect_d_spectra.txt')
    d = np.loadtxt( StringIO(prospect_d_spectraf))
    nr = d[:, 1]
    kab = d[:, 2]
    kcar = d[:, 3]
    kbrown = d[:, 5]
    kw = d[:, 6]
    km = d[:, 7]
    kant = d[:, 4]
    prospect_d_spectra = ProspectDSpectra(nr, kab, kcar, kbrown, kw, km, kant)
    # PROSPECT 5
    prospect_5_spectraf = pkgutil.get_data('prosail', 'prospect5_spectra.txt')
    nr, kab, kcar, kbrown, kw, km =  np.loadtxt(StringIO(prospect_5_spectraf), 
                                                unpack=True)
    prospect_5_spectra = Prospect5Spectra(nr, kab, kcar, kbrown, kw, km)
    # SOIL
    soil_spectraf = pkgutil.get_data('prosail', 'soil_reflectance.txt')
    rsoil1, rsoil2 =  np.loadtxt(StringIO(soil_spectraf), 
                                                unpack=True)
    soil_spectra = SoilSpectra(rsoil1, rsoil2)    
    # LIGHT
    light_spectraf = pkgutil.get_data('prosail', 'light_spectra.txt')
    es, ed =  np.loadtxt(StringIO(light_spectraf), 
                                                unpack=True)
    light_spectra = LightSpectra(es, ed)
    spectra = Spectra(prospect_d_spectra, prospect_5_spectra, 
                      soil_spectra, light_spectra)
    return spectra
    
    
    


