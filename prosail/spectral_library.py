#!/usr/bin/env python
"""Spectral libraries for PROSPECT + SAIL
"""
import pkgutil
from collections import namedtuple
from io import BytesIO

import numpy as np

Spectra = namedtuple('Spectra', 'prospect5 prospectd soil light')
Prospect5Spectra = namedtuple('Prospect5Spectra', 
                                'nr kab kcar kbrown kw km')
ProspectDSpectra = namedtuple('ProspectDSpectra', 
                                'nr kab kcar kbrown kw km kant')
SoilSpectra = namedtuple("SoilSpectra", "rsoil1 rsoil2")
LightSpectra = namedtuple("LightSpectra", "es ed")

def get_spectra():
    """Reads the spectral information and stores is for future use."""

    # PROSPECT-D
    prospect_d_spectraf = pkgutil.get_data('prosail', 'prospect_d_spectra.txt')
    _, nr, kab, kcar, kant, kbrown, kw, km= np.loadtxt(
        BytesIO(prospect_d_spectraf), unpack=True)
    prospect_d_spectra = ProspectDSpectra(nr, kab, kcar, kbrown, kw, km, kant)
    # PROSPECT 5
    prospect_5_spectraf = pkgutil.get_data('prosail', 'prospect5_spectra.txt')
    nr, kab, kcar, kbrown, kw, km =  np.loadtxt(BytesIO(prospect_5_spectraf),
                                                unpack=True)
    prospect_5_spectra = Prospect5Spectra(nr, kab, kcar, kbrown, kw, km)
    # SOIL
    soil_spectraf = pkgutil.get_data('prosail', 'soil_reflectance.txt')
    rsoil1, rsoil2 =  np.loadtxt(BytesIO(soil_spectraf),
                                                unpack=True)
    soil_spectra = SoilSpectra(rsoil1, rsoil2)    
    # LIGHT
    light_spectraf = pkgutil.get_data('prosail', 'light_spectra.txt')
    es, ed =  np.loadtxt(BytesIO(light_spectraf),
                                                unpack=True)
    light_spectra = LightSpectra(es, ed)
    spectra = Spectra(prospect_5_spectra, prospect_d_spectra, 
                      soil_spectra, light_spectra)
    return spectra
