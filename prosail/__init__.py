#!/usr/bin/env python
from collections import namedtuple

import numpy as np

from spectral_libary import get_spectra
from prospect_d import run_prospect

Spectra = namedtuple('Spectra', 'prospect5 prospectd soil light')
Prospect5Spectra = namedtuple('Prospect5Spectra', 
                                'nr kab kcar kbrown kw km')
ProspectDSpectra = namedtuple('ProspectDSpectra', 
                                'nr kab kcar kbrown kw km kant')
SoilSpectra = namedtuple("SoilSpectra", "rsoil1 rsoil2")
LightSpectra = namedtuple("LightSpectra", "es ed")


    
spectral_lib = get_spectra()
    
    


