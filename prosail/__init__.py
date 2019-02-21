#!/usr/bin/env python
__author__ = "J Gomez-Dans"
__copyright__ = "Copyright 2017, 2018 J Gomez-Dans"
__version__ = "2.0.4"
__license__ = "GPLv3"
__email__ = "j.gomez-dans@ucl.ac.uk"

from .spectral_library import get_spectra
spectral_lib = get_spectra()
from .prospect_d import run_prospect
from .sail_model import run_sail, run_prosail, run_thermal_sail
