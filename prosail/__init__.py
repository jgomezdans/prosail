import numpy as np
from prosail_fortran import run_sail, run_prosail, prospect_5b


def trans_prosail ( N, cab, car, cbrown, cw, cm, lai, lidfa, lidfb, rsoil, psoil, \
        hspot, tts, tto, psi, typelidf):
    """A version of PROSAIL that uses transformed parameters to quasi-linearise
    the   model. See http://dx.doi.org/10.1016/j.rse.2011.12.027"""
    # Define the constants
    slai = -2.0
    skab = -100.0
    skar = -100.0
    skw =  -1./50.
    skm =  -1./100.
    # Transform the parameters to real units
    xlai = slai * np.log ( lai )
    xkab = skab * np.log ( cab )
    xkar = skar * np.log ( car )
    xkw = skw * np.log ( cw )
    xdm = skm * np.log ( cm )
    # Run the PROSAIL model
    retval = run_prosail ( N, xkab, xkar, cbrown, xkw, xdm, xlai, \
            lidfa, lidfb, rsoil, psoil, hspot, tts, tto, psi, typelidf )
    return retval
