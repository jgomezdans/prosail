import numpy as np
from prosail_fortran import run_sail as sail, run_prosail as prosail
from prosail_fortran import prospect_5b
from prosail_fortran import mod_dataspec_p5b as spectral_libs

def run_prosail (n,cab,car,cbrown,cw,cm,lai,lidfa,lidfb,rsoil,psoil,hspot,
                 tts,tto,psi,typelidf, 
                 soil_spectrum1=None,soil_spectrum2=None ):
    """Run the PROSPECT_5B and SAILh radiative transfer models. The soil
    model is a linear mixture model, where two spectra are combined together as
    
         rho_soil = rsoil*(psoil*soil_spectrum1+(1-psoil)*soil_spectrum2)
    By default, ``soil_spectrum1`` is a dry soil, and ``soil_spectrum2`` is a
    wet soil, so in that case, ``psoil`` is a surface soil moisture parameter.
    ``rsoil`` is a  soil brightness term. You can provide one or the two
    soil spectra if you want.  The soil spectra must be defined
    between 400 and 2500 nm with 1nm spacing.
    
    Parameters
    ----------
    n: float
        Leaf layers
    cab: float
        leaf chlorophyll concentration
    car: float
        leaf carotenoid concentration
    cbrown: float
        senescent pigment
    cw: float
        equivalent leaf water
    cm: float
        leaf dry matter
    lai: float
        leaf area index
    lidfa: float
        a parameter for leaf angle distribution. If ``typliedf``=2, average
        leaf inclination angle.
    lidfb: float
        b parameter for leaf angle distribution. If ``typelidf``=2, ignored
    rsoil: float
        Soil scalar
    psoil: float
        Soil scalar
    tts: float
        Solar zenith angle
    tto: float
        Sensor zenith angle
    psi: float
        Relative sensor-solar azimuth angle ( saa - vaa )
    soil_spectrum1: 2101-element array
        First component of the soil spectrum
    soil_spectrum2: 2101-element array
        Second component of the soil spectrum
    
    Returns
    --------
    Directional surface reflectance between 400 and 2500 nm

        
    """

    if soil_spectrum1 is not None:
        assert ( len(soil_spectrum1) == 2101 )
    else:
        soil_spectrum1 = spectral_libs.rsoil1

    if soil_spectrum2 is not None:
        assert ( len(soil_spectrum1) == 2101 )
    else:
        soil_spectrum2 = spectral_libs.rsoil2

    rho = prosail (n,cab,car,cbrown,cw,cm,lai,lidfa,lidfb,rsoil,psoil,hspot,
                 tts,tto,psi,typelidf, soil_spectrum1, soil_spectrum2 )
    return rho

def run_sail (refl,trans,lai,lidfa,lidfb,rsoil,psoil,hspot,tts,tto,psi,typelidf,
                 soil_spectrum1=None,soil_spectrum2=None ):
    """Run the SAILh radiative transfer model. The soil model is a linear 
    mixture model, where two spectra are combined together as
    
         rho_soil = rsoil*(psoil*soil_spectrum1+(1-psoil)*soil_spectrum2)
         
    By default, ``soil_spectrum1`` is a dry soil, and ``soil_spectrum2`` is a
    wet soil, so in that case, ``psoil`` is a surface soil moisture parameter.
    ``rsoil`` is a  soil brightness term. You can provide one or the two
    soil spectra if you want. The soil spectra, and leaf spectra must be defined
    between 400 and 2500 nm with 1nm spacing.
    
    Parameters
    ----------
    refl: 2101-element array
        Leaf reflectance
    trans: 2101-element array
        leaf transmittance
    lai: float
        leaf area index
    lidfa: float
        a parameter for leaf angle distribution. If ``typliedf``=2, average
        leaf inclination angle.
    lidfb: float
        b parameter for leaf angle distribution. If ``typelidf``=2, ignored
    rsoil: float
        Soil scalar
    psoil: float
        Soil scalar
    tts: float
        Solar zenith angle
    tto: float
        Sensor zenith angle
    psi: float
        Relative sensor-solar azimuth angle ( saa - vaa )
    soil_spectrum1: 2101-element array
        First component of the soil spectrum
    soil_spectrum2: 2101-element array
        Second component of the soil spectrum
    
    Returns
    --------
    Directional surface reflectance between 400 and 2500 nm

        
    """


    if soil_spectrum1 is not None:
        assert ( len(soil_spectrum1) == 2101 )
    else:
        soil_spectrum1 = spectral_libs.rsoil1

    if soil_spectrum2 is not None:
        assert ( len(soil_spectrum1) == 2101 )
    else:
        soil_spectrum2 = spectral_libs.rsoil2

    rho = sail (refl,trans,lai,lidfa,lidfb,rsoil,psoil,hspot,tts,tto,psi,typelidf, 
                soil_spectrum1, soil_spectrum2 )
    return rho



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
