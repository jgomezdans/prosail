import numpy as np

def hapke ( omega, back, c, b0, tts, tto, psi ):
    """The Hapke soil model"""
    
    # Define angles, lots of them
    # Remember to convert to radians
    cts = np.cos (np.deg2rad(tts))
    cto = np.cos (np.deg2rad(tto))
    sts = np.sin (np.deg2rad(tts))
    sto = np.sin (np.deg2rad(tto))
    cps = np.cos (np.deg2rad(psi))
    csd = cts*cto + sts*sto*cps
    if csd > 1.0:
        csd = 1.
    adel = np.acos (csd)
    # The phase function
    p = 1.0 + back*csd + 0.5*c*(3.*csd*csd-1)
    
    # Extinction, scattering
    ks = 1./cts
    ko = 1./cto
    sigb = omega*(1. + back/4.)
    sigf = omega*(1. - back/4.)
    att = 2.0 - sigf
    m = np.sqrt (att**2 + sigb**2)
    rinf = 0.
    if omega > 0:
        rinf = (att-m)/sigb
    
    sb = omega*(ks/2. + back/4.)
    sf = omega*(ks/2. - back/4.)
    vb = omega*(ko/2. + back/4.)
    vf = omega*(ko/2. - back/4.)
    
    rddsoil = rinf
    rsdsoil = (sb + sf*rinf)/(ks+m)
    rdosoil = (vb + vf*rinf)/(ko+m)
    
    w = omega*p/(4.*cts*cto)
    
    cor = b0/(1. + np.tan ( 0.5*adel)/h)
    rsosoils = w*(1+cor)/(ks+ko)
    
    rsosoild = ((rdosoil*(sf+sb*rinf)+rsdsoil*(vf+vb*rinf))/(ks+ko) 
    rsosoild = rsosoild - rdosoil*rsdsoil*rinf)/(1.-rinf**2)
    
    rsosoil = rsosoils + rsosoild
    # Returns:
    # 1. bi-hemispherical reflectance
    # 2. HDR
    # 3. DHR
    # 4. BRF
    
    return rddsoil, rsdsoil, rdosoil, rsosoil
                
