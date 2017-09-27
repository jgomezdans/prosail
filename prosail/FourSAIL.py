#!/usr/bin/env python
from math import  exp, radians


try:
    from functools import lru_cache
except ImportError:
    from backports.functools_lru_cache import lru_cache

import numpy as np
import numba


@numba.jit('Tuple((f8, f8, f8, f8))(f8,f8,f8,f8)',
           nopython=True, cache=True)
def volscatt(tts,tto,psi,ttl) :
    '''Compute volume scattering functions and interception coefficients
    for given solar zenith, viewing zenith, azimuth and leaf inclination angle.

    Parameters
    ----------
    tts : float
        Solar Zenith Angle (degrees).
    tto : float
        View Zenight Angle (degrees).
    psi : float
        View-Sun reliative azimuth angle (degrees).
    ttl : float
        leaf inclination angle (degrees).

    Returns
    -------
    chi_s : float
        Interception function  in the solar path.
    chi_o : float
        Interception function  in the view path.
    frho : float
        Function to be multiplied by leaf reflectance to obtain the volume scattering.
    ftau : float
        Function to be multiplied by leaf transmittance to obtain the volume scattering.

    References
    ----------
    Wout Verhoef, april 2001, for CROMA.
    '''


    cts = np.cos(np.radians(tts))
    cto = np.cos(np.radians(tto))
    sts = np.sin(np.radians(tts))
    sto = np.sin(np.radians(tto))
    cospsi = np.cos(np.radians(psi))
    psir = np.radians(psi)
    cttl = np.cos(np.radians(ttl))
    sttl = np.sin(np.radians(ttl))
    cs = cttl*cts
    co = cttl*cto
    ss = sttl*sts
    so = sttl*sto
    cosbts = 5.
    if np.abs(ss) > 1e-6 :
        cosbts = -cs/ss
    cosbto = 5.
    if np.abs(so) > 1e-6 :
        cosbto = -co/so
    if np.abs(cosbts) < 1.0:
        bts = np.arccos(cosbts)
        ds = ss
    else:
        bts = np.pi
        ds = cs
    chi_s = 2./np.pi*((bts-np.pi*0.5)*cs+np.sin(bts)*ss)
    if abs(cosbto) < 1.0:
        bto = np.arccos(cosbto)
        do_ = so
    else:
        if tto < 90.:
            bto = np.pi
            do_ = co
        else:
            bto = 0.0
            do_ = -co
    chi_o = 2.0/np.pi*((bto-np.pi*0.5)*co+np.sin(bto)*so)
    btran1 = np.abs(bts-bto)
    btran2 = np.pi-np.abs(bts+bto-np.pi)
    if psir <= btran1:
        bt1 = psir
        bt2 = btran1
        bt3 = btran2
    else:
        bt1 = btran1
        if psir <= btran2:
            bt2 = psir
            bt3 = btran2
        else:
            bt2 = btran2
            bt3 = psir
    t1 = 2.*cs*co+ss*so*cospsi
    t2 = 0.
    if bt2 > 0.:
        t2 = np.sin(bt2)*(2.*ds*do_+ss*so*np.cos(bt1)*np.cos(bt3))
    denom=2.*np.pi**2
    frho = ((np.pi-bt2)*t1+t2)/denom
    ftau = (-bt2*t1+t2)/denom
    if frho < 0. :
        frho=0.
    if ftau < 0. :
        ftau=0.

    return (chi_s,chi_o,frho,ftau)



@numba.jit('Tuple((f8, f8, f8, f8, f8))(f8[:], f8, f8, f8)',
           nopython=True, cache=True)
def weighted_sum_over_lidf (lidf, tts, tto, psi):
    ks = 0.
    ko = 0.
    bf = 0.
    sob = 0.
    sof = 0.
    cts   = np.cos(np.radians(tts))
    cto   = np.cos(np.radians(tto))
    ctscto  = cts*cto

    n_angles=len(lidf)
    angle_step=float(90.0/n_angles)
    litab = np.arange(n_angles)*angle_step + (angle_step*0.5)

    for i,ili in enumerate(litab):
        ttl = 1.*ili
        cttl=np.cos(np.radians(ttl))
        # SAIL volume scattering phase function gives interception and portions to be multiplied by rho and tau
        [chi_s,chi_o,frho,ftau]=volscatt(tts,tto,psi,ttl)
        # Extinction coefficients
        ksli=chi_s/cts
        koli=chi_o/cto
        # Area scattering coefficient fractions
        sobli=frho*np.pi/ctscto
        sofli=ftau*np.pi/ctscto
        bfli=cttl**2.
        ks+=ksli*float(lidf[i])
        ko+=koli*float(lidf[i])
        bf+=bfli*float(lidf[i])
        sob+=sobli*float(lidf[i])
        sof+=sofli*float(lidf[i])
    return ks, ko, bf, sob, sof



@lru_cache(maxsize=16)
def define_geometric_constants ( tts, tto, psi):
    cts   = np.cos(np.radians(tts))
    cto   = np.cos(np.radians(tto))
    ctscto  = cts*cto
    tants = np.tan(np.radians(tts))
    tanto = np.tan(np.radians(tto))
    cospsi  = np.cos(np.radians(psi))
    dso = np.sqrt(tants**2.+tanto**2.-2.*tants*tanto*cospsi)
    return cts, cto, ctscto, tants, tanto, cospsi, dso

@numba.jit("Tuple((f8,f8))(f8,f8,f8,f8)",
           nopython=True, cache=True)
def hotspot_calculations(alf, lai, ko, ks ):
    fhot=lai*np.sqrt(ko*ks)
    # Integrate by exponential Simpson method in 20 steps the steps are arranged according to equal partitioning of the slope of the joint probability function
    x1=0.
    y1=0.
    f1=1.
    fint=(1.-np.exp(-alf))*.05
    sumint=0.
    for istep in range(1,21):
        if istep < 20 :
            x2 = -np.log(1.-istep*fint)/alf
        else:
            x2=1.
        y2 = -(ko+ks)*lai*x2+fhot*(1.-np.exp(-alf*x2))/alf
        f2 = np.exp(y2)
        sumint = sumint+(f2-f1)*(x2-x1)/(y2-y1)
        x1 = x2
        y1 = y2
        f1 = f2
    tsstoo = f1
    if np.isnan(sumint) :
        sumint=0.
    return tsstoo, sumint

def Jfunc1(k,l,t) :
    ''' J1 function with avoidance of singularity problem.'''
    try:
        nb=len(l)
    except TypeError:
        nb = 1
    del_=(k-l)*t
    if nb > 1:
        result=np.zeros(nb)
        result[np.abs(del_) > 1e-3] = (np.exp(-l[np.abs(del_)> 1e-3]*t)-
                                  np.exp(-k*t))/(k-l[np.abs(del_)> 1e-3])
        result[np.abs(del_) <= 1e-3] = 0.5*t*(np.exp(-k*t) +
                                    np.exp(-l[np.abs(del_)<= 1e-3]*t))* \
                                    (1.-(del_[np.abs(del_)<= 1e-3]**2.)/12.)
    else:
        if np.abs(del_) > 1e-3 :
            result = (np.exp(-l*t)-np.exp(-k*t))/(k-l)
        else:
            result = 0.5*t*(np.exp(-k*t)+np.exp(-l*t))*(1.-(del_**2.)/12.)
    return result

def Jfunc2(k,l,t) :
    '''J2 function.'''
    return (1.-np.exp(-(k+l)*t))/(k+l)

@numba.jit('f8[:](f8,f8,i8)',
           nopython=True, cache=True)
def verhoef_bimodal(a,b,n_elements=18):
    '''Calculate the Leaf Inclination Distribution Function based on the
    Verhoef's bimodal LIDF distribution.
    Parameters
    ----------
    a : float
        controls the average leaf slope.
    b : float
        controls the distribution's bimodality.

            * LIDF type     [a,b].
            * Planophile    [1,0].
            * Erectophile   [-1,0].
            * Plagiophile   [0,-1].
            * Extremophile  [0,1].
            * Spherical     [-0.35,-0.15].
            * Uniform       [0,0].
            * requirement: |LIDFa| + |LIDFb| < 1.
    n_elements : int
        Total number of equally spaced inclination angles.

    Returns
    -------
    lidf : list
        Leaf Inclination Distribution Function at equally spaced angles.

    References
    ----------
    .. [Verhoef1998] Verhoef, Wout. Theory of radiative transfer models applied
        in optical remote sensing of vegetation canopies.
        Nationaal Lucht en Ruimtevaartlaboratorium, 1998.
        http://library.wur.nl/WebQuery/clc/945481.
        '''

    freq=1.0
    step=90.0/n_elements
    lidf=np.zeros(n_elements)*1.
    angles = (np.arange(n_elements)*step)[::-1]
    i = 0
    for angle in angles:
        tl1=np.radians(angle)
        if a>1.0:
            f = 1.0-np.cos(tl1)
        else:
            eps=1e-8
            delx=1.0
            x=2.0*tl1
            p=float(x)
            while delx >= eps:
                y = a*np.sin(x)+.5*b*np.sin(2.*x)
                dx=.5*(y-x+p)
                x=x+dx
                delx=abs(dx)
            f = (2.*y+p)/np.pi
        freq=freq-f
        lidf[i] = freq
        freq=float(f)
        i += 1
    lidf=lidf[::-1]
    return  lidf



@numba.jit('f8[:](f8,i8)',
           nopython=True, cache=True)
def campbell(alpha,n_elements=18):
    '''Calculate the Leaf Inclination Distribution Function based on the
    mean angle of [Campbell1990] ellipsoidal LIDF distribution.
    Parameters
    ----------
    alpha : float
        Mean leaf angle (degrees) use 57 for a spherical LIDF.
    n_elements : int
        Total number of equally spaced inclination angles .

    Returns
    -------
    lidf : list
        Leaf Inclination Distribution Function for 18 equally spaced angles.

    References
    ----------
    .. [Campbell1986] G.S. Campbell, Extinction coefficients for radiation in
        plant canopies calculated using an ellipsoidal inclination angle distribution,
        Agricultural and Forest Meteorology, Volume 36, Issue 4, 1986, Pages 317-321,
        ISSN 0168-1923, http://dx.doi.org/10.1016/0168-1923(86)90010-9.
    .. [Campbell1990] G.S Campbell, Derivation of an angle density function for
        canopies with ellipsoidal leaf angle distributions,
        Agricultural and Forest Meteorology, Volume 49, Issue 3, 1990, Pages 173-176,
        ISSN 0168-1923, http://dx.doi.org/10.1016/0168-1923(90)90030-A.
    '''


    alpha=float(alpha)
    excent=exp(-1.6184e-5*alpha**3.+2.1145e-3*alpha**2.-1.2390e-1*alpha+3.2491)
    sum0 = 0.
    freq = np.zeros(n_elements)
    step=90.0/n_elements
    for  i in range (n_elements):
        tl1=radians(i*step)
        tl2=radians((i+1.)*step)
        x1  = excent/(np.sqrt(1.+excent**2.*np.tan(tl1)**2.))
        x2  = excent/(np.sqrt(1.+excent**2.*np.tan(tl2)**2.))
        if excent == 1. :
            freq[i] = abs(np.cos(tl1)-np.cos(tl2))
        else :
            alph  = excent/np.sqrt(abs(1.-excent**2.))
            alph2 = alph**2.
            x12 = x1**2.
            x22 = x2**2.
            if excent > 1. :
                alpx1 = np.sqrt(alph2+x12)
                alpx2 = np.sqrt(alph2+x22)
                dum   = x1*alpx1+alph2*np.log(x1+alpx1)
                freq[i]=abs(dum-(x2*alpx2+alph2*np.log(x2+alpx2)))
            else :
                almx1 = np.sqrt(alph2-x12)
                almx2 = np.sqrt(alph2-x22)
                dum   = x1*almx1+alph2*np.arcsin(x1/alph)
                freq[i]=abs(dum-(x2*almx2+alph2*np.arcsin(x2/alph)))
    sum0 = np.sum(freq)
    lidf=np.zeros(n_elements)
    for i in range(n_elements):
        lidf[i] = freq[i]/sum0

    return lidf


def foursail (rho, tau, lidfa, lidfb, lidftype, lai, hotspot,
    tts, tto, psi, rsoil):
    """
    Parameters
    ----------
    rho : array_like
        leaf lambertian reflectance.
    tau : array_like
        leaf transmittance.
    lidfa : float
        Leaf Inclination Distribution at regular angle steps.
    lidfb : float
        Leaf Inclination Distribution at regular angle steps.
    lidftype : float
        Leaf Inclination Distribution at regular angle steps.
    lai : float
        Leaf Area Index.
    hotspot : float
        Hotspot parameter.
    tts : float
        Sun Zenith Angle (degrees).
    tto : float
        View(sensor) Zenith Angle (degrees).
    psi : float
        Relative Sensor-Sun Azimuth Angle (degrees).
    rsoil : array_like
        soil lambertian reflectance.

    Returns
    -------
    tss : array_like
        beam transmittance in the sun-target path.
    too : array_like
        beam transmittance in the target-view path.
    tsstoo : array_like
        beam tranmittance in the sur-target-view path.
    rdd : array_like
        canopy bihemisperical reflectance factor.
    tdd : array_like
        canopy bihemishperical transmittance factor.
    rsd : array_like
        canopy directional-hemispherical reflectance factor.
    tsd : array_like
        canopy directional-hemispherical transmittance factor.
    rdo : array_like
        canopy hemispherical-directional reflectance factor.
    tdo : array_like
        canopy hemispherical-directional transmittance factor.
    rso : array_like
        canopy bidirectional reflectance factor.
    rsos : array_like
        single scattering contribution to rso.
    rsod : array_like
        multiple scattering contribution to rso.
    rddt : array_like
        surface bihemispherical reflectance factor.
    rsdt : array_like
        surface directional-hemispherical reflectance factor.
    rdot : array_like
        surface hemispherical-directional reflectance factor.
    rsodt : array_like
        reflectance factor.
    rsost : array_like
        reflectance factor.
    rsot : array_like
        surface bidirectional reflectance factor.
    gammasdf : array_like
        'Thermal gamma factor'.
    gammasdb : array_like
        'Thermal gamma factor'.
    gammaso : array_like
        'Thermal gamma factor'.

    References
    ----------
    .. [Verhoef2007] Verhoef, W.; Jia, Li; Qing Xiao; Su, Z., (2007) Unified Optical-Thermal
        Four-Stream Radiative Transfer Theory for Homogeneous Vegetation Canopies,
        IEEE Transactions on Geoscience and Remote Sensing, vol.45, no.6, pp.1808-1822,
        http://dx.doi.org/10.1109/TGRS.2007.895844 based on  in Verhoef et al. (2007).
    """


    # Define some geometric constants.
    cts, cto, ctscto, tants, tanto, cospsi, dso = \
        define_geometric_constants(tts, tto, psi)

    # Calcualte leaf angle distribution
    if lidftype == 1:
        lidf = verhoef_bimodal(lidfa, lidfb, n_elements=18)
    elif lidftype == 2:
        lidf = campbell(lidfa, n_elements=18)
    else:
        raise ValueError(
            "lidftype can only be 1 (Campbell) or 2 (ellipsoidal)"
        )
    #Calculate geometric factors associated with extinction and scattering
    #Initialise sums
    ks=0.
    ko=0.
    bf=0.
    sob=0.
    sof=0.

    ks, ko, bf, sob, sof = weighted_sum_over_lidf(lidf, tts, tto, psi)


    # Geometric factors to be used later with rho and tau
    sdb = 0.5*(ks + bf)
    sdf = 0.5*(ks - bf)
    dob = 0.5*(ko + bf)
    dof = 0.5*(ko - bf)
    ddb = 0.5*(1.0 + bf)
    ddf = 0.5*(1.0 - bf)


    sigb=ddb*rho+ddf*tau
    sigf=ddf*rho+ddb*tau
    try:
        sigf[sigf == 0.0] = 1.e-36
        sigb[sigb == 0.0] = 1.0e-36
    except TypeError:
        sigf = max(1e-36, sigf)
        sigb = max(1e-36, sigb)
    att = 1.-sigf
    m = np.sqrt(att**2.- sigb**2.)
    sb = sdb*rho + sdf*tau
    sf = sdf*rho+sdb*tau
    vb = dob*rho+dof*tau
    vf = dof*rho+dob*tau
    w = sob*rho+sof*tau

    if lai<=0:
        # No canopy...
        tss = 1
        too= 1
        tsstoo= 1
        rdd= 0
        tdd=1
        rsd=0
        tsd=0
        rdo=0
        tdo=0
        rso=0
        rsos=0
        rsod=0
        rddt= rsoil
        rsdt= rsoil
        rdot= rsoil
        rsodt= 0
        rsost= rsoil
        rsot= rsoil
        gammasdf=0
        gammaso=0
        gammasdb=0

        return [tss,too,tsstoo,rdd,tdd,rsd,tsd,rdo,tdo,
            rso,rsos,rsod,rddt,rsdt,rdot,rsodt,rsost,rsot,gammasdf,gammasdb,gammaso]

    e1 = np.exp(-m*lai)
    e2 = e1**2.
    rinf = (att-m)/sigb
    rinf2 = rinf**2.
    re = rinf*e1
    denom = 1.-rinf2*e2
    J1ks = Jfunc1(ks,m,lai)
    J2ks = Jfunc2(ks,m,lai)
    J1ko = Jfunc1(ko,m,lai)
    J2ko = Jfunc2(ko,m,lai)
    Pss = (sf+sb*rinf)*J1ks
    Qss = (sf*rinf+sb)*J2ks
    Pv = (vf+vb*rinf)*J1ko
    Qv = (vf*rinf+vb)*J2ko
    tdd = (1.-rinf2)*e1/denom
    rdd = rinf*(1.-e2)/denom
    tsd = (Pss-re*Qss)/denom
    rsd = (Qss-re*Pss)/denom
    tdo = (Pv-re*Qv)/denom
    rdo = (Qv-re*Pv)/denom
    # Thermal "sd" quantities
    gammasdf = (1.+rinf)*(J1ks-re*J2ks)/denom
    gammasdb = (1.+rinf)*(-re*J1ks+J2ks)/denom
    tss = np.exp(-ks*lai)
    too = np.exp(-ko*lai)
    z = Jfunc2(ks,ko,lai)
    g1 = (z-J1ks*too)/(ko+m)
    g2 = (z-J1ko*tss)/(ks+m)
    Tv1 = (vf*rinf+vb)*g1
    Tv2 = (vf+vb*rinf)*g2
    T1 = Tv1*(sf+sb*rinf)
    T2 = Tv2*(sf*rinf+sb)
    T3 = (rdo*Qss+tdo*Pss)*rinf
    # Multiple scattering contribution to bidirectional canopy reflectance
    rsod = (T1+T2-T3)/(1.-rinf2)
    # Thermal "sod" quantity
    T4 = Tv1*(1.+rinf)
    T5 = Tv2*(1.+rinf)
    T6 = (rdo*J2ks+tdo*J1ks)*(1.+rinf)*rinf
    gammasod = (T4+T5-T6)/(1.-rinf2)
    #Treatment of the hotspot-effect
    alf=1e36
    # Apply correction 2/(K+k) suggested by F.-M. Breon
    if hotspot > 0. :
        alf=(dso/hotspot)*2./(ks+ko)
    if alf == 0. :
        # The pure hotspot
        tsstoo = tss
        sumint=(1.-tss)/(ks*lai)
    else :
        # Outside the hotspot
        tsstoo, sumint = hotspot_calculations(alf, lai, ko, ks)

    # Bidirectional reflectance
    # Single scattering contribution
    rsos=w*lai*sumint
    gammasos=ko*lai*sumint
    # Total canopy contribution
    rso=rsos+rsod
    gammaso=gammasos+gammasod
    #Interaction with the soil
    dn=1.-rsoil*rdd
    try:
        dn[ dn < 1e-36] = 1e-36
    except TypeError:
        dn = max ( 1e-36, dn)
    rddt=rdd+tdd*rsoil*tdd/dn
    rsdt=rsd+(tsd+tss)*rsoil*tdd/dn
    rdot=rdo+tdd*rsoil*(tdo+too)/dn
    rsodt=((tss+tsd)*tdo+(tsd+tss*rsoil*rdd)*too)*rsoil/dn
    rsost=rso+tsstoo*rsoil
    rsot=rsost+rsodt

    return [tss,too,tsstoo,rdd,tdd,rsd,tsd,rdo,tdo,
          rso,rsos,rsod,rddt,rsdt,rdot,rsodt,rsost,rsot,gammasdf,gammasdb,gammaso]

