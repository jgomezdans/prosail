#!/usr/bin/env python
from functools import lru_cache
import numpy as np

import numba


def weighted_sum_over_lidf (lidf, tts, tto, psi, cts, cto, ctscto):
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
    litab=[float(angle)*angle_step+(angle_step/2.0) for angle in range(n_angles)]
    for i,ili in enumerate(litab):
        ttl=float(ili)
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

cweighted_sum_over_lidf = numba.jit(weighted_sum_over_lidf)

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
        x1 = np.array(x2)
        y1 = np.array(y2)
        f1 = np.array(f2)
    tsstoo = np.array(f1)
    if np.isnan(sumint) : 
        sumint=0.

    return tsstoo, sumint

chotspot_calculations = numba.jit(hotspot_calculations)

def run_sail (rho, tau, lidfa, lidfb, lidftype, lai, hotspot, 
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
        lidf = verhoef_bimodal(lidfa, lidfb)
    elif lidftype == 2:
        lidf = ellipsoidal(lidfa)
    else:
        raise ValueError, \
            "lidftype can only be 1 (Campbell) or 2 (ellipsoidal)"
    #Calculate geometric factors associated with extinction and scattering 
    #Initialise sums
    ks=0.
    ko=0.
    bf=0.
    sob=0.
    sof=0.

    try:
        ks, ko, bf, sob, sof = cweighted_sum_over_lidf(lidf, tts, tto, psi)
    except:
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
    if len(sigf)>1:
        sigf[sigf == 0.0] = 1e-36
        sigb[sigb == 0.0] = 1e-36
    else:
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
        try:
            tsstoo, sumint = chotspot_calculations(alf, lai, ko, ks)
        except: 
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
    if len(dn)>1:
        dn[dn < 1e-36]=1e-36
    else:
        dn=max(1e-36,dn)
    rddt=rdd+tdd*rsoil*tdd/dn
    rsdt=rsd+(tsd+tss)*rsoil*tdd/dn
    rdot=rdo+tdd*rsoil*(tdo+too)/dn
    rsodt=((tss+tsd)*tdo+(tsd+tss*rsoil*rdd)*too)*rsoil/dn
    rsost=rso+tsstoo*rsoil
    rsot=rsost+rsodt
    
    return [tss,too,tsstoo,rdd,tdd,rsd,tsd,rdo,tdo,
          rso,rsos,rsod,rddt,rsdt,rdot,rsodt,rsost,rsot,gammasdf,gammasdb,gammaso]


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

    from math import sin, cos, acos, radians, pi 
    cts=cos(radians(tts))
    cto=cos(radians(tto))
    sts=sin(radians(tts))
    sto=sin(radians(tto))
    cospsi=cos(radians(psi))
    psir=radians(psi)
    cttl=cos(radians(ttl))
    sttl=sin(radians(ttl))
    cs=cttl*cts
    co=cttl*cto
    ss=sttl*sts
    so=sttl*sto  
    cosbts=5.
    if abs(ss) > 1e-6 : cosbts=-cs/ss
    cosbto=5.
    if abs(so) > 1e-6 : cosbto=-co/so
    if abs(cosbts) < 1.0:
        bts=acos(cosbts)
        ds=ss
    else:
        bts=pi
        ds=cs
    chi_s=2./pi*((bts-pi*0.5)*cs+sin(bts)*ss)
    if abs(cosbto) < 1.0:
        bto=acos(cosbto)
        do_=so
    else:
        if tto < 90.:
            bto=pi
            do_=co
        else:
            bto=0.0
            do_=-co
    chi_o=2.0/pi*((bto-pi*0.5)*co+sin(bto)*so)
    btran1=abs(bts-bto)
    btran2=pi-abs(bts+bto-pi)
    if psir <= btran1:
        bt1=psir
        bt2=btran1
        bt3=btran2
    else:
        bt1=btran1
        if psir <= btran2:
            bt2=psir
            bt3=btran2
        else:
            bt2=btran2
            bt3=psir
    t1=2.*cs*co+ss*so*cospsi
    t2=0.
    if bt2 > 0.: t2=sin(bt2)*(2.*ds*do_+ss*so*cos(bt1)*cos(bt3))
    denom=2.*pi**2
    frho=((pi-bt2)*t1+t2)/denom
    ftau=(-bt2*t1+t2)/denom
    if frho < 0. : frho=0.
    if ftau < 0. : ftau=0.
   
    return [chi_s,chi_o,frho,ftau]    

def Jfunc1(k,l,t) :
    ''' J1 function with avoidance of singularity problem.'''
    from numpy import exp,zeros,size
    nb=size(l)
    del_=(k-l)*t
    if nb > 1:
        result=zeros(nb)
        result[abs(del_) > 1e-3]=(exp(-l[abs(del_)> 1e-3]*t)-exp(-k*t))/(k-l[abs(del_)> 1e-3])
        result[abs(del_)<= 1e-3]=0.5*t*(exp(-k*t)+exp(-l[abs(del_)<= 1e-3]*t))*(1.-(del_[abs(del_)<= 1e-3]**2.)/12.)
    else:
        if abs(del_) > 1e-3 :
            result=(exp(-l*t)-exp(-k*t))/(k-l)
        else:
            result=0.5*t*(exp(-k*t)+exp(-l*t))*(1.-(del_**2.)/12.)
    return result

def Jfunc2(k,l,t) :
    '''J2 function.'''
    from numpy import exp
    return (1.-exp(-(k+l)*t))/(k+l)

def Jfunc1_wl(k,l,t) :
    '''J1 function with avoidance of singularity problem.'''
    from math import exp
    del_=(k-l)*t
    if abs(del_) > 1e-3 :
      result=(exp(-l*t)-exp(-k*t))/(k-l)
    else:
      result=0.5*t*(exp(-k*t)+exp(-l*t))*(1.-(del_**2.)/12.)
    return result

def Jfunc2_wl(k,l,t) :
    '''J2 function.'''
    from math import exp
    return (1.-exp(-(k+l)*t))/(k+l)


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

    import math as m
    freq=1.0
    step=90.0/n_elements
    lidf=[]
    angles=[i*step for i in reversed(range(n_elements))]
    for angle in angles:
        tl1=m.radians(angle)
        if a>1.0:
            f = 1.0-m.cos(tl1)
        else:
            eps=1e-8
            delx=1.0
            x=2.0*tl1
            p=float(x)
            while delx >= eps:
                y = a*m.sin(x)+.5*b*m.sin(2.*x)
                dx=.5*(y-x+p)
                x=x+dx
                delx=abs(dx)
            f = (2.*y+p)/m.pi
        freq=freq-f
        lidf.append(freq)
        freq=float(f)
    lidf=list(reversed(lidf))
    return  lidf
    
def ellipsoidal(alpha,n_elements=18):
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
    
    from math import cos, asin,tan, log, exp, sqrt, radians
    
    alpha=float(alpha)
    excent=exp(-1.6184e-5*alpha**3.+2.1145e-3*alpha**2.-1.2390e-1*alpha+3.2491)
    sum0 = 0.
    freq=[]
    step=90.0/n_elements
    for  i in range (n_elements):
        tl1=radians(i*step)
        tl2=radians((i+1.)*step)
        x1  = excent/(sqrt(1.+excent**2.*tan(tl1)**2.))
        x2  = excent/(sqrt(1.+excent**2.*tan(tl2)**2.))
        if excent == 1. :
            freq.append(abs(cos(tl1)-cos(tl2)))
        else :
            alph  = excent/sqrt(abs(1.-excent**2.))
            alph2 = alph**2.
            x12 = x1**2.
            x22 = x2**2.
            if excent > 1. :
                alpx1 = sqrt(alph2+x12)
                alpx2 = sqrt(alph2+x22)
                dum   = x1*alpx1+alph2*log(x1+alpx1)
                freq.append(abs(dum-(x2*alpx2+alph2*log(x2+alpx2))))
            else :
                almx1 = sqrt(alph2-x12)
                almx2 = sqrt(alph2-x22)
                dum   = x1*almx1+alph2*asin(x1/alph)
                freq.append(abs(dum-(x2*almx2+alph2*asin(x2/alph))))
    sum0 = sum(freq)
    lidf=[]
    for i in range(n_elements):
        lidf.append(float(freq[i])/sum0)
    
    return lidf
