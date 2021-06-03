#!/usr/bin/env python
"""The PROSPECT leaf optical properties model
Versions 5, D and PRO

Thanks for @jajberni for ProspectPRO implementation!

"""
import numpy as np
from scipy.special import expi

from prosail import spectral_lib


def run_prospect(
    n,
    cab,
    car,
    cbrown,
    cw,
    cm,
    ant=0.0,
    prot=0.0,
    cbc=0.0,
    prospect_version="D",
    nr=None,
    kab=None,
    kcar=None,
    kbrown=None,
    kw=None,
    km=None,
    kant=None,
    kprot=None,
    kcbc=None,
    alpha=40.0,
):
    """The PROSPECT model, versions 5, D and PRO.
    This function runs PROSPECT. You can select the version using the
    `prospect_version` argument, and you can also set some of the spectra
    used for model calculations.

    Parameters
    -----------
    n: float
        The number of leaf layers. Unitless [-].
    cab: float
        The chlorophyll a+b concentration. [ug cm^{-2}].
    car: float
        Carotenoid concentration.  [ug cm^{-2}].
    cbrown: float
        The brown/senescent pigment. Unitless [-], often between 0 and 1
        but the literature on it is wide ranging!
    cw: float
        Equivalent leaf water. [cm]
    cm: float
        Dry matter [g cm^{-2}]
    ant: float, optional
        Anthocyanins content. Used in Prospect-D and Prospect-PRO [ug cm^{-2}]
    prot: float, optional
        Protein content. Used in Prospect-PRO. [g cm^{-2}]
    cbc: float, optional
        Carbon based constituents. Used in Prospect-PRO. [ug cm^{-2}]
    prospect_version: string, optiona, default "D".
        The version of PROSPECT, "5", "D" or "PRO".
    nr: array, optional
        The refractive index of the leaf. If `None` (default), will use the
        values for the selected PROPSECT version. [-].
    kab: 2101-element array, optional
        The specific absorption coefficient of chlorophyll (a+b) [cm^2 ug^{-1}].
    kcar: 2101-element array, optional
        The specific absorption coefficient of carotenoids [cm^2 ug^{-1}].
    kbrown:  2101-element array, optional
        The specific absorption coefficient of brown pigments (arbitrary units).
    kw:  2101-element array, optional
        The specific absorption coefficient of water (cm^{-1}).
    km: 2101-element array, optional
        The specific absorption coefficient of dry matter [cm^2 g^{-1}].
    kant: 2101-element array, optional
        The specific absorption coefficient of Anthocyanins [cm^2 nmol^{-1}].
    kprot: 2101-element array, optional
        The specific absorption coefficient of proteins [cm^2 g^{-1}].
    kcbc: 2101-element array, optional
        The specific absorption coefficient of carbon based constituents [cm^2 ug^{-1}].
    alpha: float, optional, default 40..
        Maximum incident angle relative to the normal of the leaf plane. [deg]


    Returns
    -------

    3 arrays of the size 2101: the wavelengths in [nm], the leaf reflectance
    and transmittance.

    """

    if prospect_version == "5":
        # Call the original PROSPECT-5. In case the user has supplied
        # spectra, use them.
        wv, refl, trans = prospect_d(
            n,
            cab,
            car,
            cbrown,
            cw,
            cm,
            0.0,
            0.0,
            0.0,
            spectral_lib.prospect5.nr if nr is None else nr,
            spectral_lib.prospect5.kab if kab is None else kab,
            spectral_lib.prospect5.kcar if kcar is None else kcar,
            spectral_lib.prospect5.kbrown if kbrown is None else kbrown,
            spectral_lib.prospect5.kw if kw is None else kw,
            spectral_lib.prospect5.km if km is None else km,
            np.zeros_like(spectral_lib.prospect5.km),
            np.zeros_like(spectral_lib.prospect5.km),
            np.zeros_like(spectral_lib.prospect5.km),
            alpha=alpha,
        )
    elif prospect_version.upper() == "D":
        wv, refl, trans = prospect_d(
            n,
            cab,
            car,
            cbrown,
            cw,
            cm,
            ant,
            0.0,
            0.0,
            spectral_lib.prospectd.nr if nr is None else nr,
            spectral_lib.prospectd.kab if kab is None else kab,
            spectral_lib.prospectd.kcar if kcar is None else kcar,
            spectral_lib.prospectd.kbrown if kbrown is None else kbrown,
            spectral_lib.prospectd.kw if kw is None else kw,
            spectral_lib.prospectd.km if km is None else km,
            spectral_lib.prospectd.kant if kant is None else kant,
            np.zeros_like(spectral_lib.prospect5.km),
            np.zeros_like(spectral_lib.prospect5.km),
            alpha=alpha,
        )
    elif prospect_version.upper() == "PRO":
        wv, refl, trans = prospect_d(
            n,
            cab,
            car,
            cbrown,
            cw,
            cm,
            ant,
            prot,
            cbc,
            spectral_lib.prospectpro.nr if nr is None else nr,
            spectral_lib.prospectpro.kab if kab is None else kab,
            spectral_lib.prospectpro.kcar if kcar is None else kcar,
            spectral_lib.prospectpro.kbrown if kbrown is None else kbrown,
            spectral_lib.prospectpro.kw if kw is None else kw,
            spectral_lib.prospectpro.km if km is None else km,
            spectral_lib.prospectpro.kant if kant is None else kant,
            spectral_lib.prospectpro.kprot if kprot is None else kprot,
            spectral_lib.prospectpro.kcbc if kcbc is None else kcbc,
            alpha=alpha,
        )
    else:
        raise ValueError("prospect_version can only be 5 or D!")

    return wv, refl, trans


def calctav(alpha, nr):
    # ***********************************************************************
    # calctav
    # ***********************************************************************
    # Stern F. (1964), Transmission of isotropic radiation across an
    # interface between two dielectrics, Appl. Opt., 3(1):111-113.
    # Allen W.A. (1973), Transmission of isotropic light across a
    # dielectric surface in two and three dimensions, J. Opt. Soc. Am.,
    # 63(6):664-666.
    # ***********************************************************************

    # rd  = pi/180 np.deg2rad
    n2 = nr * nr
    npx = n2 + 1
    nm = n2 - 1
    a = (nr + 1) * (nr + 1) / 2.0
    k = -(n2 - 1) * (n2 - 1) / 4.0
    sa = np.sin(np.deg2rad(alpha))

    if alpha != 90:
        b1 = np.sqrt((sa * sa - npx / 2) * (sa * sa - npx / 2) + k)
    else:
        b1 = 0.0
    b2 = sa * sa - npx / 2
    b = b1 - b2
    b3 = b ** 3
    a3 = a ** 3
    ts = (k ** 2 / (6 * b3) + k / b - b / 2) - (
        k ** 2.0 / (6 * a3) + k / a - a / 2
    )

    tp1 = -2 * n2 * (b - a) / (npx ** 2)
    tp2 = -2 * n2 * npx * np.log(b / a) / (nm ** 2)
    tp3 = n2 * (1 / b - 1 / a) / 2
    tp4 = (
        16
        * n2 ** 2
        * (n2 ** 2 + 1)
        * np.log((2 * npx * b - nm ** 2) / (2 * npx * a - nm ** 2))
        / (npx ** 3 * nm ** 2)
    )
    tp5 = (
        16
        * n2 ** 3
        * (1.0 / (2 * npx * b - nm ** 2) - 1 / (2 * npx * a - nm ** 2))
        / (npx ** 3)
    )
    tp = tp1 + tp2 + tp3 + tp4 + tp5
    tav = (ts + tp) / (2 * sa ** 2)

    return tav


def refl_trans_one_layer(alpha, nr, tau):
    # ***********************************************************************
    # reflectance and transmittance of one layer
    # ***********************************************************************
    # Allen W.A., Gausman H.W., Richardson A.J., Thomas J.R. (1969),
    # Interaction of isotropic ligth with a compact plant leaf, J. Opt.
    # Soc. Am., 59(10):1376-1379.
    # ***********************************************************************
    # reflectivity and transmissivity at the interface
    # -------------------------------------------------
    talf = calctav(alpha, nr)
    ralf = 1.0 - talf
    t12 = calctav(90, nr)
    r12 = 1.0 - t12
    t21 = t12 / (nr * nr)
    r21 = 1 - t21

    # top surface side
    denom = 1.0 - r21 * r21 * tau * tau
    Ta = talf * tau * t21 / denom
    Ra = ralf + r21 * tau * Ta

    # bottom surface side
    t = t12 * tau * t21 / denom
    r = r12 + r21 * tau * t

    return r, t, Ra, Ta, denom


def prospect_d(
    N,
    cab,
    car,
    cbrown,
    cw,
    cm,
    ant,
    prot,
    cbc,
    nr,
    kab,
    kcar,
    kbrown,
    kw,
    km,
    kant,
    kprot,
    kcbc,
    alpha=40.0,
):

    lambdas = np.arange(400, 2501)  # wavelengths
    n_lambdas = len(lambdas)
    n_elems_list = [
        len(spectrum)
        for spectrum in [nr, kab, kcar, kbrown, kw, km, kant, kprot, kcbc]
    ]
    if not all(n_elems == n_lambdas for n_elems in n_elems_list):
        raise ValueError("Leaf spectra don't have the right shape!")

    kall = (
        cab * kab
        + car * kcar
        + ant * kant
        + cbrown * kbrown
        + cw * kw
        + cm * km
        + prot * kprot
        + cbc * kcbc
    ) / N
    j = kall > 0
    t1 = (1 - kall) * np.exp(-kall)
    t2 = kall ** 2 * (-expi(-kall))
    tau = np.ones_like(t1)
    tau[j] = t1[j] + t2[j]

    r, t, Ra, Ta, denom = refl_trans_one_layer(alpha, nr, tau)

    # ***********************************************************************
    # reflectance and transmittance of N layers
    # Stokes equations to compute properties of next N-1 layers (N real)
    # Normal case
    # ***********************************************************************
    # Stokes G.G. (1862), On the intensity of the light reflected from
    # or transmitted through a pile of plates, Proc. Roy. Soc. Lond.,
    # 11:545-556.
    # ***********************************************************************
    D = np.sqrt((1 + r + t) * (1 + r - t) * (1.0 - r + t) * (1.0 - r - t))
    rq = r * r
    tq = t * t
    a = (1 + rq - tq + D) / (2 * r)
    b = (1 - rq + tq + D) / (2 * t)

    bNm1 = np.power(b, N - 1)
    bN2 = bNm1 * bNm1
    a2 = a * a
    denom = a2 * bN2 - 1
    Rsub = a * (bN2 - 1) / denom
    Tsub = bNm1 * (a2 - 1) / denom

    # Case of zero absorption
    j = r + t >= 1.0
    Tsub[j] = t[j] / (t[j] + (1 - t[j]) * (N - 1))
    Rsub[j] = 1 - Tsub[j]

    # Reflectance and transmittance of the leaf: combine top layer with next N-1 layers
    denom = 1 - Rsub * r
    tran = Ta * Tsub / denom
    refl = Ra + Ta * Rsub * t / denom

    return lambdas, refl, tran


###if __name__ == "__main__":

###k_cab = prosail.spectral_libs.k_cab
###k_w = prosail.spectral_libs.k_cw
###k_m = prosail.spectral_libs.k_cm
###k_car = prosail.spectral_libs.k_car
###k_brown = prosail.spectral_libs.k_brown
###nr = prosail.spectral_libs.refractive

###wv, r, t = prospect_d (2.1, 60., 10., 0.1, 0.013, 0.016, 0,
###nr, k_cab, k_car, k_brown, k_w, k_m, k_m*0.,
###alpha=40.)

###rt = prosail.prospect_5b(2.1, 60., 10., 0.1, 0.013, 0.016)
###plt.plot(wv, r-rt[:,0], '--')
###plt.plot(wv, t-rt[:,1], '--')

###wv, r, t = prospect_d (2.1, 10., 10., 0.1, 0.013, 0.016, 0,
###nr, k_cab, k_car, k_brown, k_w, k_m, k_m*0.,
###alpha=40.)
####    plt.plot(wv, r)
###rt = prosail.prospect_5b(2.1, 10., 10., 0.1, 0.013, 0.016)
###plt.plot(wv, r-rt[:,0], '-')
###plt.plot(wv, t-rt[:,1], '-')
