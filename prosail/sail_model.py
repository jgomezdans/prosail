#!/usr/bin/env python
import numpy as np

from prosail import spectral_lib

from .FourSAIL import foursail
from .prospect_d import run_prospect


def calculate_fcover_fapar_albedo(
    n,
    cab,
    car,
    cbrown,
    cw,
    cm,
    lai,
    lidfa,
    hspot,
    tts,
    tto,
    psi,
    ant=0.0,
    prot=0.0,
    cbc=0.0,
    alpha=40.0,
    prospect_version="5",
    typelidf=2,
    lidfb=0.0,
    rsoil0=None,
    rsoil=None,
    psoil=None,
    soil_spectrum1=None,
    soil_spectrum2=None,
    fapar_slice=slice(0, 300),
    albedo_slice=slice(0, -100),
):
    """A function to calculate fCover, fAPAR and albedo.
    fCover is just `1-too`. You can select the spectral
    range of albedo and fAPAR using the slice notation. By
    default fAPAR is calculated between 400 and 700 nm,
    and albedo is calculated between 400 and 2100 nm.

    Parameters
    -----------
    n: float
        The number of leaf layers. Unitless [-].
    cab: float
        The chlorophyll a+b concentration. [g cm^{-2}].
    car: float
        Carotenoid concentration.  [g cm^{-2}].
    cbrown: float
        The brown/senescent pigment. Unitless [-], often between 0 and 1
        but the literature on it is wide ranging!
    cw: float
        Equivalent leaf water. [cm]
    cm: float
        Dry matter [g cm^{-2}]
    lai: float
        leaf area index
    lidfa: float
        a parameter for leaf angle distribution. If ``typliedf``=2, average
        leaf inclination angle.
    tts: float
        Solar zenith angle
    tto: float
        Sensor zenith angle
    psi: float
        Relative sensor-solar azimuth angle ( saa - vaa )
    ant: float, optional
        Anthocyanins content. Used in Prospect-D and Prospect-PRO [g cm^{-2}]
    prot: float, optional
        Protein content. Used in Prospect-PRO. [g cm^{-2}]
    cbc: float, optional
        Carbon based constituents. Used in Prospect-PRO. [g cm^{-2}]
    alpha: float
        The alpha angle (in degrees) used in the surface scattering
        calculations. By default it's set to 40 degrees.
    prospect_version: str
        Which PROSPECT version to use. We have "5", "D" and "PRO"
    typelidf: int, optional
        The type of leaf angle distribution function to use. By default, is set
        to 2.
    lidfb: float, optional
        b parameter for leaf angle distribution. If ``typelidf``=2, ignored
    rsoil0: float, optional
        The soil reflectance spectrum
    rsoil: float, optional
        Soil scalar 1 (brightness)
    psoil: float, optional
        Soil scalar 2 (moisture)
    soil_spectrum1: 2101-element array
        First component of the soil spectrum
    soil_spectrum2: 2101-element array
        Second component of the soil spectrum
    fapar_slice: slice
        Wavelengths to use for fAPAR calculation. Selects slice positions for
        the interval 400-2500 (inclusive). By default, use `slice(0,300)`,
        equivalent to 400 to 700 nm.
    albedo_slice: slice
        Wavelengths to use for albedo calculation. Selects slice positions for
        the interval 400-2500 (inclusive). By default, use `slice(0,-100)`,
        equivalent to 400 to 2100 nm.
    Returns
    --------
    fCover, fAPAR, albedo


    """
    (
        tss,
        too,
        tsstoo,
        rdd,
        tdd,
        rsd,
        tsd,
        rdo,
        tdo,
        rso,
        rsos,
        rsod,
        rddt,
        rsdt,
        rdot,
        rsodt,
        rsost,
        rsot,
        gammasdf,
        gammasdb,
        gammaso,
    ) = run_prosail(
        n,
        cab,
        car,
        cbrown,
        cw,
        cm,
        lai,
        lidfa,
        hspot,
        tts,
        tto,
        psi,
        ant=ant,
        prot=prot,
        cbc=cbc,
        alpha=alpha,
        prospect_version=prospect_version,
        typelidf=typelidf,
        lidfb=lidfb,
        factor="ALLALL",
        rsoil0=rsoil0,
        rsoil=rsoil,
        psoil=psoil,
        soil_spectrum1=soil_spectrum1,
        soil_spectrum2=soil_spectrum2,
    )
    # fCover
    fCover = 1 - too

    # direct diffuse illumination
    es = spectral_lib.light.es
    ed = spectral_lib.light.ed
    skyl = (
        0.847
        - 1.61 * np.sin(np.deg2rad(90 - tts))
        + 1.04 * np.sin(np.deg2rad(90 - tts)) * np.sin(np.deg2rad(90 - tts))
    )
    edir = (1 - skyl) * es
    edif = skyl * ed
    # Interaction with the soil
    dn = 1.0 - rsoil * rdd
    try:
        dn[dn < 1e-36] = 1e-36
    except TypeError:
        dn = max(1e-36, dn)
    # fAPAR calculations
    # Caclulate absorbances of the isolated canopy
    alfas = 1.0 - tss - tsd - rsd  # direct flux
    alfad = 1.0 - tdd - rdd  # diffuse
    alfasx = alfas + (rsoil * (tss + tsd) / dn) * alfad
    alfadx = alfad + ((tdd * rsoil) / dn) * alfas
    top = alfasx * edir + alfadx * edif
    fAPAR = np.sum((top[fapar_slice])) / np.sum((edir + edif)[fapar_slice])

    # Albedo calculatins
    rsdstar = rsd + (tss + tsd) * rsoil * tdd / dn
    rddstar = rdd + (tdd * tdd * rsoil) / dn

    num = ((rsdstar * edir + rddstar * edif)[albedo_slice]).sum()
    den = ((edir + edif)[albedo_slice]).sum()
    albedo = num / den

    return fCover, fAPAR, albedo


def run_prosail(
    n,
    cab,
    car,
    cbrown,
    cw,
    cm,
    lai,
    lidfa,
    hspot,
    tts,
    tto,
    psi,
    ant=0.0,
    prot=0.0,
    cbc=0.0,
    alpha=40.0,
    prospect_version="5",
    typelidf=2,
    lidfb=0.0,
    factor="SDR",
    rsoil0=None,
    rsoil=None,
    psoil=None,
    soil_spectrum1=None,
    soil_spectrum2=None,
):
    """Run the PROSPECT_5B and SAILh radiative transfer models. The soil
    model is a linear mixture model, where two spectra are combined together as

         rho_soil = rsoil*(psoil*soil_spectrum1+(1-psoil)*soil_spectrum2)
    By default, ``soil_spectrum1`` is a dry soil, and ``soil_spectrum2`` is a
    wet soil, so in that case, ``psoil`` is a surface soil moisture parameter.
    ``rsoil`` is a  soil brightness term. You can provide one or the two
    soil spectra if you want.  The soil spectra must be defined
    between 400 and 2500 nm with 1nm spacing.

    Parameters
    -----------
    n: float
        The number of leaf layers. Unitless [-].
    cab: float
        The chlorophyll a+b concentration. [g cm^{-2}].
    car: float
        Carotenoid concentration.  [g cm^{-2}].
    cbrown: float
        The brown/senescent pigment. Unitless [-], often between 0 and 1
        but the literature on it is wide ranging!
    cw: float
        Equivalent leaf water. [cm]
    cm: float
        Dry matter [g cm^{-2}]
    lai: float
        leaf area index
    lidfa: float
        a parameter for leaf angle distribution. If ``typliedf``=2, average
        leaf inclination angle.
    tts: float
        Solar zenith angle
    tto: float
        Sensor zenith angle
    psi: float
        Relative sensor-solar azimuth angle ( saa - vaa )
    ant: float, optional
        Anthocyanins content. Used in Prospect-D and Prospect-PRO [g cm^{-2}]
    prot: float, optional
        Protein content. Used in Prospect-PRO. [g cm^{-2}]
    cbc: float, optional
        Carbon based constituents. Used in Prospect-PRO. [g cm^{-2}]
    alpha: float
        The alpha angle (in degrees) used in the surface scattering
        calculations. By default it's set to 40 degrees.
    prospect_version: str
        Which PROSPECT version to use. We have "5", "D" and "PRO"
    typelidf: int, optional
        The type of leaf angle distribution function to use. By default, is set
        to 2.
    lidfb: float, optional
        b parameter for leaf angle distribution. If ``typelidf``=2, ignored
    factor: str, optional
        What reflectance factor to return:
        * "SDR": directional reflectance factor (default)
        * "BHR": bi-hemispherical r. f.
        * "DHR": Directional-Hemispherical r. f. (directional illumination)
        * "HDR": Hemispherical-Directional r. f. (directional view)
        * "ALL": All of them
        * "ALLALL": All of the terms calculated by SAIL, including the above
    rsoil0: float, optional
        The soil reflectance spectrum
    rsoil: float, optional
        Soil scalar 1 (brightness)
    psoil: float, optional
        Soil scalar 2 (moisture)
    soil_spectrum1: 2101-element array
        First component of the soil spectrum
    soil_spectrum2: 2101-element array
        Second component of the soil spectrum
    Returns
    --------
    A reflectance factor between 400 and 2500 nm


    """

    if soil_spectrum1 is not None:
        assert len(soil_spectrum1) == 2101
    else:
        soil_spectrum1 = spectral_lib.soil.rsoil1

    if soil_spectrum2 is not None:
        assert len(soil_spectrum1) == 2101
    else:
        soil_spectrum2 = spectral_lib.soil.rsoil2

    if rsoil0 is None:
        if (rsoil is None) or (psoil is None):
            raise ValueError(
                "If rsoil0 isn't define, then rsoil and psoil"
                " need to be defined!"
            )
        else:
            rsoil0 = rsoil * (
                psoil * soil_spectrum1 + (1.0 - psoil) * soil_spectrum2
            )

    _, refl, trans = run_prospect(
        n,
        cab,
        car,
        cbrown,
        cw,
        cm,
        ant=ant,
        prot=prot,
        cbc=cbc,
        prospect_version=prospect_version,
        alpha=alpha,
    )

    return run_sail(
        refl,
        trans,
        lai,
        lidfa,
        hspot,
        tts,
        tto,
        psi,
        typelidf=typelidf,
        lidfb=lidfb,
        factor=factor,
        rsoil0=rsoil0,
    )


def run_sail(
    refl,
    trans,
    lai,
    lidfa,
    hspot,
    tts,
    tto,
    psi,
    typelidf=2,
    lidfb=0.0,
    factor="SDR",
    rsoil0=None,
    rsoil=None,
    psoil=None,
    soil_spectrum1=None,
    soil_spectrum2=None,
):
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
    hspot: float
        The hotspot parameter
    tts: float
        Solar zenith angle
    tto: float
        Sensor zenith angle
    psi: float
        Relative sensor-solar azimuth angle ( saa - vaa )
    typelidf: int, optional
        The type of leaf angle distribution function to use. By default, is set
        to 2.
    lidfb: float, optional
        b parameter for leaf angle distribution. If ``typelidf``=2, ignored
    factor: str, optional
        What reflectance factor to return:
        * "SDR": directional reflectance factor (default)
        * "BHR": bi-hemispherical r. f.
        * "DHR": Directional-Hemispherical r. f. (directional illumination)
        * "HDR": Hemispherical-Directional r. f. (directional view)
        * "ALL": All of them
        * "ALLALL": All of the terms calculated by SAIL, including the above
    rsoil0: float, optional
        The soil reflectance spectrum
    rsoil: float, optional
        Soil scalar 1 (brightness)
    psoil: float, optional
        Soil scalar 2 (moisture)
    soil_spectrum1: 2101-element array
        First component of the soil spectrum
    soil_spectrum2: 2101-element array
        Second component of the soil spectrum

    Returns
    --------
    Directional surface reflectance between 400 and 2500 nm


    """

    factor = factor.upper()
    if factor not in ["SDR", "BHR", "DHR", "HDR", "ALL", "BRF", "ALLALL"]:
        raise ValueError(
            "'factor' must be one of BRF, SDR, BHR, DHR, HDR, ALL or ALLALL"
        )
    if soil_spectrum1 is not None:
        assert len(soil_spectrum1) == 2101
    else:
        soil_spectrum1 = spectral_lib.soil.rsoil1

    if soil_spectrum2 is not None:
        assert len(soil_spectrum1) == 2101
    else:
        soil_spectrum2 = spectral_lib.soil.rsoil2

    if rsoil0 is None:
        if (rsoil is None) or (psoil is None):
            raise ValueError(
                "If rsoil0 isn't define, then rsoil and psoil"
                " need to be defined!"
            )
        else:
            rsoil0 = rsoil * (
                psoil * soil_spectrum1 + (1.0 - psoil) * soil_spectrum2
            )

    [
        tss,
        too,
        tsstoo,
        rdd,
        tdd,
        rsd,
        tsd,
        rdo,
        tdo,
        rso,
        rsos,
        rsod,
        rddt,
        rsdt,
        rdot,
        rsodt,
        rsost,
        rsot,
        gammasdf,
        gammasdb,
        gammaso,
    ] = foursail(
        refl, trans, lidfa, lidfb, typelidf, lai, hspot, tts, tto, psi, rsoil0
    )
    es = spectral_lib.light.es
    ed = spectral_lib.light.ed
    skyl = (
        0.847
        - 1.61 * np.sin(np.deg2rad(90 - tts))
        + 1.04 * np.sin(np.deg2rad(90 - tts)) * np.sin(np.deg2rad(90 - tts))
    )
    edir = (1 - skyl) * es
    edif = skyl * ed
    BRF = (rdot * edif + rsot * edir) / (edir + edif)
    if factor == "SDR":
        return rsot
    elif factor == "BRF":
        return BRF
    elif factor == "BHR":
        return rddt
    elif factor == "DHR":
        return rsdt
    elif factor == "HDR":
        return rdot
    elif factor == "ALL":
        return [rsot, rddt, rsdt, rdot]
    elif factor == "ALLALL":
        return [
            tss,
            too,
            tsstoo,
            rdd,
            tdd,
            rsd,
            tsd,
            rdo,
            tdo,
            rso,
            rsos,
            rsod,
            rddt,
            rsdt,
            rdot,
            rsodt,
            rsost,
            rsot,
            gammasdf,
            gammasdb,
            gammaso,
        ]


def run_thermal_sail(
    lam,
    tveg,
    tsoil,
    tveg_sunlit,
    tsoil_sunlit,
    t_atm,
    lai,
    lidfa,
    hspot,
    tts,
    tto,
    psi,
    rsoil=None,
    refl=None,
    emv=None,
    ems=None,
    typelidf=2,
    lidfb=0,
):
    c1 = 3.741856e-16
    c2 = 14388.0
    # Calculate the thermal emission from the different
    # components using Planck's Law
    top = (1.0e-6) * c1 * (lam * 1e-6) ** (-5.0)
    Hc = top / (np.exp(c2 / (lam * tveg)) - 1.0)  # Shade leaves
    Hh = top / (np.exp(c2 / (lam * tveg_sunlit)) - 1.0)  # Sunlit leaves
    Hd = top / (np.exp(c2 / (lam * tsoil)) - 1.0)  # shade soil
    Hs = top / (np.exp(c2 / (lam * tsoil_sunlit)) - 1.0)  # Sunlit soil
    Hsky = top / (np.exp(c2 / (lam * t_atm)) - 1.0)  # Sky emission

    # Emissivity calculations
    if refl is not None and emv is None:
        emv = 1.0 - refl  # Assuming absorption is 1

    if rsoil is not None and ems is None:
        ems = 1.0 - rsoil

    if rsoil is None and ems is not None:
        rsoil = 1.0 - ems
    if refl is None and emv is not None:
        refl = 1.0 - emv

    [
        tss,
        too,
        tsstoo,
        rdd,
        tdd,
        rsd,
        tsd,
        rdo,
        tdo,
        rso,
        rsos,
        rsod,
        rddt,
        rsdt,
        rdot,
        rsodt,
        rsost,
        rsot,
        gammasdf,
        gammasdb,
        gammaso,
    ] = foursail(
        refl,
        np.zeros_like(refl),
        lidfa,
        lidfb,
        typelidf,
        lai,
        hspot,
        tts,
        tto,
        psi,
        rsoil,
    )

    gammad = 1.0 - rdd - tdd
    gammao = 1.0 - rdo - tdo - too

    # tso = tss * too + tss * (tdo + rsoil * rdd * too) / (1.0 - rsoil * rdd)
    tso = tsstoo + tss * (tdo + rsoil * rdd * too) / (1.0 - rsoil * rdd)
    ttot = (too + tdo) / (1.0 - rsoil * rdd)
    gammaot = gammao + ttot * rsoil * gammad
    gammasot = gammaso + ttot * rsoil * gammasdf

    aeev = gammaot
    aees = ttot * ems

    Lw = (
        (rdot * Hsky) / np.pi
        + (
            aeev * Hc
            + gammasot * emv * (Hh - Hc)
            + aees * Hd
            + tso * ems * (Hs - Hd)
        )
    ) / np.pi

    dnoem1 = top / (Lw * np.pi)
    Tbright = c2 / (lam * np.log(dnoem1 + 1.0))
    dir_em = 1.0 - rdot
    return Lw, Tbright, dir_em
