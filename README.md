<p><img src="https://www.nceo.ac.uk/wp-content/themes/nceo/assets/images/logos/img_logo_purple.svg" align="left" />

<img src="http://www.esa.int/esalogo/images/logotype/img_colorlogo_darkblue.gif" scale="20%" align="right" />
</p>

<br/>
<br/>

---



# PROSAIL Python Bindings

#### J Gomez-Dans (NCEO & UCL) ``j.gomez-dans@ucl.ac.uk``

[![DOI](https://zenodo.org/badge/19469/jgomezdans/prosail.svg)](https://zenodo.org/badge/latestdoi/19469/jgomezdans/prosail)

[![Build Status](https://travis-ci.org/jgomezdans/prosail.png)](https://travis-ci.org/jgomezdans/prosail)

[![Coverage Status](https://coveralls.io/repos/github/jgomezdans/prosail/badge.svg?branch=master)](https://coveralls.io/github/jgomezdans/prosail?branch=master)
[![codecov](https://codecov.io/gh/jgomezdans/prosail/branch/master/graph/badge.svg?longCache=true&style=flat)](https://codecov.io/gh/jgomezdans/prosail)
[![Anaconda-Server Badge](https://anaconda.org/jgomezdans/prosail/badges/version.svg)](https://anaconda.org/jgomezdans/prosail)
[![PyPI version](https://badge.fury.io/py/prosail.svg)](https://badge.fury.io/py/prsoail)


## Install using Anaconda

You should be able to easily install this using Anaconda (only tested on Linux!) with

`conda install -c jgomezdans prosail`

I **think** it might work on both Python 2.7 and 3.6. But I'm only a scientist, so expect car crashes!


## Description

This repository contains the Python bindings to the PROSPECT and SAIL leaf and 
canopy reflectance models, respectively. Both models have been rewritten and
coupled in Python, with some changes to improve on efficiency. The bindings implement
the following models:

* **PROSPECT**: versions 5 and D. Flexibility to add/modify leaf absorption profiles.
* **SAIL**: FourSAIL version. The thermal extension of the model is also implemented, although this hasn't been widely tested.
* Simple Lambertian soil reflectance model

I have used as a benchmark the codes available from [Jussieu](http://teledetection.ipgp.jussieu.fr/prosail/). 

A recent(ish) review on the use of both RT
models is availabe in [this paper](http://webdocs.dow.wur.nl/internet/grs/Workshops/Environmental_Applications_Imaging_Spectroscopy/12_Jacquemoud_Prospect/IEEE_Jacquemoud_PROSPECT.pdf)_.


## Installing the bindings

The installation of the bindings is quite straightforward: unpack the distribution
and run the following command   

    python setup.py install
    
This assumes that you have the following things installed:

* [Numpy](http://www.numpy.org/)
* [Scipy](http://www.scipy.org/)
* The [LRU_Cache library for Python 2.7 ](https://pypi.python.org/pypi/backports.functools_lru_cache/1.0.1)
* [Numba](http://numba.pydata.org)

Most of these things can be installed quite easily using [Anaconda Python](https://www.continuum.io/downloads). 
In this case, you can probably just install everything you need with

      conda install python=2.7 numpy numba scipy
      pip install -U backports.functools_lru_cache
      
The bindings should then install without any issue.


## Using the bindings

Once you import the bindings into the namespace with

    import prosail
    
you can then run SAIL (using prescribed leaf reflectance and transmittance spectra, as well as canopy structure/soil parameters), PROSPECT and both (e.g. use PROSPECT to provide the spectral leaf optical properties).

### `run_sail`

To run SAIl with two element arrays of leaf reflectance and transmittance sampled at 1nm between 400 and 2500 nm `rho` and `tau`, using a black soil (e.g. zero reflectance), you can just do 

    rho_canopy = prosail.run_sail(rho, tau, lai, lidfa, hspot, sza, vza, raa, rsoil0=np.zeros(2101))

Here, `lai` is the LAI, `lidfa` is the mean leaf angle in degrees, `hspot` is the hotspot parameter, `sza`, `vza` and `raa` are the solar zenith, sensor zenith and relative azimuth angles, and `rsoil0` is set to an array of 0s to define the soil reflectance.

You have quite a few other options:

* You can use a different way of specifying the leaf angle distribution (by default we use a Campbell distribution with one single parameter, but you might want to use the Verhoef distribution). The Verhoef distribution is selected by adding the extra keyword `typelidf=1` and the two parameters are given by `lidfa` and the additional optional parameter `lidfb`.
* You can use the internal soil spectrum model. This model is basically `rho_soil = rsoil*(psoil*soil_spectrum1+(1-psoil)*soil_spectrum2)`. The first spectrum is a dry soil, the second one a wet one. You can also set the spectra using the `soil_spectrum1` and `soil_spectrum2` keywords.
* By default, we return the surface directional reflectance, but you can choose other reflectance factors (e.g. BHR, DHR, HDR).


### `run_prospect`

To calculate leaf reflectance and transmittance using the PROSPECT model, you can use the `run_prospect` function. You can select either the PROSPECT-5 or PROSPECT-D versions (by default, version 'D' is used). A call to this would look like:
   
    lam, rho, tau = prosail.run_prospect(n, cab, car, cbrown, cw, cm, ant=8.0)
    
Where the parameters are all scalars, and have their usual PROSPECT meanings (see table below). `ant` stands for anthocyannins, which isn't present in PROSPECT-5.

To do the same for PROSPECT-5...

    lam, rho, tau = prosail.run_prospect(n, cab, car, cbrown, cw, cm, prospect_version='5')
    
You can change a number of things when calling PROSPECT, but I can't be arsed documenting it now.

### `run_prosail`

The marriage of heaven and hell, PROSPECT being fed into SAIL in one go! Same options as the two other functions put together:

    rho_canopy = prosail.run_prosail(n, cab, car, cbrown, cw, cm, lai, lidfa, hspot, tts, tto, psi, \
                        ant=0.0, alpha=40.0, prospect_version='5', typelidf=2, lidfb=0.0, \
                        factor='SDR', rsoil0=None, rsoil=None, psoil=None, \
                        soil_spectrum1=None, soil_spectrum2=None)
    

   


## The parameters

The parameters used by the models and their units are introduced below:

| Parameter   | Description of parameter        | Units        |Typical min | Typical max |
|-------------|---------------------------------|--------------|------------|-------------|
|   N         | Leaf structure parameter        | N/A          | 0.8        | 2.5         |
|  cab        | Chlorophyll a+b concentration   | ug/cm2       | 0          | 80          |
|  caw        | Equivalent water thickiness     | cm           | 0          | 200         |
|  car        | Carotenoid concentration        | ug/cm2       | 0          | 20          |
|  cbrown     | Brown pigment                   | NA           | 0          | 1           |
|  cm         | Dry matter content              | g/cm2        | 0          | 200         |
|  lai        | Leaf Area Index                 | N/A          | 0          | 10          |
|  lidfa      | Leaf angle distribution         | N/A          | -          | -           |
|  lidfb      | Leaf angle distribution         | N/A          | -          | -           |
|  psoil      | Dry/Wet soil factor             | N/A          | 0          | 1           |
|  rsoil      | Soil brigthness factor          | N/A          | -          | -           |
|  hspot      | Hotspot parameter               | N/A          | -          | -           |
|  tts        | Solar zenith angle              | deg          | 0          | 90          |
|  tto        | Observer zenith angle           | deg          | 0          | 90          |
|  phi        | Relative azimuth angle          | deg          | 0          | 360         |
| typelidf    | Leaf angle distribution type    | Integer      | -          | -           |

### Specifying the leaf angle distribution

The parameter ``typelidf`` regulates the leaf angle distribution family being used. The following options are understood:

* ``typelidf = 1``: use the two parameter LAD parameterisation, where ``a`` and ``b`` control the average leaf slope and the distribution bimodality, respectively. Typical distributions
are given by the following parameter  choices:

| LIDF type    | ``LIDFa`` |  ``LIDFb``       |
|--------------|-----------|------------------|
| Planophile   |    1      |  0               |
|   Erectophile|    -1     |   0              |
|   Plagiophile|     0     |  -1              |
|  Extremophile|    0      |  1               |
|   Spherical  |    -0.35  |  -0.15           |
|   Uniform    |     0     |   0              |

* ``typelidf = 2`` Ellipsoidal distribution, where ``LIDFa`` parameter stands for mean leaf angle (0 degrees is planophile, 90 degrees is erectophile). ``LIDFb`` parameter is ignored.
   
### The soil model

The soil model is a fairly simple linear mixture model, where two spectra are mixed and then a brightness term added:

    rho_soil = rsoil*(psoil*soil_spectrum1+(1-psoil)*soil_spectrum2)


The idea is that one of the spectra is a dry soil and the other a wet soil, so soil moisture is then contorlled by ``psoil``. ``rsoil`` is just a brightness scaling term.


