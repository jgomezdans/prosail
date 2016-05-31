<img src="https://www.nceo.ac.uk/wp-content/themes/nceo/assets/images//logos/img_logo_white.svg" scale=50% alt="NCEO logo" align="right" />
<img src="http://www.esa.int/esalogo/images/logotype/img_colorlogo_darkblue.gif" scale=20% alt="ESA logo" align="left" />

<br/>
<br/>
<br/>

---

# PROSAIL Python Bindings

#### J Gomez-Dans (NCEO & UCL) ``j.gomez-dans@ucl.ac.uk``

[![DOI](https://zenodo.org/badge/19469/jgomezdans/prosail.svg)](https://zenodo.org/badge/latestdoi/19469/jgomezdans/prosail)

## Description

This repository contains the Python bindings to the PROSPECT and SAIL leaf and 
canopy reflectance models. The code is written in FORTRAN. The original fortran
code was downloaded from [Jussieu](http://teledetection.ipgp.jussieu.fr/prosail/). 
I only added a function to simplify writing the wrappers, and you should go to
that page to get newer versions of the code. A recent review of the use of both
models is availabe in [this paper](http://webdocs.dow.wur.nl/internet/grs/Workshops/Environmental_Applications_Imaging_Spectroscopy/12_Jacquemoud_Prospect/IEEE_Jacquemoud_PROSPECT.pdf)_.


## Installing the bindings

The installation of the bindings is quite straightforward: unpack the distribution
and run the following command   

    python setup.py install
    
You can usually install it to your user's directory (if you haven't got superuser
privileges) by 

    python setup.py install --user
    
#### **Note**

    
You will need a working FORTRAN compiler. I have only tested this with GCC on Linux, but it should work on other systems. You can also pass optimisation flags to the compiler: 
    
    python setup.py config_fc  --fcompiler=gnu95   --arch=-march=native --opt=-O3  install --user
    
## Using the bindings

The bindings offer several functions, which will be described in detail below:.

* ``run_prospect``: This function runs the PROSPECT 5B model in Feret et al 2008. The input parameters are the usual ``(n,cab,car,cbrown,cw,cm)`` (e.g. leaf layers, leaf Chlorophyll concentration, leaf Carotenoid concentration, leaf senescent fraction, Equivalent leaf water, leaf dry matter). It returns a spectrum between 400 and 2500 nm.
* ``run_sail``:  The SAILh model, which in this case requires leaf reflectance and transmittance to be fed to the model (e.g. you have already measured these spectra in the field). The rest of the parameters are ``(refl,trans,lai,lidfa,lidfb,rsoil,psoil,hspot,tts,tto,psi,typelidf)``. Additionally, there are two optional parameters, ``soil_spectrum1``, ``soil_spectrum2``, which allow you to set the soil spectra (otherwise, some default spectra get used). Output is a spectrum in the 400-2500 nm range.
* ``run_prosail``: PROSPECT_5B and SAILh coupled together, with input parameters given by ``(n,cab,car,cbrown,cw,cm,lai,lidfa,lidfb,rsoil,psoil,hspot,tts,tto,psi,typelidf)``. Output is a spectrum in the 400-2500 nm range.


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


