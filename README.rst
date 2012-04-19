PROSAIL
==========

Description
--------------

This repository contains the Python bindings to the PROSPECT and SAIL leaf and 
canopy reflectance models. The code is written in FORTRAN. The original fortran
code was downloaded from `Jussieu <http://teledetection.ipgp.jussieu.fr/prosail/>`_. 
I only added a function to simplify writing the wrappers, and you should go to
that page to get newer versions of the code. A recent review of the use of both
models is availabe in `this paper <http://webdocs.dow.wur.nl/internet/grs/Workshops/Environmental_Applications_Imaging_Spectroscopy/12_Jacquemoud_Prospect/IEEE_Jacquemoud_PROSPECT.pdf>`_.


Installing the bindings
-------------------------

The installation of the bindings is quite straightforward: unpack the distribution
and run the following command   

    python setup.py install
    
You can usually install it to your user's directory (if you haven't got superuser
privileges) by 

    python setup.py install --user
    
Note
*******
    
You will need a working FORTRAN compiler. I have only tested this with GCC on Linux, but it should work on other systems. You can also pass optimisation flags to the compiler: 
    
    python setup.py config_fc  --fcompiler=gnu95   --arch=-march=native --opt=-O3  install --user
    
    
Using the bindings
---------------------

The bindings offer a single function, ``run_prosail``.

    retval = run_prosail(n,cab,car,cbrown,cw,cm,lai,lidfa,lidfb,psoil,hspot,tts,tto,psi)
    
The parameters that are fed into ``run_prosail`` are

+-------------+---------------------------------+--------------+------------+-------------+
| Parameter   | Description of parameter        | Units        |Typical min | Typical max |
+=============+=================================+==============+============+=============+
|   N         | Leaf structure parameter        | N/A          | 0.8        | 2.5         |
+-------------+---------------------------------+--------------+------------+-------------+
|  cab        | Chlorophyll a+b concentration   | ug/cm2       | 0          | 200         |
+-------------+---------------------------------+--------------+------------+-------------+
|  caw        | Equivalent water thickiness     | cm           | 0          | 200         |
+-------------+---------------------------------+--------------+------------+-------------+
|  car        | Carotenoid concentration        | ug/cm2       | 0          | 200         |
+-------------+---------------------------------+--------------+------------+-------------+
|  cbrown     | Brown pigment                   | NA           | 0          | 1           |
+-------------+---------------------------------+--------------+------------+-------------+
|  cm         | Dry matter content              | g/cm2        | 0          | 200         |
+-------------+---------------------------------+--------------+------------+-------------+
|  lai        | Leaf Area Index                 | N/A          | 0.1        | 10          |
+-------------+---------------------------------+--------------+------------+-------------+
|  lidfa      | Leaf angle distribution         | N/A          | -          | -           |
+-------------+---------------------------------+--------------+------------+-------------+
|  lidfb      | Leaf angle distribution         | N/A          | -          | -           |
+-------------+---------------------------------+--------------+------------+-------------+
|  psoil      | Dry/Wet soil factor             | N/A          | 0          | 1           |
+-------------+---------------------------------+--------------+------------+-------------+
|  hspot      | Hotspot parameter               | N/A          | 0          | 0.          |
+-------------+---------------------------------+--------------+------------+-------------+
|  tts        | Solar zenith angle              | deg          | 0          | 90          |
+-------------+---------------------------------+--------------+------------+-------------+
|  tto        | Observer zenith angle           | deg          | 0          | 90          |
+-------------+---------------------------------+--------------+------------+-------------+
|  phi        | Relative azimuth angle          | deg          | 0          | 360         |
+-------------+---------------------------------+--------------+------------+-------------+

``lidfa`` and ``lidfb`` parameters control the leaf angle distribution. Typical distributions
are given by the following parameter  choices:

+--------------+-----------+------------------+
|LIDF type     |  lidfa    |    lidfb         |
+==============+===========+==================+
|Planophile    |    1      |  b               |
+--------------+-----------+------------------+
|   Erectophile|    -1     |   0              |
+--------------+-----------+------------------+
|   Plagiophile|     0     |  -1              |
+--------------+-----------+------------------+
|  Extremophile|    0      |  1               |
+--------------+-----------+------------------+
|   Spherical  |    -0.35  |  -0.15           |
+--------------+-----------+------------------+
|   Uniform    |     0     |   0              |
+--------------+-----------+------------------+
   
   
    

A simple example of using bindings would be

 cm=0.009
 cab = 80
 car = 15
 cbrown = 0
 cw = 0.001
 lidfa = 0
 lidfb = 0
 psoil = 0
 hspot = 0.01
 tts = 30
 tto = 10
 phi = 0
 lai = np.arange ( 0, 5, 0.2 )
 for l in lai:
    plt.plot ( run_prosail(1.5,cab,car,cbrown,cw,cm,l,lidfa,lidfb,psoil,hspot,tts,tto,psi))
    
This results in a simulation of surface reflectance (from 400 to 2500 nm) as a function of LAI and under the other parameters' prescribed values. This yields the following

.. image:: http://i.imgur.com/2Hh0z.png