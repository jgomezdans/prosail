PROSAIL
==========

Description
--------------

This repository contains the Python bindings to the PROSPECT and SAIL leaf and 
canopy reflectance models. The code is written in FORTRAN. The original fortran
code was downloaded from `Jussieu <http://teledetection.ipgp.jussieu.fr/prosail/>`_. 
I only added a function to simplify writing the wrappers, and you should go to
that page to get newer versions of the code. A recent review of the use of both
models is availabe in `this paper <http://webdocs.dow.wur.nl/internet/grs/Workshops/Environmental_Applications_Imaging_Spectroscopy/12_Jacquemoud_Prospect/IEEE_Jacquemoud_PROSPECT.pdf>`.


Installing the bindings
-------------------------

The installation of the bindings is quite straightforward: unpack the distribution
and run the following command:..

    python setup.py install
    
You can usually install it to your user's directory (if you haven't got superuser
privileges) by ..

    python setup.py install --user
    
.. note::
    
    You will need a working FORTRAN compiler. I have only tested this with GCC on Linux, but it should work on other systems. You can also pass optimisation flags to the compiler: ..
    
    python setup.py config_fc  --fcompiler=gnu95   --arch=-march=native --opt=-O3  install --user
    
    
Using the bindings
---------------------

The bindings offer a single function, `run_prosail`.