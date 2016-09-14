#!/usr/bin/env python
"""Setup script for building prosail's python bindings"""

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration(parent_package,top_path)
    config.add_extension('prosail_fortran', ['prosail/dataSpec_P5B.f90', \
        "prosail/dladgen.f", "prosail/LIDF.f90", "prosail/MODULE_PRO4SAIL.f90", \
        "prosail/PRO4SAIL.f90", "prosail/prospect_5B.f90", \
        "prosail/run_prosail.f90","prosail/tav_abs.f90", \
        "prosail/volscatt.f90","prosail/prosail_fortran.pyf"] )
    return config

if __name__ == "__main__":
    import os
    from numpy.distutils.core import setup
    # Global variables for this extension:
    name         = "prosail"  # name of the generated python extension (.so)
    description  = "PROSPECT, SAIL and PROSAIL Python wrappers"
    long_description = "The PROSPECT + SAILh radiative transfer models from Python."
    if os.path.exists ( "README.txt" ):
        long_description = open( "README.txt", 'r').read()

    author       = "J Gomez-Dans/NCEO & University College London"
    author_email = "j.gomez-dans@ucl.ac.uk"
    url = "http://github.com/jgomezdans/prosail"
    classifiers=[
	'Development Status :: 5 - Production/Stable',
        'Natural Language :: English',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 2',
        'Topic :: Scientific/Engineering',
        'Topic :: Software Development :: Libraries :: Python Modules',
        "Topic :: Scientific/Engineering :: Atmospheric Science",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: GIS",
        'Intended Audience :: Science/Research',
        'Intended Audience :: End Users/Desktop',
        'Intended Audience :: Developers',
        'Environment :: Console']
    setup( name=name,
        description=description,
        long_description=long_description,
        author=author,
	url=url,
        author_email = author_email,
        classifiers = classifiers,
        configuration = configuration,
        version="1.3.0",
        packages=["prosail"])
