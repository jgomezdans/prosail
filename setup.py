#!/usr/bin/env python
"""Setup script for building prosail's python bindings"""
from distutils.core import setup
import os
# Global variables for this extension:
name = "prosail"  # name of the generated python extension (.so)
description = "PROSPECT, SAIL and PROSAIL Python wrappers"
long_description = "The PROSPECT + SAILh radiative transfer models from Python"


def read(filename):
    with open(os.path.join(os.path.dirname(__file__), filename)) as f:
        return f.read()


if os.path.exists("README.txt"):
    long_description = read("README.txt")

author = "J Gomez-Dans/NCEO & University College London"
author_email = "j.gomez-dans@ucl.ac.uk"
url = "http://github.com/jgomezdans/prosail"
classifiers = [
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
    'Environment :: Console'
]

setup(
    name=name,
    description=description,
    long_description=long_description,
    author=author,
    url=url,
    author_email=author_email,
    classifiers=classifiers,
    package_data={"prosail": ["*.txt"]},
    include_package_data=True,
    install_requires=[
        "numpy",
        "numba",
        "scipy",
        "backports.functools_lru_cache;python_version<\"3.2\"",
    ],
    version="2.0.0alpha",
    packages=["prosail"]
)
