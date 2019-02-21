#!/usr/bin/env python
"""Setup script for building prosail's python bindings"""
import os
import codecs
import re
from os import path
from setuptools import setup

# Global variables for this extension:
name = "prosail"  # name of the generated python extension (.so)
description = "PROSPECT, SAIL and PROSAIL Python wrappers"
long_description = "The PROSPECT + SAILh radiative transfer models from Python"

this_directory = path.abspath(path.dirname(__file__))                                                              

def read(filename):
    with open(os.path.join(this_directory, filename), "rb") as f:
        return f.read().decode("utf-8")


if os.path.exists("README.md"):
    long_description = read("README.md")

def read(*parts):                                                                                                  
    with codecs.open(os.path.join(this_directory, *parts), 'r') as fp:                                             
        return fp.read()                                                                                           
                                                                                                                   
def find_version(*file_paths):                                                                                     
    version_file = read(*file_paths)                                                                               
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]",                                               
                              version_file, re.M)                                                                  
    if version_match:                                                                                              
        return version_match.group(1)                                                                              
    raise RuntimeError("Unable to find version string.")


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
    long_description_content_type='text/markdown',
    author=author,
    url=url,
    author_email=author_email,
    classifiers=classifiers,
    package_data={"prosail": ["*.txt"]},
    install_requires=[
        "numpy",
        "numba",
        "scipy",
        "pytest",
    ],
    version=find_version("prosail", "__init__.py"),
    packages=["prosail"],
    zip_safe=False # Apparently needed for conda
)
