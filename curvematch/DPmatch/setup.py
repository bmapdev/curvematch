#!/usr/bin/env python
"""Setup file"""

__author__ = "Shantanu H. Joshi"
__copyright__ = "Copyright 2013, Shantanu H. Joshi, Brandon Ayers, \
                 Ahmanson-Lovelace Brain Mapping Center, University of California Los Angeles"
__email__ = "s.joshi@g.ucla.edu"
try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup
from distutils.extension import Extension
import numpy as np

USE_CYTHON = False

if USE_CYTHON:
    try:
        from Cython.Distutils import build_ext
    except ImportError:
        if USE_CYTHON == 'auto':
            USE_CYTHON = False
        else:
            raise

cmdclass = {}
ext_modules = []
base_dir = '.'

if USE_CYTHON:
    ext_modules += \
        [
        Extension("DPmatchcy", sources=["DPmatchcy.pyx",
                                                           "DPmatch.cpp",
                                                           "shape.cpp"],
        include_dirs = [np.get_include(), '.'], language='c++', pyrex_gdb=True)
        ]
    cmdclass.update({'build_ext': build_ext })
else:
    ext_modules += \
        [
         Extension("DPmatchcy", sources=["DPmatchcy.cpp",
                                                            "DPmatch.cpp",
                                                            "shape.cpp"],
         include_dirs = [np.get_include(), '.'], language='c++')]

setup(
    name='DPmatchcy',
    version='0.1_dev',
    cmdclass=cmdclass,
    ext_modules=ext_modules,
    packages=['DPmatchcy'],
    package_dir={
        'DPmatchcy': base_dir,
        },
    test_suite='nose.collector',
    url='',
    license='TBD',
    author='Shantanu H. Joshi, Brandon Ayers',
    author_email='s.joshi@ucla.edu, ayersb@ucla.edu',
    description='Cython: Dynamic programming for n-dimensional functions',
    requires=['numpy (>=1.3.0)','cython (>=0.15.1)'],
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'License :: TBD',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.4',
        'Programming Language :: Python :: 2.5',
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Cython',
        'Topic :: Scientific/Engineering :: Mathematics',
    ],
    keywords='dynamic-programming optimization',
    )
