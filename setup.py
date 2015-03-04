#!/usr/bin/env python
"""Setup file"""

__author__ = "Shantanu H. Joshi"
__copyright__ = "Copyright 2013, Shantanu H. Joshi, Brandon Ayers, \
                 Ahmanson-Lovelace Brain Mapping Center, University of California Los Angeles"
__email__ = "s.joshi@ucla.edu"

from distutils.core import setup
from distutils.extension import Extension
import numpy as np

USE_CYTHON = True

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
base_dir = 'curvematch'

if USE_CYTHON:
    ext_modules += \
        [
        Extension("curvematch.DPmatch.DPmatchcy", sources=["curvematch/DPmatch/DPmatchcy.pyx",
                                                           "curvematch/DPmatch/DPmatch.cpp",
                                                           "curvematch/DPmatch/shape.cpp"],
                  include_dirs=[np.get_include(), '.'], language='c++', pyrex_gdb=True)
        ]
    cmdclass.update({'build_ext': build_ext })
else:
    ext_modules += \
        [
         Extension("curvematch.DPmatch.DPmatchcy", sources=["curvematch/DPmatch/DPmatchcy.cpp",
                                                            "curvematch/DPmatch/DPmatch.cpp",
                                                            "curvematch/DPmatch/shape.cpp"],
                   include_dirs=[np.get_include(), '.'], language='c++')]


setup(
    name='curvematch',
    version='0.1_dev',
    packages=['curvematch', 'curvematch.test', 'curvematch.DPmatch'],
    package_dir={
            'curvematch': base_dir,
            'curvematch.DPmatch': base_dir + '/DPmatch',
            'curvematch.test': base_dir + '/test'
            },
    package_data={'curvematch':['test/data/*.ucf', 'test/data/group/*']},
    scripts=['bin/curve_match.py',
             'bin/curve_match_group.py',
             ],
    test_suite='nose.collector',
    url='',
    license='TBD',
    author='Shantanu H. Joshi, Brandon Ayers',
    author_email='s.joshi@ucla.edu, ayersb@ucla.edu',
    description='Shape Matching and Analysis of Curves',
    requires=['numpy (>=1.3.0)'],
    cmdclass=cmdclass,
    ext_modules=ext_modules,
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
    keywords='shape curves matching dynamic programming',
    )
