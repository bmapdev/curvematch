from distutils.core import setup
from distutils.extension import Extension
setup(
    name='DPmatchcy',
    version='0.1_dev',
    cmdclass=cmdclass,
    ext_modules=ext_modules,
    packages=['test'],
    package_dir={
        'test': base_dir,
        },
    package_data={'test': ['/data/*.ucf', '/data/group/*']},
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
