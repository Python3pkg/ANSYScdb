# -*- coding: utf-8 -*-
"""
Setup.py for ANSYScdb

"""
from setuptools import setup, Extension
from Cython.Distutils import build_ext

import numpy


setup(
    name='ANSYScdb',
    packages = ['ANSYScdb', 'ANSYScdb.Tests'],

    # Version
    version='0.10.3',

    description='Loads ANSYS cdb files',
    long_description=open('README.rst').read(),

    # Author details
    author='Alex Kaszynski',
    author_email='akascap@gmail.com',

    license='MIT',
    classifiers=[
        'Development Status :: 4 - Beta',

        # Target audience
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Information Analysis',

        # MIT License
        'License :: OSI Approved :: MIT License',

        # Tested only on Python 2.7 (untested with 3)
        'Programming Language :: Python :: 2.7',
    ],

    # Website
    url = 'https://github.com/akaszynski/ANSYScdb',

    # Build cython modules
    cmdclass={'build_ext': build_ext},
    ext_modules=[Extension("ANSYScdb.CDBparser", 
                           ["ANSYScdb/cython/CDBparser.pyx"],
                           language='c'),

                 Extension("ANSYScdb.CDBreader", 
                           ["ANSYScdb/cython/CDBreader.pyx"],
                           language='c'),
                ],
                           
    keywords='vtk ANSYS cdb',                           
                           
    include_dirs=[numpy.get_include()],
                  
    package_data={'ANSYScdb.Tests': ['HexBeam.cdb', 'TetBeam.cdb']},

    # Might work with earlier versions (untested)
    install_requires=['numpy>1.9.3', 'cython>0.23.1']

)

