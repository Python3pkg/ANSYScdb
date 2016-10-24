"""
Setup.py for ANSYScdb
"""
#import os
from setuptools import setup, Extension
from Cython.Distutils import build_ext

import numpy
import sys

# Check compiler and assign compile arguments accordingly
def compilerName():
  import re
  import distutils.ccompiler
  comp = distutils.ccompiler.get_default_compiler()
  getnext = False

  for a in sys.argv[2:]:
    if getnext:
      comp = a
      getnext = False
      continue
    #separated by space
    if a == '--compiler'  or  re.search('^-[a-z]*c$', a):
      getnext = True
      continue
    #without space
    m = re.search('^--compiler=(.+)', a)
    if m == None:
      m = re.search('^-[a-z]*c(.+)', a)
    if m:
      comp = m.group(1)

  return comp

# Assign arguments based on compiler
compiler = compilerName()
if compiler == 'unix' or compiler == 'msvc':
    cmp_arg = ['-O3']
else:
    cmp_arg = ['/Ox']


setup(
    name='ANSYScdb',
    packages = ['ANSYScdb', 'ANSYScdb.Tests'],

    # Version
    version='0.13.1',

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

        # Tested with on Python 2.7 and 3.5
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.5',
    ],

    # Website
    url = 'https://github.com/akaszynski/ANSYScdb',

    # Build cython modules
    cmdclass={'build_ext': build_ext},
    ext_modules=[Extension("ANSYScdb.CDBparser", 
                           ["ANSYScdb/cython/CDBparser.pyx"],
                           language='c'),

                 Extension('ANSYScdb._reader', 
                           ['ANSYScdb/cython/_reader.pyx',
                            'ANSYScdb/cython/reader.c'],
                           extra_compile_args=cmp_arg,
                           language='c',),

                ],
                           
    keywords='vtk ANSYS cdb',                           
                           
    include_dirs=[numpy.get_include()],
                  
    package_data={'ANSYScdb.Tests': ['HexBeam.cdb', 'TetBeam.cdb']},

    # Might work with earlier versions (untested)
    install_requires=['numpy>1.9.3', 'cython>0.23.1']

)

