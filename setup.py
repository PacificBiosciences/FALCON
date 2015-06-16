#!/usr/bin/env python

#from setuptools import setup

from distutils.core import Extension, setup
#from distutils.command.build_clib import build_clib

import glob

#libfalcon = ('falcon', {'sources': ['src/c/DW_banded.c', 'src/c/kmer_lookup.c', 'src/c/falcon.c']})

#install_requires=[ "pbcore >= 0.6.3", "networkx >= 1.7" ]
install_requires=[ "networkx >= 1.7" ]

scripts = glob.glob("src/py_scripts/*.py")

setup(name='falcon_kit',
      version='0.3.0',
      description='a small toolkit for DNA seqeucne alignment, overlapping, and assembly',
      author='Jason Chin',
      author_email='jchin@pacificbiosciences.com',
      packages=['falcon_kit'],
      package_dir={'falcon_kit':'src/py/'},
      #libraries=[libfalcon],
      #cmdclass = {'build_clib': build_clib},
      ext_modules=[
                   Extension('falcon_kit.ext_falcon', ['src/c/ext_falcon.c', 'src/c/DW_banded.c', 'src/c/kmer_lookup.c', 'src/c/falcon.c'],
                    extra_link_args=["-fPIC",  "-O3"]),
                   ],
      scripts = scripts,
      zip_safe = False,
      install_requires=install_requires
     )

