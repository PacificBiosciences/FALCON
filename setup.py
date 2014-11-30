#!/usr/bin/env python

from setuptools import setup

from distutils.core import Extension

import glob

scripts = ["src/py_scripts/falcon_asm.py", 
         "src/py_scripts/falcon_overlap.py",
         "src/py_scripts/falcon_overlap2.py",
         "src/py_scripts/falcon_qrm.py",
         "src/py_scripts/falcon_fixasm.py",
         "src/py_scripts/falcon_dedup.py",
         "src/py_scripts/falcon_sense.py",
         "src/py_scripts/falcon_ucns_data.py",
         "src/py_scripts/falcon_utgcns.py",
         "src/py_scripts/get_rdata.py",
         "src/py_scripts/remove_dup_ctg.py"]

scripts_v2 = glob.glob("src/py_scripts_v2/*.py")

#install_requires=[ "pbcore >= 0.6.3", "networkx >= 1.7" ]
install_requires=[ "networkx >= 1.7" ]

scripts.extend( scripts_v2 )

setup(name='falcon_kit',
      version='0.1.3',
      description='a small toolkit for DNA seqeucne alignment, overlapping, and assembly',
      author='Jason Chin',
      author_email='jchin@pacificbiosciences.com',
      packages=['falcon_kit'],
      package_dir={'falcon_kit':'src/py/'},
      ext_modules=[Extension('falcon_kit.DW_align', ['src/c/DW_banded.c'], 
                   extra_link_args=["-fPIC",  "-O3"]),
                   Extension('falcon_kit.kmer_lookup', ['src/c/kmer_lookup.c'],
                   extra_link_args=["-fPIC",  "-O3"]),
                   Extension('falcon_kit.falcon', ['src/c/DW_banded.c', 'src/c/kmer_lookup.c', 'src/c/falcon.c'],
                   extra_link_args=["-fPIC",  "-O3"]),
                   ],
      scripts = scripts,
      zip_safe = False,
      install_requires=install_requires
     )

