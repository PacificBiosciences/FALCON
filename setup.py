#!/usr/bin/env python

from setuptools import setup

from distutils.core import Extension

setup(name='falcon_kit',
      version='0.1',
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
      scripts = ["src/py_scripts/falcon_asm.py", 
                 "src/py_scripts/falcon_asm_dev.py",
                 "src/py_scripts/falcon_overlap.py",
                 "src/py_scripts/falcon_overlap2.py",
                 "src/py_scripts/falcon_fixasm.py",
                 "src/py_scripts/falcon_wrap.py",
                 "src/py_scripts/get_rdata.py",
                 "src/py_scripts/remove_dup_ctg.py"],
      zip_safe = False,
      install_requires=[ "pbcore >= 0.6.3", "networkx >= 1.7" ]
     )

