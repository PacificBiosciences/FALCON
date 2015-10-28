#!/usr/bin/env python

from setuptools import setup, Extension

import glob

#install_requires=[ "pbcore >= 0.6.3", "networkx >= 1.7" ]
install_requires=[ "networkx >= 1.7" ]

scripts = glob.glob("src/py_scripts/*.py")

setup(name='falcon_kit',
      version='0.4.0',
      description='a small toolkit for DNA seqeucne alignment, overlapping, and assembly',
      author='Jason Chin',
      author_email='jchin@pacificbiosciences.com',
      packages=['falcon_kit',
          'falcon_kit.mains',
          'falcon_kit.util',
          ],
      package_dir={'falcon_kit':'src/py/'},
      ext_modules=[
                   Extension('falcon_kit.ext_falcon', ['src/c/ext_falcon.c', 'src/c/DW_banded.c', 'src/c/kmer_lookup.c', 'src/c/falcon.c'],
                    extra_link_args=[],
                    extra_compile_args=['-fPIC', '-O3', '-fno-omit-frame-pointer'],
                    # '-fno-omit-frame-pointer' can help with gperftools.
                    #libraries=['profiler'],
                    #include_dirs=['/home/cdunn/local/include'],
                    #library_dirs=['/home/cdunn/local/lib'],
                    #language="c++", # c for now
                    #export_symbols=['generate_consensus'], # for windows?
                   ),
                  ],
      entry_points = {'console_scripts': [
          'falcon-task=falcon_kit.mains.tasks:main',
          'fc_actg_coordinate=falcon_kit.mains.actg_coordinate:main',
          'fc_consensus=falcon_kit.mains.consensus:main',
          'fc_contig_annotate=falcon_kit.mains.contig_annotate:main',
          'fc_ctg_link_analysis=falcon_kit.mains.ctg_link_analysis:main',
          'fc_dedup_a_tigs=falcon_kit.mains.dedup_a_tigs:main',
          'fc_graph_to_contig=falcon_kit.mains.graph_to_contig:main',
          'fc_graph_to_utgs=falcon_kit.mains.graph_to_utgs:main',
          'fc_ovlp_filter=falcon_kit.mains.ovlp_filter:main',
          'fc_ovlp_stats=falcon_kit.mains.ovlp_stats:main',
          'fc_ovlp_to_graph=falcon_kit.mains.ovlp_to_graph:main',
          'fc_run=falcon_kit.mains.run:main',
          ],
      },
      extras_require = {
          'falcon-task':  ['falcon_kit'],
      },
      scripts = scripts,
      zip_safe = False,
      setup_requires=install_requires,
      install_requires=install_requires
     )
