#!/usr/bin/env python2.7

from setuptools import setup, Extension
import subprocess
import sys
import glob

install_requires = [
    "networkx >=1.7, <=1.11",
    "msgpack",
    "future >= 0.16.0",
    #"logging_tree",
    #"pbcore >= 0.6.3",
    "pypeFLOW >= 2.0.0",
]

scripts = []

try:
    local_version = '+git.{}'.format(
        subprocess.check_output('git rev-parse HEAD', shell=True))
except Exception:
    local_version = ''

setup(name='falcon_kit',
      version='1.1.1' + local_version,
      description='a small toolkit for DNA seqeucne alignment, overlapping, and assembly',
      author='Jason Chin',
      author_email='jchin@pacificbiosciences.com',
      packages=['falcon_kit',
                'falcon_kit.mains',
                'falcon_kit.util',
                ],
      package_dir={'falcon_kit': 'falcon_kit/'},
      ext_modules=[
          Extension('ext_falcon', ['src/c/ext_falcon.c', 'src/c/DW_banded.c', 'src/c/kmer_lookup.c', 'src/c/falcon.c'],
                    extra_link_args=[],
                    extra_compile_args=['-fPIC', '-O3',
                                        '-fno-omit-frame-pointer'],
                    # '-fno-omit-frame-pointer' can help with gperftools.
                    # libraries=['profiler'],
                    # include_dirs=['/home/cdunn/local/include'],
                    # library_dirs=['/home/cdunn/local/lib'],
                    # language="c++", # c for now
                    # export_symbols=['generate_consensus'], # for windows?
                    ),
      ],
      entry_points={'console_scripts': [
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
          'fc_calc_cutoff=falcon_kit.mains.calc_cutoff:main',
          'fc_run=falcon_kit.mains.run1:main',
          'fc_run1=falcon_kit.mains.run1:main',
          'fc_fasta2fasta=falcon_kit.mains.fasta2fasta:main',
          'fc_fetch_reads=falcon_kit.mains.fetch_reads:main',
          'fc_get_read_ctg_map=falcon_kit.mains.get_read_ctg_map:main',
          'fc_pr_ctg_track=falcon_kit.mains.pr_ctg_track:main',
          'fc_rr_ctg_track=falcon_kit.mains.rr_ctg_track:main',
          'fc_gen_gfa_v1 = falcon_kit.mains.gen_gfa_v1:main',

          # Some scripts still use old .py scripts, which need shebang lines.
          # (I guess pip does not add shebang to scripts.)
          'fc_run.py=falcon_kit.mains.run1:main',
      ],
      },
      extras_require={
          'falcon-task':  ['falcon_kit'],
      },
      scripts=scripts,
      zip_safe=False,
      install_requires=install_requires,
)
