from __future__ import absolute_import
from .. import pype_tasks
# pylint: disable=no-name-in-module, import-error, fixme, line-too-long
from pypeflow.simple_pwatcher_bridge import (PypeProcWatcherWorkflow, MyFakePypeThreadTaskBase,
                                             makePypeLocalFile, fn, PypeTask)
PypeThreadTaskBase = MyFakePypeThreadTaskBase
import argparse
import glob
import logging
import sys
import subprocess as sp
import shlex
import os

LOG = logging.getLogger(__name__)


def make_dirs(d):
    if not os.path.isdir(d):
        LOG.debug('mkdirs {}'.format(d))
        os.makedirs(d)


def get_read_ctg_map(rawread_dir, pread_dir, asm_dir):
    read_map_dir = os.path.abspath(os.path.join(asm_dir, 'read_maps'))
    make_dirs(read_map_dir)

    wf = PypeProcWatcherWorkflow(
        max_jobs=12,
    )
    """
            job_type=config['job_type'],
            job_queue=config['job_queue'],
            sge_option=config.get('sge_option', ''),
            watcher_type=config['pwatcher_type'],
            watcher_directory=config['pwatcher_directory'])
    """

    rawread_db = makePypeLocalFile(os.path.join(rawread_dir, 'raw_reads.db'))
    rawread_id_file = makePypeLocalFile(os.path.join(
        read_map_dir, 'dump_rawread_ids', 'rawread_ids'))

    task = PypeTask(
        inputs={'rawread_db': rawread_db},
        outputs={'rawread_id_file': rawread_id_file},
    )
    wf.addTask(task(pype_tasks.task_dump_rawread_ids))

    pread_db = makePypeLocalFile(os.path.join(pread_dir, 'preads.db'))
    pread_id_file = makePypeLocalFile(os.path.join(
        read_map_dir, 'dump_pread_ids', 'pread_ids'))

    task = PypeTask(
        inputs={'pread_db': pread_db},
        outputs={'pread_id_file': pread_id_file},
    )
    wf.addTask(task(pype_tasks.task_dump_pread_ids))

    wf.refreshTargets()  # block

    sg_edges_list = makePypeLocalFile(os.path.join(asm_dir, 'sg_edges_list'))
    utg_data = makePypeLocalFile(os.path.join(asm_dir, 'utg_data'))
    ctg_paths = makePypeLocalFile(os.path.join(asm_dir, 'ctg_paths'))

    inputs = {'rawread_id_file': rawread_id_file,
              'pread_id_file': pread_id_file,
              'sg_edges_list': sg_edges_list,
              'utg_data': utg_data,
              'ctg_paths': ctg_paths}

    read_to_contig_map = makePypeLocalFile(os.path.join(
        read_map_dir, 'get_ctg_read_map', 'read_to_contig_map'))

    task = PypeTask(
        inputs=inputs,
        outputs={'read_to_contig_map': read_to_contig_map},
    )
    wf.addTask(task(pype_tasks.task_generate_read_to_ctg_map))

    wf.refreshTargets()  # block


def parse_args(argv):
    parser = argparse.ArgumentParser(description='generate `2-asm-falcon/read_maps/read_to_contig_map` that contains the \
information from the chain of mapping: (contig id) -> (internal p-read id) -> (internal raw-read id) -> (original read id)',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--basedir', type=str, default='./',
                        help='the base working dir of a FALCON assembly')
    args = parser.parse_args(argv[1:])
    return args


def main(argv=sys.argv):
    logging.basicConfig()
    args = parse_args(argv)
    basedir = args.basedir
    rawread_dir = os.path.abspath(os.path.join(basedir, '0-rawreads'))
    pread_dir = os.path.abspath(os.path.join(basedir, '1-preads_ovl'))
    asm_dir = os.path.abspath(os.path.join(basedir, '2-asm-falcon'))

    get_read_ctg_map(rawread_dir=rawread_dir,
                     pread_dir=pread_dir, asm_dir=asm_dir)


if __name__ == '__main__':
    main()
