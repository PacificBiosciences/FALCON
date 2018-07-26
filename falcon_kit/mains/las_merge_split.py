from __future__ import absolute_import

import argparse
import collections
import logging
import os
import re
import sys
from .. import io
from .. import bash
from .. import pype_tasks

LOG = logging.getLogger()


# Could be L1.* or preads.* or raw_reads.*
re_las = re.compile(r'\.(\d*)(?:\.[^\.]+\.\d*)?\.las$') # see daligner_gather_las() in run_support.py

def get_block(las_fn):
    """
    >>> get_block('foo.42.las')
    42
    >>> get_block('foo.12.foo.34.las')
    12
    """
    mo = re_las.search(las_fn)
    if not mo:
        msg = 'No match of las file {!r} with {}'.format(las_fn, re_las.pattern)
        raise Exception(msg)
    block = int(mo.group(1))
    return block

def read_gathered_las(las_fns):
    """Return dict of block->[las_paths].
    The input is one .las file per line.
    """
    result = collections.defaultdict(list)
    for las_fn in las_fns:
        block = get_block(las_fn)
        result[block].append(las_fn)
    return result


def run(run_jobs_fn, gathered_las_fn, db_prefix, wildcards,
        split_fn, bash_template_fn):
    with open(bash_template_fn, 'w') as stream:
        stream.write(pype_tasks.TASK_LAS_MERGE_SCRIPT)
    LOG.info('Splitting las file-list from {!r} (based on {!r}) into {!r}.'.format(
        gathered_las_fn, run_jobs_fn, split_fn))
    config = dict() # not used anyway
    merge_scripts = list(bash.scripts_merge(config, db_prefix, run_jobs_fn))
    gathered_las = io.deserialize(gathered_las_fn)
    gathered_dict = read_gathered_las(gathered_las)
    LOG.info('Gathered {} las files.'.format(len(gathered_dict)))
    gathered_dict_dir = os.path.normpath(os.path.dirname(gathered_las_fn))

    if 0 == len(gathered_dict):
        raise Exception('No .las files. There must have been some kind of error.\n gathered_las_fn={!r},\n run_jobs_fn={!r}'.format(
            gathered_las_fn, run_jobs_fn))
    if 1 == len(gathered_dict) and 1 == len(gathered_dict.values()[0]):
        LOG.warning('There is only 1 .las file in {!r}.\n The las-merge step will be a no-op.'.format(gathered_las_fn))
        assert len(merge_scripts) == 0, 'Expected more than 1 .las file. merge_scripts={}. See {!r} and {!r}'.format(
                merge_scripts, run_jobs_fn, gathered_las_fn)
        solo_block, solo_las_fns = gathered_dict.items()[0]
        merge_scripts = [(solo_block, '# no-op', os.path.basename(solo_las_fns[0]))]

    cwd = os.getcwd()
    jobs = list()
    for p_id, merge_script, merged_las_fn in merge_scripts:
        job_id = 'm_%05d' %p_id
        # Write the split inputs.
        merge_script_fn = '{cwd}/las-merge-scripts/{job_id}/merge-script.sh'.format(**locals())
        las_paths_fn = '{cwd}/las-merge-scripts/{job_id}/las_paths.json'.format(**locals())
        merged_las_json_fn = '{cwd}/las-merge-scripts/{job_id}/merged_las.json'.format(**locals())
        io.mkdirs(os.path.dirname(merge_script_fn), os.path.dirname(las_paths_fn), os.path.dirname(merged_las_json_fn))
        with open(merge_script_fn, 'w') as stream:
            stream.write(merge_script)
        # las_paths must be updated for new relative paths
        las_paths_dir = os.path.normpath(os.path.dirname(las_paths_fn))
        las_paths = list()
        for p in gathered_dict[p_id]:
            updated_relative_path = os.path.relpath(os.path.join(gathered_dict_dir, p), las_paths_dir)
            las_paths.append(updated_relative_path)
        io.serialize(las_paths_fn, las_paths)
        io.serialize(merged_las_json_fn, merged_las_fn)

        # Record in a job-dict.
        job = dict()
        job['input'] = dict(
                las_paths = las_paths_fn,
                merge_script = merge_script_fn,
                merged_las_json = merged_las_json_fn,
        )
        job['output'] = dict(
                merged_las = 'merged.las',
                job_done = 'merge.done',
                p_id = 'p_id'
        )
        job['params'] = dict(
                actual_merged_las = merged_las_fn,
                p_id_num = p_id,
        )
        job['wildcards'] = {wildcards: job_id}
        jobs.append(job)

    io.serialize(split_fn, jobs)


class HelpF(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass


def parse_args(argv):
    description = 'Prepare for parallel las-merge jobs.'
    epilog = ''
    parser = argparse.ArgumentParser(
        description=description,
        epilog=epilog,
        formatter_class=HelpF,
    )
    parser.add_argument(
        '--run-jobs-fn',
        help='Input. Result of HPC.daligner.')
    parser.add_argument(
        '--gathered-las-fn',
        help='Input. (Not sure of content yet.)')
    parser.add_argument(
        '--db-prefix', default='raw_reads',
        help='Either preads or raw_reads.')
    parser.add_argument(
        '--wildcards',
        help='Input. Comma-separated wildcard names. Might be needed downstream.')
    parser.add_argument(
        '--split-fn',
        help='Output. JSON list of jobs, where each is a dict of input/output/params/wildcards.')
    parser.add_argument(
        '--bash-template-fn',
        help='Output. Known template of bash for running las-merge.')
    args = parser.parse_args(argv[1:])
    return args


def main(argv=sys.argv):
    args = parse_args(argv)
    logging.basicConfig(level=logging.INFO)
    run(**vars(args))


if __name__ == '__main__':  # pragma: no cover
    main()
