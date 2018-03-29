from __future__ import absolute_import

import argparse
import collections
import logging
import os
import sys
from .. import io
from .. import bash
from .. import run_support
from .. import pype_tasks

LOG = logging.getLogger()

def run(bash_template_fn, db_prefix, pread_aln, skip_checks, run_jobs_fn, db_fn, wildcards, nproc, split_fn):
    with open(bash_template_fn, 'w') as stream:
        stream.write(pype_tasks.TASK_DALIGNER_SCRIPT)

    nblock = run_support.get_nblock(db_fn)

    # TODO: nproc can restrict -T => 1 or 2 (never 3)
    if 4 != nproc:
        LOG.warning('Currently ignoring nproc={} in daligner_split.'.format(nproc))

    db_build_done_fn = None
    assert isinstance(pread_aln, int)
    daligner_scripts = bash.scripts_daligner(run_jobs_fn, db_prefix, db_build_done_fn, nblock=nblock, pread_aln=bool(pread_aln), skip_check=skip_checks)
    #basedir = os.path.dirname(os.path.abspath(split_fn))
    #rootdir = os.path.dirname(os.path.dirname(basedir)) # for now
    cwd = os.getcwd()
    jobs = list()
    for job_uid, script in daligner_scripts:
        job_id = 'j_{}'.format(job_uid)
        daligner_settings = dict(db_prefix=db_prefix)

        # Write script/settings for a unit-of-work.
        daligner_script_fn = '{cwd}/daligner-scripts/{job_id}/daligner-script.sh'.format(**locals())
        daligner_settings_fn = '{cwd}/daligner-scripts/{job_id}/settings.json'.format(**locals())
        io.mkdirs(os.path.dirname(daligner_script_fn))
        with open(daligner_script_fn, 'w') as stream:
            stream.write(script)
        io.serialize(daligner_settings_fn, daligner_settings)

        # Record in a job-dict.
        job = dict()
        job['input'] = dict(
                daligner_script = daligner_script_fn,
                daligner_settings = daligner_settings_fn, # not used today, but maybe someday
        )
        job['output'] = dict(
                job_done = 'daligner.done'
        )
        job['params'] = dict(
        )
        job['wildcards'] = {wildcards: job_id}
        jobs.append(job)

    io.serialize(split_fn, jobs)

class HelpF(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass


def parse_args(argv):
    description = 'Split HPC.daligner output into multiple daligner jobs (possible in parallel).'
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
        '--db-prefix', default='raw_reads',
        help='Either preads or raw_reads.')
    parser.add_argument(
        '--skip-checks', default=0, type=int,
        help='Skip LAcheck calls after daligner. (0 => do not skip)')
    parser.add_argument(
        '--db-fn',
        help='Dazzler DB of reads. (Used to calculate number of blocks.)')
    parser.add_argument(
        '--pread-aln', default=0, type=int,
        help='If non-zero, use pread alignment mode. (Run daligner_p instead of daligner.)')
    parser.add_argument(
        '--wildcards',
        help='Input. Comma-separated wildcard names. Might be needed downstream.')
    parser.add_argument(
        '--nproc', type=int, default=4,
        help='(Ignored for now.) Number of processors available to this run. Should restrict -T if < 4.',
    )
    parser.add_argument(
        '--split-fn',
        help='Output. JSON list of units of work.')
    parser.add_argument(
        '--bash-template-fn',
        help='Output. Copy of known daligner bash template, for use later.')
    args = parser.parse_args(argv[1:])
    return args


def main(argv=sys.argv):
    args = parse_args(argv)
    logging.basicConfig(level=logging.INFO)
    run(**vars(args))


if __name__ == '__main__':  # pragma: no cover
    main()
