from __future__ import absolute_import
from __future__ import unicode_literals
import argparse
import logging
import os
import sys
from .. import io
from .. import bash  # for write_script

LOG = logging.getLogger()

def run(config_fn, input_fofn_fn, job_done_fn, run_jobs_fn):
    LOG.info('Building pdb from {!r} ({!r}), to write {!r}'.format(
        input_fofn_fn, config_fn, run_jobs_fn))

    config = io.deserialize(config_fn)
    run_jobs_fn = os.path.basename(run_jobs_fn)
    script = bash.script_build_pdb(config, input_fofn_fn, run_jobs_fn)
    script_fn = 'build_pdb.sh'
    bash.write_script(script, script_fn, job_done_fn)
    io.syscall('bash -vex {}'.format(script_fn))


class HelpF(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass


def parse_args(argv):
    description = 'Build Dazzler DB for preads, and run HPC.daligner to generate run-script.'
    epilog = ''
    parser = argparse.ArgumentParser(
        description=description,
        epilog=epilog,
        formatter_class=HelpF,
    )
    parser.add_argument(
        '--input-fofn-fn',
        help='Input. (Basename only.) User-provided file of fasta filename.',
    )
    parser.add_argument(
        '--config-fn',
        help='Input. JSON of relevant configuration (currently from General section of full-prog config).',
    )
    parser.add_argument(
        '--job-done-fn',
        help='Output. (Basename only.) Sentinel (to be dropped someday).',
    )
    parser.add_argument(
        '--run-jobs-fn',
        help='Output. (Basename only.) Produced by HPC.daligner.',
    )
    args = parser.parse_args(argv[1:])
    return args


def main(argv=sys.argv):
    args = parse_args(argv)
    logging.basicConfig(level=logging.INFO)
    run(**vars(args))


if __name__ == '__main__':  # pragma: no cover
    main()
