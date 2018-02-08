from __future__ import unicode_literals
import argparse
import logging
import os
import sys
from .. import io
from .. import bash
from .. import run_support

LOG = logging.getLogger()


def run(config_fn, length_cutoff_fn, raw_reads_db_fn, preads_fofn_fn, pre_assembly_report_fn):
    config = io.deserialize(config_fn)
    genome_length = int(config['genome_size'])
    length_cutoff_user = int(config['length_cutoff'])
    # Update length_cutoff if auto-calc (when length_cutoff is negative).
    # length_cutoff_fn was created long ago, so no filesystem issues.
    length_cutoff = run_support.get_length_cutoff(
        length_cutoff_user, length_cutoff_fn)
    # Hmmm. Actually, I think we now write the user length_cutoff into the length_cutoff file,
    # if not -1. TODO(CD): Check on that, and simplify here if so.

    script = bash.script_run_report_pre_assembly(
        raw_reads_db_fn, preads_fofn_fn, genome_length, length_cutoff, pre_assembly_report_fn)
    script_fn = 'run-report-pre-assembly.sh'
    job_done_fn = 'job.done'
    bash.write_script(script, script_fn, job_done_fn)
    io.syscall('bash -vex {}'.format(script_fn))


class HelpF(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass


def parse_args(argv):
    description = 'Prepare to run the pre-assembly report generator, and run it.'
    epilog = 'length_cutoff might be cleaned up someday. For now, yeah, it is confusing.'
    parser = argparse.ArgumentParser(
        description=description,
        epilog=epilog,
        formatter_class=HelpF,
    )
    parser.add_argument(
        '--config-fn',
        help='Input. JSON configuration. We use "length_cutoff" (if positive) and "genome_size".',
    )
    parser.add_argument(
        '--length-cutoff-fn',
        help='Input. File of a single number: the length-cutoff for raw reads.',
    )
    parser.add_argument(
        '--raw-reads-db-fn',
        help='Input. Dazzler DB of raw reads.',
    )
    parser.add_argument(
        '--preads-fofn-fn',
        help='Input. FOFN of preads las files.',
    )
    parser.add_argument(
        '--pre-assembly-report-fn',
        help='Output. In JSON format.',
    )
    args = parser.parse_args(argv[1:])
    return args


def main(argv=sys.argv):
    args = parse_args(argv)
    logging.basicConfig(level=logging.INFO)
    run(**vars(args))


if __name__ == '__main__':  # pragma: no cover
    main()
