from __future__ import absolute_import
from __future__ import unicode_literals
import argparse
import logging
import os
import sys
from .. import io
from .. import bash

LOG = logging.getLogger()

# This function was copied from bash.py and modified.
def script_run_consensus(config, db_fn, las_fn, out_file_fn):
    """config: dazcon, falcon_sense_greedy, falcon_sense_skip_contained, LA4Falcon_preload
    """
    io.rm(out_file_fn) # in case of resume
    out_file_bfn = out_file_fn + '.tmp'
    params = dict(config)
    length_cutoff = params['length_cutoff']
    bash_cutoff = '{}'.format(length_cutoff)
    params.update(locals())
    LA4Falcon_flags = 'P' if params.get('LA4Falcon_preload') else ''
    if config["falcon_sense_skip_contained"]:
        LA4Falcon_flags += 'fso'
    elif config["falcon_sense_greedy"]:
        LA4Falcon_flags += 'fog'
    else:
        LA4Falcon_flags += 'fo'
    if LA4Falcon_flags:
        LA4Falcon_flags = '-' + ''.join(set(LA4Falcon_flags))
    run_consensus = "LA4Falcon -H$CUTOFF %s {db_fn} {las_fn} | python -m falcon_kit.mains.consensus {falcon_sense_option} >| {out_file_bfn}" % LA4Falcon_flags

    if config.get('dazcon', False):
        run_consensus = """
which dazcon
dazcon {pa_dazcon_option} -s {db_fn} -a {las_fn} >| {out_file_bfn}
"""

    script = """
set -o pipefail
CUTOFF=%(bash_cutoff)s
%(run_consensus)s
mv -f {out_file_bfn} {out_file_fn}
""" % (locals())
    return script.format(**params)


def run(config_fn, length_cutoff_fn, las_fn, db_fn, fasta_fn):
    job_done_fn = 'job.done'
    length_cutoff = int(open(length_cutoff_fn).read())
    config = io.deserialize(config_fn)
    config['length_cutoff'] = length_cutoff
    script = script_run_consensus(
        config, db_fn, las_fn, os.path.basename(fasta_fn)) # not sure basename is really needed here
    script_fn = 'run_consensus.sh'
    bash.write_script(script, script_fn, job_done_fn)
    io.syscall('bash -vex {}'.format(script_fn))


class HelpF(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass


def parse_args(argv):
    description = 'Run consensus on a merged .las file, to produce a fasta file of preads.'
    epilog = ''
    parser = argparse.ArgumentParser(
        description=description,
        epilog=epilog,
        formatter_class=HelpF,
    )
    parser.add_argument(
        '--las-fn',
        help='Input. Merged .las file.',
    )
    parser.add_argument(
        '--db-fn',
        help='Input. Dazzler DB of raw-reads.',
    )
    parser.add_argument(
        '--length-cutoff-fn',
        help='Input. Contains a single integer, the length-cutoff.',
    )
    parser.add_argument(
        '--config-fn',
        help='Input. JSON of relevant configuration (currently from General section of full-prog config).',
    )
    parser.add_argument(
        '--fasta-fn',
        help='Output. Consensus fasta file.',
    )
    args = parser.parse_args(argv[1:])
    return args


def main(argv=sys.argv):
    args = parse_args(argv)
    logging.basicConfig(level=logging.INFO)
    run(**vars(args))


if __name__ == '__main__':  # pragma: no cover
    main()
