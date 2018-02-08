"""
"""
from __future__ import print_function
import argparse
import logging
import os
import sys
from .. import io

LOG = logging.getLogger()

def convert_job_id_to_p_id(job_id):
    """
    >>> convert_job_id_to_p_id('cns_0011')
    11
    """
    return int(job_id[4:], base=10)

def run(gathered_fn, preads_fofn_fn):
    gathered = io.deserialize(gathered_fn)
    fasta_fns = dict()
    for key, desc in gathered.items():
        job_id = key.split(',')[0].split('=')[1]
        p_id = convert_job_id_to_p_id(job_id)
        fasta_fns[p_id] = desc['fns']['fasta']
    with open(preads_fofn_fn,  'w') as f:
        for filename in sorted(fasta_fns.values()):
            print(filename, file=f)


class HelpF(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass


def parse_args(argv):
    description = 'Turn gathered file into FOFN of fasta files.'
    epilog = ''
    parser = argparse.ArgumentParser(
        description=description,
        epilog=epilog,
        formatter_class=HelpF,
    )
    parser.add_argument(
        '--gathered-fn',
        help='Input. (Not sure of content yet.)',
    )
    parser.add_argument(
        '--preads-fofn-fn',
        help='Output. FOFN of preads (fasta files).',
    )
    args = parser.parse_args(argv[1:])
    return args


def main(argv=sys.argv):
    args = parse_args(argv)
    logging.basicConfig(level=logging.INFO)
    run(**vars(args))


if __name__ == '__main__':  # pragma: no cover
    main()
