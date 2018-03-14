"""
"""
from __future__ import absolute_import
from __future__ import print_function
from __future__ import unicode_literals

from future.utils import viewitems
import argparse
import logging
import os
import sys
from .. import io

LOG = logging.getLogger()


def run(gathered_fn, preads_fofn_fn):
    gathered = io.deserialize(gathered_fn)
    fasta_fns = list()
    for desc in gathered:
        fasta_fns.append(desc['fasta'])
    with open(preads_fofn_fn,  'w') as f:
        for filename in sorted(fasta_fns):
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
        help='Input. JSON list of output dicts.')
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
