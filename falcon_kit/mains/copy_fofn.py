from __future__ import absolute_import
from __future__ import unicode_literals
import argparse
import logging
import os
import sys
from .. import io

LOG = logging.getLogger()


def run(abs, in_fn, out_fn):
    out_dir = os.path.normpath(os.path.dirname(out_fn))
    io.mkdirs(out_dir)
    def identity(fn): return fn
    def relative(fn): return os.path.relpath(fn, out_dir)
    adjusted_fn = identity if abs else relative
    with open(out_fn, 'w') as stream:
        for abs_fn in io.yield_abspath_from_fofn(in_fn):
            fn = adjusted_fn(abs_fn)
            stream.write('{}\n'.format(fn))


class HelpF(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass


def parse_args(argv):
    description = 'Copy FOFN. If directory changes, then relative paths must change too.'
    epilog = ''
    parser = argparse.ArgumentParser(
        description=description,
        epilog=epilog,
        formatter_class=HelpF,
    )
    parser.add_argument(
        '--in-fn',
        help='Input. FOFN of paths relative to its own directory.'
    )
    parser.add_argument(
        '--abs', action='store_true',
        help='Store absolute paths. (Otherwise, paths will be relative to directory of output FOFN.)'
    )
    parser.add_argument(
        '--out-fn',
        help='Output. FOFN of paths relative to its own directory.'
    )
    args = parser.parse_args(argv[1:])
    return args


def main(argv=sys.argv):
    args = parse_args(argv)
    logging.basicConfig(level=logging.INFO)
    run(**vars(args))


if __name__ == '__main__':  # pragma: no cover
    main()
