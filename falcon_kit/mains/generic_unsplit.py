from __future__ import absolute_import
from __future__ import unicode_literals
import argparse
import logging
import os
import sys
from .. import io

LOG = logging.getLogger()


def run(result_fn_list_fn, gathered_fn):
    thatdir = os.path.dirname(result_fn_list_fn)
    thisdir = os.path.dirname(gathered_fn)
    result_fn_list = io.deserialize(result_fn_list_fn)
    io.serialize(gathered_fn, result_fn_list)
    gathered = list()
    for result_fn in result_fn_list:
        some_results = io.deserialize(result_fn)
        gathered.extend(some_results)
    io.serialize(gathered_fn, gathered)


class HelpF(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass


def parse_args(argv):
    description = 'Gather the contents of contents of result-lists into a single gathered-list.'
    epilog = 'results-list is known already, so that is a pseudo output. Its filenames point to the actual, unknown results.'
    # Question: Do we need to know the wildcards for each result?
    parser = argparse.ArgumentParser(
        description=description,
        epilog=epilog,
        formatter_class=HelpF,
    )
    parser.add_argument(
        '--result-fn-list-fn',
        help='Input: Combined list of filenames of results (pseudo output, expected to exist already in our run-dir)')
    parser.add_argument(
        '--gathered-fn',
        help='Output: serialized something-or-other')
    args = parser.parse_args(argv[1:])
    return args


def main(argv=sys.argv):
    args = parse_args(argv)
    logging.basicConfig(level=logging.INFO)
    run(**vars(args))


if __name__ == '__main__':  # pragma: no cover
    main()
