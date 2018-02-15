"""Not sure anything uses the fopfn anymore.
"""
from __future__ import absolute_import
from __future__ import unicode_literals

from future.utils import viewitems
from future.utils import itervalues
import argparse
import logging
import os
import sys
from .. import io, run_support

LOG = logging.getLogger()

def convert_job_id_to_num(job_id):
    """
    >>> convert_job_id_to_num('m_001d')
    29
    """
    return int(job_id[2:], base=16)

def run(gathered_fn, las_paths_fn):
    gathered = io.deserialize(gathered_fn)
    job_done_fns = dict()
    for (key, desc) in viewitems(gathered):
        job_id = key.split('=')[1]
        job_num = convert_job_id_to_num(job_id)
        job_done_fns[job_num] = desc['fns']['job_done']
    import pprint
    LOG.error(pprint.pformat(job_done_fns))
    job_rundirs = sorted(os.path.dirname(fn) for fn in itervalues(job_done_fns))
    # Find all .las leaves so far.
    with open(las_paths_fn, 'w') as stream:
        for block, las_path in run_support.daligner_gather_las(job_rundirs):
            stream.write('{} {}\n'.format(block, las_path))


class HelpF(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass


def parse_args(argv):
    description = 'Turn gathered file into .las paths.'
    epilog = ''
    parser = argparse.ArgumentParser(
        description=description,
        epilog=epilog,
        formatter_class=HelpF,
    )
    parser.add_argument(
        '--gathered-fn',
        help='Input. (From generic pypeflow gen_parallel.)',
    )
    parser.add_argument(
        '--las-paths-fn',
        help='Output. JSON of las files.',
    )
    args = parser.parse_args(argv[1:])
    return args


def main(argv=sys.argv):
    args = parse_args(argv)
    logging.basicConfig(level=logging.INFO)
    run(**vars(args))


if __name__ == '__main__':  # pragma: no cover
    main()
