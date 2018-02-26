"""Not sure anything uses the fopfn anymore.
"""
from __future__ import absolute_import
from __future__ import unicode_literals

#from future.utils import viewitems
#from future.utils import itervalues
import argparse
import glob
import logging
import os
import sys
from .. import io #, run_support

LOG = logging.getLogger()


def run(gathered_fn, las_paths_fn):
    gathered = io.deserialize(gathered_fn)
    job_done_fns = list()
    for job_output in gathered:
        for fn in job_output.values():
            job_done_fns.append(fn)
    import pprint
    LOG.info('job_done_fns: {}'.format(pprint.pformat(job_done_fns)))
    job_rundirs = sorted(os.path.dirname(fn) for fn in job_done_fns)
    # Find all .las leaves so far.
    #[block, las_path in run_support.daligner_gather_las(job_rundirs)]
    las_paths = list()
    for uow_dir in job_rundirs:
        # We could assert the existence of a job_done file here.
        d = os.path.abspath(uow_dir)
        las_path_glob = glob.glob(os.path.join(d, '*.las'))
        LOG.debug('dir={!r}, glob={!r}'.format(d, las_path_glob))
        las_paths.extend(las_path_glob)
    io.serialize(las_paths_fn, las_paths)


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
