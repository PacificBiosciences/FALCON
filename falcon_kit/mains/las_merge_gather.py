"""Not sure anything uses the fopfn anymore.
"""
import argparse
import logging
import os
import sys
from .. import io

LOG = logging.getLogger()

def convert_job_id_to_p_id(job_id):
    """
    >>> convert_job_id_to_p_id('m_0011')
    11
    """
    return int(job_id[2:], base=10)

def run(gathered_fn, las_fofn_fn, las_fopfn_fn):
    gathered = io.deserialize(gathered_fn)
    las_fns = dict()
    for key, desc in gathered.items():
        job_id = key.split('=')[1]
        p_id = convert_job_id_to_p_id(job_id)
        las_fns[p_id] = desc['fns']['merged_las']
    with open(las_fofn_fn,  'w') as f:
        for filename in sorted(las_fns.values()):
            print >>f, filename
    with open(las_fopfn_fn,  'w') as f:
        # The keys are p_ids.
        for p_id, filename in sorted(las_fns.items()):
            print >>f, p_id, filename
    #wdir = os.path.dirname(las_fofn_fn)
    # pread_dir = os.path.dirname(wdir) # by convention, for now
    # Generate las.fofn in run-dir. # No longer needed!
    #system('find {}/m_*/ -name "preads.*.las" >| {}'.format(pread_dir, las_fofn_fn))


class HelpF(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass


def parse_args(argv):
    description = 'Turn gathered file into .las FOFN (and FOPFN).'
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
        '--las-fofn-fn',
        help='Output. FOFN of las files.',
    )
    parser.add_argument(
        '--las-fopfn-fn',
        help='Output. FOFN of las files, but with a p_id first, then a space.',
    )
    args = parser.parse_args(argv[1:])
    return args


def main(argv=sys.argv):
    args = parse_args(argv)
    logging.basicConfig(level=logging.INFO)
    run(**vars(args))


if __name__ == '__main__':  # pragma: no cover
    main()
