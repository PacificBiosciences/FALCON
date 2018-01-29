import argparse
import logging
import os
import sys
from .. import io

LOG = logging.getLogger()


def run(gathered_fn, scattered_fn):
    thatdir = os.path.dirname(scattered_fn)
    thisdir = os.path.dirname(gathered_fn)
    scattered = io.deserialize(scattered_fn)
    gathered = dict()
    for job in scattered:
        job_output = dict()
        #job_output['wildcards'] = dict()
        fn_dict = dict(job['output'])
        for key in fn_dict.keys():
            # Fix path to be relative to gathered_fn.
            fn = fn_dict[key]
            if not os.path.isabs(fn):
                thatfn = os.path.join(thatdir, fn)
            else:
                thatfn = fn
            thisfn = os.path.relpath(thatfn, thisdir)
            fn_dict[key] = thisfn
        job_output['fns'] = fn_dict
        wildcards = job['wildcards']
        wildkey = ','.join('{}={}'.format(k,v) for k,v in sorted(wildcards.items()))
        gathered[wildkey] = job_output
    io.serialize(gathered_fn, gathered)


class HelpF(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass


def parse_args(argv):
    description = 'Gather generic filenames into ... something. For now, just serialize.'
    epilog = 'We expect the scattered file to have a specific format.'
    parser = argparse.ArgumentParser(
        description=description,
        epilog=epilog,
        formatter_class=HelpF,
    )
    parser.add_argument(
        '--scattered-fn',
        help='Input: result of scattering',
    )
    parser.add_argument(
        '--gathered-fn',
        help='Output: serialized something-or-other',
    )
    args = parser.parse_args(argv[1:])
    return args


def main(argv=sys.argv):
    args = parse_args(argv)
    logging.basicConfig(level=logging.INFO)
    run(**vars(args))


if __name__ == '__main__':  # pragma: no cover
    main()
