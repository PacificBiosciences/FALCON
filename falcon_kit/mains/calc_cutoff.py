from .. import functional as f
import argparse
import os
import sys
import traceback


def main(argv=sys.argv):
    import argparse

    description = """
Given the result of 'DBstats -u -b1' on stdin,
print the lowest read-length required for sufficient coverage of the genome
(i.e. 'length_cutoff').
"""
    epilog = """
This is useful when length_cutoff is not provided but the genome-size
can be estimated. The purpose is to *reduce* the amount of data seen by
DALIGNER, since otherwise it will miss many alignments when it
encounters resource limits.

Note: If PBFALCON_ERRFILE is defined (and its directory is writable),
we will write errors there in addition to stderr.
"""
    parser = argparse.ArgumentParser(description=description, epilog=epilog,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--coverage', type=float, default=20,
                        help='Desired coverage ratio (i.e. over-sampling)')
    parser.add_argument('genome_size', type=int,
                        help='Estimated number of bases in genome. (haploid?)')
    parser.add_argument('capture',  # default='-', # I guess default is not allowed for required args.
                        help='File with captured output of DBstats. (Otherwise, stdin.)')
    args = parser.parse_args(argv[1:])

    target = int(args.genome_size * args.coverage)
    capture = open(args.capture) if args.capture != '-' else sys.stdin
    stats = capture.read()
    try:
        cutoff = f.calc_cutoff(target, stats)
    except Exception:
        msg = traceback.format_exc()
        msg += 'User-provided genome_size: {}\nDesired coverage: {}\n'.format(
            args.genome_size, args.coverage)
        # pbfalcon wants us to write errs here.
        errfile = os.environ.get('PBFALCON_ERRFILE')
        if errfile:
            with open(errfile, 'w') as ofs:
                ofs.write(msg)
        raise Exception(msg)
    sys.stdout.write(str(cutoff))


if __name__ == "__main__":
    main(sys.argv)
