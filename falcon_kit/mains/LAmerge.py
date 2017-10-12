#!/usr/bin/env python
"""Usage:

    LAmerge.py DB <args>

Run LAcheck on each input in args. Exclude any failures from
the arglist. Then run LAmerge on the remaining arglist.

This differs from LAsort.py in that the first las file is actually
an *explicit* output, whereas LAsort relies on *implicit* outputs.
"""
import sys
import os


def log(msg):
    sys.stderr.write(msg + '\n')


def system(call, checked=False):
    log('!{}'.format(call))
    rc = os.system(call)
    if rc:
        msg = '{} <- {!r}'.format(rc, call)
        if checked:
            raise Exception(msg)
        log(msg)
    return rc


def main(argv=sys.argv):
    db = argv[1]
    args = argv[2:]  # Skip program name
    lass = list()
    new_args = list()
    new_args.append('LAmerge')
    for arg in args:
        if arg.startswith('-'):
            new_args.append(arg)
        else:
            lass.append(arg)
    outlas = lass[0]
    new_args.append(outlas)  # This is the output las.
    for las in lass[1:]:
        rc = system('LAcheck -vS {} {}.las'.format(db, las))  # Assume sorted.
        if rc:
            log('Skipping {}.las'.format(las))
        else:
            new_args.append(las)
    system(' '.join(new_args))
    system('LAcheck -vS {} {}.las'.format(db, outlas))  # Assume sorted.


if __name__ == "__main__":
    main()
