#!/usr/bin/env python
"""Usage:

    LAsort.py DB <args>

Run LAcheck on each input in args. Exclude any failures from
the arglist. Then run LAsort on the remaining arglist.
"""
from __future__ import unicode_literals
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
    log('argv:{!r}'.format(argv))
    db = argv[1]
    args = argv[2:]  # Skip program name
    lass = list()
    new_args = list()
    new_args.append('LAsort')
    for arg in args:
        if arg.startswith('-'):
            new_args.append(arg)
        else:
            lass.append(arg)
    for las in lass:
        rc = system('LAcheck -v {} {}.las'.format(db, las))
        if rc:
            log('Skipping {}.las'.format(las))
        else:
            new_args.append(las)
    system(' '.join(new_args))


if __name__ == "__main__":
    main()
