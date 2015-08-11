"""Executable tasks.

To be called by pbsmrtpipe.

pypeFLOW uses its own adaptors instead.
"""
from .. import run_support as support
import sys


def help():
    print("""
Usage:
    falcon-task [task] <[task-args]>

tasks:
    make-fofn-abs
""")
    sys.exit(2)

def main_make_fofn_abs(i_fofn_fn, o_fofn_fn):
    support.make_fofn_abs(i_fofn_fn, o_fofn_fn)

def main(argv=sys.argv):
    if len(argv) < 2 or argv[1].startswith('-'):
        help()
    task = argv[1]
    tasks = {
        'make-fofn-abs': main_make_fofn_abs,
    }
    return tasks[task](*argv[2:])
