#!/bin/env python2.7
import argparse
import json
import os
import sys

def deserialize(fn):
    with open(fn) as ifs:
        return json.loads(ifs.read())

def assert_exists(fn):
    if not os.path.isfile(fn):
        raise Exception('Does not exist: {!r}'.format(fn))

def mkdir(dirname):
    if not os.path.isdir(dirname):
        # Possible race-condition, so dirs must be created serially.
        os.makedirs(dirname)

def symlink(name, target):
    msg = '{} -> {}'.format(name, target)
    assert not os.path.lexists(name), msg
    #print msg
    os.symlink(target, name)

def run(special_split_fn, fn_patterns):
    """
    Symlink targets will be relative to cwd.
    For each pattern, each wildcard will be substituted everywhere, e.g.
        fn_pattern == 'top/{key}/input_{key}.txt'
    """
    fnkeypattdict = dict(fnkeypatt.split('=') for fnkeypatt in fn_patterns)
    jobs = deserialize(special_split_fn)
    mapdir = os.path.normpath(os.path.dirname(os.path.normpath(special_split_fn)))
    for job in jobs:
        inputs = job['input']
        wildcards = job['wildcards']
        for fnkey, fn_pattern in fnkeypattdict.iteritems():
            val = inputs[fnkey]
            # val should be relative to the location of the special_split_fn.
            #assert not os.path.isabs(val), 'mapped input (dynamic output) filename {!r} must be relative (to serialzed file location {!r})'.format(
            #        val, special_split_fn)
            if not os.path.isabs(val):
                mapped_input_fn = os.path.join(mapdir, val)
            else:
                mapped_input_fn = val
            assert_exists(mapped_input_fn)
            try:
                symlink_name = fn_pattern.format(**wildcards)
            except Exception as err:
                import pprint
                msg = str(err) + ': for pattern {!r} and wildcards\n{!r}'.format(
                        fn_pattern, pprint.pformat(wildcards))
                raise Exception(msg)
            outdir = os.path.normpath(os.path.dirname(symlink_name))
            mkdir(outdir)
            target_name = os.path.relpath(mapped_input_fn, outdir)
            symlink(symlink_name, target_name)

def parse_args(argv):
    description = 'Create symlinks named after "fn_pattern", targeting values in "mapped_fn".'
    parser = argparse.ArgumentParser(
            description=description,
    )
    parser.add_argument(
            '--special-split-fn', required=True,
            help='Serialized split-file (in our special format), where "mapped_inputs" has a map with key to filename, relative to the directory of this file.')
    parser.add_argument(
            'fn_patterns', nargs='+',
            help='"fnkey=pattern" Can appear multiple times. Each is a pattern for symlinks, to be substituted with keys in special_split_fn. Each fnkey=filename must appear in the input section of each job listed in special-split.')
    return parser.parse_args(argv[1:])

def main(argv=sys.argv):
    args = parse_args(argv)
    run(**vars(args))

if __name__ == "__main__":
    main()
