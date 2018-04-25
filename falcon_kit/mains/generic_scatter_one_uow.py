from __future__ import absolute_import

import argparse
import collections
import glob
import logging
import os
import sys
import pypeflow.do_task
from .. import io

LOG = logging.getLogger()


# Here is some stuff basically copied from pypeflow.sample_tasks.py.
def validate(bash_template, inputs, outputs, parameterss):
    LOG.info('bash_script_from_template({}\n\tinputs={!r},\n\toutputs={!r})'.format(
        bash_template, inputs, outputs))
    def validate_dict(mydict):
        "Python identifiers are illegal as keys."
        try:
            collections.namedtuple('validate', mydict.keys())
        except ValueError as exc:
            LOG.exception('Bad key name in task definition dict {!r}'.format(mydict))
            raise
    validate_dict(inputs)
    validate_dict(outputs)
    validate_dict(parameterss)


def run(all_uow_list_fn, split_idx, one_uow_list_fn):
    all_uows = io.deserialize(all_uow_list_fn)
    dn = os.path.abspath(os.path.dirname(all_uow_list_fn))
    one_uow = all_uows[split_idx]

    def fixpath(rel):
        try:
            if rel.startswith('./'):
                return os.path.join(dn, rel)
        except Exception:
            pass
        return rel
    if isinstance(one_uow, dict):
        input_dict = one_uow['input']
        for k, v in input_dict.items():
            input_dict[k] = fixpath(v)

    io.serialize(one_uow_list_fn, [one_uow])


class HelpF(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass


def parse_args(argv):
    description = 'Scatter a single unit-of-work from many units-of-work.'
    epilog = ''
    parser = argparse.ArgumentParser(
        description=description,
        epilog=epilog,
        formatter_class=HelpF,
    )
    parser.add_argument(
        '--all-uow-list-fn',
        help='Input. JSON list of all units of work.')
    parser.add_argument(
        '--split-idx', type=int,
        help='Input. Index into the all-uow-list for our single unit-of-work.')
    parser.add_argument(
        '--one-uow-list-fn',
        help='Output. JSON list of a single unit-of-work.')
    args = parser.parse_args(argv[1:])
    return args


def main(argv=sys.argv):
    args = parse_args(argv)
    logging.basicConfig(level=logging.INFO)
    run(**vars(args))


if __name__ == '__main__':  # pragma: no cover
    main()
