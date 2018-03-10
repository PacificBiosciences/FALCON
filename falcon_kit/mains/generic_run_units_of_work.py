from __future__ import absolute_import
from __future__ import unicode_literals
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


def run(bash_template_fn, units_of_work_fn, results_fn):
    uows = io.deserialize(units_of_work_fn)
    uow_dirs = list()
    results = list()
    for i, uow in enumerate(uows):
        job = uow
        inputs = job['input']
        outputs = job['output'] # assumed to be relative to run-dir
        params = job['params']
        #params.update({k: v for (k, v) in viewitems(job['wildcards'])}) # include expanded wildcards
        LOG.warning('INPUT:{}'.format(inputs))
        LOG.warning('OUTPUT:{}'.format(outputs))
        LOG.warning('PARAMS:{}'.format(params))
        uow_dir = 'uow-{:02d}'.format(i)
        uow_dirs.append(uow_dir)
        io.mkdirs(uow_dir)
        script = open(bash_template_fn).read()
        with io.cd(uow_dir):
            pypeflow.do_task.run_bash(script, inputs, outputs, params)
            resolved_outputs = {k: os.path.abspath(v) for k,v in outputs.items()}
            results.append(resolved_outputs)
        #wildcards_str = '_'.join(w for w in itervalues(job['wildcards']))
        #job_name = 'job{}'.format(wildcards_str)
        #for (output_name, output_fn) in viewitems(outputs):
        #    giname = '{}_{}'.format(job_name, output_name)
        #    gather_inputs[giname] = output_fn
    io.serialize(results_fn, results)


class HelpF(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass


def parse_args(argv):
    description = 'Run a bash script once for each unit-of-work, in its own sub-dir.'
    epilog = 'For now, runs will be in series, since we do not know how many processors we can use.'
    parser = argparse.ArgumentParser(
        description=description,
        epilog=epilog,
        formatter_class=HelpF,
    )
    parser.add_argument(
        '--bash-template-fn',
        help='Input. Template of bash script to run on each unit-of-work, with snakemake-style substitutions.')
    parser.add_argument(
        '--units-of-work-fn',
        help='Input. JSON list of records of unit-of-work. Each record is a dict of input, output, and params (snakemake-style).')
    parser.add_argument(
        '--results-fn',
        help='Output. JSON list of result records.')
    args = parser.parse_args(argv[1:])
    return args


def main(argv=sys.argv):
    args = parse_args(argv)
    logging.basicConfig(level=logging.INFO)
    run(**vars(args))


if __name__ == '__main__':  # pragma: no cover
    main()
