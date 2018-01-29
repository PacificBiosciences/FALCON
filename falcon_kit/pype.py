"""This was copied from falcon_unzip, but we
needed to modify the TASK SCRIPT to use our copy of
generic_gather.py
"""
from pypeflow.sample_tasks import gen_task as pype_gen_task
from . import io
import os

import logging
LOG = logging.getLogger(__name__)

TASK_GENERIC_GATHER_SCRIPT = """
python -m falcon_kit.mains.generic_gather --scattered-fn={input.scattered} --gathered-fn={output.gathered}
"""


def gen_task(rule_writer, script, inputs, outputs, parameters={}):
    rel_inputs = dict()
    rel_outputs = dict()
    # Make relative to CWD. (But better if caller does this.)
    def get_rel(maybe_abs):
        rel = dict()
        for k,v in maybe_abs.items():
            if os.path.isabs(v):
                v = os.path.relpath(v)
            rel[k] = v
        return rel
    inputs = get_rel(inputs)
    outputs = get_rel(outputs)

    first_output_dir = os.path.normpath(os.path.dirname(outputs.values()[0]))
    rel_topdir = os.path.relpath('.', first_output_dir) # redundant for rel-inputs, but fine
    params = dict(parameters)
    params['topdir'] = rel_topdir

    pt = pype_gen_task(script, inputs, outputs, params)

    # Run pype_gen_task first because it can valid some stuff.
    rule_writer(inputs, outputs, params, script)
    return pt


def gen_parallel_tasks(
        wf, rule_writer,
        scattered_fn,
        gathered_fn,
        run_dict,
):
    # run_dict['inputs'] should be patterns to match the inputs in scattered_fn, by convention.

    # Write 3 wildcard rules for snakemake, 2 with dynamic.
    rule_writer.write_dynamic_rules(
            rule_name="foo",
            input_json=scattered_fn,
            inputs=dict_rel_paths(run_dict['inputs']),
            shell_template=run_dict['script'],
            parameters=run_dict['parameters'],
            wildcard_outputs=dict_rel_paths(run_dict['outputs']),
            output_json=gathered_fn,
    )

    #outputs = {k:patt.format(**jobkv) for k,patt in output_patterns}
    #inputs =  {k:patt.format(**jobkv) for k,patt in input_patterns}
    #inputs['SCATTERED'] = scattered_fn # presumably ignored by script; might not be needed at all
    #scattered_fn = scatter_dict['outputs']['scattered'] # by convention
    wf.refreshTargets()
    scattered = io.deserialize(scattered_fn)

    gather_inputs = {}
    for job in scattered:
        inputs = job['input']
        outputs = job['output']
        params = job['params']
        params.update({k: v for k,v in job['wildcards'].items()}) # include expanded wildcards
        LOG.warning('OUT:{}'.format(outputs))
        wf.addTask(pype_gen_task(
                script=run_dict['script'],
                inputs=inputs,
                outputs=outputs,
                parameters=params,
        ))
        wildcards_str = '_'.join(w for w in job['wildcards'].values())
        job_name = 'job{}'.format(wildcards_str)
        for output_name, output_fn in outputs.iteritems():
            giname = '{}_{}'.format(job_name, output_name)
            gather_inputs[giname] = output_fn

    # An implicit "gatherer" simply takes the output filenames and writes them into a FOFN.
    assert 'scattered' not in gather_inputs
    gather_inputs['scattered'] = scattered_fn
    #gather_inputs.update(gather_dict['inputs'])
    wf.addTask(pype_gen_task(
        script=TASK_GENERIC_GATHER_SCRIPT,
        inputs=gather_inputs,
        outputs={
            'gathered': gathered_fn,
        },
        parameters={},
    ))


def dict_rel_paths(dict_paths):
    return {k: os.path.relpath(v) for k,v in dict_paths.iteritems()}
