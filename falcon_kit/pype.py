"""This was copied from falcon_unzip, but we
needed to modify the TASK SCRIPT to use our copy of
generic_gather.py
"""
from __future__ import absolute_import
from __future__ import unicode_literals

from future.utils import viewitems
from future.utils import itervalues
from pypeflow.sample_tasks import gen_task as pype_gen_task
from . import io
import os
import tempfile

import logging
LOG = logging.getLogger(__name__)

TASK_GENERIC_RUN_UNITS_SCRIPT = """\
python -m falcon_kit.mains.generic_run_units_of_work --units-of-work-fn={input.units_of_work} --bash-template-fn={input.bash_template} --results-fn={output.results}
"""
TASK_GENERIC_SCATTER_ONE_UOW_SCRIPT = """\
python -m falcon_kit.mains.generic_scatter_one_uow --all-uow-list-fn={input.all} --one-uow-list-fn={output.one} --split-idx={params.split_idx}
"""
TASK_GENERIC_UNSPLIT_SCRIPT = """
python -m falcon_kit.mains.generic_unsplit --result-fn-list-fn={output.result_fn_list} --gathered-fn={output.gathered}
"""
#TASK_GENERIC_CHUNKING_SCRIPT = """\
#python -m falcon_kit.mains.generic_chunking split-fn={input.split} --bash-template-temp-fn={input.bash_template_temp} --units-of-work-fn={output.units_of_work} --uow-template-fn={output.uow_template} --split-idx={params.split_idx}
#"""


def gen_task(rule_writer, script, inputs, outputs, parameters={}):
    rel_inputs = dict()
    rel_outputs = dict()
    # Make relative to CWD. (But better if caller does this.)
    def get_rel(maybe_abs):
        rel = dict()
        for (k, v) in viewitems(maybe_abs):
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
        split_fn,
        gathered_fn,
        run_dict,
):
    """
    By convention, the first (wildcard) output in run_dict['outputs'] must be the gatherable list,
    in the same format as the gathered_fn to be generated from them.

    For now, we require a single such output, since we do not yet test for wildcards.
    """
    # run_dict['inputs'] should be patterns to match the inputs in split_fn, by convention.

    # Write 3 wildcard rules for snakemake, 2 with dynamic.
    rule_writer.write_dynamic_rules(
            rule_name="foo",
            input_json=split_fn,
            inputs=dict_rel_paths(run_dict['inputs']),
            shell_template=run_dict['script'],
            parameters=run_dict['parameters'],
            wildcard_outputs=dict_rel_paths(run_dict['outputs']),
            output_json=gathered_fn,
    )

    #outputs = {k:patt.format(**jobkv) for k,patt in output_patterns}
    #inputs =  {k:patt.format(**jobkv) for k,patt in input_patterns}
    #inputs['SPLIT'] = split_fn # presumably ignored by script; might not be needed at all
    #split_fn = scatter_dict['outputs']['split'] # by convention
    wf.refreshTargets()
    split = io.deserialize(split_fn)
    bash_template_fn = run_dict['bash_template_fn']

    def find_wildcard_input(inputs):
        for k,v in inputs.items():
            if '{' in v:
                return v
        else:
            raise Exception('No wildcard inputs among {!r}'.format(inputs))

    LOG.warning('PARALLEL OUTPUTS:{}'.format(run_dict['outputs']))
    task_results = dict()
    for split_idx, job in enumerate(split):
        #inputs = job['input']
        #outputs = job['output']
        #params = job['params']
        wildcards = job['wildcards']
        #params.update({k: v for (k, v) in viewitems(job['wildcards'])}) # include expanded wildcards
        #LOG.warning('OUT:{}'.format(outputs))

        wildcards = job['wildcards']
        def resolved(v):
            return v.format(**wildcards)
        def resolved_dict(d):
            result = dict(d)
            LOG.warning('wildcards={!r}'.format(wildcards))
            for k,v in d.items():
                LOG.warning('v={!r}'.format(v))
                result[k] = v.format(**wildcards)
            return result
        task_inputs = resolved_dict(run_dict['inputs'])
        task_outputs = resolved_dict(run_dict['outputs'])
        task_parameters = resolved_dict(run_dict['parameters'])

        wild_input = find_wildcard_input(run_dict['inputs'])
        one_uow_fn = os.path.abspath(wild_input.format(**wildcards))

        wf.addTask(pype_gen_task(
                script=TASK_GENERIC_SCATTER_ONE_UOW_SCRIPT,
                inputs={
                    'all': split_fn,
                },
                outputs={
                    'one': one_uow_fn,
                },
                parameters={
                    'split_idx': split_idx,
                },
        ))

        wf.addTask(pype_gen_task(
                script=TASK_GENERIC_RUN_UNITS_SCRIPT,
                inputs={
                    'units_of_work': one_uow_fn,
                    'bash_template': bash_template_fn,
                },
                outputs=task_outputs,
                parameters=task_parameters,
        ))
        wildcards_str = '_'.join(w for w in itervalues(job['wildcards']))
        job_name = 'job{}'.format(wildcards_str)
        task_results[job_name] = os.path.abspath(task_outputs.values()[0])

    gather_inputs = dict(task_results)
    ## An implicit "gatherer" simply takes the output filenames and combines their contents.
    result_fn_list_fn = os.path.join(os.path.dirname(gathered_fn), 'result-fn-list.json')
    io.serialize(result_fn_list_fn, list(task_results.values())) # dump into next task-dir before next task starts
    #assert 'result_fn_list' not in gather_inputs
    #gather_inputs['result_fn_list'] = result_fn_list_fn # No! pseudo output, since it must exist in a known directory
    LOG.warning('gather_inputs:{!r}'.format(gather_inputs))
    wf.addTask(pype_gen_task(
        script=TASK_GENERIC_UNSPLIT_SCRIPT,
        inputs=gather_inputs,
        outputs={
            'gathered': gathered_fn,
            'result_fn_list': result_fn_list_fn,
        },
        parameters={},
    ))


def dict_rel_paths(dict_paths):
    return {k: os.path.relpath(v) for (k, v) in viewitems(dict_paths)}
