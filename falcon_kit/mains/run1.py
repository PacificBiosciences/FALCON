from __future__ import absolute_import

from ..pype import (wrap_gen_task as gen_task, gen_parallel_tasks, Dist)
from .. import run_support as support
from .. import bash, pype_tasks, snakemake
from ..util.system import (only_these_symlinks, lfs_setstripe_maybe)
from .. import io
# pylint: disable=no-name-in-module, import-error, fixme, line-too-long
from pypeflow.simple_pwatcher_bridge import (PypeProcWatcherWorkflow, MyFakePypeThreadTaskBase,
                                             makePypeLocalFile, fn, PypeTask)
import argparse
import glob
import json
import logging
import os
import re
import sys
import time


LOG = logging.getLogger(__name__)  # default, for remote tasks


def create_daligner_tasks(basedir, scatter_fn):
    tasks = []
    tasks_out = {}
    try:
        content = json.loads(open(scatter_fn).read())  # array of descriptions
    except Exception:
        msg = 'Failed to read JSON from {!r}'.format(scatter_fn)
        LOG.exception(msg)
        raise Exception(msg)
    for section in content:
        parameters = section['parameters']
        inputs = section['inputs']
        inputs['scatter_fn'] = os.path.abspath(scatter_fn)
        outputs = section['outputs']
        URL = section['URL']
        job_uid = parameters['job_uid']
        wdir = os.path.join(basedir, 'job_%s' % job_uid)
        make_daligner_task = PypeTask(inputs=inputs,
                                      outputs=outputs,
                                      parameters=parameters,
                                      wdir=wdir,
                                      )
        daligner_task = make_daligner_task(pype_tasks.task_run_daligner)
        tasks.append(daligner_task)
        # these are relative, so we need the PypeLocalFiles
        tasks_out['ajob_%s' % job_uid] = daligner_task.outputs['job_done']
    return tasks, tasks_out


def create_merge_tasks(basedir, scatter_fn):
    tasks = []
    p_ids_merged_las = {}  # for consensus
    content = json.loads(open(scatter_fn).read())  # array of descriptions
    for section in content:
        parameters = section['parameters']
        inputs = section['inputs']
        inputs['scatter_fn'] = os.path.abspath(scatter_fn)
        outputs = section['outputs']
        URL = section['URL']
        p_id = parameters['job_id']
        wdir = os.path.join(basedir, 'm_%05d' % p_id)
        make_task = PypeTask(inputs=inputs,
                             outputs=outputs,
                             parameters=parameters,
                             wdir=wdir,
                             )
        task = make_task(pype_tasks.task_run_las_merge)
        tasks.append(task)
        # these are relative, so we need the PypeLocalFiles
        las_fn = task.outputs['merged_las']
        p_ids_merged_las[p_id] = las_fn
    return tasks, p_ids_merged_las


def create_consensus_tasks(basedir, scatter_fn):
    consensus_tasks = []
    consensus_out = {}
    content = json.loads(open(scatter_fn).read())  # array of descriptions
    for section in content:
        parameters = section['parameters']
        inputs = section['inputs']
        inputs['scatter_fn'] = scatter_fn
        outputs = section['outputs']
        URL = section['URL']
        p_id = int(parameters['job_id'])
        cns_label = 'cns_%05d' % int(p_id)
        wdir = os.path.join(basedir, 'preads', cns_label)
        make_c_task = PypeTask(inputs=inputs,
                               outputs=outputs,
                               parameters=parameters,
                               wdir=wdir,
                               )
        c_task = make_c_task(pype_tasks.task_run_consensus)
        consensus_tasks.append(c_task)
        consensus_out['cjob_%d' % p_id] = outputs['out_file']
    return consensus_tasks, consensus_out


def create_merge_gather_task(wd, inputs):
    las_fofn_plf = makePypeLocalFile(os.path.join(wd, 'las.fofn'))
    las_fopfn_plf = makePypeLocalFile(os.path.join(wd, 'las.fopfn'))

    make_task = PypeTask(inputs=inputs,  # p_ids_merged_las
                         outputs={'las_fofn': las_fofn_plf,
                                  'las_fopfn': las_fopfn_plf,
                                  },
                         )
    task = make_task(pype_tasks.task_merge_gather)
    return task, las_fofn_plf, las_fopfn_plf


def create_consensus_gather_task(wd, inputs):
    # Happens only in stage-0.
    preads_fofn_plf = makePypeLocalFile(os.path.join(wd, 'input_preads.fofn'))

    make_cns_gather_task = PypeTask(
        inputs=inputs,  # consensus_out
        outputs={'preads_fofn': preads_fofn_plf},
    )
    task = make_cns_gather_task(pype_tasks.task_cns_gather)
    return task, preads_fofn_plf


def main1(prog_name, input_config_fn, logger_config_fn=None):
    global LOG
    LOG = support.setup_logger(logger_config_fn)
    lfs_setstripe_maybe(path='.', stripe=12)

    LOG.info('fc_run started with configuration %s', input_config_fn)
    try:
        config = support.parse_cfg_file(input_config_fn)
        import json
        dumped = json.dumps(config, indent=2, separators=(',', ': '), sort_keys=True)
        LOG.info('cfg=\n{}'.format(dumped))
    except Exception:
        LOG.exception('Failed to parse config "{}".'.format(input_config_fn))
        raise
    general_config = config['General']
    assert 'input_fofn' in general_config, 'Missing "input_fofn" in {}.'.format(input_config_fn)
    input_fofn_plf = makePypeLocalFile(general_config['input_fofn'])
    genome_size = int(general_config.get('genome_size', '0'))
    squash = True if 0 < genome_size < 1000000 else False
    wf = PypeProcWatcherWorkflow(job_defaults=config['job.defaults'],
                                 squash=squash,
    )
    general_config['ver'] = '100'
    config_fn = './config.json' # must not be in a task-dir
    io.serialize(config_fn, config)
    with open('foo.snake', 'w') as snakemake_writer:
        rule_writer = snakemake.SnakemakeRuleWriter(snakemake_writer)
        run(wf, config, rule_writer,
            os.path.abspath(config_fn),
            input_fofn_plf=input_fofn_plf,
            )


def run(wf, config, rule_writer,
        config_fn,
        input_fofn_plf,
        ):
    """
    Preconditions (for now):
    * LOG
    * run_support.logger
    """
    parsed_config = io.deserialize(config_fn)
    if parsed_config != config:
        msg = 'Config from {!r} != passed config'.format(config_fn)
        raise Exception(msg)
    general_config = config['General']
    general_config_fn = os.path.join(os.path.dirname(config_fn), 'General_config.json')
    io.serialize(general_config_fn, general_config) # Some tasks use this.
    rawread_dir = '0-rawreads'
    pread_dir = '1-preads_ovl'
    falcon_asm_dir = '2-asm-falcon'

    for d in (rawread_dir, pread_dir, falcon_asm_dir):
        support.make_dirs(d)

    # only matter for parallel jobs
    job_defaults = config['job.defaults']
    exitOnFailure = bool(job_defaults.get('stop_all_jobs_on_failure', False))
    default_njobs = int(job_defaults.get('njobs', 7))
    wf.max_jobs = default_njobs

    assert general_config['input_type'] in (
        'raw', 'preads'), 'Invalid input_type=={!r}'.format(general_config['input_type'])

    # Store config as JSON, available to many tasks.

    if general_config['input_type'] == 'raw':
        parameters = {}

        # import sequences into daligner DB
        # calculate length_cutoff (if specified as -1)
        # split DB
        # run DBdust
        r_db_dust_fn = os.path.join(rawread_dir, 'build', 'raw_reads.db')
        length_cutoff_fn = os.path.join(rawread_dir, 'build', 'length_cutoff')
        wf.addTask(gen_task(
            script=pype_tasks.TASK_DB_BUILD_SCRIPT,
            inputs={
                'config': general_config_fn,
                'input_fofn': fn(input_fofn_plf),
            },
            outputs={
                'length_cutoff': length_cutoff_fn,
                'db': r_db_dust_fn,
                # Also .raw_reads.*, of course. And dust track.
            },
            parameters=dict(
            ),
            rule_writer=rule_writer,
            dist=Dist(NPROC=1, job_dict=config['job.step.bd']),
        ))

        # run TANmask
        tan_uows_fn = os.path.join(
            rawread_dir, 'tan-split', 'tan-uows.json')
        tan_bash_template_fn = os.path.join(
            rawread_dir, 'tan-split', 'bash_template.sh')
        wf.addTask(gen_task(
            script=pype_tasks.TASK_DB_TAN_SPLIT_SCRIPT,
            inputs={
                'config': general_config_fn,
                'db': r_db_dust_fn,
            },
            outputs={
                'split': tan_uows_fn,
                'bash_template': tan_bash_template_fn,
            },
            parameters={},
            rule_writer=rule_writer,
            dist=Dist(NPROC=1),
        ))

        gathered_fn = os.path.join(rawread_dir, 'tan-gathered', 'gathered-done-files.json')
        gen_parallel_tasks(
            wf, rule_writer,
            tan_uows_fn, gathered_fn,
            run_dict=dict(
                bash_template_fn=tan_bash_template_fn,
                script='fubar-TODO', #pype_tasks.TASK_DB_TAN_APPLY_SCRIPT, # for snakemake stuff
                inputs={
                    'units_of_work': '0-rawreads/tan-chunks/{tan0_id}/some-units-of-work.json',
                },
                outputs={
                    #'job_done': '0-rawreads/{dal0_id}/daligner.done',
                    'results': '0-rawreads/tan-runs/{tan0_id}/some-done-files.json',
                },
                parameters={},

            ),
            dist=Dist(NPROC=4, MB=4000, job_dict=config['job.step.da']),
        )

        r_db_tan_fn = os.path.join(rawread_dir, 'tan-combine', 'raw_reads.db')
        wf.addTask(gen_task(
            script=pype_tasks.TASK_DB_TAN_COMBINE_SCRIPT,
            inputs={
                'config': general_config_fn,
                'db': r_db_dust_fn,
                'gathered': gathered_fn,
            },
            outputs={
                'new_db': r_db_tan_fn,
            },
            parameters={},
            rule_writer=rule_writer,
            dist=Dist(local=True),
        ))

        # run daligner
        wf.max_jobs = config['job.step.da'].get('njobs', default_njobs)
        #rawreads_db_fn = os.path.join(rawread_dir, 'raw_reads.db')
        daligner_all_units_fn = os.path.join(
            rawread_dir, 'daligner-split', 'all-units-of-work.json')
        daligner_bash_template_fn = os.path.join(
            rawread_dir, 'daligner-split', 'daligner_bash_template.sh')
        params = dict(parameters)
        #params['db_prefix'] = 'raw_reads'
        #params['pread_aln'] = 0
        params['skip_checks'] = int(general_config.get('skip_checks', 0))
        params['wildcards'] = 'dal0_id'
        wf.addTask(gen_task(
            script=pype_tasks.TASK_DB_DALIGNER_SPLIT_SCRIPT,
            inputs={
                'config': general_config_fn,
                'db': r_db_tan_fn,
                'length_cutoff': length_cutoff_fn,
            },
            outputs={
                'split': daligner_all_units_fn,
                'bash_template': daligner_bash_template_fn
            },
            parameters=params,
            rule_writer=rule_writer,
            dist=Dist(local=True, NPROC=4), # really, NPROC=1, but we need to know the max
        ))

        gathered_fn = os.path.join(rawread_dir, 'daligner-gathered', 'gathered-done-files.json')
        gen_parallel_tasks(
            wf, rule_writer,
            daligner_all_units_fn, gathered_fn,
            run_dict=dict(
                bash_template_fn=daligner_bash_template_fn,
                script=pype_tasks.TASK_DB_DALIGNER_APPLY_SCRIPT, # for snakemake stuff
                inputs={
                    'units_of_work': os.path.join(rawread_dir, 'daligner-chunks/{dal0_id}/some-units-of-work.json'),
                },
                outputs={
                    'results': os.path.join(rawread_dir, 'daligner-runs/{dal0_id}/some-done-files.json'),
                },
                parameters={},
            ),
            dist=Dist(NPROC=4, MB=4000, job_dict=config['job.step.da']),
        )

        r_gathered_las_fn = os.path.join(rawread_dir, 'daligner-combine', 'gathered-las.json')
        wf.addTask(gen_task(
            script=pype_tasks.TASK_DB_DALIGNER_COMBINE_SCRIPT,
            inputs={
                'config': general_config_fn,
                'db': r_db_tan_fn,
                'gathered': gathered_fn,
            },
            outputs={
                'las_paths': r_gathered_las_fn,
            },
            parameters={},
            rule_writer=rule_writer,
            #dist=Dist(NPROC=1, MB=4000, job_dict=config['job.step.da'])
            dist=Dist(local=True),
        ))

        # Merge .las files.
        wf.max_jobs = config['job.step.la'].get('njobs', default_njobs)
        las_merge_all_units_fn = os.path.join(rawread_dir, 'las-merge-split', 'all-units-of-work.json')
        bash_template_fn = os.path.join(rawread_dir, 'las-merge-split', 'las-merge-bash-template.sh')
        params = dict(parameters)
        params['db_prefix'] = 'raw_reads'
        params['wildcards'] = 'mer0_id'
        wf.addTask(gen_task(
            script=pype_tasks.TASK_DB_LAMERGE_SPLIT_SCRIPT,
            inputs={
                'config': general_config_fn,
                'las_paths': r_gathered_las_fn,
            },
            outputs={
                'split': las_merge_all_units_fn,
                'bash_template': bash_template_fn,
            },
            parameters=params,
            rule_writer=rule_writer,
            dist=Dist(local=True),
        ))

        gathered_fn = os.path.join(rawread_dir, 'las-merge-gathered', 'gathered.json')
        gen_parallel_tasks(
            wf, rule_writer,
            las_merge_all_units_fn, gathered_fn,
            run_dict=dict(
                bash_template_fn=bash_template_fn,
                script=pype_tasks.TASK_DB_LAMERGE_APPLY_SCRIPT, # for snakemake
                inputs={
                    #'las_paths': './0-rawreads/merge-scripts/{mer0_id}/las_paths.json',
                    #'merge_script': './0-rawreads/merge-scripts/{mer0_id}/merge-script.sh',
                    #'merged_las_json': './0-rawreads/merge-scripts/{mer0_id}/merged_las.json',
                    'units_of_work': '0-rawreads/las-merge-chunks/{mer0_id}/some-units-of-work.json',
                },
                outputs={
                    #'merged_las': './0-rawreads/{mer0_id}/merged.las',
                    #'job_done': './0-rawreads/{mer0_id}/merge.done',
                    'results': '0-rawreads/las-merge-runs/{mer0_id}/some-las-paths.json',
                },
                parameters={},
            ),
            dist=Dist(NPROC=1, job_dict=config['job.step.la']),
        )

        p_id2las_fn = os.path.join(rawread_dir, 'las-merge-combine', 'p_id2las.json')
        las_fofn_fn = os.path.join(rawread_dir, 'las-merge-combine', 'las_fofn.json')
        wf.addTask(gen_task(
            script=pype_tasks.TASK_DB_LAMERGE_COMBINE_SCRIPT,
            inputs={
                'config': general_config_fn,
                'gathered': gathered_fn,
            },
            outputs={
                'block2las': p_id2las_fn,
                'las_paths': las_fofn_fn,
            },
            parameters={},
            rule_writer=rule_writer,
            dist=Dist(local=True),
        ))

        if general_config['target'] == 'overlapping':
            sys.exit(0)

        # Produce new FOFN of preads fasta, based on consensus of overlaps.
        wf.max_jobs = config['job.step.cns'].get('njobs', default_njobs)

        split_fn = os.path.join(
            rawread_dir, 'cns-split', 'split.json')
        bash_template_fn = os.path.join(
            rawread_dir, 'cns-split', 'consensus-bash-template.sh')
        params = dict(parameters)
        params['wildcards'] = 'cns0_id,cns0_id2'
        wf.addTask(gen_task(
            script=pype_tasks.TASK_CONSENSUS_SPLIT_SCRIPT,
            inputs={
                'p_id2las': p_id2las_fn,
                'raw_reads_db': r_db_tan_fn,
                'length_cutoff': length_cutoff_fn,
                'config': general_config_fn,
            },
            outputs={
                'split': split_fn,
                'bash_template': bash_template_fn,
            },
            parameters=params,
            rule_writer=rule_writer,
            dist=Dist(local=True),
        ))

        gathered_fn = os.path.join(rawread_dir, 'cns-gather', 'gathered.json')
        gen_parallel_tasks(
            wf, rule_writer,
            split_fn, gathered_fn,
            run_dict=dict(
                bash_template_fn=bash_template_fn,
                script=pype_tasks.TASK_CONSENSUS_TASK_SCRIPT, # for snakemake only
                inputs = {
                    #'las': '0-rawreads/cns-split/{cns0_id}/merged.{cns0_id2}.las',
                    #'db': r_db_tan_fn,
                    #'length_cutoff': length_cutoff_fn,
                    #'config': general_config_fn,
                    'units_of_work': '0-rawreads/cns-chunks/{cns0_id}/some-units-of-work.json',
                },
                outputs = {
                    #'fasta': '0-rawreads/consensus/{cns0_id}/consensus.{cns0_id2}.fasta',
                    'results': '0-rawreads/cns-runs/{cns0_id}/some-done-files.json',
                },
                parameters={},
            ),
            dist=Dist(NPROC=6, job_dict=config['job.step.cns']),
        )
        preads_fofn_fn = os.path.join(rawread_dir, 'preads', 'input_preads.fofn')
        wf.addTask(gen_task(
            script=pype_tasks.TASK_CONSENSUS_GATHER_SCRIPT,
            inputs={
                'gathered': gathered_fn,
            },
            outputs={
                'preads_fofn': preads_fofn_fn,
            },
            parameters=parameters, #{},
            rule_writer=rule_writer,
            dist=Dist(local=True),
        ))

        rdir = os.path.join(rawread_dir, 'report')
        pre_assembly_report_fn = os.path.join(rdir, 'pre_assembly_stats.json')
        params = dict(parameters)
        params['length_cutoff_user'] = general_config['length_cutoff']
        params['genome_length'] = general_config['genome_size'] # note different name; historical
        wf.addTask(gen_task(
            script=pype_tasks.TASK_REPORT_PRE_ASSEMBLY_SCRIPT,
            inputs={'length_cutoff': length_cutoff_fn,
                    'raw_reads_db': r_db_tan_fn,
                    'preads_fofn': preads_fofn_fn,
                    'config': general_config_fn,
            },
            outputs={'pre_assembly_report': pre_assembly_report_fn,
            },
            parameters=params,
            rule_writer=rule_writer,
            dist=Dist(local=True),
        ))

    if general_config['target'] == 'pre-assembly':
        LOG.info('Quitting after stage-0 for "pre-assembly" target.')
        sys.exit(0)

    # build pread database
    if general_config['input_type'] == 'preads':
        """
        preads_fofn_plf = makePypeLocalFile(os.path.join(
            pread_dir, 'preads-fofn-abs', os.path.basename(general_config['input_fofn'])))
        make_fofn_abs_task = PypeTask(inputs={'i_fofn': input_fofn_plf},
                                      outputs={'o_fofn': preads_fofn_plf},
                                      parameters={},
                                      )
        fofn_abs_task = make_fofn_abs_task(
            pype_tasks.task_make_fofn_abs_preads)
        wf.addTasks([fofn_abs_task])
        wf.refreshTargets([fofn_abs_task])
        """
        raise Exception('TODO')

    pdb_build_done = os.path.join(pread_dir, 'pdb_build_done')
    run_jobs_fn = os.path.join(pread_dir, 'run_jobs.sh')
    preads_db_fn = os.path.join(pread_dir, 'build', 'preads.db')
    length_cutoff_pr_fn = os.path.join(pread_dir, 'build', 'length_cutoff')

    wf.addTask(gen_task(
        script=pype_tasks.TASK_DB_BUILD_SCRIPT,
        inputs={
            'config': general_config_fn,
            'input_fofn': preads_fofn_fn,
        },
        outputs={
            'length_cutoff': length_cutoff_pr_fn,
            'db': preads_db_fn,
            # Also .preads.*, of course.
        },
        parameters=dict(
        ),
        rule_writer=rule_writer,
        dist=Dist(NPROC=1),
    ))

    # run daligner
    wf.max_jobs = config['job.step.pda'].get('njobs', default_njobs)
    daligner_all_units_fn = os.path.join(
        pread_dir, 'daligner-split', 'all-units-of-work.json')
    daligner_bash_template_fn = os.path.join(
        pread_dir, 'daligner-split', 'daligner_bash_template.sh')
    params = dict(parameters)
    params['skip_checks'] = int(general_config.get('skip_checks', 0))
    params['wildcards'] = 'dal1_id'
    wf.addTask(gen_task(
        script=pype_tasks.TASK_DB_DALIGNER_SPLIT_SCRIPT,
        inputs={
            'config': general_config_fn,
            'db': preads_db_fn, #not tan, yet
            'length_cutoff': length_cutoff_pr_fn,
        },
        outputs={
            'split': daligner_all_units_fn,
            'bash_template': daligner_bash_template_fn
        },
        parameters=params,
        rule_writer=rule_writer,
        dist=Dist(local=True, NPROC=4), # really, NPROC=1, but we need to know the max
    ))

    gathered_fn = os.path.join(pread_dir, 'daligner-gathered', 'gathered-done-files.json')
    gen_parallel_tasks(
        wf, rule_writer,
        daligner_all_units_fn, gathered_fn,
        run_dict=dict(
            bash_template_fn=daligner_bash_template_fn,
            script=pype_tasks.TASK_DB_DALIGNER_APPLY_SCRIPT, # for snakemake stuff
            inputs={
                'units_of_work': os.path.join(pread_dir, 'daligner-chunks/{dal1_id}/some-units-of-work.json'),
            },
            outputs={
                'results': os.path.join(pread_dir, 'daligner-runs/{dal1_id}/some-done-files.json'),
            },
            parameters={},
        ),
        dist=Dist(NPROC=4, MB=4000, job_dict=config['job.step.pda']),
    )

    gathered_las_fn = os.path.join(pread_dir, 'daligner-combine', 'gathered-las.json')
    wf.addTask(gen_task(
        script=pype_tasks.TASK_DB_DALIGNER_COMBINE_SCRIPT,
        inputs={
            'config': general_config_fn,
            'db': preads_db_fn, #r_db_tan_fn,
            'gathered': gathered_fn,
        },
        outputs={
            'las_paths': gathered_las_fn,
        },
        parameters={},
        rule_writer=rule_writer,
        #dist=Dist(NPROC=1, MB=4000, job_dict=config['job.step.pda'])
        dist=Dist(local=True),
    ))

    # Merge .las files.
    wf.max_jobs = config['job.step.pla'].get('njobs', default_njobs)
    las_merge_all_units_fn = os.path.join(pread_dir, 'las-merge-split', 'all-units-of-work.json')
    bash_template_fn = os.path.join(pread_dir, 'las-merge-split', 'las-merge-bash-template.sh')
    params = dict(parameters)
    params['db_prefix'] = 'preads'
    params['wildcards'] = 'mer1_id'
    wf.addTask(gen_task(
        script=pype_tasks.TASK_DB_LAMERGE_SPLIT_SCRIPT,
        inputs={
            'config': general_config_fn,
            'las_paths': gathered_las_fn,
        },
        outputs={
            'split': las_merge_all_units_fn,
            'bash_template': bash_template_fn,
        },
        parameters=params,
        rule_writer=rule_writer,
        dist=Dist(local=True),
    ))

    gathered_fn = os.path.join(pread_dir, 'las-merge-gathered', 'gathered.json')
    gen_parallel_tasks(
        wf, rule_writer,
        las_merge_all_units_fn, gathered_fn,
        run_dict=dict(
            bash_template_fn=bash_template_fn,
            script=pype_tasks.TASK_DB_LAMERGE_APPLY_SCRIPT, # for snakemake
            inputs={
                'units_of_work': os.path.join(pread_dir, 'las-merge-chunks/{mer0_id}/some-units-of-work.json'),
            },
            outputs={
                'results': os.path.join(pread_dir, 'las-merge-runs/{mer0_id}/some-las-paths.json'),
            },
            parameters={},
        ),
        dist=Dist(NPROC=1, job_dict=config['job.step.la']),
    )

    p_id2las_fn = os.path.join(pread_dir, 'las-merge-combine', 'block2las.json')
    las_fofn_fn = os.path.join(pread_dir, 'las-merge-combine', 'las_fofn.json')
    wf.addTask(gen_task(
        script=pype_tasks.TASK_DB_LAMERGE_COMBINE_SCRIPT,
        inputs={
            'config': general_config_fn,
            'gathered': gathered_fn,
        },
        outputs={
            'block2las': p_id2las_fn,
            'las_paths': las_fofn_fn,
        },
        parameters={},
        rule_writer=rule_writer,
        dist=Dist(local=True),
    ))

    wf.max_jobs = config['job.step.asm'].get('njobs', default_njobs)
    db2falcon_dir = os.path.join(pread_dir, 'db2falcon')
    db2falcon_done_fn = os.path.join(db2falcon_dir, 'db2falcon_done')
    preads4falcon_fn = os.path.join(db2falcon_dir, 'preads4falcon.fasta')
    wf.addTask(gen_task(
        script=pype_tasks.TASK_RUN_DB_TO_FALCON_SCRIPT,
        inputs={'p_id2las': p_id2las_fn,
                'preads_db': preads_db_fn,
                },
        outputs={'job_done': db2falcon_done_fn,
                 'preads4falcon': preads4falcon_fn,
                 },
        parameters={},
        rule_writer=rule_writer,
        dist=Dist(NPROC=4, job_dict=config['job.step.asm']),
    ))

    falcon_asm_done_fn = os.path.join(falcon_asm_dir, 'falcon_asm_done')
    for key in ('overlap_filtering_setting', 'length_cutoff_pr', 'fc_ovlp_to_graph_option'):
        parameters[key] = general_config[key]
    wf.addTask(gen_task(
        script=pype_tasks.TASK_RUN_FALCON_ASM_SCRIPT,
        inputs={'db2falcon_done': db2falcon_done_fn, 'db_file': preads_db_fn,
                'preads4falcon_fasta': preads4falcon_fn,
                'las_fofn': las_fofn_fn,
                'config': general_config_fn,
                },
        outputs={'falcon_asm_done': falcon_asm_done_fn},
        parameters=parameters,
        rule_writer=rule_writer,
        dist=Dist(NPROC=4, job_dict=config['job.step.asm']),
    ))
    wf.refreshTargets()

    with io.cd('0-rawreads'):
        # for backwards-compatibility
        io.symlink('las-merge-combine', 'las-gather')

    #return falcon_asm_done


def main(argv=sys.argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('config',
                        help='.cfg/.ini/.json')
    parser.add_argument('logger',
                        nargs='?',
                        help='(Optional)JSON config for standard Python logging module')
    args = parser.parse_args(argv[1:])
    main1(argv[0], args.config, args.logger)


if __name__ == '__main__':
    main()
