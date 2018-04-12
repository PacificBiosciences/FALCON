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
    # Copy General section to top, for now.
    for key, val in config['General'].items():
        config[key] = val
    input_fofn_plf = makePypeLocalFile(config['input_fofn'])
    genome_size = config.get('genome_size')
    squash = True if 0 < genome_size < 1000000 else False
    wf = PypeProcWatcherWorkflow(job_defaults=config['job.defaults'],
                                 squash=squash,
    )
    general_config_fn = './config.json' # must not be in a task-dir
    config['ver'] = '100'
    io.serialize(general_config_fn, config)
    with open('foo.snake', 'w') as snakemake_writer:
        rule_writer = snakemake.SnakemakeRuleWriter(snakemake_writer)
        run(wf, config, rule_writer,
            os.path.abspath(general_config_fn),
            input_fofn_plf=input_fofn_plf,
            )


def run(wf, config, rule_writer,
        general_config_fn,
        input_fofn_plf,
        ):
    """
    Preconditions (for now):
    * LOG
    * run_support.logger
    """
    general_config = io.deserialize(general_config_fn)
    if general_config != config:
        msg = 'Config from {!r} != passed config'.format(general_config_fn)
        LOG.error(msg)
        raise Exception(msg)
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

    assert config['input_type'] in (
        'raw', 'preads'), 'Invalid input_type=={!r}'.format(config['input_type'])

    # Store config as JSON, available to many tasks.

    if config['input_type'] == 'raw':
        parameters = {}
        rdb_build_done = os.path.join(rawread_dir, 'rdb_build_done')
        run_jobs_fn = os.path.join(rawread_dir, 'run_jobs.sh')
        length_cutoff_fn = os.path.join(rawread_dir, 'length_cutoff')
        raw_reads_db_fn = os.path.join(rawread_dir, 'raw_reads.db')
        # Also .raw_reads.*, of course.

        # import sequences into daligner DB
        # and calculate length_cutoff (if specified as -1)
        wf.addTask(gen_task(
            script=pype_tasks.TASK_BUILD_RDB_SCRIPT,
            inputs={
                'config': general_config_fn,
                'raw_reads_fofn': fn(input_fofn_plf),
            },
            outputs={
                'run_jobs': run_jobs_fn,
                'raw_reads_db': raw_reads_db_fn,
                'length_cutoff': length_cutoff_fn,
                'db_build_done': rdb_build_done, # only for ordering
            },
            parameters={},
            rule_writer=rule_writer,
            dist=Dist(NPROC=1, job_dict=config['job.step.da']),
        ))

        # run daligner
        wf.max_jobs = config['job.step.da'].get('njobs', default_njobs)
        rawreads_db_fn = os.path.join(rawread_dir, 'raw_reads.db')
        daligner_all_units_fn = os.path.join(
            rawread_dir, 'daligner-split', 'all-units-of-work.json')
        daligner_bash_template_fn = os.path.join(
            rawread_dir, 'daligner-split', 'daligner_bash_template.sh')
        params = dict(parameters)
        params['db_prefix'] = 'raw_reads'
        params['pread_aln'] = 0
        params['skip_checks'] = int(config.get('skip_checks', 0))
        params['wildcards'] = 'dal0_id'
        wf.addTask(gen_task(
            script=pype_tasks.TASK_DALIGNER_SPLIT_SCRIPT,
            inputs={
                'run_jobs': run_jobs_fn,
                'db': rawreads_db_fn,
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
                script=pype_tasks.TASK_DALIGNER_SCRIPT, # for snakemake stuff
                inputs={
                    #'daligner_script': '0-rawreads/daligner-scripts/{dal0_id}/daligner-script.sh',
                    #'daligner_settings': '0-rawreads/daligner-scripts/{dal0_id}/settings.json',
                    'units_of_work': '0-rawreads/daligner-chunks/{dal0_id}/some-units-of-work.json',
                },
                outputs={
                    #'job_done': '0-rawreads/{dal0_id}/daligner.done',
                    'results': '0-rawreads/daligner-runs/{dal0_id}/some-done-files.json',
                },
                parameters={},
            ),
            dist=Dist(NPROC=4, MB=4000, job_dict=config['job.step.da']),
        )

        r_gathered_las_fn = os.path.join(rawread_dir, 'daligner-intermediate-gathered-las', 'gathered-las.json')
        wf.addTask(gen_task(
            script=pype_tasks.TASK_DALIGNER_FIND_LAS_SCRIPT,
            inputs={'gathered': gathered_fn,
            },
            outputs={'las_paths': r_gathered_las_fn,
            },
            parameters={},
            rule_writer=rule_writer,
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
            script=pype_tasks.TASK_LAS_MERGE_SPLIT_SCRIPT,
            inputs={
                'run_jobs': run_jobs_fn,
                'gathered_las': r_gathered_las_fn,
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
                script=pype_tasks.TASK_LAS_MERGE_SCRIPT, # for snakemake
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

        p_id2las_fn = os.path.join(rawread_dir, 'las-gather', 'p_id2las.json')
        las_fofn_fn = os.path.join(rawread_dir, 'las-gather', 'las_fofn.json')
        wf.addTask(gen_task(
            script=pype_tasks.TASK_LAS_MERGE_GATHER_SCRIPT,
            inputs={'gathered': gathered_fn,
            },
            outputs={'p_id2las': p_id2las_fn,
                     'las': las_fofn_fn,
            },
            parameters={},
            rule_writer=rule_writer,
            dist=Dist(local=True),
        ))

        if config['target'] == 'overlapping':
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
                'raw_reads_db': raw_reads_db_fn,
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
                    #'db': raw_reads_db_fn,
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
        params['length_cutoff_user'] = config['length_cutoff']
        params['genome_length'] = config['genome_size'] # note different name; historical
        wf.addTask(gen_task(
            script=pype_tasks.TASK_REPORT_PRE_ASSEMBLY_SCRIPT,
            inputs={'length_cutoff': length_cutoff_fn,
                    'raw_reads_db': raw_reads_db_fn,
                    'preads_fofn': preads_fofn_fn,
                    'config': general_config_fn,
            },
            outputs={'pre_assembly_report': pre_assembly_report_fn,
            },
            parameters=params,
            rule_writer=rule_writer,
            dist=Dist(local=True),
        ))

    if config['target'] == 'pre-assembly':
        LOG.info('Quitting after stage-0 for "pre-assembly" target.')
        sys.exit(0)

    # build pread database
    if config['input_type'] == 'preads':
        """
        preads_fofn_plf = makePypeLocalFile(os.path.join(
            pread_dir, 'preads-fofn-abs', os.path.basename(config['input_fofn'])))
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
    preads_db_fn = os.path.join(pread_dir, 'preads.db')
    # Also .preads.*, of course.

    wf.addTask(gen_task(
        script=pype_tasks.TASK_BUILD_PDB_SCRIPT,
        inputs={
            'config': general_config_fn,
            'preads_fofn': preads_fofn_fn,
        },
        outputs={
            'run_jobs': run_jobs_fn,
            'preads_db': preads_db_fn,
            'db_build_done': pdb_build_done, # only for ordering
        },
        parameters=parameters, #{},
        rule_writer=rule_writer,
        dist=Dist(NPROC=1, job_dict=config['job.step.pda']),
    ))

    # run daligner
    wf.max_jobs = config['job.step.pda'].get('njobs', default_njobs)
    daligner_all_units_fn = os.path.join(
        pread_dir, 'daligner-split', 'all-units-of-work.json')
    daligner_bash_template_fn = os.path.join(
        pread_dir, 'daligner-split', 'daligner_bash_template.sh')
    params = dict(parameters)
    params['db_prefix'] = 'preads'
    #params['stage'] = os.path.basename(pread_dir)
    params['pread_aln'] = 1
    #params['nblock'] = preads_nblock
    params['skip_checks'] = int(config.get('skip_checks', 0))
    params['wildcards'] = 'dal1_id'
    wf.addTask(gen_task(
        script=pype_tasks.TASK_DALIGNER_SPLIT_SCRIPT,
        inputs={
            'run_jobs': run_jobs_fn,
            'db': preads_db_fn,
        },
        outputs={
            'split': daligner_all_units_fn,
            'bash_template': daligner_bash_template_fn
        },
        parameters=params,
        rule_writer=rule_writer,
        dist=Dist(local=True, NPROC=4),
    ))

    gathered_fn = os.path.join(pread_dir, 'daligner-gathered', 'gathered-done-files.json')
    gen_parallel_tasks(
        wf, rule_writer,
        daligner_all_units_fn, gathered_fn,
        run_dict=dict(
            bash_template_fn=daligner_bash_template_fn,
            script=pype_tasks.TASK_DALIGNER_SCRIPT, # for snakemake stuff
            inputs={
                #'daligner_script': '1-preads_ovl/daligner-scripts/{dal1_id}/daligner-script.sh',
                #'daligner_settings': '1-preads_ovl/daligner-scripts/{dal1_id}/settings.json',
                'units_of_work': '1-preads_ovl/daligner-chunks/{dal1_id}/some-units-of-work.json',
            },
            outputs={
                #'job_done': '1-preads_ovl/{dal1_id}/daligner.done',
                'results': '1-preads_ovl/daligner-runs/{dal1_id}/some-done-files.json',
            },
            parameters={},
        ),
        dist=Dist(NPROC=4, MB=4000, job_dict=config['job.step.pda']),
    )

    p_gathered_las_fn = os.path.join(pread_dir, 'daligner-intermediate-gathered-las', 'gathered-las.json')
    wf.addTask(gen_task(
        script=pype_tasks.TASK_DALIGNER_FIND_LAS_SCRIPT,
        inputs={'gathered': gathered_fn,
        },
        outputs={'las_paths': p_gathered_las_fn,
        },
        parameters={},
        rule_writer=rule_writer,
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
        script=pype_tasks.TASK_LAS_MERGE_SPLIT_SCRIPT,
        inputs={
            'run_jobs': run_jobs_fn,
            'gathered_las': p_gathered_las_fn,
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
            script=pype_tasks.TASK_LAS_MERGE_SCRIPT, # for snakemake
            inputs={
                #'las_paths': './1-preads_ovl/merge-scripts/{mer1_id}/las_paths.json',
                #'merge_script': './1-preads_ovl/merge-scripts/{mer1_id}/merge-script.sh',
                #'merged_las_json': './1-preads_ovl/merge-scripts/{mer1_id}/merged_las.json',
                'units_of_work': '1-preads_ovl/las-merge-chunks/{mer1_id}/some-units-of-work.json',
            },
            outputs={
                #'merged_las': './1-preads_ovl/{mer1_id}/merged.las',
                #'job_done': './1-preads_ovl/{mer1_id}/merge.done',
                'results': '1-preads_ovl/las-merge-runs/{mer1_id}/some-las-paths.json',
            },
            parameters={},
        ),
        dist=Dist(NPROC=1, job_dict=config['job.step.pla']),
    )

    p_id2las_fn = os.path.join(pread_dir, 'las-gather', 'p_id2las.json')
    las_fofn_fn = os.path.join(pread_dir, 'las-gather', 'las_fofn.json')
    wf.addTask(gen_task(
        script=pype_tasks.TASK_LAS_MERGE_GATHER_SCRIPT,
        inputs={'gathered': gathered_fn,
        },
        outputs={'p_id2las': p_id2las_fn,
                 'las': las_fofn_fn,
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
        parameters[key] = config[key]
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
