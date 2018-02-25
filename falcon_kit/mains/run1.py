from __future__ import absolute_import
from __future__ import unicode_literals
from ..pype import (gen_task, gen_parallel_tasks) # copied verbatim from falcon_unzip
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
        #merge_script = parameters['merge_script']
        #sge_option = parameters['sge_option']
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

    LOG.info('fc_run started with configuration %s', input_config_fn)
    try:
        config = support.get_dict_from_old_falcon_cfg(
            support.parse_config(input_config_fn))
    except Exception:
        LOG.exception('Failed to parse config "{}".'.format(input_config_fn))
        raise
    input_fofn_plf = makePypeLocalFile(config['input_fofn'])
    genome_size = config.get('genome_size')
    squash = True if 0 < genome_size < 1000000 else False
    wf = PypeProcWatcherWorkflow(job_type=config['job_type'],
                                 job_queue=config['job_queue'],
                                 job_name_style=config['job_name_style'],
                                 sge_option=config.get('sge_option', ''),
                                 watcher_type=config['pwatcher_type'],
                                 watcher_directory=config['pwatcher_directory'],
                                 use_tmpdir=config.get('use_tmpdir'),
                                 squash=squash
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
    exitOnFailure = config['stop_all_jobs_on_failure']
    wf.max_jobs = config['default_concurrent_jobs']

    assert config['input_type'] in (
        'raw', 'preads'), 'Invalid input_type=={!r}'.format(config['input_type'])

    # Store config as JSON, available to many tasks.

    if config['input_type'] == 'raw':
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
        ))

        # run daligner
        wf.max_jobs = config['da_concurrent_jobs']
        rawreads_db_fn = os.path.join(rawread_dir, 'raw_reads.db')
        #config['sge_option_da'] = config['sge_option_pda']
        scattered_fn = os.path.join(
            rawread_dir, 'daligner-scatter', 'scattered.json')
        params = dict() #dict(parameters)
        params['db_prefix'] = 'raw_reads'
        params['stage'] = os.path.basename(rawread_dir)
        params['pread_aln'] = 0
        params['skip_checks'] = int(config.get('skip_checks', 0))
        params['wildcards'] = 'dal0_id'
        wf.addTask(gen_task(
            script=pype_tasks.TASK_DALIGNER_SCATTER_SCRIPT,
            inputs={
                'run_jobs': run_jobs_fn,
                'db': rawreads_db_fn,
            },
            outputs={
                'scattered': scattered_fn,
            },
            parameters=params,
            rule_writer=rule_writer,
        ))

        gathered_fn = os.path.join(rawread_dir, 'daligner-gathered', 'gathered.json')
        gen_parallel_tasks(
            wf, rule_writer,
            scattered_fn, gathered_fn,
            run_dict=dict(
                script=pype_tasks.TASK_DALIGNER_SCRIPT,
                inputs={
                    'daligner_script': '0-rawreads/daligner-scripts/{dal0_id}/daligner-script.sh',
                    'daligner_settings': '0-rawreads/daligner-scripts/{dal0_id}/settings.json',
                },
                outputs={
                    'job_done': '0-rawreads/{dal0_id}/daligner.done',
                },
                parameters={},
            ),
        )

        r_gathered_las_fn = os.path.join(rawread_dir, 'daligner-intermediate-gathered-las', 'gathered-las.txt')
        wf.addTask(gen_task(
            script=pype_tasks.TASK_DALIGNER_GATHER_SCRIPT,
            inputs={'gathered': gathered_fn,
            },
            outputs={'las_paths': r_gathered_las_fn,
            },
            parameters={},
            rule_writer=rule_writer,
        ))

        # Merge .las files.
        wf.max_jobs = config['la_concurrent_jobs']
        scattered_fn = os.path.join(rawread_dir, 'merge-scatter', 'scattered.json')
        params = dict() #(parameters)
        params['db_prefix'] = 'raw_reads'
        params['stage'] = os.path.basename(rawread_dir) # TODO(CD): Make this more clearly constant.
        params['wildcards'] = 'mer0_id'
        wf.addTask(gen_task(
            script=pype_tasks.TASK_LAS_MERGE_SCATTER_SCRIPT,
            inputs={
                'run_jobs': run_jobs_fn,
                'gathered_las': r_gathered_las_fn,
            },
            outputs={
                'scattered': scattered_fn,
            },
            parameters=params,
            rule_writer=rule_writer,
        ))

        gathered_fn = os.path.join(rawread_dir, 'merge-gathered', 'gathered.json')
        gen_parallel_tasks(
            wf, rule_writer,
            scattered_fn, gathered_fn,
            run_dict=dict(
                script=pype_tasks.TASK_LAS_MERGE_SCRIPT,
                inputs={
                    'las_paths': './0-rawreads/merge-scripts/{mer0_id}/las_paths.json',
                    'merge_script': './0-rawreads/merge-scripts/{mer0_id}/merge-script.sh',
                    'merged_las_json': './0-rawreads/merge-scripts/{mer0_id}/merged_las.json',
                },
                outputs={
                    'merged_las': './0-rawreads/{mer0_id}/merged.las',
                    'job_done': './0-rawreads/{mer0_id}/merge.done',
                },
                parameters={},
            ),
        )

        las_fofn_fn = os.path.join(rawread_dir, 'merged-las-fofn', 'las.fofn')
        las_fopfn_fn = os.path.join(rawread_dir, 'merged-las-fofn', 'las.fopfn')
        wf.addTask(gen_task(
            script=pype_tasks.TASK_LAS_MERGE_GATHER_SCRIPT,
            inputs={'gathered': gathered_fn,
            },
            outputs={'las_fofn': las_fofn_fn,
                     'las_fopfn': las_fopfn_fn,
            },
            parameters={},
            rule_writer=rule_writer,
        ))

        if config['target'] == 'overlapping':
            sys.exit(0)

        # Produce new FOFN of preads fasta, based on consensus of overlaps.
        wf.max_jobs = config['cns_concurrent_jobs']

        scattered_fn = os.path.join(
            rawread_dir, 'cns-scatter', 'scattered.json')
        params = dict()
        params['wildcards'] = 'cns0_id,cns0_id2'
        wf.addTask(gen_task(
            script=pype_tasks.TASK_CONSENSUS_SCATTER_SCRIPT,
            inputs={
                'las_fopfn': las_fopfn_fn,
                'raw_reads_db': raw_reads_db_fn,
                'length_cutoff': length_cutoff_fn,
                'config': general_config_fn,
            },
            outputs={
                'scattered': scattered_fn,
            },
            parameters=params,
            rule_writer=rule_writer,
        ))

        gathered_fn = os.path.join(rawread_dir, 'cns-gather', 'gathered.json')
        gen_parallel_tasks(
            wf, rule_writer,
            scattered_fn, gathered_fn,
            run_dict=dict(
                script=pype_tasks.TASK_CONSENSUS_TASK_SCRIPT,
                inputs = {
                    'las': '0-rawreads/cns-scatter/{cns0_id}/merged.{cns0_id2}.las',
                    'db': raw_reads_db_fn,
                    'length_cutoff': length_cutoff_fn,
                    'config': general_config_fn,
                },
                outputs = {
                    'fasta': '0-rawreads/consensus/{cns0_id}/consensus.{cns0_id2}.fasta',
                },
                parameters={},
            ),
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
            parameters={}, #=parameters,
            rule_writer=rule_writer,
        ))

        rdir = os.path.join(rawread_dir, 'report')
        pre_assembly_report_fn = os.path.join(rdir, 'pre_assembly_stats.json')
        params = dict() #dict(parameters)
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

    #parameters = {'work_dir': pread_dir,
    #              'sge_option': config['sge_option_pda'],
    #              #'config_fn': input_config_fn,
    #              'config': config}
    parameters = dict()

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
        parameters={},
        rule_writer=rule_writer,
    ))

    # run daligner
    wf.max_jobs = config['pda_concurrent_jobs']
    #config['sge_option_da'] = config['sge_option_pda']
    scattered_fn = os.path.join(
        pread_dir, 'daligner-scatter', 'scattered.json')
    params = dict() #dict(parameters)
    params['db_prefix'] = 'preads'
    params['stage'] = os.path.basename(pread_dir)
    params['pread_aln'] = 1
    #params['nblock'] = preads_nblock
    params['skip_checks'] = int(config.get('skip_checks', 0))
    params['wildcards'] = 'dal1_id'
    wf.addTask(gen_task(
        script=pype_tasks.TASK_DALIGNER_SCATTER_SCRIPT,
        inputs={
            'run_jobs': run_jobs_fn,
            'db': preads_db_fn,
        },
        outputs={
            'scattered': scattered_fn,
        },
        parameters=params,
        rule_writer=rule_writer,
    ))

    gathered_fn = os.path.join(pread_dir, 'daligner-gathered', 'gathered.json')
    gen_parallel_tasks(
        wf, rule_writer,
        scattered_fn, gathered_fn,
        run_dict=dict(
            script=pype_tasks.TASK_DALIGNER_SCRIPT,
            inputs={
                'daligner_script': '1-preads_ovl/daligner-scripts/{dal1_id}/daligner-script.sh',
                'daligner_settings': '1-preads_ovl/daligner-scripts/{dal1_id}/settings.json',
            },
            outputs={
                'job_done': '1-preads_ovl/{dal1_id}/daligner.done',
            },
            parameters={},
        ),
    )

    p_gathered_las_fn = os.path.join(pread_dir, 'daligner-intermediate-gathered-las', 'gathered-las.txt')
    wf.addTask(gen_task(
        script=pype_tasks.TASK_DALIGNER_GATHER_SCRIPT,
        inputs={'gathered': gathered_fn,
        },
        outputs={'las_paths': p_gathered_las_fn,
        },
        parameters={},
        rule_writer=rule_writer,
    ))

    # Merge .las files.
    #wf.max_jobs = config['pla_concurrent_jobs']
    #config['sge_option_la'] = config['sge_option_pla']
    scattered_fn = os.path.join(pread_dir, 'merge-scatter', 'scattered.json')
    #task = make_task(pype_tasks.task_merge_scatter)
    #wf.refreshTargets(exitOnFailure=exitOnFailure)
    params = dict() #(parameters)
    params['db_prefix'] = 'preads'
    params['stage'] = os.path.basename(pread_dir) # TODO(CD): Make this more clearly constant.
    params['wildcards'] = 'mer1_id'
    wf.addTask(gen_task(
        script=pype_tasks.TASK_LAS_MERGE_SCATTER_SCRIPT,
        inputs={
            'run_jobs': run_jobs_fn,
            'gathered_las': p_gathered_las_fn,
        },
        outputs={
            'scattered': scattered_fn,
        },
        parameters=params,
        rule_writer=rule_writer,
    ))

    gathered_fn = os.path.join(pread_dir, 'merge-gathered', 'gathered.json')
    gen_parallel_tasks(
        wf, rule_writer,
        scattered_fn, gathered_fn,
        run_dict=dict(
            script=pype_tasks.TASK_LAS_MERGE_SCRIPT,
            inputs={
                'las_paths': './1-preads_ovl/merge-scripts/{mer1_id}/las_paths.json',
                'merge_script': './1-preads_ovl/merge-scripts/{mer1_id}/merge-script.sh',
                'merged_las_json': './1-preads_ovl/merge-scripts/{mer1_id}/merged_las.json',
            },
            outputs={
                'merged_las': './1-preads_ovl/{mer1_id}/merged.las',
                'job_done': './1-preads_ovl/{mer1_id}/merge.done',
            },
            parameters={},
        ),
    )

    las_fofn_fn = os.path.join(pread_dir, 'merged-las-fofn', 'las.fofn')
    las_fopfn_fn = os.path.join(pread_dir, 'merged-las-fofn', 'las.fopfn')
    wf.addTask(gen_task(
        script=pype_tasks.TASK_LAS_MERGE_GATHER_SCRIPT,
        inputs={'gathered': gathered_fn,
        },
        outputs={'las_fofn': las_fofn_fn,
                 'las_fopfn': las_fopfn_fn,
        },
        parameters={},
        rule_writer=rule_writer,
    ))

    # Draft assembly (called 'fc_' for now)
    wf.max_jobs = config['fc_concurrent_jobs']
    db2falcon_dir = os.path.join(pread_dir, 'db2falcon')
    db2falcon_done_fn = os.path.join(db2falcon_dir, 'db2falcon_done')
    preads4falcon_fn = os.path.join(db2falcon_dir, 'preads4falcon.fasta')
    wf.addTask(gen_task(
        script=pype_tasks.TASK_RUN_DB_TO_FALCON_SCRIPT,
        inputs={'las_fofn': las_fofn_fn,
                'preads_db': preads_db_fn,
                },
        outputs={'job_done': db2falcon_done_fn,
                 'preads4falcon': preads4falcon_fn,
                 },
        parameters={},
        rule_writer=rule_writer,
    ))

    falcon_asm_done_fn = os.path.join(falcon_asm_dir, 'falcon_asm_done')
    parameters = {
        #'sge_option': config['sge_option_fc'], # if we want this, we must use double curlies
    }
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
    ))
    wf.refreshTargets()

    #return falcon_asm_done


def main(argv=sys.argv):
    lfs_setstripe_maybe(path='.', stripe=12)
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
