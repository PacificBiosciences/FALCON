from .. import run_support as support
from .. import bash, pype_tasks
from ..util.system import only_these_symlinks
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


fc_run_logger = logging.getLogger(__name__) # default, for remote tasks


def create_daligner_tasks(basedir, scatter_fn):
    tasks = []
    tasks_out = {}
    content = json.loads(open(scatter_fn).read()) # array of descriptions
    for section in content:
        parameters = section['parameters']
        inputs = section['inputs']
        inputs['scatter_fn'] = scatter_fn
        outputs = section['outputs']
        URL = section['URL']
        job_uid = parameters['job_uid']
        wdir = os.path.join(basedir, 'job_%s' %job_uid)
        make_daligner_task = PypeTask(inputs = inputs,
                                      outputs = outputs,
                                      parameters = parameters,
                                      wdir = wdir,
        )
        daligner_task = make_daligner_task(pype_tasks.task_run_daligner)
        tasks.append(daligner_task)
        tasks_out['ajob_%s' % job_uid] = daligner_task.outputs['job_done'] # these are relative, so we need the PypeLocalFiles
    return tasks, tasks_out

def create_merge_tasks(basedir, scatter_fn):
    tasks = []
    p_ids_merged_las = {} # for consensus
    content = json.loads(open(scatter_fn).read()) # array of descriptions
    for section in content:
        parameters = section['parameters']
        inputs = section['inputs']
        inputs['scatter_fn'] = scatter_fn
        outputs = section['outputs']
        URL = section['URL']
        p_id = parameters['job_id']
        #merge_script = parameters['merge_script']
        #sge_option = parameters['sge_option']
        wdir = os.path.join(basedir, 'm_%05d' %p_id)
        make_task = PypeTask(inputs = inputs,
                             outputs = outputs,
                             parameters = parameters,
                             wdir = wdir,
        )
        task = make_task(pype_tasks.task_run_las_merge)
        tasks.append(task)
        las_fn = task.outputs['merged_las'] # these are relative, so we need the PypeLocalFiles
        p_ids_merged_las[p_id] = las_fn
    return tasks, p_ids_merged_las

def create_consensus_tasks(basedir, scatter_fn):
    consensus_tasks = []
    consensus_out ={}
    content = json.loads(open(scatter_fn).read()) # array of descriptions
    for section in content:
        parameters = section['parameters']
        inputs = section['inputs']
        inputs['scatter_fn'] = scatter_fn
        outputs = section['outputs']
        URL = section['URL']
        p_id = int(parameters['job_id'])
        cns_label = 'cns_%05d' %int(p_id)
        wdir = os.path.join(basedir, 'preads', cns_label)
        make_c_task = PypeTask(inputs = inputs,
                               outputs = outputs,
                               parameters = parameters,
                               wdir = wdir,
        )
        c_task = make_c_task(pype_tasks.task_run_consensus)
        consensus_tasks.append(c_task)
        consensus_out['cjob_%d' % p_id] = outputs['out_file']
    return consensus_tasks, consensus_out

def create_merge_gather_task(wd, inputs):
    las_fofn_plf = makePypeLocalFile(os.path.join(wd, 'las.fofn'))
    las_fopfn_plf = makePypeLocalFile(os.path.join(wd, 'las.fopfn'))

    make_task = PypeTask(inputs = inputs, # p_ids_merged_las
                         outputs =  {'las_fofn': las_fofn_plf,
                                     'las_fopfn': las_fopfn_plf,
                         },
    )
    task = make_task(pype_tasks.task_merge_gather)
    return task, las_fofn_plf, las_fopfn_plf

def create_consensus_gather_task(wd, inputs):
    # Happens only in stage-0.
    preads_fofn_plf = makePypeLocalFile(os.path.join(wd, 'input_preads.fofn'))

    make_cns_gather_task = PypeTask(
                inputs = inputs, # consensus_out
                outputs =  {'preads_fofn': preads_fofn_plf},
    )
    task = make_cns_gather_task(pype_tasks.task_cns_gather)
    return task, preads_fofn_plf


def main1(prog_name, input_config_fn, logger_config_fn=None):
    global fc_run_logger
    fc_run_logger = support.setup_logger(logger_config_fn)

    fc_run_logger.info('fc_run started with configuration %s', input_config_fn)
    try:
        config = support.get_dict_from_old_falcon_cfg(support.parse_config(input_config_fn))
    except Exception:
        fc_run_logger.exception('Failed to parse config "{}".'.format(input_config_fn))
        raise
    input_fofn_plf = makePypeLocalFile(config['input_fofn'])
    genome_size = config.get('genome_size')
    squash = True if 0 < genome_size < 1000000 else False
    wf = PypeProcWatcherWorkflow(job_type=config['job_type'],
            job_queue=config['job_queue'],
            sge_option=config.get('sge_option', ''),
            watcher_type=config['pwatcher_type'],
            watcher_directory=config['pwatcher_directory'],
            use_tmpdir=config.get('use_tmpdir'),
            squash=squash
    )
    run(wf, config,
            os.path.abspath(input_config_fn),
            input_fofn_plf=input_fofn_plf,
    )


def run(wf, config,
        input_config_fn,
        input_fofn_plf,
        ):
    """
    Preconditions (for now):
    * fc_run_logger
    * run_support.logger
    """
    rawread_dir = os.path.abspath('./0-rawreads')
    pread_dir = os.path.abspath('./1-preads_ovl')
    falcon_asm_dir  = os.path.abspath('./2-asm-falcon')
    script_dir = os.path.abspath('./scripts')
    sge_log_dir = os.path.abspath('./sge_log')

    for d in (rawread_dir, pread_dir, falcon_asm_dir, script_dir, sge_log_dir):
        support.make_dirs(d)

    exitOnFailure=config['stop_all_jobs_on_failure'] # only matter for parallel jobs
    concurrent_jobs = config['pa_concurrent_jobs']
    wf.max_jobs = concurrent_jobs

    rawread_fofn_plf = makePypeLocalFile(os.path.join(rawread_dir, 'raw-fofn-abs', os.path.basename(config['input_fofn'])))
    make_fofn_abs_task = PypeTask(inputs = {'i_fofn': input_fofn_plf},
                                  outputs = {'o_fofn': rawread_fofn_plf},
                                  parameters = {},
    )
    fofn_abs_task = make_fofn_abs_task(pype_tasks.task_make_fofn_abs_raw)

    wf.addTasks([fofn_abs_task])
    wf.refreshTargets([fofn_abs_task])

    if config['input_type'] == 'raw':
        #### import sequences into daligner DB
        sleep_done = makePypeLocalFile( os.path.join( rawread_dir, 'sleep_done') )
        rdb_build_done = makePypeLocalFile( os.path.join( rawread_dir, 'rdb_build_done') )
        run_jobs = makePypeLocalFile( os.path.join( rawread_dir, 'run_jobs.sh') )
        parameters = {'work_dir': rawread_dir,
                      'sge_option': config['sge_option_da'],
                      'config_fn': input_config_fn,
                      'config': config}

        length_cutoff_plf = makePypeLocalFile(os.path.join(rawread_dir, 'length_cutoff'))
        raw_reads_db_plf = makePypeLocalFile(os.path.join(rawread_dir, '%s.db' % 'raw_reads'))
        make_build_rdb_task = PypeTask(inputs = {'input_fofn': rawread_fofn_plf},
                                      outputs = {'rdb_build_done': rdb_build_done,
                                                 'raw_reads_db': raw_reads_db_plf,
                                                 'length_cutoff': length_cutoff_plf,
                                                 'run_jobs': run_jobs,
                                      },
                                      parameters = parameters,
        )
        build_rdb_task = make_build_rdb_task(pype_tasks.task_build_rdb)

        wf.addTasks([build_rdb_task])
        wf.refreshTargets([rdb_build_done])

        raw_reads_nblock = support.get_nblock(fn(raw_reads_db_plf))
        #### run daligner
        scattered_plf = os.path.join(rawread_dir, 'daligner-scatter', 'scattered.json')
        make_daligner_scatter = PypeTask(
                inputs = {
                    'run_jobs_fn': run_jobs,
                    'db_build_done': rdb_build_done,
                },
                outputs = {
                    'scatter_fn': scattered_plf,
                },
                parameters = {
                    'db_prefix': 'raw_reads',
                    'nblock': raw_reads_nblock,
                    'pread_aln': False,
                    'config': config,
                },
        )
        task = make_daligner_scatter(pype_tasks.task_daligner_scatter)
        wf.addTask(task)
        wf.refreshTargets(exitOnFailure=exitOnFailure)

        daligner_tasks, daligner_out = create_daligner_tasks(rawread_dir, scattered_plf)

        wf.addTasks(daligner_tasks)
        r_gathered_las_plf = makePypeLocalFile(os.path.join(rawread_dir, 'raw-gather', 'gathered_las.txt'))

        parameters =  {
                'nblock': raw_reads_nblock,
        }
        make_daligner_gather = PypeTask(
                   inputs = daligner_out,
                   outputs =  {'gathered': r_gathered_las_plf},
                   parameters = parameters,
        )
        check_r_da_task = make_daligner_gather(pype_tasks.task_daligner_gather)
        wf.addTask(check_r_da_task)
        wf.refreshTargets(exitOnFailure=exitOnFailure)

        # Merge .las files.
        scattered_plf = os.path.join(rawread_dir, 'merge-scatter', 'scattered.json')
        make_task = PypeTask(
                inputs = {
                    'run_jobs': run_jobs,
                    'gathered_las': r_gathered_las_plf,
                },
                outputs = {
                    'scattered': scattered_plf,
                },
                parameters = {
                    'db_prefix': 'raw_reads',
                    'config': config,
                },
        )
        task = make_task(pype_tasks.task_merge_scatter)
        wf.addTask(task)
        wf.refreshTargets(exitOnFailure=exitOnFailure)

        merge_tasks, p_ids_merged_las = create_merge_tasks(rawread_dir, scattered_plf)
        wf.addTasks(merge_tasks)
        task, _, las_fopfn_plf = create_merge_gather_task(os.path.join(rawread_dir, 'merge-gather'), p_ids_merged_las)
        wf.addTask(task)
        wf.refreshTargets(exitOnFailure=exitOnFailure)

        if config['target'] == 'overlapping':
            sys.exit(0)

        # Produce new FOFN of preads fasta, based on consensus of overlaps.
        scattered_plf = os.path.join(rawread_dir, 'cns-scatter', 'scattered.json')
        make_task = PypeTask(
                inputs = {
                    'gathered': las_fopfn_plf,
                    'db': raw_reads_db_plf,
                },
                outputs = {
                    'scattered': scattered_plf,
                },
                parameters = {
                    'db_prefix': 'raw_reads',
                    'config': config,
                },
        )
        task = make_task(pype_tasks.task_consensus_scatter)
        wf.addTask(task)
        wf.refreshTargets(exitOnFailure=exitOnFailure)

        tasks, consensus_out = create_consensus_tasks(rawread_dir, scattered_plf)
        wf.addTasks(tasks)
        wf.refreshTargets(exitOnFailure=exitOnFailure)

        task, preads_fofn_plf = create_consensus_gather_task(os.path.join(rawread_dir, 'preads'), consensus_out)
        wf.addTask(task)

        rdir = os.path.join(rawread_dir, 'report')
        pre_assembly_report_plf = makePypeLocalFile(os.path.join(rdir, 'pre_assembly_stats.json'))
        parameters = dict(config)
        parameters['cwd'] = rdir
        make_task = PypeTask(
                inputs = {'length_cutoff_fn': length_cutoff_plf,
                          'raw_reads_db': raw_reads_db_plf,
                          'preads_fofn': preads_fofn_plf, },
                outputs = {'pre_assembly_report': pre_assembly_report_plf, },
                parameters = parameters,
        )
        task = make_task(pype_tasks.task_report_pre_assembly)
        wf.addTask(task)

        concurrent_jobs = config['cns_concurrent_jobs']
        wf.max_jobs = concurrent_jobs
        wf.refreshTargets(exitOnFailure=exitOnFailure)


    if config['target'] == 'pre-assembly':
        log.info('Quitting after stage-0 for "pre-assembly" target.')
        sys.exit(0)

    # build pread database
    if config['input_type'] == 'preads':
        preads_fofn_plf = makePypeLocalFile(os.path.join(pread_dir, 'preads-fofn-abs', os.path.basename(config['input_fofn'])))
        make_fofn_abs_task = PypeTask(inputs = {'i_fofn': rawread_fofn_plf},
                                     outputs = {'o_fofn': preads_fofn_plf},
                                     parameters = {},
        )
        fofn_abs_task = make_fofn_abs_task(pype_tasks.task_make_fofn_abs_preads)
        wf.addTasks([fofn_abs_task])
        wf.refreshTargets([fofn_abs_task])

    pdb_build_done = makePypeLocalFile( os.path.join( pread_dir, 'pdb_build_done') )
    parameters = {'work_dir': pread_dir,
                  'sge_option': config['sge_option_pda'],
                  'config_fn': input_config_fn,
                  'config': config}

    run_jobs = makePypeLocalFile(os.path.join(pread_dir, 'run_jobs.sh'))
    preads_db = makePypeLocalFile(os.path.join(pread_dir, 'preads.db')) # Also .preads.*, of course.
    make_build_pdb_task  = PypeTask(inputs = {'preads_fofn': preads_fofn_plf },
                                    outputs = {'pdb_build_done': pdb_build_done,
                                               'preads_db': preads_db,
                                               'run_jobs': run_jobs,
                                    },
                                    parameters = parameters,
    )
    build_pdb_task = make_build_pdb_task(pype_tasks.task_build_pdb)

    wf.addTasks([build_pdb_task])
    wf.refreshTargets([pdb_build_done])


    preads_nblock = support.get_nblock(fn(preads_db))
    #### run daligner
    config['sge_option_da'] = config['sge_option_pda']

    scattered_plf = os.path.join(pread_dir, 'daligner-scatter', 'scattered.json')
    make_daligner_scatter = PypeTask(
            inputs = {
                'run_jobs_fn': run_jobs,
                'db_build_done': pdb_build_done,
            },
            outputs = {
                'scatter_fn': scattered_plf,
            },
            parameters = {
                'db_prefix': 'preads',
                'nblock': preads_nblock,
                'pread_aln': True,
                'config': config,
            },
    )
    task = make_daligner_scatter(pype_tasks.task_daligner_scatter)
    wf.addTask(task)
    wf.refreshTargets(exitOnFailure=exitOnFailure)

    daligner_tasks, daligner_out = create_daligner_tasks(pread_dir, scattered_plf)
    wf.addTasks(daligner_tasks)

    p_gathered_las_plf = makePypeLocalFile(os.path.join(pread_dir, 'gathered-las', 'gathered-las.txt'))
    parameters =  {
            'nblock': preads_nblock,
    }
    make_daligner_gather = PypeTask(
                inputs = daligner_out,
                outputs =  {'gathered': p_gathered_las_plf},
                parameters = parameters,
    )
    check_p_da_task = make_daligner_gather(pype_tasks.task_daligner_gather)
    wf.addTask(check_p_da_task)
    wf.refreshTargets(exitOnFailure=exitOnFailure)

    # Merge .las files.
    config['sge_option_la'] = config['sge_option_pla']
    scattered_plf = os.path.join(pread_dir, 'merge-scatter', 'scattered.json')
    make_task = PypeTask(
            inputs = {
                'run_jobs': run_jobs,
                'gathered_las': p_gathered_las_plf,
            },
            outputs = {
                'scattered': scattered_plf,
            },
            parameters = {
                'db_prefix': 'preads',
                'config': config,
            },
    )
    task = make_task(pype_tasks.task_merge_scatter)
    wf.addTask(task)
    wf.refreshTargets(exitOnFailure=exitOnFailure)

    merge_tasks, p_ids_merged_las = create_merge_tasks(pread_dir, scattered_plf)
    wf.addTasks(merge_tasks)
    task, las_fofn_plf, las_fopfn_plf = create_merge_gather_task(os.path.join(pread_dir, 'merge-gather'), p_ids_merged_las)
    wf.addTask(task)

    concurrent_jobs = config['ovlp_concurrent_jobs']
    wf.max_jobs = concurrent_jobs

    wf.refreshTargets(exitOnFailure=exitOnFailure)


    db2falcon_dir = os.path.join(pread_dir, 'db2falcon')
    db2falcon_done = makePypeLocalFile(os.path.join(db2falcon_dir, 'db2falcon_done'))
    preads4falcon_plf = makePypeLocalFile(os.path.join(db2falcon_dir, 'preads4falcon.fasta'))
    make_run_db2falcon = PypeTask(
               inputs = {'las_fofn_plf': las_fofn_plf,
                         'preads_db': preads_db,
                        },
               outputs =  {'db2falcon_done': db2falcon_done,
                           'preads4falcon': preads4falcon_plf,
                          },
               parameters = {'wd': db2falcon_dir,
                             'config': config,
                             'sge_option': config['sge_option_fc'],
                            },
    )
    wf.addTask(make_run_db2falcon(pype_tasks.task_run_db2falcon))

    falcon_asm_done = makePypeLocalFile( os.path.join(falcon_asm_dir, 'falcon_asm_done'))
    make_run_falcon_asm = PypeTask(
               inputs = {'db2falcon_done': db2falcon_done, 'db_file': preads_db,
                         'preads4falcon': preads4falcon_plf,
                         'las_fofn': las_fofn_plf,
                        },
               outputs =  {'falcon_asm_done': falcon_asm_done},
               parameters = {'wd': falcon_asm_dir,
                             'config': config,
                             'pread_dir': pread_dir,
                             'sge_option': config['sge_option_fc'],
               },
    )
    wf.addTask(make_run_falcon_asm(pype_tasks.task_run_falcon_asm))
    wf.refreshTargets()

    return falcon_asm_done


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
