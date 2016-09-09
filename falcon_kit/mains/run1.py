from .. import run_support as support
from .. import bash
from ..util.system import only_these_symlinks
from falcon_kit import stats_preassembly
from pypeflow.pwatcher_bridge import PypeProcWatcherWorkflow, MyFakePypeThreadTaskBase
from pypeflow.data import PypeLocalFile, makePypeLocalFile, fn
from pypeflow.task import PypeTask
import argparse
import collections
import glob
import json
import os
import re
import sys
import time


fc_run_logger = None

def remove(*fns):
    for fn in fns:
        if os.path.exists(fn):
            os.remove(fn)
        assert not os.path.exists(fn)

def system(call, check=False):
    fc_run_logger.debug('$(%s)' %repr(call))
    rc = os.system(call)
    msg = "Call %r returned %d." % (call, rc)
    if rc:
        fc_run_logger.warning(msg)
        if check:
            raise Exception(msg)
    else:
        fc_run_logger.debug(msg)
    return rc

def task_make_fofn_abs_raw(self):
    #script_fn = 'noop.sh'
    #open(script_fn, 'w').write('echo NOOP raw')
    #self.generated_script_fn = script_fn
    support.make_fofn_abs(self.i_fofn.path, self.o_fofn.path)

def task_make_fofn_abs_preads(self):
    #script_fn = 'noop.sh'
    #open(script_fn, 'w').write('echo NOOP preads')
    #self.generated_script_fn = script_fn
    support.make_fofn_abs(self.i_fofn.path, self.o_fofn.path)

def task_build_rdb(self):
    input_fofn_fn = fn(self.input_fofn)
    job_done = fn(self.rdb_build_done)
    db = fn(self.raw_reads_db)
    run_jobs = fn(self.run_jobs)
    remove(job_done, db, run_jobs)
    work_dir = self.parameters["work_dir"]
    config = self.parameters["config"]

    script_fn = os.path.join( work_dir, "prepare_rdb.sh" )
    args = {
        'input_fofn_fn': input_fofn_fn,
        'config': config,
        'job_done': job_done,
        'script_fn': script_fn,
        'run_jobs_fn': run_jobs,
    }
    support.build_rdb(**args)
    self.generated_script_fn = script_fn

def task_build_pdb(self):  #essential the same as build_rdb() but the subtle differences are tricky to consolidate to one function
    input_fofn_fn = fn(self.pread_fofn)
    job_done = fn(self.pdb_build_done)
    db = fn(self.preads_db)
    run_jobs = fn(self.run_jobs)
    remove(job_done, db, run_jobs)
    work_dir = self.parameters["work_dir"]
    config = self.parameters["config"]

    script_fn = os.path.join( work_dir, "prepare_pdb.sh" )
    args = {
        'input_fofn_fn': input_fofn_fn,
        'config': config,
        'job_done': job_done,
        'script_fn': script_fn,
        'run_jobs_fn': run_jobs,
    }
    support.build_pdb(**args)
    self.generated_script_fn = script_fn

def task_run_db2falcon(self):
    wd = self.parameters["wd"]
    #self.p_merge_done
    job_done = fn(self.db2falcon_done)
    config = self.parameters["config"]
    script_dir = os.path.join(wd)
    script_fn = os.path.join(script_dir ,"run_db2falcon.sh")
    args = {
        'config': config,
        'job_done': job_done,
        'script_fn': script_fn,
    }
    support.run_db2falcon(**args)
    self.generated_script_fn = script_fn

def task_run_falcon_asm(self):
    wd = self.parameters["wd"]
    #self.db2falcon_done
    db_file = fn(self.db_file)
    job_done = fn(self.falcon_asm_done)
    config = self.parameters["config"]
    pread_dir = self.parameters["pread_dir"]
    script_dir = os.path.join( wd )
    script_fn =  os.path.join( script_dir ,"run_falcon_asm.sh" )
    # Generate las.fofn in run-dir.
    system('cd {}; find {}/m_*/ -name "preads.*.las" >| las.fofn'.format(wd, pread_dir))
    las_fofn_fn = 'las.fofn'
    args = {
        'las_fofn_fn': las_fofn_fn,
        'preads4falcon_fasta_fn': os.path.join(pread_dir, 'preads4falcon.fasta'),
        'db_file_fn': db_file,
        'config': config,
        'job_done': job_done,
        'script_fn': script_fn,
    }
    support.run_falcon_asm(**args)
    self.generated_script_fn = script_fn

def task_report_pre_assembly(self):
    # TODO(CD): Bashify this, in case it is slow.
    i_raw_reads_db_fn = fn(self.raw_reads_db)
    i_preads_fofn_fn = fn(self.preads_fofn)
    i_length_cutoff_fn = fn(self.length_cutoff_fn)
    o_json_fn = fn(self.pre_assembly_report)
    cfg = self.parameters
    genome_length = int(cfg.get('genome_size', 0)) # different name in falcon
    length_cutoff = int(cfg['length_cutoff'])
    length_cutoff = support.get_length_cutoff(length_cutoff, i_length_cutoff_fn)
    kwds = {
        'i_raw_reads_db_fn': i_raw_reads_db_fn,
        'i_preads_fofn_fn': i_preads_fofn_fn,
        'genome_length': genome_length,
        'length_cutoff': length_cutoff,
    }
    fc_run_logger.info('Report inputs: {}'.format(repr(kwds)))
    report_dict = stats_preassembly.calc_dict(**kwds)
    content = json.dumps(report_dict, sort_keys=True, indent=4, separators=(',', ': '))
    fc_run_logger.info('Report stats:\n{}'.format(content))
    open(o_json_fn, 'w').write(content)

def task_run_daligner(self):
    job_done = fn(self.job_done)
    daligner_script = self.parameters["daligner_script"]
    job_uid = self.parameters["job_uid"]
    cwd = self.parameters["cwd"]
    mkdir(cwd)
    db_prefix = self.parameters["db_prefix"]
    config = self.parameters["config"]
    script_dir = os.path.join(cwd)
    script_fn =  os.path.join(script_dir , "rj_%s.sh" % (job_uid))
    args = {
        'daligner_script': daligner_script,
        'db_prefix': db_prefix,
        'config': config,
        'job_done': job_done,
        'script_fn': script_fn,
    }
    support.run_daligner(**args)
    self.generated_script_fn = script_fn

def task_run_las_merge(self):
    script = self.parameters["merge_script"]
    job_id = self.parameters["job_id"]
    cwd = self.parameters["cwd"]
    job_done = fn(self.job_done)
    config = self.parameters["config"]

    script_dir = os.path.join( cwd )
    script_fn =  os.path.join( script_dir , "rp_%05d.sh" % (job_id))
    args = {
        'script': script,
        'config': config,
        'job_done': job_done,
        'script_fn': script_fn,
    }
    support.run_las_merge(**args)
    self.generated_script_fn = script_fn

def task_run_consensus(self):
    out_file_fn = fn(self.out_file)
    job_id = self.parameters["job_id"]
    cwd = self.parameters["cwd"]
    config = self.parameters["config"]
    prefix = self.parameters["prefix"]
    script_dir = os.path.join( cwd )
    job_done = os.path.join( cwd, "c_%05d_done" % job_id )
    script_fn =  os.path.join( script_dir , "c_%05d.sh" % (job_id))
    db_fn = os.path.abspath('{cwd}/../{prefix}'.format(**locals()))
    las_fn = os.path.abspath('{cwd}/../m_{job_id:05d}/{prefix}.{job_id}.las'.format(**locals()))
    args = {
        'db_fn': db_fn,
        'las_fn': las_fn,
        'out_file_fn': out_file_fn,
        'config': config,
        'job_done': job_done,
        'script_fn': script_fn,
    }
    support.run_consensus(**args)
    self.generated_script_fn = script_fn

def mkdir(d):
    if not os.path.isdir(d):
        os.makedirs(d)

def task_daligner_gather(self):
    da_done = fn(self.da_done)
    main_dir = os.path.dirname(da_done)
    out_dict = self.inputDataObjs
    nblock = self.parameters['nblock']
    fc_run_logger.debug('nblock=%d, out_dir:\n%s'%(nblock, out_dict))

    # Create m_* dirs.
    for block in xrange(1, nblock+1):
        mdir = os.path.join(main_dir, 'm_%05d' %block) # By convention. pbsmrtpipe works differently.
        mkdir(mdir)
        # TODO: Remove existing symlinks?
    job_rundirs = [os.path.dirname(fn(dal_done)) for dal_done in out_dict.values()]

    # Symlink all daligner *.las.
    links = collections.defaultdict(list)
    for block, las_path in support.daligner_gather_las(job_rundirs):
            mdir = os.path.join(main_dir, 'm_%05d' %block) # By convention. pbsmrtpipe works differently.
            #las_path = os.path.relpath(las_path, mdir)
            links[mdir].append(las_path)
    only_these_symlinks(links)
    system("touch %s" %da_done)

def create_daligner_tasks(run_jobs_fn, wd, db_prefix, rdb_build_done, config, pread_aln=False):
    tasks = []
    tasks_out = {}
    skip_checks = config.get('skip_checks')
    fc_run_logger.info('Skip LAcheck after daligner? {}'.format(skip_checks))
    for job_uid, script in bash.scripts_daligner(run_jobs_fn, db_prefix, rdb_build_done, pread_aln, skip_checks):
        run_dir = "job_%s" %job_uid
        cwd = os.path.join(wd, run_dir)
        job_done_fn = os.path.abspath(os.path.join(cwd, "job_%s_done" %job_uid))
        job_done = makePypeLocalFile(job_done_fn)
        parameters =  {"daligner_script": script,
                       "cwd": cwd,
                       "job_uid": job_uid,
                       "config": config,
                       "sge_option": config["sge_option_da"],
                       "db_prefix": db_prefix}
        make_daligner_task = PypeTask(inputs = {"rdb_build_done": rdb_build_done},
                                      outputs = {"job_done": job_done},
                                      parameters = parameters,
                                      TaskType = MyFakePypeThreadTaskBase,
                                      URL = "task://localhost/d_%s_%s" %(job_uid, db_prefix))
        daligner_task = make_daligner_task(task_run_daligner)
        tasks.append(daligner_task)
        tasks_out[ "ajob_%s" % job_uid ] = job_done
    return tasks, tasks_out

def create_merge_tasks(run_jobs_fn, wd, db_prefix, input_dep, config):
    merge_tasks = []
    merge_out = {}
    p_ids_merge_job_done = [] # for consensus

    merge_scripts = bash.scripts_merge(config, db_prefix, run_jobs_fn)
    for p_id, merge_script in merge_scripts:
        job_done = makePypeLocalFile(os.path.abspath("%s/m_%05d/m_%05d_done" % (wd, p_id, p_id)))
        parameters =  {"merge_script": merge_script,
                       "cwd": os.path.join(wd, "m_%05d" % p_id),
                       "job_id": p_id,
                       "sge_option": config["sge_option_la"],
                       "config": config}
        make_merge_task = PypeTask(inputs = {"input_dep": input_dep},
                                   outputs = {"job_done": job_done},
                                   parameters = parameters,
                                   TaskType = MyFakePypeThreadTaskBase,
                                   URL = "task://localhost/m_%05d_%s" % (p_id, db_prefix))
        merge_task = make_merge_task(task_run_las_merge)
        merge_out["mjob_%d" % p_id] = job_done
        merge_tasks.append(merge_task)
        p_ids_merge_job_done.append((p_id, job_done))
    return merge_tasks, merge_out, p_ids_merge_job_done

def create_consensus_tasks(wd, db_prefix, config, p_ids_merge_job_done):
    consensus_tasks = []
    consensus_out ={}
    # Unlike the merge tasks, consensus occurs in a single directory.
    rdir = os.path.join(wd, 'preads')
    mkdir(rdir)
    for p_id, job_done in p_ids_merge_job_done:
        out_file = makePypeLocalFile(os.path.abspath("%s/preads/out.%05d.fasta" % (wd, p_id)))
        out_done = makePypeLocalFile(os.path.abspath("%s/preads/c_%05d_done" % (wd, p_id)))
        parameters =  {"cwd": rdir,
                       "job_id": p_id,
                       "prefix": db_prefix,
                       "sge_option": config["sge_option_cns"],
                       "config": config}
        make_c_task = PypeTask(inputs = {"job_done": job_done},
                               outputs = {"out_file": out_file, "out_done": out_done},
                               parameters = parameters,
                               TaskType = MyFakePypeThreadTaskBase,
                               URL = "task://localhost/ct_%05d" % p_id)
        c_task = make_c_task(task_run_consensus)
        consensus_tasks.append(c_task)
        consensus_out["cjob_%d" % p_id] = out_done
    return consensus_tasks, consensus_out


def main1(prog_name, input_config_fn, logger_config_fn=None):
    global fc_run_logger
    fc_run_logger = support.setup_logger(logger_config_fn)

    fc_run_logger.info("fc_run started with configuration %s", input_config_fn)
    try:
        config = support.get_dict_from_old_falcon_cfg(support.parse_config(input_config_fn))
    except Exception:
        fc_run_logger.exception('Failed to parse config "{}".'.format(input_config_fn))
        raise
    input_fofn_plf = makePypeLocalFile(config["input_fofn"])
    #Workflow = PypeProcWatcherWorkflow
    wf = PypeProcWatcherWorkflow(job_type=config['job_type'], watcher_type=config['watcher_type'], watcher_directory=config['watcher_directory'])
    run(wf, config,
            input_fofn_plf=input_fofn_plf,
            setNumThreadAllowed=PypeProcWatcherWorkflow.setNumThreadAllowed)


def run(wf, config,
        input_fofn_plf,
        setNumThreadAllowed,
        ):
    """
    Preconditions (for now):
    * fc_run_logger
    * run_support.logger
    """
    rawread_dir = os.path.abspath("./0-rawreads")
    pread_dir = os.path.abspath("./1-preads_ovl")
    falcon_asm_dir  = os.path.abspath("./2-asm-falcon")
    script_dir = os.path.abspath("./scripts")
    sge_log_dir = os.path.abspath("./sge_log")

    for d in (rawread_dir, pread_dir, falcon_asm_dir, script_dir, sge_log_dir):
        support.make_dirs(d)

    exitOnFailure=config['stop_all_jobs_on_failure'] # only matter for parallel jobs
    concurrent_jobs = config["pa_concurrent_jobs"]
    setNumThreadAllowed(concurrent_jobs, concurrent_jobs)

    rawread_fofn_plf = makePypeLocalFile(os.path.join(rawread_dir, os.path.basename(config["input_fofn"])))
    make_fofn_abs_task = PypeTask(inputs = {"i_fofn": input_fofn_plf},
                                  outputs = {"o_fofn": rawread_fofn_plf},
                                  parameters = {},
                                  TaskType = MyFakePypeThreadTaskBase)
    fofn_abs_task = make_fofn_abs_task(task_make_fofn_abs_raw)

    wf.addTasks([fofn_abs_task])
    wf.refreshTargets([fofn_abs_task])

    if config["input_type"] == "raw":
        #### import sequences into daligner DB
        sleep_done = makePypeLocalFile( os.path.join( rawread_dir, "sleep_done") )
        rdb_build_done = makePypeLocalFile( os.path.join( rawread_dir, "rdb_build_done") )
        run_jobs = makePypeLocalFile( os.path.join( rawread_dir, "run_jobs.sh") )
        parameters = {"work_dir": rawread_dir,
                      "sge_option": config["sge_option_da"],
                      "config": config}

        length_cutoff_plf = makePypeLocalFile(os.path.join(rawread_dir, "length_cutoff"))
        raw_reads_db_plf = makePypeLocalFile(os.path.join(rawread_dir, "%s.db" % "raw_reads"))
        make_build_rdb_task = PypeTask(inputs = {"input_fofn": rawread_fofn_plf},
                                      outputs = {"rdb_build_done": rdb_build_done,
                                                 "raw_reads_db": raw_reads_db_plf,
                                                 "length_cutoff": length_cutoff_plf,
                                                 "run_jobs": run_jobs,
                                      },
                                      parameters = parameters,
                                      TaskType = MyFakePypeThreadTaskBase)
        build_rdb_task = make_build_rdb_task(task_build_rdb)

        wf.addTasks([build_rdb_task])
        wf.refreshTargets([rdb_build_done])

        raw_reads_nblock = support.get_nblock(fn(raw_reads_db_plf))
        #### run daligner
        daligner_tasks, daligner_out = create_daligner_tasks(fn(run_jobs), rawread_dir, "raw_reads", rdb_build_done, config)

        wf.addTasks(daligner_tasks)
        r_da_done = makePypeLocalFile( os.path.join( rawread_dir, "da_done") )

        parameters =  {
                "nblock": raw_reads_nblock,
        }
        make_daligner_gather = PypeTask(
                   inputs = daligner_out,
                   outputs =  {"da_done":r_da_done},
                   parameters = parameters,
                   TaskType = MyFakePypeThreadTaskBase,
                   URL = "task://localhost/rda_check" )
        check_r_da_task = make_daligner_gather(task_daligner_gather)
        wf.addTask(check_r_da_task)
        wf.refreshTargets(exitOnFailure=exitOnFailure)

        merge_tasks, merge_out, p_ids_merge_job_done = create_merge_tasks(fn(run_jobs), rawread_dir, "raw_reads", r_da_done, config)
        wf.addTasks( merge_tasks )
        wf.refreshTargets(exitOnFailure=exitOnFailure)

        if config["target"] == "overlapping":
            sys.exit(0)
        consensus_tasks, consensus_out = create_consensus_tasks(rawread_dir, "raw_reads", config, p_ids_merge_job_done)
        wf.addTasks( consensus_tasks )

        r_cns_done = makePypeLocalFile( os.path.join( rawread_dir, "cns_done") )
        pread_fofn = makePypeLocalFile( os.path.join( pread_dir,  "input_preads.fofn" ) )

        @PypeTask( inputs = consensus_out,
                   outputs =  {"cns_done":r_cns_done, "pread_fofn": pread_fofn},
                   TaskType = MyFakePypeThreadTaskBase,
                   URL = "task://localhost/cns_check" )
        def check_r_cns_task(self):
            with open(fn(self.pread_fofn),  "w") as f:
                fn_list =  glob.glob("%s/preads/out*.fasta" % rawread_dir)
                fn_list.sort()
                for fa_fn in fn_list:
                    print >>f, fa_fn
            system("touch %s" % fn(self.cns_done))
        wf.addTask(check_r_cns_task)

        pre_assembly_report_plf = makePypeLocalFile(os.path.join(rawread_dir, "pre_assembly_stats.json")) #tho technically it needs pread_fofn
        make_task = PypeTask(
                inputs = {"length_cutoff_fn": length_cutoff_plf,
                          "raw_reads_db": raw_reads_db_plf,
                          "preads_fofn": pread_fofn, },
                outputs = {"pre_assembly_report": pre_assembly_report_plf, },
                parameters = config,
                TaskType = MyFakePypeThreadTaskBase,
                URL = "task://localhost/report_pre_assembly")
        task = make_task(task_report_pre_assembly)
        wf.addTask(task)

        concurrent_jobs = config["cns_concurrent_jobs"]
        setNumThreadAllowed(concurrent_jobs, concurrent_jobs)
        wf.refreshTargets(exitOnFailure=exitOnFailure)


    if config["target"] == "pre-assembly":
        log.info("Quitting after stage-0 for 'pre-assembly' target.")
        sys.exit(0)

    # build pread database
    if config["input_type"] == "preads":
        pread_fofn = makePypeLocalFile(os.path.join(pread_dir, os.path.basename(config["input_fofn"])))
        make_fofn_abs_task = PypeTask(inputs = {"i_fofn": rawread_fofn_plf},
                                     outputs = {"o_fofn": pread_fofn},
                                     parameters = {},
                                     TaskType = MyFakePypeThreadTaskBase)
        fofn_abs_task = make_fofn_abs_task(task_make_fofn_abs_preads)
        wf.addTasks([fofn_abs_task])
        wf.refreshTargets([fofn_abs_task])

    pdb_build_done = makePypeLocalFile( os.path.join( pread_dir, "pdb_build_done") )
    parameters = {"work_dir": pread_dir,
                  "sge_option": config["sge_option_pda"],
                  "config": config}

    run_jobs = makePypeLocalFile(os.path.join(pread_dir, 'run_jobs.sh'))
    preads_db = makePypeLocalFile(os.path.join(pread_dir, 'preads.db')) # Also .preads.*, of course.
    make_build_pdb_task  = PypeTask(inputs = {"pread_fofn": pread_fofn },
                                    outputs = {"pdb_build_done": pdb_build_done,
                                               "preads_db": preads_db,
                                               "run_jobs": run_jobs,
                                    },
                                    parameters = parameters,
                                    TaskType = MyFakePypeThreadTaskBase,
                                    URL = "task://localhost/build_pdb")
    build_pdb_task = make_build_pdb_task(task_build_pdb)

    wf.addTasks([build_pdb_task])
    wf.refreshTargets([pdb_build_done])


    preads_nblock = support.get_nblock(fn(preads_db))
    #### run daligner
    config["sge_option_da"] = config["sge_option_pda"]
    daligner_tasks, daligner_out = create_daligner_tasks(fn(run_jobs), pread_dir, "preads", pdb_build_done, config, pread_aln=True)
    wf.addTasks(daligner_tasks)

    p_da_done = makePypeLocalFile(os.path.join( pread_dir, "da_done"))
    parameters =  {
            "nblock": preads_nblock,
    }
    make_daligner_gather = PypeTask(
                inputs = daligner_out,
                outputs =  {"da_done":p_da_done},
                parameters = parameters,
                TaskType = MyFakePypeThreadTaskBase,
                URL = "task://localhost/pda_check" )
    check_p_da_task = make_daligner_gather(task_daligner_gather)
    wf.addTask(check_p_da_task)
    wf.refreshTargets(exitOnFailure=exitOnFailure)

    config["sge_option_la"] = config["sge_option_pla"]
    merge_tasks, merge_out, _ = create_merge_tasks(fn(run_jobs), pread_dir, "preads", p_da_done, config)
    wf.addTasks( merge_tasks )

    p_merge_done = makePypeLocalFile(os.path.join( pread_dir, "p_merge_done"))

    @PypeTask( inputs = merge_out,
               outputs =  {"p_merge_done": p_merge_done},
               TaskType = MyFakePypeThreadTaskBase,
               URL = "task://localhost/pmerge_check" )
    def check_p_merge_check_task(self):
        system("touch %s" % fn(self.p_merge_done))
    wf.addTask(check_p_merge_check_task)

    concurrent_jobs = config["ovlp_concurrent_jobs"]
    setNumThreadAllowed(concurrent_jobs, concurrent_jobs)

    wf.refreshTargets(exitOnFailure=exitOnFailure)


    db2falcon_done = makePypeLocalFile( os.path.join(pread_dir, "db2falcon_done"))
    make_run_db2falcon = PypeTask(
               inputs = {"p_merge_done": p_merge_done,},
               outputs =  {"db2falcon_done": db2falcon_done},
               parameters = {"wd": pread_dir,
                             "config": config,
                             "sge_option": config["sge_option_fc"],
                            },
               TaskType = MyFakePypeThreadTaskBase,
               URL = "task://localhost/db2falcon" )
    wf.addTask(make_run_db2falcon(task_run_db2falcon))

    falcon_asm_done = makePypeLocalFile( os.path.join( falcon_asm_dir, "falcon_asm_done") )
    make_run_falcon_asm = PypeTask(
               inputs = {"db2falcon_done": db2falcon_done, "db_file": preads_db},
               outputs =  {"falcon_asm_done": falcon_asm_done},
               parameters = {"wd": falcon_asm_dir,
                             "config": config,
                             "pread_dir": pread_dir,
                             "sge_option": config["sge_option_fc"],
               },
               TaskType = MyFakePypeThreadTaskBase,
               URL = "task://localhost/falcon_asm" )
    wf.addTask(make_run_falcon_asm(task_run_falcon_asm))
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
