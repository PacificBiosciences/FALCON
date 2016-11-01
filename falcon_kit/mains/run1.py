from .. import run_support as support
from .. import bash
from ..util.system import only_these_symlinks
from pypeflow.pwatcher_bridge import PypeProcWatcherWorkflow, MyFakePypeThreadTaskBase
from pypeflow.data import makePypeLocalFile, fn
from pypeflow.task import PypeTask
#from pypeflow.simple_pwatcher_bridge import (PypeProcWatcherWorkflow, MyFakePypeThreadTaskBase,
#        makePypeLocalFile, fn, PypeTask)
import argparse
import collections
import glob
import json
import logging
import os
import re
import sys
import time


fc_run_logger = logging.getLogger(__name__) # default, for remote tasks

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
    script_fn = 'noop.sh'
    open(script_fn, 'w').write('echo NOOP raw')
    self.generated_script_fn = script_fn
    mkdir(os.path.dirname(fn(self.o_fofn)))
    support.make_fofn_abs(fn(self.i_fofn), fn(self.o_fofn))

def task_make_fofn_abs_preads(self):
    script_fn = 'noop.sh'
    open(script_fn, 'w').write('echo NOOP preads')
    self.generated_script_fn = script_fn
    mkdir(os.path.dirname(fn(self.o_fofn)))
    support.make_fofn_abs(fn(self.i_fofn), fn(self.o_fofn))

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
    input_fofn_fn = fn(self.preads_fofn)
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
    mkdir(wd)
    #self.p_merge_done
    job_done = fn(self.db2falcon_done)
    preads4falcon_fn = fn(self.preads4falcon)
    preads_db = fn(self.preads_db)
    config = self.parameters["config"]
    script_dir = os.path.join(wd)
    script_fn = os.path.join(script_dir ,"run_db2falcon.sh")
    args = {
        'config': config,
        'job_done': job_done,
        'script_fn': script_fn,
        'preads4falcon_fn': preads4falcon_fn,
        'preads_db': preads_db,
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
    preads4falcon_fn = fn(self.preads4falcon)
    script_dir = os.path.join( wd )
    script_fn =  os.path.join( script_dir ,"run_falcon_asm.sh" )
    # Generate las.fofn in run-dir.
    system('cd {}; find {}/m_*/ -name "preads.*.las" >| las.fofn'.format(wd, pread_dir))
    las_fofn_fn = 'las.fofn'
    args = {
        'las_fofn_fn': las_fofn_fn,
        'preads4falcon_fasta_fn': preads4falcon_fn,
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
    # Update length_cutoff if auto-calc (when length_cutoff is negative).
    # i_length_cutoff_fn was created long ago, so no filesystem issues.
    length_cutoff = support.get_length_cutoff(length_cutoff, i_length_cutoff_fn)
    cwd = self.parameters['cwd']
    mkdir(cwd)
    script_fn = os.path.join(cwd , 'run_report_pre_assembly.sh')
    job_done = os.path.join(cwd, 'report_pa_done')
    kwds = {
        'i_raw_reads_db_fn': i_raw_reads_db_fn,
        'i_preads_fofn_fn': i_preads_fofn_fn,
        'genome_length': genome_length,
        'length_cutoff': length_cutoff,
        'o_json_fn': o_json_fn,
        'job_done': job_done,
        'script_fn': script_fn,
    }
    fc_run_logger.info('Report inputs: {}'.format(repr(kwds)))
    support.run_report_pre_assembly(**kwds)
    self.generated_script_fn = script_fn

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

def read_gathered_las(path):
    """Return dict of block->[las_paths].
    For now, these are ws separated on each line of input.
    """
    result = collections.defaultdict(list)
    with open(path) as ifs:
        for line in ifs:
            block, las_path = line.split()
            result[int(block)].append(las_path)
    #fc_run_logger.warning('path={!r}, result={}'.format(
    #    path, pprint.pformat(result)))
    return result

def task_run_las_merge(self):
    gathered_las_fn = fn(self.gathered_las)
    script = self.parameters["merge_script"]
    job_id = self.parameters["job_id"] # aka "block"
    cwd = self.parameters["cwd"]
    mkdir(cwd)

    gathered_dict = read_gathered_las(gathered_las_fn)
    las_paths = gathered_dict[job_id]
    for las_path in las_paths:
        src = os.path.relpath(las_path, cwd)
        tgt = os.path.join(cwd, os.path.basename(las_path))
        fc_run_logger.debug('symlink {!r} -> {!r}'.format(src, tgt))
        os.symlink(src, tgt)

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
    merge_job_done = fn(self.job_done)
    out_file_fn = fn(self.out_file)
    out_done = fn(self.out_done)
    job_id = self.parameters["job_id"]
    cwd = self.parameters["cwd"]
    config = self.parameters["config"]
    prefix = self.parameters["prefix"]
    script_dir = os.path.join( cwd )
    script_fn =  os.path.join( script_dir , "c_%05d.sh" % (job_id))
    db_fn = os.path.abspath('{cwd}/../../{prefix}'.format(**locals())) # ASSUMING 2-levels deep
    merge_job_dir = os.path.dirname(merge_job_done)
    # by convention, we assume the name of the .las file
    las_fn = os.path.abspath('{merge_job_dir}/{prefix}.{job_id}.las'.format(**locals()))
    args = {
        'db_fn': db_fn,
        'las_fn': las_fn,
        'out_file_fn': out_file_fn,
        'config': config,
        'job_done': out_done,
        'script_fn': script_fn,
    }
    support.run_consensus(**args)
    self.generated_script_fn = script_fn

def mkdir(d):
    if not os.path.isdir(d):
        os.makedirs(d)

def task_daligner_gather(self):
    """Find all .las leaves so far.
    """
    out_dict = self.inputs
    gathered_fn = fn(self.gathered)
    nblock = self.parameters['nblock']
    fc_run_logger.debug('nblock=%d, out_dir:\n%s'%(nblock, out_dict))
    job_rundirs = [os.path.dirname(fn(dal_done)) for dal_done in out_dict.values()]
    wdir = os.path.dirname(gathered_fn)
    mkdir(wdir)
    with open(gathered_fn, 'w') as ofs:
        for block, las_path in support.daligner_gather_las(job_rundirs):
            ofs.write('{} {}\n'.format(block, las_path))

    # Because we need a script always, for now.
    script_fn = os.path.join(wdir, 'noop.sh')
    open(script_fn, 'w').write('echo NOOP raw')
    self.generated_script_fn = script_fn
def old_task_daligner_gather(self):
    """Create symlinks in the m_* directories.
    Kinda messy, and it must run in the base-dir for now.
    """
    da_done = fn(self.da_done)
    main_dir = os.path.dirname(da_done)
    out_dict = self.inputs
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

    # Because we need a script always, for now.
    script_fn = 'noop.sh'
    open(script_fn, 'w').write('echo NOOP raw')
    self.generated_script_fn = script_fn

def create_daligner_tasks(run_jobs_fn, wd, db_prefix, rdb_build_done, nblock, config, pread_aln=False):
    tasks = []
    tasks_out = {}
    skip_checks = config.get('skip_checks')
    fc_run_logger.info('Skip LAcheck after daligner? {}'.format(skip_checks))
    for job_uid, script in bash.scripts_daligner(run_jobs_fn, db_prefix, rdb_build_done, nblock, pread_aln, skip_checks):
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

def create_merge_tasks(run_jobs_fn, wd, db_prefix, gathered_las_plf, config):
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
        make_merge_task = PypeTask(inputs = {"gathered_las": gathered_las_plf,
                                   },
                                   outputs = {"job_done": job_done,
                                   },
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
    fasta_plfs = []
    for p_id, job_done in p_ids_merge_job_done:
        cns_label = 'cns_%05d' %p_id
        rdir = os.path.join(wd, 'preads', cns_label)
        mkdir(rdir)
        out_done = makePypeLocalFile(os.path.abspath("%s/%s_done" % (rdir, cns_label)))
        out_file = makePypeLocalFile(os.path.abspath("%s/%s.fasta" % (rdir, cns_label)))
        fasta_plfs.append(out_file)
        parameters =  {"cwd": rdir,
                       "job_id": p_id,
                       "prefix": db_prefix,
                       "sge_option": config["sge_option_cns"],
                       "config": config}
        make_c_task = PypeTask(inputs = {"job_done": job_done},
                               outputs = {"out_file": out_file, "out_done": out_done},
                               parameters = parameters,
                               TaskType = MyFakePypeThreadTaskBase,
                               URL = "task://localhost/%s" % cns_label)
        c_task = make_c_task(task_run_consensus)
        consensus_tasks.append(c_task)
        #consensus_out["cjob_%d" % p_id] = out_done
        consensus_out["cjob_%d" % p_id] = out_file

    r_cns_done_plf = makePypeLocalFile(os.path.join(wd, 'preads', "cns_done"))
    preads_fofn_plf = makePypeLocalFile(os.path.join(wd, 'preads', "input_preads.fofn"))

    make_check_r_cns_task = PypeTask(
                inputs = consensus_out,
                outputs =  {"cns_done":r_cns_done_plf, "preads_fofn": preads_fofn_plf},
                TaskType = MyFakePypeThreadTaskBase,
                URL = "task://localhost/cns_check" )
    consensus_tasks.append(make_check_r_cns_task(check_r_cns_task))
    return consensus_tasks, preads_fofn_plf

# PypeTask functions now need to be module-level.
def check_r_cns_task(self):
    with open(fn(self.preads_fofn),  "w") as f:
        for fa_fn in sorted(fn(plf) for plf in self.inputs.values()):
            print >>f, fa_fn
    wdir = os.path.dirname(fn(self.cns_done))
    #mkdir(wdir) We SHOULD need this! TODO
    system("touch %s" % fn(self.cns_done))
    script_fn = os.path.join(wdir, 'noop.sh')
    open(script_fn, 'w').write('echo NOOP raw')
    self.generated_script_fn = script_fn
def check_p_merge_check_task(self):
    wdir = os.path.dirname(fn(self.p_merge_done))
    mkdir(wdir)
    system("touch %s" % fn(self.p_merge_done))
    script_fn = os.path.join(wdir, 'noop.sh')
    open(script_fn, 'w').write('echo NOOP raw')
    self.generated_script_fn = script_fn


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
    wf = PypeProcWatcherWorkflow(job_type=config['job_type'],
            job_queue=config['job_queue'],
            watcher_type=config['pwatcher_type'],
            watcher_directory=config['pwatcher_directory'])
    run(wf, config,
            os.path.abspath(input_config_fn),
            input_fofn_plf=input_fofn_plf,
            setNumThreadAllowed=PypeProcWatcherWorkflow.setNumThreadAllowed)


def run(wf, config,
        input_config_fn,
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

    rawread_fofn_plf = makePypeLocalFile(os.path.join(rawread_dir, 'raw-fofn-abs', os.path.basename(config["input_fofn"])))
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
                      "config_fn": input_config_fn,
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
        daligner_tasks, daligner_out = create_daligner_tasks(fn(run_jobs), rawread_dir, "raw_reads", rdb_build_done,
                nblock=raw_reads_nblock, config=config)

        wf.addTasks(daligner_tasks)
        r_gathered_las_plf = makePypeLocalFile( os.path.join( rawread_dir, 'raw_gather', 'gathered_las.txt') )

        parameters =  {
                "nblock": raw_reads_nblock,
        }
        make_daligner_gather = PypeTask(
                   inputs = daligner_out,
                   outputs =  {"gathered": r_gathered_las_plf},
                   parameters = parameters,
                   TaskType = MyFakePypeThreadTaskBase,
                   URL = "task://localhost/rda_check" )
        check_r_da_task = make_daligner_gather(task_daligner_gather)
        wf.addTask(check_r_da_task)
        wf.refreshTargets(exitOnFailure=exitOnFailure)

        merge_tasks, merge_out, p_ids_merge_job_done = create_merge_tasks(fn(run_jobs), rawread_dir, "raw_reads", r_gathered_las_plf, config)
        wf.addTasks( merge_tasks )
        wf.refreshTargets(exitOnFailure=exitOnFailure)

        if config["target"] == "overlapping":
            sys.exit(0)
        consensus_tasks, preads_fofn_plf = create_consensus_tasks(rawread_dir, "raw_reads", config, p_ids_merge_job_done)
        wf.addTasks( consensus_tasks )

        rdir = os.path.join(rawread_dir, 'report')
        pre_assembly_report_plf = makePypeLocalFile(os.path.join(rdir, "pre_assembly_stats.json"))
        parameters = dict(config)
        parameters['cwd'] = rdir
        make_task = PypeTask(
                inputs = {"length_cutoff_fn": length_cutoff_plf,
                          "raw_reads_db": raw_reads_db_plf,
                          "preads_fofn": preads_fofn_plf, },
                outputs = {"pre_assembly_report": pre_assembly_report_plf, },
                parameters = parameters,
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
        preads_fofn_plf = makePypeLocalFile(os.path.join(pread_dir, 'preads-fofn-abs', os.path.basename(config["input_fofn"])))
        make_fofn_abs_task = PypeTask(inputs = {"i_fofn": rawread_fofn_plf},
                                     outputs = {"o_fofn": preads_fofn_plf},
                                     parameters = {},
                                     TaskType = MyFakePypeThreadTaskBase)
        fofn_abs_task = make_fofn_abs_task(task_make_fofn_abs_preads)
        wf.addTasks([fofn_abs_task])
        wf.refreshTargets([fofn_abs_task])

    pdb_build_done = makePypeLocalFile( os.path.join( pread_dir, "pdb_build_done") )
    parameters = {"work_dir": pread_dir,
                  "sge_option": config["sge_option_pda"],
                  "config_fn": input_config_fn,
                  "config": config}

    run_jobs = makePypeLocalFile(os.path.join(pread_dir, 'run_jobs.sh'))
    preads_db = makePypeLocalFile(os.path.join(pread_dir, 'preads.db')) # Also .preads.*, of course.
    make_build_pdb_task  = PypeTask(inputs = {"preads_fofn": preads_fofn_plf },
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
    daligner_tasks, daligner_out = create_daligner_tasks(fn(run_jobs), pread_dir, "preads", pdb_build_done,
                nblock=preads_nblock, config=config, pread_aln=True)
    wf.addTasks(daligner_tasks)

    p_gathered_las_plf = makePypeLocalFile(os.path.join(pread_dir, 'gathered-las', 'gathered-las.txt'))
    parameters =  {
            "nblock": preads_nblock,
    }
    make_daligner_gather = PypeTask(
                inputs = daligner_out,
                outputs =  {"gathered": p_gathered_las_plf},
                parameters = parameters,
                TaskType = MyFakePypeThreadTaskBase,
                URL = "task://localhost/pda_check" )
    check_p_da_task = make_daligner_gather(task_daligner_gather)
    wf.addTask(check_p_da_task)
    wf.refreshTargets(exitOnFailure=exitOnFailure)

    config["sge_option_la"] = config["sge_option_pla"]
    merge_tasks, merge_out, _ = create_merge_tasks(fn(run_jobs), pread_dir, 'preads-merge', p_gathered_las_plf, config)
    wf.addTasks( merge_tasks )

    p_merge_done = makePypeLocalFile(os.path.join( pread_dir, 'preads-merge', 'p_merge_done'))

    make_check_p_merge_task = PypeTask( inputs = merge_out,
               outputs =  {"p_merge_done": p_merge_done},
               TaskType = MyFakePypeThreadTaskBase,
               URL = "task://localhost/pmerge_check" )
    wf.addTask(make_check_p_merge_task(check_p_merge_check_task))

    concurrent_jobs = config["ovlp_concurrent_jobs"]
    setNumThreadAllowed(concurrent_jobs, concurrent_jobs)

    wf.refreshTargets(exitOnFailure=exitOnFailure)


    db2falcon_dir = os.path.join(pread_dir, 'db2falcon')
    db2falcon_done = makePypeLocalFile(os.path.join(db2falcon_dir, 'db2falcon_done'))
    preads4falcon_plf = makePypeLocalFile(os.path.join(db2falcon_dir, 'preads4falcon.fasta'))
    make_run_db2falcon = PypeTask(
               inputs = {"p_merge_done": p_merge_done,
                         "preads_db": preads_db,
                        },
               outputs =  {"db2falcon_done": db2falcon_done,
                           "preads4falcon": preads4falcon_plf,
                          },
               parameters = {"wd": db2falcon_dir,
                             "config": config,
                             "sge_option": config["sge_option_fc"],
                            },
               TaskType = MyFakePypeThreadTaskBase,
               URL = "task://localhost/db2falcon" )
    wf.addTask(make_run_db2falcon(task_run_db2falcon))

    falcon_asm_done = makePypeLocalFile( os.path.join(falcon_asm_dir, 'falcon_asm_done'))
    make_run_falcon_asm = PypeTask(
               inputs = {"db2falcon_done": db2falcon_done, "db_file": preads_db,
                         "preads4falcon": preads4falcon_plf,
                        },
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
