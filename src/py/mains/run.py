from .. import run_support as support
from ..functional import get_daligner_job_descriptions, get_script_xformer
from pypeflow.data import PypeLocalFile, makePypeLocalFile, fn
from pypeflow.task import PypeTask, PypeThreadTaskBase, PypeTaskBase
from pypeflow.controller import PypeWorkflow, PypeThreadWorkflow
from falcon_kit.FastaReader import FastaReader
import glob
import os
import re
import sys
import time


wait_time = 5
fc_run_logger = None

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

def _qsub_script(job_data, specific):
        script_fn = job_data["script_fn"]
        job_name = job_data["job_name"]
        cwd = job_data["cwd"]
        sge_option = job_data["sge_option"]
        sge_cmd="qsub -N {job_name} {sge_option} -o {cwd}/sge_log {specific}\
                 -S /bin/bash {script}".format(job_name=job_name,
                                               cwd=os.getcwd(),
                                               specific=specific,
                                               sge_option=sge_option,
                                               script=script_fn)
        system(sge_cmd, check=True)

def _run_script_sge(job_data):
    specific = '-j y'
    _qsub_script(job_data, specific)

def _run_script_torque(job_data):
    # See https://github.com/PacificBiosciences/FALCON/pull/227
    specific = '-j oe'
    _qsub_script(job_data, specific)

def _run_script_slurm(job_data):
        script_fn = job_data["script_fn"]
        job_name = job_data["job_name"]
        cwd = job_data["cwd"]
        sge_option = job_data["sge_option"]
        with open(script_fn, 'r') as original: data = original.read()
        with open(script_fn, 'w') as modified: modified.write("#!/bin/sh" + "\n" + data)
        sge_cmd="sbatch -J {job_name} {sge_option} {script}".format(job_name=job_name, cwd=os.getcwd(),sge_option=sge_option, script=script_fn)
        system(sge_cmd, check=True)

def _run_script_local(job_data):
        script_fn = job_data["script_fn"]
        job_name = job_data["job_name"]
        log_fn = '{0}.log'.format(script_fn)
        cmd = "bash {0} 1> {1} 2>&1".format(script_fn, log_fn)
        try:
            system(cmd, check=True)
        except Exception:
            out = open(log_fn).read()
            fc_run_logger.exception('Contents of %r:\n%s' %(log_fn, out))
            raise

_run_scripts = {
        'SGE': _run_script_sge,
        'TORQUE': _run_script_torque,
        'SLURM': _run_script_slurm,
        'LOCAL': _run_script_local,
}

def run_script(job_data, job_type = "SGE" ):
    """For now, we actually modify the script before running it.
    This assume a simple bash script.
    We will have a better solution eventually.
    """
    try:
        _run_script = _run_scripts[job_type.upper()]
    except LookupError as e:
        msg = 'Unknown job_type=%s' %repr(job_type)
        fc_run_logger.exception(msg)
        raise
    job_name = job_data["job_name"]
    script_fn = job_data["script_fn"]
    support.update_env_in_script(script_fn,
        ['PATH', 'PYTHONPATH', 'LD_LIBRARY_PATH'])
    fc_run_logger.info('(%s) %r' %(job_type, script_fn))
    fc_run_logger.debug('%s (job %r)' %(_run_script.__name__, job_name))
    rc = _run_script(job_data)
    # Someday, we might trap exceptions here, as a failure would be caught later anyway.

def wait_for_file(filename, task, job_name = ""):
    """We could be in the thread or sub-process which spawned a qsub job,
    so we must check for the shutdown_event.
    """
    while 1:
        time.sleep(wait_time)
        # We prefer all jobs to rely on `*done.exit`, but not all do yet. So we check that 1st.
        exit_fn = filename + '.exit'
        if os.path.exists(exit_fn):
            fc_run_logger.info( "%r found." % (exit_fn) )
            fc_run_logger.debug( " job: %r exited." % (job_name) )
            os.unlink(exit_fn) # to allow a restart later, if not done
            if not os.path.exists(filename):
                fc_run_logger.warning( "%r is missing. job: %r failed!" % (filename, job_name) )
            break
        if os.path.exists(filename) and not os.path.exists(exit_fn):
            # (rechecked exit_fn to avoid race condition)
            fc_run_logger.info( "%r not found, but job is done." % (exit_fn) )
            fc_run_logger.debug( " job: %r exited." % (job_name) )
            break
        if task.shutdown_event is not None and task.shutdown_event.is_set():
            fc_run_logger.warning( "shutdown_event received (Keyboard Interrupt maybe?), %r not finished."
                % (job_name) )
            if support.job_type == "SGE":
                fc_run_logger.info( "deleting the job by `qdel` now..." )
                system("qdel %s" % job_name) # Failure is ok.
            if support.job_type == "SLURM":
                fc_run_logger.info( "Deleting the job by 'scancel' now...")
                system("scancel -n %s" % job_name)
            break

def task_make_fofn_abs_raw(self):
    return support.make_fofn_abs(self.i_fofn.path, self.o_fofn.path)

def task_make_fofn_abs_preads(self):
    return support.make_fofn_abs(self.i_fofn.path, self.o_fofn.path)

def task_build_rdb(self):
    input_fofn_fn = fn(self.input_fofn)
    job_done = fn(self.rdb_build_done)
    work_dir = self.parameters["work_dir"]
    config = self.parameters["config"]
    sge_option_da = config["sge_option_da"]

    script_fn = os.path.join( work_dir, "prepare_rdb.sh" )
    args = {
        'input_fofn_fn': input_fofn_fn,
        'work_dir': work_dir,
        'config': config,
        'job_done': job_done,
        'script_fn': script_fn,
        'run_jobs_fn': fn(self.run_jobs),
    }
    support.build_rdb(**args)

    job_data = support.make_job_data(self.URL, script_fn)
    job_data["sge_option"] = sge_option_da
    run_script(job_data, job_type = config["job_type"])
    wait_for_file(job_done, task=self, job_name=job_data['job_name'])

def task_build_pdb(self):  #essential the same as build_rdb() but the subtle differences are tricky to consolidate to one function
    input_fofn_fn = fn(self.pread_fofn)
    job_done = fn(self.pdb_build_done)
    work_dir = self.parameters["work_dir"]
    config = self.parameters["config"]
    sge_option_pda = config["sge_option_pda"]

    script_fn = os.path.join( work_dir, "prepare_pdb.sh" )
    args = {
        'input_fofn_fn': input_fofn_fn,
        'work_dir': work_dir,
        'config': config,
        'job_done': job_done,
        'script_fn': script_fn,
        'run_jobs_fn': fn(self.run_jobs),
    }
    support.build_pdb(**args)

    job_data = support.make_job_data(self.URL, script_fn)
    job_data["sge_option"] = sge_option_pda
    run_script(job_data, job_type = config["job_type"])
    wait_for_file(job_done, task=self, job_name=job_data['job_name'])

def task_run_falcon_asm(self):
    wd = self.parameters["wd"]
    #p_merge_done = self.p_merge_done
    db_file = fn(self.db_file)
    job_done = fn(self.falcon_asm_done)
    config = self.parameters["config"]
    pread_dir = self.parameters["pread_dir"]
    script_dir = os.path.join( wd )
    script_fn =  os.path.join( script_dir ,"run_falcon_asm.sh" )
    args = {
        'pread_dir': pread_dir,
        'db_file': db_file,
        'config': config,
        'job_done': job_done,
        'script_fn': script_fn,
    }
    support.run_falcon_asm(**args)

    job_data = support.make_job_data(self.URL, script_fn)
    job_data["sge_option"] = config["sge_option_fc"]
    run_script(job_data, job_type = config["job_type"])
    wait_for_file(job_done, task=self, job_name=job_data['job_name'])

def task_run_daligner(self):
    job_done = fn(self.job_done)
    daligner_cmd = self.parameters["daligner_cmd"]
    job_uid = self.parameters["job_uid"]
    cwd = self.parameters["cwd"]
    db_prefix = self.parameters["db_prefix"]
    nblock = self.parameters["nblock"]
    config = self.parameters["config"]
    sge_option_da = config["sge_option_da"]
    script_dir = os.path.join( cwd )
    script_fn =  os.path.join( script_dir , "rj_%s.sh" % (job_uid))
    args = {
        'daligner_cmd': daligner_cmd,
        'db_prefix': db_prefix,
        'nblock': nblock,
        'config': config,
        'job_done': job_done,
        'script_fn': script_fn,
    }
    support.run_daligner(**args)

    job_data = support.make_job_data(self.URL, script_fn)
    job_data["sge_option"] = sge_option_da
    run_script(job_data, job_type = config["job_type"])
    wait_for_file(job_done, task=self, job_name=job_data['job_name'])

def task_run_las_merge(self):
    p_script_fn = self.parameters["merge_script"]
    job_id = self.parameters["job_id"]
    cwd = self.parameters["cwd"]
    job_done = fn(self.job_done)
    config = self.parameters["config"]
    sge_option_la = config["sge_option_la"]

    script_dir = os.path.join( cwd )
    script_fn =  os.path.join( script_dir , "rp_%05d.sh" % (job_id))
    args = {
        'p_script_fn': p_script_fn,
        'config': config,
        'job_done': job_done,
        'script_fn': script_fn,
    }
    support.run_las_merge(**args)

    job_data = support.make_job_data(self.URL, script_fn)
    job_data["sge_option"] = sge_option_la
    run_script(job_data, job_type = config["job_type"])
    wait_for_file(job_done, task=self, job_name=job_data['job_name'])

def task_run_consensus(self):
    out_file_fn = fn(self.out_file)
    job_id = self.parameters["job_id"]
    cwd = self.parameters["cwd"]
    config = self.parameters["config"]
    prefix = self.parameters["prefix"]
    sge_option_cns = config["sge_option_cns"]
    script_dir = os.path.join( cwd )
    job_done = os.path.join( cwd, "c_%05d_done" % job_id )
    script_fn =  os.path.join( script_dir , "c_%05d.sh" % (job_id))
    args = {
        'job_id': job_id,
        'out_file_fn': out_file_fn,
        'prefix': prefix,
        'config': config,
        'job_done': job_done,
        'script_fn': script_fn,
    }
    support.run_consensus(**args)

    job_data = support.make_job_data(self.URL, script_fn)
    job_data["sge_option"] = sge_option_cns
    run_script(job_data, job_type = config["job_type"])
    wait_for_file(job_done, task=self, job_name=job_data['job_name'])

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

    # Symlink all daligner *.las.
    # Could be L1.* or preads.*
    re_las = re.compile(r'\.(\d*)(\.\d*)?\.las$')
    for dal_done in out_dict.values():
        job_rundir = os.path.dirname(fn(dal_done))
        for las_fn in os.listdir(job_rundir):
            mo = re_las.search(las_fn)
            if not mo:
                continue
            block = int(mo.group(1)) # We will merge in the m_* dir of the left block.
            mdir = os.path.join(main_dir, 'm_%05d' %block) # By convention. pbsmrtpipe works differently.
            las_path = os.path.join('..', os.path.basename(job_rundir), las_fn)
            cmd = 'ln -sf {} {}'.format(las_path, mdir)
            system(cmd)
    system("touch %s" %da_done)

def get_nblock(db_file):
    nblock = 1
    new_db = True
    if os.path.exists(db_file):
        with open(db_file) as f:
            for l in f:
                l = l.strip().split()
                if l[0] == "blocks" and l[1] == "=":
                    nblock = int(l[2])
                    new_db = False
                    break
    # Ignore new_db for now.
    return nblock

def create_daligner_tasks(run_jobs_fn, wd, db_prefix, nblock, rdb_build_done, config, pread_aln = False):
    job_id = 0
    tasks = []
    tasks_out = {}

    xform_script = get_script_xformer(pread_aln)

    line_count = 0
    job_descs = get_daligner_job_descriptions(open(run_jobs_fn), db_prefix)
    for desc, bash in job_descs.iteritems():
        job_uid = '%04x' %line_count
        line_count += 1

        support.make_dirs(os.path.join( wd, "./job_%s" % job_uid))
        call = "cd %s/job_%s;ln -sf ../.%s.bps .; ln -sf ../.%s.idx .; ln -sf ../%s.db ." % (wd, job_uid, db_prefix, db_prefix, db_prefix)
        rc = system(call)
        if rc:
            raise Exception("Failure in system call: %r -> %d" %(call, rc))
        job_done = makePypeLocalFile(os.path.abspath( "%s/job_%s/job_%s_done" % (wd, job_uid, job_uid)  ))
        bash = xform_script(bash)
        parameters =  {"daligner_cmd": bash,
                        "cwd": os.path.join(wd, "job_%s" % job_uid),
                        "job_uid": job_uid,
                        "config": config,
                        "nblock": nblock,
                        "db_prefix": db_prefix}
        make_daligner_task = PypeTask( inputs = {"rdb_build_done": rdb_build_done},
                                        outputs = {"job_done": job_done},
                                        parameters = parameters,
                                        TaskType = PypeThreadTaskBase,
                                        URL = "task://localhost/d_%s_%s" % (job_uid, db_prefix) )
        daligner_task = make_daligner_task( task_run_daligner )
        tasks.append( daligner_task )
        tasks_out[ "ajob_%s" % job_uid ] = job_done
        job_id += 1
    return tasks, tasks_out

def create_merge_tasks(run_jobs_fn, wd, db_prefix, input_dep, config):
    merge_tasks = []
    consensus_tasks = []
    merge_out = {}
    consensus_out ={}
    mjob_data = {}

    with open(run_jobs_fn) as f :
        for l in f:
            l = l.strip().split()
            if l[0] not in ( "LAsort", "LAmerge", "mv" ):
                continue
            if l[0] == "LAsort":
                # We now run this part w/ daligner, but we still need
                # a small script for some book-keeping.
                p_id = int( l[2].split(".")[1] )
                mjob_data.setdefault( p_id, [] )
                #mjob_data[p_id].append(  " ".join(l) ) # Already done w/ daligner!
            if l[0] == "LAmerge":
                l2 = l[2].split(".")
                if l2[1][0] == "L":
                    p_id = int(  l[2].split(".")[2] )
                    mjob_data.setdefault( p_id, [] )
                    mjob_data[p_id].append(  " ".join(l) )
                else:
                    p_id = int( l[2].split(".")[1] )
                    mjob_data.setdefault( p_id, [] )
                    mjob_data[p_id].append(  " ".join(l) )
            if l[0] == "mv":
                l2 = l[1].split(".")
                if l2[1][0] == "L":
                    p_id = int(  l[1].split(".")[2] )
                    mjob_data.setdefault( p_id, [] )
                    mjob_data[p_id].append(  " ".join(l) )
                else:
                    p_id = int( l[1].split(".")[1] )
                    mjob_data.setdefault( p_id, [] )
                    mjob_data[p_id].append(  " ".join(l) )

    for p_id in mjob_data:
        s_data = mjob_data[p_id]

        support.make_dirs("%s/m_%05d" % (wd, p_id))
        support.make_dirs("%s/preads" % (wd) )
        support.make_dirs("%s/las_files" % (wd) )

        merge_script_file = os.path.abspath( "%s/m_%05d/m_%05d.sh" % (wd, p_id, p_id) )
        with open(merge_script_file, "w") as merge_script:
            #print >> merge_script, """for f in `find .. -wholename "*job*/%s.%d.%s.*.*.las"`; do ln -sf $f .; done""" % (db_prefix, p_id, db_prefix)
            for l in s_data:
                print >> merge_script, l
            print >> merge_script, "ln -sf ../m_%05d/%s.%d.las ../las_files" % (p_id, db_prefix, p_id) 
            print >> merge_script, "ln -sf ./m_%05d/%s.%d.las .. " % (p_id, db_prefix, p_id) 
            
        job_done = makePypeLocalFile(os.path.abspath( "%s/m_%05d/m_%05d_done" % (wd, p_id, p_id)  ))
        parameters =  {"merge_script": merge_script_file, 
                       "cwd": os.path.join(wd, "m_%05d" % p_id),
                       "job_id": p_id,
                       "config": config}

        make_merge_task = PypeTask( inputs = {"input_dep": input_dep},
                                       outputs = {"job_done": job_done},
                                       parameters = parameters,
                                       TaskType = PypeThreadTaskBase,
                                       URL = "task://localhost/m_%05d_%s" % (p_id, db_prefix) )
        merge_task = make_merge_task ( task_run_las_merge)

        merge_out["mjob_%d" % p_id] = job_done
        merge_tasks.append(merge_task)


        out_file = makePypeLocalFile(os.path.abspath( "%s/preads/out.%05d.fasta" % (wd, p_id)  ))
        out_done = makePypeLocalFile(os.path.abspath( "%s/preads/c_%05d_done" % (wd, p_id)  ))
        parameters =  {"cwd": os.path.join(wd, "preads" ),
                       "job_id": p_id, 
                       "prefix": db_prefix,
                       "config": config}
        make_c_task = PypeTask( inputs = {"job_done": job_done},
                                outputs = {"out_file": out_file, "out_done": out_done },
                                parameters = parameters,
                                TaskType = PypeThreadTaskBase,
                                URL = "task://localhost/ct_%05d" % p_id )
        
        c_task = make_c_task( task_run_consensus)
        consensus_tasks.append(c_task)
        consensus_out["cjob_%d" % p_id] = out_done 

    return merge_tasks, merge_out, consensus_tasks, consensus_out



def main1(prog_name, input_config_fn, logger_config_fn=None):
    global fc_run_logger
    fc_run_logger = support.setup_logger(logger_config_fn)

    fc_run_logger.info( "fc_run started with configuration %s", input_config_fn ) 
    config = support.get_dict_from_old_falcon_cfg(support.parse_config(input_config_fn))
    rawread_dir = os.path.abspath("./0-rawreads")
    pread_dir = os.path.abspath("./1-preads_ovl")
    falcon_asm_dir  = os.path.abspath("./2-asm-falcon")
    script_dir = os.path.abspath("./scripts")
    sge_log_dir = os.path.abspath("./sge_log")

    for d in (rawread_dir, pread_dir, falcon_asm_dir, script_dir, sge_log_dir):
        support.make_dirs(d)

    concurrent_jobs = config["pa_concurrent_jobs"]
    PypeThreadWorkflow.setNumThreadAllowed(concurrent_jobs, concurrent_jobs)
    wf = PypeThreadWorkflow()

    input_fofn_plf = makePypeLocalFile(os.path.basename(config["input_fofn_fn"]))
    rawread_fofn_plf = makePypeLocalFile(os.path.join(rawread_dir, os.path.basename(config["input_fofn_fn"])))
    make_fofn_abs_task = PypeTask(inputs = {"i_fofn": input_fofn_plf},
                                  outputs = {"o_fofn": rawread_fofn_plf},
                                  parameters = {},
                                  TaskType = PypeThreadTaskBase)
    fofn_abs_task = make_fofn_abs_task(task_make_fofn_abs_raw)

    wf.addTasks([fofn_abs_task])
    wf.refreshTargets([fofn_abs_task])

    if config["input_type"] == "raw":
        #### import sequences into daligner DB
        sleep_done = makePypeLocalFile( os.path.join( rawread_dir, "sleep_done") )
        rdb_build_done = makePypeLocalFile( os.path.join( rawread_dir, "rdb_build_done") ) 
        run_jobs = makePypeLocalFile( os.path.join( rawread_dir, "run_jobs.sh") ) 
        parameters = {"work_dir": rawread_dir,
                      "config": config}

        make_build_rdb_task = PypeTask(inputs = {"input_fofn": rawread_fofn_plf},
                                      outputs = {"rdb_build_done": rdb_build_done,
                                                 "run_jobs": run_jobs}, 
                                      parameters = parameters,
                                      TaskType = PypeThreadTaskBase)
        build_rdb_task = make_build_rdb_task(task_build_rdb)

        wf.addTasks([build_rdb_task])
        wf.refreshTargets([rdb_build_done]) 

        db_file = makePypeLocalFile(os.path.join( rawread_dir, "%s.db" % "raw_reads" ))
        raw_reads_nblock = get_nblock(fn(db_file))
        #### run daligner
        daligner_tasks, daligner_out = create_daligner_tasks(fn(run_jobs), rawread_dir, "raw_reads", raw_reads_nblock, rdb_build_done, config) 

        wf.addTasks(daligner_tasks)
        #wf.refreshTargets(updateFreq = 60) # larger number better for more jobs
        r_da_done = makePypeLocalFile( os.path.join( rawread_dir, "da_done") )

        parameters =  {
                "nblock": raw_reads_nblock,
        }
        make_daligner_gather = PypeTask(
                   inputs = daligner_out, 
                   outputs =  {"da_done":r_da_done},
                   parameters = parameters,
                   TaskType = PypeThreadTaskBase,
                   URL = "task://localhost/rda_check" )
        check_r_da_task = make_daligner_gather(task_daligner_gather)
        wf.addTask(check_r_da_task)
        wf.refreshTargets(updateFreq = wait_time) # larger number better for more jobs, need to call to run jobs here or the # of concurrency is changed
        
        concurrent_jobs = config["cns_concurrent_jobs"]
        PypeThreadWorkflow.setNumThreadAllowed(concurrent_jobs, concurrent_jobs)
        merge_tasks, merge_out, consensus_tasks, consensus_out = create_merge_tasks(fn(run_jobs), rawread_dir, "raw_reads", r_da_done, config)
        wf.addTasks( merge_tasks )
        if config["target"] == "overlapping":
            wf.refreshTargets(updateFreq = wait_time) # larger number better for more jobs, need to call to run jobs here or the # of concurrency is changed
            sys.exit(0)
        wf.addTasks( consensus_tasks )

        r_cns_done = makePypeLocalFile( os.path.join( rawread_dir, "cns_done") )
        pread_fofn = makePypeLocalFile( os.path.join( pread_dir,  "input_preads.fofn" ) )

        @PypeTask( inputs = consensus_out, 
                   outputs =  {"cns_done":r_cns_done, "pread_fofn": pread_fofn},
                   TaskType = PypeThreadTaskBase,
                   URL = "task://localhost/cns_check" )
        def check_r_cns_task(self):
            with open(fn(self.pread_fofn),  "w") as f:
                fn_list =  glob.glob("%s/preads/out*.fasta" % rawread_dir)
                fn_list.sort()
                for fa_fn in fn_list:
                    print >>f, fa_fn
            system("touch %s" % fn(self.cns_done))

        wf.addTask(check_r_cns_task)
        wf.refreshTargets(updateFreq = wait_time) # larger number better for more jobs

    if config["target"] == "pre-assembly":
        sys.exit(0)

    # build pread database
    if config["input_type"] == "preads":
        pread_fofn = makePypeLocalFile(os.path.join(pread_dir, os.path.basename(config["input_fofn_fn"])))
        make_fofn_abs_task = PypeTask(inputs = {"i_fofn": rawread_fofn_plf},
                                     outputs = {"o_fofn": pread_fofn},
                                     parameters = {},
                                     TaskType = PypeThreadTaskBase)
        fofn_abs_task = make_fofn_abs_task(task_make_fofn_abs_preads)
        wf.addTasks([fofn_abs_task])
        wf.refreshTargets([fofn_abs_task])

    pdb_build_done = makePypeLocalFile( os.path.join( pread_dir, "pdb_build_done") ) 
    parameters = {"work_dir": pread_dir,
                  "config": config}

    run_jobs = makePypeLocalFile(os.path.join(pread_dir, 'run_jobs.sh'))
    make_build_pdb_task  = PypeTask(inputs = { "pread_fofn": pread_fofn },
                                    outputs = { "pdb_build_done": pdb_build_done,
                                                "run_jobs": run_jobs},
                                    parameters = parameters,
                                    TaskType = PypeThreadTaskBase,
                                    URL = "task://localhost/build_pdb")
    build_pdb_task = make_build_pdb_task(task_build_pdb)

    wf.addTasks([build_pdb_task])
    wf.refreshTargets([pdb_build_done]) 



    db_file = makePypeLocalFile(os.path.join( pread_dir, "%s.db" % "preads" ))
    preads_nblock = get_nblock(fn(db_file))
    #### run daligner
    concurrent_jobs = config["ovlp_concurrent_jobs"]
    PypeThreadWorkflow.setNumThreadAllowed(concurrent_jobs, concurrent_jobs)
    config["sge_option_da"] = config["sge_option_pda"]
    config["sge_option_la"] = config["sge_option_pla"]
    daligner_tasks, daligner_out = create_daligner_tasks(fn(run_jobs), pread_dir, "preads", preads_nblock, pdb_build_done, config, pread_aln= True) 
    wf.addTasks(daligner_tasks)
    #wf.refreshTargets(updateFreq = 30) # larger number better for more jobs

    p_da_done = makePypeLocalFile( os.path.join( pread_dir, "da_done") )
    parameters =  {
            "nblock": preads_nblock,
    }
    make_daligner_gather = PypeTask(
                inputs = daligner_out, 
                outputs =  {"da_done":p_da_done},
                parameters = parameters,
                TaskType = PypeThreadTaskBase,
                URL = "task://localhost/pda_check" )
    check_p_da_task = make_daligner_gather(task_daligner_gather)
    wf.addTask(check_p_da_task)

    merge_tasks, merge_out, consensus_tasks, consensus_out = create_merge_tasks(fn(run_jobs), pread_dir, "preads", p_da_done, config)
    wf.addTasks( merge_tasks )
    #wf.refreshTargets(updateFreq = 30) #all

    p_merge_done = makePypeLocalFile( os.path.join( pread_dir, "p_merge_done") )

    @PypeTask( inputs = merge_out, 
               outputs =  {"p_merge_done":p_merge_done},
               TaskType = PypeThreadTaskBase,
               URL = "task://localhost/pmerge_check" )
    def check_p_merge_check_task(self):
        system("touch %s" % fn(self.p_merge_done))
    
    wf.addTask(check_p_merge_check_task)
    wf.refreshTargets(updateFreq = wait_time) #all

    
    falcon_asm_done = makePypeLocalFile( os.path.join( falcon_asm_dir, "falcon_asm_done") )
    make_run_falcon_asm = PypeTask(
               inputs = {"p_merge_done": p_merge_done, "db_file":db_file},
               outputs =  {"falcon_asm_done":falcon_asm_done},
               parameters = {"wd": falcon_asm_dir,
                             "config": config,
                             "pread_dir": pread_dir},
               TaskType = PypeThreadTaskBase,
               URL = "task://localhost/falcon" )
    wf.addTask(make_run_falcon_asm(task_run_falcon_asm))
    wf.refreshTargets(updateFreq = wait_time) #all


def main(argv=sys.argv):
    print(argv)
    if len(argv) < 2:
        sys.stderr.write( """
you need to specify a configuration file"
usage: fc_run fc_run.cfg [logging.cfg]
""")
        sys.exit(2)
    main1(*argv)
