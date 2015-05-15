#!/usr/bin/env python

#################################################################################$$
# Copyright (c) 2011-2014, Pacific Biosciences of California, Inc.
#
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted (subject to the limitations in the
# disclaimer below) provided that the following conditions are met:
#
#  * Redistributions of source code must retain the above copyright
#  notice, this list of conditions and the following disclaimer.
#
#  * Redistributions in binary form must reproduce the above
#  copyright notice, this list of conditions and the following
#  disclaimer in the documentation and/or other materials provided
#  with the distribution.
#
#  * Neither the name of Pacific Biosciences nor the names of its
#  contributors may be used to endorse or promote products derived
#  from this software without specific prior written permission.
#
# NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
# GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
# BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
# WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
# OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
# USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
# OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
# SUCH DAMAGE.
#################################################################################$$


from pypeflow.common import * 
from pypeflow.data import PypeLocalFile, makePypeLocalFile, fn
from pypeflow.task import PypeTask, PypeThreadTaskBase, PypeTaskBase
from pypeflow.controller import PypeWorkflow, PypeThreadWorkflow
from falcon_kit.FastaReader import FastaReader
import glob
import sys
import os
import re
import time
import logging
import uuid


wait_time = 5
pypeflow_log = 0
if pypeflow_log:
    pypeflow_logger = logging.getLogger("pypeflow.common")
    #pypeflow_logger.setLevel(logging.INFO)
    pypeflow_logger.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh = logging.FileHandler('pypeflow.log')
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(formatter)
    pypeflow_logger.addHandler(fh)

fc_run_log = 1
if fc_run_log:
    fc_run_logger = logging.getLogger("fc_run")
    #fc_run_logger.setLevel(logging.INFO)
    fc_run_logger.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh = logging.FileHandler('fc_run.log')
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(formatter)
    fc_run_logger.addHandler(fh)

def run_script(job_data, job_type = "SGE" ):
    if job_type == "SGE":
        job_name = job_data["job_name"]
        cwd = job_data["cwd"]
        sge_option = job_data["sge_option"]
        script_fn = job_data["script_fn"]
        sge_cmd="qsub -N {job_name} {sge_option} -o {cwd}/sge_log -j y\
                 -S /bin/bash {script}".format(job_name=job_name,  
                                               cwd=os.getcwd(), 
                                               sge_option=sge_option, 
                                               script=script_fn)

        fc_run_logger.info( "submitting %s for SGE, start job: %s " % (script_fn, job_name) )
        cmd = sge_cmd
        rc = os.system(cmd)
    elif job_type == "local":
        script_fn = job_data["script_fn"]
        job_name = job_data["job_name"]
        fc_run_logger.info( "executing %r locally, start job: %r " % (script_fn, job_name) )
        cmd = "bash %s" % script_fn
        rc = os.system(cmd)
    if rc:
        msg = "Cmd %r (job %r) returned %d." % (cmd, job_name, rc)
        fc_run_logger.info(msg)
        # For non-qsub, this might still help with debugging. But technically
        # we should not raise here, as a failure should be noticed later.
        # When we are confident that script failures are handled well,
        # we can make this optional.
        raise Exception(msg)
    else:
        msg = "Cmd %r (job %r) returned %d" % (cmd, job_name, rc)
        fc_run_logger.debug(msg)

def wait_for_file(filename, task = None, job_name = ""):
    while 1:
        time.sleep(wait_time)
        if os.path.exists(filename):
            fc_run_logger.info( "%s generated. job: %s finished." % (filename, job_name) )
            if os.path.exists(filename + '.incomplete'):
                os.unlink(filename) # Tell __call__() that we actually failed.
                fc_run_logger.info( "%s generated. job: %s is incomplete!" % (filename, job_name) )
            break

        if task != None:
            if task.shutdown_event != None and task.shutdown_event.is_set(): 
                fc_run_logger.info( "Keyborad Interrupt Detect, %s not finished, deleting the job by `qdel` now " % (job_name) )
                os.system("qdel %s" % job_name)
                break



def build_rdb(self):  #essential the same as build_rdb() but the subtle differences are tricky to consolidate to one function

    input_fofn = self.input_fofn
    input_fofn_fn = fn(input_fofn)
    rdb_build_done = self.rdb_build_done
    work_dir = self.parameters["work_dir"]
    config = self.parameters["config"]
    sge_option_da = config["sge_option_da"]
    install_prefix = config["install_prefix"]
    length_cutoff = config["length_cutoff"]
    pa_HPCdaligner_option = config["pa_HPCdaligner_option"]
    pa_DBsplit_option = config["pa_DBsplit_option"]
    openending = config["openending"]



    script_fn = os.path.join( work_dir, "prepare_db.sh" )
    
    last_block = 1
    new_db = True
    if os.path.exists( os.path.join(work_dir, "raw_reads.db") ):
        with open(  os.path.join(work_dir, "raw_reads.db") ) as f:
            for l in f:
                l = l.strip().split()
                if l[0] == "blocks" and l[1] == "=":
                    last_block = int(l[2])
                    new_db = False
                    break


    with open(script_fn,"w") as script_file:
        script_file.write("set -e\n")
        script_file.write("touch {rdb_build_done}.incomplete\n".format(rdb_build_done = fn(rdb_build_done)))
        script_file.write("trap 'touch {rdb_build_done}' EXIT\n".format(rdb_build_done = fn(rdb_build_done)))
        script_file.write("source {install_prefix}/bin/activate\n".format(install_prefix = install_prefix))
        script_file.write("cd {work_dir}\n".format(work_dir = work_dir))
        script_file.write("hostname >> db_build.log\n")
        script_file.write("date >> db_build.log\n")
        #script_file.write("for f in `cat {input_fofn_fn}`; do fasta2DB raw_reads $f; done >> db_build.log \n".format(input_fofn_fn = input_fofn_fn))
        script_file.write("fasta2DB -v raw_reads -f{input_fofn_fn} >> db_build.log \n".format(input_fofn_fn = input_fofn_fn))
        if new_db  == True:
            script_file.write("DBsplit %s raw_reads\n" % pa_DBsplit_option)
        if openending == True:
            script_file.write("""LB=$(cat raw_reads.db | awk '$1 == "blocks" {print $3-1}')\n""")
        else:
            script_file.write("""LB=$(cat raw_reads.db | awk '$1 == "blocks" {print $3}')\n""")
        script_file.write("HPCdaligner %s -H%d raw_reads %d-$LB > run_jobs.sh\n" % (pa_HPCdaligner_option, length_cutoff, last_block))
        
        script_file.write("\\rm -f {rdb_build_done}.incomplete\n".format(rdb_build_done = fn(rdb_build_done)))

    job_name = self.URL.split("/")[-1]
    job_name += "-"+str(uuid.uuid4())[:8]
    job_data = {"job_name": job_name,
                "cwd": os.getcwd(),
                "sge_option": sge_option_da,
                "script_fn": script_fn }
    run_script(job_data, job_type = config["job_type"])
    wait_for_file( fn(rdb_build_done), task=self, job_name=job_name )

def build_pdb(self):

    input_fofn = self.pread_fofn
    input_fofn_fn = fn(input_fofn)
    pdb_build_done = self.pdb_build_done
    work_dir = self.parameters["work_dir"]
    config = self.parameters["config"]
    sge_option_pda = config["sge_option_pda"]
    install_prefix = config["install_prefix"]
    length_cutoff = config["length_cutoff_pr"]
    ovlp_HPCdaligner_option = config["ovlp_HPCdaligner_option"]
    ovlp_DBsplit_option = config["ovlp_DBsplit_option"]


    script_fn = os.path.join( work_dir, "prepare_db.sh" )

    with open(script_fn,"w") as script_file:
        script_file.write("set -e\n")
        script_file.write("source {install_prefix}/bin/activate\n".format(install_prefix = install_prefix))
        script_file.write("cd {work_dir}\n".format(work_dir = work_dir))
        script_file.write("hostname >> db_build.log\n")
        script_file.write("date >> db_build.log\n")
        script_file.write("fasta2DB -v preads -f{input_fofn_fn} >> db_build.log \n".format(input_fofn_fn = input_fofn_fn))
        script_file.write("DBsplit -x%d %s preads\n" % (length_cutoff, ovlp_DBsplit_option))
        script_file.write("HPCdaligner %s -H%d preads > run_jobs.sh\n" % (ovlp_HPCdaligner_option, length_cutoff))
        script_file.write("touch {pdb_build_done}\n".format(pdb_build_done = fn(pdb_build_done)))

    job_name = self.URL.split("/")[-1]
    job_name += "-"+str(uuid.uuid4())[:8]
    job_data = {"job_name": job_name,
                "cwd": os.getcwd(),
                "sge_option": sge_option_pda,
                "script_fn": script_fn }
    run_script(job_data, job_type = config["job_type"])
    wait_for_file( fn(pdb_build_done), task=self, job_name=job_name )
    

def run_daligner(self):
    daligner_cmd = self.parameters["daligner_cmd"]
    job_uid = self.parameters["job_uid"]
    cwd = self.parameters["cwd"]
    config = self.parameters["config"]
    sge_option_da = config["sge_option_da"]
    install_prefix = config["install_prefix"]
    db_prefix = self.parameters["db_prefix"]
    nblock = self.parameters["nblock"]

    script_dir = os.path.join( cwd )
    script_fn =  os.path.join( script_dir , "rj_%s.sh" % (job_uid))
    log_path = os.path.join( script_dir, "rj_%s.log" % (job_uid))

    script = []
    script.append( "set -e" )
    script.append( "source {install_prefix}/bin/activate".format(install_prefix = install_prefix) )
    script.append( "cd %s" % cwd )
    script.append( "hostname >> %s" % log_path )
    script.append( "date >> %s" % log_path )
    script.append( "/usr/bin/time "+ daligner_cmd + ( " >> %s 2>&1 " % log_path ) + ( " && touch %s" % fn( self.job_done ) ) )

    for p_id in xrange( 1, nblock+1 ):
        script.append( """ for f in `find $PWD -wholename "*%s.%d.%s.*.*.las"`; do ln -sf $f ../m_%05d; done """  % (db_prefix, p_id, db_prefix, p_id) )

    with open(script_fn,"w") as script_file:
        script_file.write("\n".join(script))

    job_name = self.URL.split("/")[-1]
    job_name += "-"+str(uuid.uuid4())[:8]
    job_data = {"job_name": job_name,
                "cwd": cwd,
                "sge_option": sge_option_da,
                "script_fn": script_fn }
    run_script(job_data, job_type = config["job_type"])
    wait_for_file( fn( self.job_done ), task=self, job_name=job_name )

def run_merge_task(self):
    p_script_fn = self.parameters["merge_script"]
    job_id = self.parameters["job_id"]
    cwd = self.parameters["cwd"]
    config = self.parameters["config"]
    sge_option_la = config["sge_option_la"]
    install_prefix = config["install_prefix"]

    script_dir = os.path.join( cwd )
    script_fn =  os.path.join( script_dir , "rp_%05d.sh" % (job_id))
    log_path = os.path.join( script_dir, "rp_%05d.log" % (job_id))

    script = []
    script.append( "set -e" )
    script.append( "source {install_prefix}/bin/activate".format(install_prefix = install_prefix) )
    script.append( "cd %s" % cwd )
    script.append( "hostname >> %s" % log_path )
    script.append( "date >> %s" % log_path )
    script.append( ("/usr/bin/time bash %s " % p_script_fn)  + ( " >> %s 2>&1" % log_path ) + ( " && touch %s" % fn( self.job_done ) ) )

    with open(script_fn,"w") as script_file:
        script_file.write("\n".join(script))

    job_name = self.URL.split("/")[-1]
    job_name += "-"+str(uuid.uuid4())[:8]
    job_data = {"job_name": job_name,
                "cwd": cwd,
                "sge_option": sge_option_la,
                "script_fn": script_fn }
    run_script(job_data, job_type = config["job_type"])
    wait_for_file( fn( self.job_done ), task=self, job_name=job_name )

def run_consensus_task(self):
    job_id = self.parameters["job_id"]
    cwd = self.parameters["cwd"]
    config = self.parameters["config"]
    sge_option_cns = config["sge_option_cns"]
    install_prefix = config["install_prefix"]
    script_dir = os.path.join( cwd )
    script_fn =  os.path.join( script_dir , "c_%05d.sh" % (job_id))
    log_path = os.path.join( script_dir, "c_%05d.log" % (job_id))
    prefix = self.parameters["prefix"]
    falcon_sense_option = config["falcon_sense_option"]
    length_cutoff = config["length_cutoff"]

    with open( os.path.join(cwd, "cp_%05d.sh" % job_id), "w") as c_script:
        print >> c_script, "source {install_prefix}/bin/activate\n".format(install_prefix = install_prefix)
        print >> c_script, "cd .."
        if config["falcon_sense_skip_contained"] == True:
            print >> c_script, """LA4Falcon -H%d -fso %s las_files/%s.%d.las | """ % (length_cutoff, prefix, prefix, job_id),
        else:
            print >> c_script, """LA4Falcon -H%d -fo %s las_files/%s.%d.las | """ % (length_cutoff, prefix, prefix, job_id),
        print >> c_script, """fc_consensus.py %s > %s""" % (falcon_sense_option, fn(self.out_file))

    script = []
    script.append( "set -e" )
    script.append( "source {install_prefix}/bin/activate".format(install_prefix = install_prefix) )
    script.append( "cd %s" % cwd )
    script.append( "hostname >> %s" % log_path )
    script.append( "date >> %s" % log_path )
    script.append( ("/usr/bin/time bash cp_%05d.sh " % job_id )  + ( " >> %s 2>&1 " % log_path ) + ( " && touch c_%05d_done" % job_id  ) )

    with open(script_fn,"w") as script_file:
        script_file.write("\n".join(script))

    job_name = self.URL.split("/")[-1]
    job_name += "-"+str(uuid.uuid4())[:8]
    job_data = {"job_name": job_name,
                "cwd": cwd,
                "sge_option": sge_option_cns,
                "script_fn": script_fn }
    run_script(job_data, job_type = config["job_type"])
    wait_for_file( os.path.join(cwd,"c_%05d_done" % job_id) , task=self, job_name=job_name )


def create_daligner_tasks(wd, db_prefix, db_file, rdb_build_done, config, pread_aln = False):

    import hashlib
    job_id = 0
    tasks = []
    tasks_out = {}

    nblock = 1
    new_db = True
    if os.path.exists( fn(db_file) ):
        with open( fn(db_file) ) as f:
            for l in f:
                l = l.strip().split()
                if l[0] == "blocks" and l[1] == "=":
                    nblock = int(l[2])
                    new_db = False
                    break

    for pid in xrange(1, nblock + 1):
        try:
            os.makedirs("%s/m_%05d" % (wd, pid))
        except OSError:
            pass


    with open(os.path.join(wd,  "run_jobs.sh")) as f :
        for l in f :
            l = l.strip()
            job_uid = hashlib.md5(l).hexdigest()
            job_uid = job_uid[:8]
            l = l.split()
            if l[0] == "daligner":
                try:
                    os.makedirs(os.path.join( wd, "./job_%s" % job_uid))
                except OSError:
                    pass
                os.system("cd %s/job_%s;ln -sf ../.%s.bps .; ln -sf ../.%s.idx .; ln -sf ../%s.db ." % (wd, job_uid, db_prefix, db_prefix, db_prefix) )
                job_done = makePypeLocalFile(os.path.abspath( "%s/job_%s/job_%s_done" % (wd, job_uid, job_uid)  ))
                if pread_aln == True:
                    l[0] = "daligner_p"
                parameters =  {"daligner_cmd": " ".join(l),
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
                daligner_task = make_daligner_task ( run_daligner )
                tasks.append( daligner_task )
                tasks_out[ "ajob_%s" % job_uid ] = job_done
                job_id += 1
    return tasks, tasks_out

def create_merge_tasks(wd, db_prefix, input_dep, config):
    merge_tasks = []
    consensus_tasks = []
    merge_out = {}
    consensus_out ={}
    mjob_data = {}

    with open(os.path.join(wd,  "run_jobs.sh")) as f :
        for l in f:
            l = l.strip().split()
            if l[0] not in ( "LAsort", "LAmerge", "mv" ):
                continue
            if l[0] == "LAsort":
                p_id = int( l[2].split(".")[1] )
                mjob_data.setdefault( p_id, [] )
                mjob_data[p_id].append(  " ".join(l) )
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

        try:
            os.makedirs("%s/m_%05d" % (wd, p_id))
        except OSError:
            pass
        try:
            os.makedirs("%s/preads" % (wd) )
        except OSError:
            pass
        try:
            os.makedirs("%s/las_files" % (wd) )
        except OSError:
            pass

        with open("%s/m_%05d/m_%05d.sh" % (wd, p_id, p_id), "w") as merge_script:
            #print >> merge_script, """for f in `find .. -wholename "*job*/%s.%d.%s.*.*.las"`; do ln -sf $f .; done""" % (db_prefix, p_id, db_prefix)
            for l in s_data:
                print >> merge_script, l
            print >> merge_script, "ln -sf ../m_%05d/%s.%d.las ../las_files" % (p_id, db_prefix, p_id) 
            print >> merge_script, "ln -sf ./m_%05d/%s.%d.las .. " % (p_id, db_prefix, p_id) 
            
        merge_script_file = os.path.abspath( "%s/m_%05d/m_%05d.sh" % (wd, p_id, p_id) )
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
        merge_task = make_merge_task ( run_merge_task )

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
        
        c_task = make_c_task( run_consensus_task )
        consensus_tasks.append(c_task)
        consensus_out["cjob_%d" % p_id] = out_done 

    return merge_tasks, merge_out, consensus_tasks, consensus_out



def get_config(config_fn):

    import ConfigParser

    config = ConfigParser.ConfigParser()

    config.read(config_fn)
    
    job_type = "SGE"
    if config.has_option('General', 'job_type'):
        job_type = config.get('General', 'job_type')
    
    pa_concurrent_jobs = 8
    if config.has_option('General', 'pa_concurrent_jobs'):
        pa_concurrent_jobs = config.getint('General', 'pa_concurrent_jobs')
    
    cns_concurrent_jobs = 8
    if config.has_option('General', 'cns_concurrent_jobs'):
        cns_concurrent_jobs = config.getint('General', 'cns_concurrent_jobs')

    ovlp_concurrent_jobs = 8
    if config.has_option('General', 'ovlp_concurrent_jobs'):
        ovlp_concurrent_jobs = config.getint('General', 'ovlp_concurrent_jobs')

    #appending = False
    #if config.has_option('General', 'appending'):
    #    appending = config.get('General', 'appending')
    #    if appending == "True":
    #        appending = True

    openending = False
    if config.has_option('General', 'openending'):
        openending = config.get('General', 'openending')
        if openending == "True":
            openending = True

    input_type = "raw"
    if config.has_option('General', 'input_type'):
        input_type = config.get('General', 'input_type')

    overlap_filtering_setting =  """--max_diff 1000 --max_cov 1000 --min_cov 2"""
    if config.has_option('General', 'overlap_filtering_setting'):
        overlap_filtering_setting = config.get('General', 'overlap_filtering_setting')

    pa_HPCdaligner_option = """-v -dal4 -t16 -e.70 -l1000 -s100"""
    if config.has_option('General', 'pa_HPCdaligner_option'):
        pa_HPCdaligner_option = config.get('General', 'pa_HPCdaligner_option')

    ovlp_HPCdaligner_option = """ -v -dal24 -t32 -h60 -e.96 -l500 -s1000"""
    if config.has_option('General', 'ovlp_HPCdaligner_option'):
        ovlp_HPCdaligner_option = config.get('General', 'ovlp_HPCdaligner_option')

    pa_DBsplit_option = """ -x500 -s200"""
    if config.has_option('General', 'pa_DBsplit_option'):
        pa_DBsplit_option = config.get('General', 'pa_DBsplit_option')

    ovlp_DBsplit_option = """ -x500 -s200"""
    if config.has_option('General', 'ovlp_DBsplit_option'):
        ovlp_DBsplit_option = config.get('General', 'ovlp_DBsplit_option')

    falcon_sense_option = """ --output_multi --min_idt 0.70 --min_cov 2 --local_match_count_threshold 0 --max_n_read 1800 --n_core 6"""
    if config.has_option('General', 'falcon_sense_option'):
        falcon_sense_option = config.get('General', 'falcon_sense_option')

    falcon_sense_skip_contained = "False"
    if config.has_option('General', 'falcon_sense_skip_contained'):
        falcon_sense_skip_contained = config.get('General', 'falcon_sense_skip_contained')
        if falcon_sense_skip_contained in ["True", "true", "1"]:
            falcon_sense_skip_contained = True
        else:
            falcon_sense_skip_contained = False

    length_cutoff = config.getint('General', 'length_cutoff')
    input_fofn_fn = config.get('General', 'input_fofn')
    
    length_cutoff_pr = config.getint('General', 'length_cutoff_pr')
    

    bestn = 12
    if config.has_option('General', 'bestn'):
        bestn = config.getint('General', 'bestn')


    if config.has_option('General', 'target'):
        target = config.get('General', 'target')
        if target not in ["overlapping", "pre-assembly", "assembly"]:
            print """ Target has to be "overlapping", "pre-assembly" or "assembly" in this verison. You have an unknown target %s in the configuration file.  """ % target
            sys.exit(1)
    else:
        print """ No target specified, assuming "assembly" as target """
        target = "assembly"


    hgap_config = {"input_fofn_fn" : input_fofn_fn,
                   "target" : target,
                   "job_type" : job_type,
                   "input_type": input_type,
                   "openending": openending,
                   "pa_concurrent_jobs" : pa_concurrent_jobs,
                   "ovlp_concurrent_jobs" : ovlp_concurrent_jobs,
                   "cns_concurrent_jobs" : cns_concurrent_jobs,
                   "overlap_filtering_setting": overlap_filtering_setting,
                   "length_cutoff" : length_cutoff,
                   "length_cutoff_pr" : length_cutoff_pr,
                   "sge_option_da": config.get('General', 'sge_option_da'),
                   "sge_option_la": config.get('General', 'sge_option_la'),
                   "sge_option_pda": config.get('General', 'sge_option_pda'),
                   "sge_option_pla": config.get('General', 'sge_option_pla'),
                   "sge_option_fc": config.get('General', 'sge_option_fc'),
                   "sge_option_cns": config.get('General', 'sge_option_cns'),
                   "pa_HPCdaligner_option": pa_HPCdaligner_option,
                   "ovlp_HPCdaligner_option": ovlp_HPCdaligner_option,
                   "pa_DBsplit_option": pa_DBsplit_option,
                   "ovlp_DBsplit_option": ovlp_DBsplit_option,
                   "falcon_sense_option": falcon_sense_option,
                   "falcon_sense_skip_contained": falcon_sense_skip_contained
                   }

    hgap_config["install_prefix"] = sys.prefix
    
    return hgap_config


if __name__ == '__main__':

    if len(sys.argv) < 2:
        print "you need to specify a configuration file"
        print "usage: fc_run.py fc_run.cfg"
        sys.exit(1)
    
    fc_run_logger.info( "fc_run started with configuration %s", sys.argv[1] ) 
    rawread_dir = os.path.abspath("./0-rawreads")
    pread_dir = os.path.abspath("./1-preads_ovl")
    falcon_asm_dir  = os.path.abspath("./2-asm-falcon")
    script_dir = os.path.abspath("./scripts")
    sge_log_dir = os.path.abspath("./sge_log")

    for d in (rawread_dir, pread_dir, falcon_asm_dir, script_dir, sge_log_dir):
        try:
            os.makedirs(d)
        except:
            pass

    config = get_config(sys.argv[1])
    concurrent_jobs = config["pa_concurrent_jobs"]
    PypeThreadWorkflow.setNumThreadAllowed(concurrent_jobs, concurrent_jobs)
    wf = PypeThreadWorkflow()

    if config["input_type"] == "raw":
        #### import sequences into daligner DB
        input_h5_fofn = makePypeLocalFile( os.path.abspath( config["input_fofn_fn"] ) )
        rdb_build_done = makePypeLocalFile( os.path.join( rawread_dir, "rdb_build_done") ) 
        parameters = {"work_dir": rawread_dir,
                      "config": config}

        make_build_rdb_task = PypeTask(inputs = {"input_fofn": input_h5_fofn},
                                      outputs = {"rdb_build_done": rdb_build_done}, 
                                      parameters = parameters,
                                      TaskType = PypeThreadTaskBase)

        build_rdb_task = make_build_rdb_task(build_rdb)

        wf.addTasks([build_rdb_task])
        wf.refreshTargets([rdb_build_done]) 
        

        db_file = makePypeLocalFile(os.path.join( rawread_dir, "%s.db" % "raw_reads" ))
        #### run daligner
        daligner_tasks, daligner_out = create_daligner_tasks( rawread_dir, "raw_reads", db_file, rdb_build_done, config) 

        wf.addTasks(daligner_tasks)
        #wf.refreshTargets(updateFreq = 60) # larger number better for more jobs

        r_da_done = makePypeLocalFile( os.path.join( rawread_dir, "da_done") )

        @PypeTask( inputs = daligner_out, 
                   outputs =  {"da_done":r_da_done},
                   TaskType = PypeThreadTaskBase,
                   URL = "task://localhost/rda_check" )
        def check_r_da_task(self):
            os.system("touch %s" % fn(self.da_done))
        
        wf.addTask(check_r_da_task)
        wf.refreshTargets(updateFreq = wait_time) # larger number better for more jobs, need to call to run jobs here or the # of concurrency is changed
        
        concurrent_jobs = config["cns_concurrent_jobs"]
        PypeThreadWorkflow.setNumThreadAllowed(concurrent_jobs, concurrent_jobs)
        merge_tasks, merge_out, consensus_tasks, consensus_out = create_merge_tasks( rawread_dir, "raw_reads", r_da_done, config )
        wf.addTasks( merge_tasks )
        if config["target"] == "overlapping":
            wf.refreshTargets(updateFreq = wait_time) # larger number better for more jobs, need to call to run jobs here or the # of concurrency is changed
            exit(0)
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
            os.system("touch %s" % fn(self.cns_done))

        wf.addTask(check_r_cns_task)
        wf.refreshTargets(updateFreq = wait_time) # larger number better for more jobs

    if config["target"] == "pre-assembly":
        exit(0)

    # build pread database
    if config["input_type"] == "preads":
        if not os.path.exists( "%s/input_preads.fofn" % pread_dir):
            os.system( "cp %s %s/input_preads.fofn" % (os.path.abspath( config["input_fofn_fn"] ), pread_dir) )
        pread_fofn = makePypeLocalFile( os.path.join( pread_dir,  "input_preads.fofn" ) )

    pdb_build_done = makePypeLocalFile( os.path.join( pread_dir, "pdb_build_done") ) 
    parameters = {"work_dir": pread_dir,
                  "config": config}

    make_build_pdb_task  = PypeTask( inputs = { "pread_fofn": pread_fofn },
                                    outputs = { "pdb_build_done": pdb_build_done },
                                    parameters = parameters,
                                    TaskType = PypeThreadTaskBase,
                                    URL = "task://localhost/build_pdb")
    build_pdb_task = make_build_pdb_task(build_pdb)

    wf.addTasks([build_pdb_task])
    wf.refreshTargets([pdb_build_done]) 



    db_file = makePypeLocalFile(os.path.join( pread_dir, "%s.db" % "preads" ))
    #### run daligner
    concurrent_jobs = config["ovlp_concurrent_jobs"]
    PypeThreadWorkflow.setNumThreadAllowed(concurrent_jobs, concurrent_jobs)
    config["sge_option_da"] = config["sge_option_pda"]
    config["sge_option_la"] = config["sge_option_pla"]
    daligner_tasks, daligner_out = create_daligner_tasks( pread_dir, "preads", db_file, pdb_build_done, config, pread_aln= True) 
    wf.addTasks(daligner_tasks)
    #wf.refreshTargets(updateFreq = 30) # larger number better for more jobs

    p_da_done = makePypeLocalFile( os.path.join( pread_dir, "da_done") )

    @PypeTask( inputs = daligner_out, 
               outputs =  {"da_done":p_da_done},
               TaskType = PypeThreadTaskBase,
               URL = "task://localhost/pda_check" )
    def check_p_da_task(self):
        os.system("touch %s" % fn(self.da_done))
    
    wf.addTask(check_p_da_task)

    merge_tasks, merge_out, consensus_tasks, consensus_out = create_merge_tasks( pread_dir, "preads", p_da_done, config )
    wf.addTasks( merge_tasks )
    #wf.refreshTargets(updateFreq = 30) #all            

    p_merge_done = makePypeLocalFile( os.path.join( pread_dir, "p_merge_done") )

    @PypeTask( inputs = merge_out, 
               outputs =  {"p_merge_done":p_merge_done},
               TaskType = PypeThreadTaskBase,
               URL = "task://localhost/pmerge_check" )
    def check_p_merge_check_task(self):
        os.system("touch %s" % fn(self.p_merge_done))
    
    wf.addTask(check_p_merge_check_task)
    wf.refreshTargets(updateFreq = wait_time) #all            

    
    falcon_asm_done = makePypeLocalFile( os.path.join( falcon_asm_dir, "falcon_asm_done") )
    @PypeTask( inputs = {"p_merge_done": p_merge_done, "db_file":db_file}, 
               outputs =  {"falcon_asm_done":falcon_asm_done},
               parameters = {"wd": falcon_asm_dir,
                             "config": config,
                             "pread_dir": pread_dir},
               TaskType = PypeThreadTaskBase,
               URL = "task://localhost/falcon" )

    def run_falcon_asm_task(self):
        wd = self.parameters["wd"]
        config = self.parameters["config"]
        install_prefix = config["install_prefix"]
        pread_dir = self.parameters["pread_dir"]
        script_dir = os.path.join( wd )
        script_fn =  os.path.join( script_dir ,"run_falcon_asm.sh" )
        
        script = []
        script.append( "set -e" )
        script.append( "source {install_prefix}/bin/activate".format(install_prefix = install_prefix) )
        script.append( "cd %s" % pread_dir )
        script.append( "DB2Falcon preads")
        script.append( "cd %s" % wd )
        script.append( """find %s/las_files -name "*.las" > las.fofn """ % pread_dir )
        overlap_filtering_setting = config["overlap_filtering_setting"]
        length_cutoff_pr = config["length_cutoff_pr"]
        script.append( """fc_ovlp_filter.py --db %s --fofn las.fofn %s --min_len %d > preads.ovl""" %\
                (fn(db_file), overlap_filtering_setting, length_cutoff_pr) )

        script.append( "ln -sf %s/preads4falcon.fasta ." % pread_dir)
        script.append( """fc_ovlp_to_graph.py preads.ovl --min_len %d > fc_ovlp_to_graph.log""" % length_cutoff_pr)
        script.append( """fc_graph_to_contig.py""" )
        script.append( """touch %s\n""" % fn(self.falcon_asm_done))

        with open(script_fn, "w") as script_file:
            script_file.write("\n".join(script))

        job_name = self.URL.split("/")[-1]
        job_name += "-"+str(uuid.uuid4())[:8]
        job_data = {"job_name": job_name,
                    "cwd": wd,
                    "sge_option": config["sge_option_fc"],
                    "script_fn": script_fn }
        run_script(job_data, job_type = config["job_type"])
        wait_for_file( fn(self.falcon_asm_done), task=self, job_name=job_name )
    
    wf.addTask( run_falcon_asm_task )
    wf.refreshTargets(updateFreq = wait_time) #all            
