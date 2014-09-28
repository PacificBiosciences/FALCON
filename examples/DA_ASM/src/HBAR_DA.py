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
from pbcore.io import FastaReader
import glob
import sys
import os
import re
import time
import logging
import uuid

log = 0
if log:
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    logger.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    ch.setFormatter(formatter)
    logger.addHandler(ch)

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

        os.system( sge_cmd )
    elif job_type == "local":
        os.system( "bash %s" % job_data["script_fn"] )

def wait_for_file(filename, task = None, job_name = ""):
    while 1:
        time.sleep(30)
        if os.path.exists(filename):
            break

        if task != None:
            if task.shutdown_event != None and task.shutdown_event.is_set(): 
                os.system("qdel %s" % job_name)
                break


def build_rdb(self):

    input_fofn = self.input_fofn
    input_fofn_fn = fn(input_fofn)
    rdb_build_done = self.rdb_build_done
    work_dir = self.parameters["work_dir"]
    config = self.parameters["config"]
    sge_option_da = config["sge_option_da"]
    install_prefix = config["install_prefix"]


    script_fn = os.path.join( work_dir, "prepare_db.sh" )
    with open(script_fn,"w") as script_file:
        script_file.write("source {install_prefix}/bin/activate\n".format(install_prefix = install_prefix))
        script_file.write("cd {work_dir}\n".format(work_dir = work_dir))
        script_file.write("for f in `cat {input_fofn_fn}`; do fasta2DB raw_reads $f; done\n".format(input_fofn_fn = input_fofn_fn))
        script_file.write("DBsplit -x500 -s400 raw_reads\n")
        script_file.write("HPCdaligner -v -dal4 -t16 -e.70 -H1000 -l1000 -s1000 raw_reads > run_jobs.sh\n")
        script_file.write("touch {rdb_build_done}\n".format(rdb_build_done = fn(rdb_build_done)))

    job_name = self.URL.split("/")[-1]
    job_name += "-"+str(uuid.uuid1())[:8]
    job_data = {"job_name": job_name,
                "cwd": os.getcwd(),
                "sge_option": sge_option_da,
                "script_fn": script_fn }
    run_script(job_data, job_type = config["job_type"])
    wait_for_file( fn(rdb_build_done), task=self, job_name=job_name )

def run_daligner(self):
    daligner_cmd = self.parameters["daligner_cmd"]
    job_id = self.parameters["job_id"]
    cwd = self.parameters["cwd"]
    config = self.parameters["config"]
    sge_option_da = config["sge_option_da"]
    install_prefix = config["install_prefix"]

    script_dir = os.path.join( cwd )
    script_fn =  os.path.join( script_dir , "rj_%05d.sh" % (job_id))
    log_path = os.path.join( script_dir, "rj_%05d.log" % (job_id))

    script = []
    script.append( "source {install_prefix}/bin/activate\n".format(install_prefix = install_prefix) )
    script.append( "cd %s" % cwd )
    script.append( "/usr/bin/time "+ daligner_cmd + ( " >& %s " % log_path ) + ( " && touch %s" % fn( self.job_done ) ) )

    with open(script_fn,"w") as script_file:
        script_file.write("\n".join(script))

    job_name = self.URL.split("/")[-1]
    job_name += "-"+str(uuid.uuid1())[:8]
    job_data = {"job_name": job_name,
                "cwd": cwd,
                "sge_option": sge_option_da,
                "script_fn": script_fn }
    run_script(job_data, job_type = "SGE")
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
    script.append( "source {install_prefix}/bin/activate\n".format(install_prefix = install_prefix) )
    script.append( "cd %s" % cwd )
    script.append( ("/usr/bin/time bash %s " % p_script_fn)  + ( " >& %s " % log_path ) + ( " && touch %s" % fn( self.job_done ) ) )

    with open(script_fn,"w") as script_file:
        script_file.write("\n".join(script))



    job_name = self.URL.split("/")[-1]
    job_name += "-"+str(uuid.uuid1())[:8]
    job_data = {"job_name": job_name,
                "cwd": cwd,
                "sge_option": sge_option_la,
                "script_fn": script_fn }
    run_script(job_data, job_type = "SGE")
    wait_for_file( fn( self.job_done ), task=self, job_name=job_name )

def run_consensus_task(self):
    job_id = self.parameters["job_id"]
    cwd = self.parameters["cwd"]
    config = self.parameters["config"]
    sge_option_la = config["sge_option_da"]
    install_prefix = config["install_prefix"]
    script_dir = os.path.join( cwd )
    script_fn =  os.path.join( script_dir , "c_%05d.sh" % (job_id))
    log_path = os.path.join( script_dir, "c_%05d.log" % (job_id))
    prefix = self.parameters["prefix"]

    with open( os.path.join(cwd, "cp_%05d.sh" % job_id), "w") as c_script:
        print >> c_script, "source {install_prefix}/bin/activate\n".format(install_prefix = install_prefix)
        print >> c_script, "cd .."
        print >> c_script, """LA4Falcon -o -H1000 -f:%s las_files/%s.%d.las | """ % (prefix, prefix, job_id),
        print >> c_script, """ falcon_sense.py --output_multi  --min_idt 0.70 --min_cov 2 --local_match_count_threshold 0 --max_n_read 1800 --n_core 6 > %s""" % fn(self.out_file)

    script = []
    script.append( "source {install_prefix}/bin/activate\n".format(install_prefix = install_prefix) )
    script.append( "cd %s" % cwd )
    script.append( ("/usr/bin/time bash cp_%05d.sh " % job_id )  + ( " >& %s " % log_path ) + ( " && touch c_%05d_done" % job_id  ) )

    with open(script_fn,"w") as script_file:
        script_file.write("\n".join(script))

    job_name = self.URL.split("/")[-1]
    job_name += "-"+str(uuid.uuid1())[:8]
    job_data = {"job_name": job_name,
                "cwd": cwd,
                "sge_option": sge_option_la,
                "script_fn": script_fn }
    run_script(job_data, job_type = "SGE")
    wait_for_file( os.path.join(cwd,"c_%05d_done" % job_id) , task=self, job_name=job_name )


def create_daligner_tasks(wd, db_prefix, db_file, rdb_build_done, config, pread_aln = False):

    job_id = 0
    tasks = []
    with open(os.path.join(wd,  "run_jobs.sh")) as f :
        for l in f :
            l = l.strip().split()
            if l[0] == "daligner":
                try:
                    os.makedirs(os.path.join( wd, "./job_%05d" % job_id))
                except OSError:
                    pass
                os.system("cd %s/job_%05d;ln -sf ../.%s.bps .; ln -sf ../.%s.idx .; ln -sf ../%s.db ." % (wd, job_id, db_prefix, db_prefix, db_prefix) )
                job_done = makePypeLocalFile(os.path.abspath( "%s/job_%05d/job_%05d_done" % (wd, job_id, job_id)  ))
                if pread_aln == True:
                    l[0] = "daligner_p"
                parameters =  {"daligner_cmd": " ".join(l),
                               "cwd": os.path.join(wd, "job_%05d" % job_id),
                               "job_id": job_id,
                               "config": config}
                make_daligner_task = PypeTask( inputs = {"rdb_build_done": rdb_build_done},
                                               outputs = {"job_done": job_done},
                                               parameters = parameters,
                                               TaskType = PypeThreadTaskBase,
                                               URL = "task://localhost/d_%05d_%s" % (job_id, db_prefix) )
                daligner_task = make_daligner_task ( run_daligner )
                tasks.append( daligner_task )
                job_id += 1
    return tasks

def create_merge_tasks(wd, db_prefix, db_file, config):
    merge_tasks = []
    consensus_tasks = []
    mjob_data = {}

    with open(os.path.join(wd,  "run_jobs.sh")) as f :
        for l in f:
            l = l.strip().split()
            if l[0] not in ( "LAsort", "LAmerge" ):
                continue
            if l[0] == "LAsort":
                p_id = int( l[2].split(".")[1] )
                mjob_data.setdefault( p_id, [] )
                mjob_data[p_id].append(  " ".join(l) )
            if l[0] == "LAmerge":
                l2 = l[2].split(".")
                if l2[1] == "L2":
                    p_id = int(  l[2].split(".")[2] )
                    mjob_data.setdefault( p_id, [] )
                    mjob_data[p_id].append(  " ".join(l) )
                else:
                    p_id = int( l[2].split(".")[1] )
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
            print >> merge_script, """for f in `find .. -wholename "*job*/%s.%d.%s.*.*.las"`; do ln -sf $f .; done""" % (db_prefix, p_id, db_prefix)
            for l in s_data:
                print >> merge_script, l
                print >> merge_script, "mv %s.%d.las ../las_files" % (db_prefix, p_id) 
            
        merge_script_file = os.path.abspath( "%s/m_%05d/m_%05d.sh" % (wd, p_id, p_id) )
        job_done = makePypeLocalFile(os.path.abspath( "%s/m_%05d/m_%05d_done" % (wd, p_id, p_id)  ))
        parameters =  {"merge_script": merge_script_file, 
                       "cwd": os.path.join(wd, "m_%05d" % p_id),
                       "job_id": p_id,
                       "config": config}

        make_merge_task = PypeTask( inputs = {"db_file": db_file},
                                       outputs = {"job_done": job_done},
                                       parameters = parameters,
                                       TaskType = PypeThreadTaskBase,
                                       URL = "task://localhost/m_%05d_%s" % (p_id, db_prefix) )
        merge_task = make_merge_task ( run_merge_task )

        merge_tasks.append(merge_task)


        out_file = makePypeLocalFile(os.path.abspath( "%s/preads/out.%04d.fa" % (wd, p_id)  ))
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

    return merge_tasks, consensus_tasks



def get_config(config_fn):

    import ConfigParser

    config = ConfigParser.RawConfigParser()

    config.read(config_fn)
    
    job_type = "SGE"
    if config.has_option('General', 'job_type'):
        job_type = config.get('General', 'job_type')
    
    concurrent_jobs = 8
    if config.has_option('General', 'concurrent_jobs'):
        concurrent_jobs = config.getint('General', 'concurrent_jobs')

    length_cutoff = config.getint('General', 'length_cutoff')
    input_fofn_fn = config.get('General', 'input_fofn')
    
    length_cutoff_pr = config.getint('General', 'length_cutoff_pr')
    

    bestn = 12
    if config.has_option('General', 'bestn'):
        bestn = config.getint('General', 'bestn')


    if config.has_option('General', 'target'):
        target = config.get('General', 'target')
        if target not in ["mapping", "pre_assembly", "falcon_asm"]:
            print """ Target has to be "mapping", "pre_assembly" or "falcon_asm" in this verison. You have an unknown target %s in the configuration file.  """ % target
            sys.exit(1)
    else:
        print """ No target specified, assuming a "pre-assembly" only target """
        target = "pre_assembly"


    hgap_config = {"input_fofn_fn" : input_fofn_fn,
                   "job_type" : job_type,
                   "concurrent_jobs" : concurrent_jobs,
                   "length_cutoff" : length_cutoff,
                   "length_cutoff_pr" : length_cutoff_pr,
                   "sge_option_da": config.get('General', 'sge_option_da'),
                   "sge_option_la": config.get('General', 'sge_option_la'),
                   "sge_option_pda": config.get('General', 'sge_option_pda'),
                   "sge_option_pla": config.get('General', 'sge_option_pla')
                   }

    hgap_config["install_prefix"] = sys.prefix
    
    return hgap_config


if __name__ == '__main__':

    if len(sys.argv) < 2:
        print "you need to specify a configuration file"
        print "example: HGAP.py HGAP_run.cfg"
        sys.exit(1)
    
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
    concurrent_jobs = config["concurrent_jobs"]
    PypeThreadWorkflow.setNumThreadAllowed(concurrent_jobs, concurrent_jobs)
    wf = PypeThreadWorkflow()


    #### import sequences into daligner DB
    input_h5_fofn = makePypeLocalFile( os.path.abspath( config["input_fofn_fn"] ) )
    rdb_build_done = makePypeLocalFile( os.path.join( rawread_dir, "rdb_build_done") ) 
    parameters = {"work_dir": rawread_dir,
                  "config": config} 

    make_buid_rdb_task = PypeTask(inputs = {"input_fofn": input_h5_fofn},
                                  outputs = {"rdb_build_done": rdb_build_done}, 
                          parameters = parameters,
                          TaskType = PypeThreadTaskBase)
    buid_rdb_task = make_buid_rdb_task(build_rdb)

    wf.addTasks([buid_rdb_task])
    wf.refreshTargets([rdb_build_done]) 

    

    db_file = makePypeLocalFile(os.path.join( rawread_dir, "%s.db" % "raw_reads" ))
    #### run daligner
    daligner_tasks = create_daligner_tasks( rawread_dir, "raw_reads", db_file, rdb_build_done, config) 

    wf.addTasks(daligner_tasks)
    #wf.refreshTargets(updateFreq = 45) # larger number better for more jobs
    wf.refreshTargets(updateFreq = 1) # larger number better for more jobs
    
    merge_tasks, consensus_tasks = create_merge_tasks( rawread_dir, "raw_reads", db_file, config )
    wf.addTasks( merge_tasks )
    wf.addTasks( consensus_tasks )
    wf.refreshTargets(updateFreq = 15) #all            


    if not os.path.exists("%s/rdb_build_done"  % pread_dir):
        with open("%s/preads_norm.fasta" % pread_dir, "w") as p_norm:
            for fa_fn in glob.glob("%s/preads/out*.fa" % rawread_dir):
                f = FastaReader(fa_fn)
                for r in f:
                    name = r.name
                    name = name.replace("_","")
                    print >> p_norm, ">prolog/%s/%d_%d" % ( name, 0, len(r.sequence) )
                    for i in range(0, len(r.sequence)/80):
                        print >> p_norm, r.sequence[ i *80 : (i + 1) * 80]
                    print >> p_norm, r.sequence[(i+1)*80:]

        os.system("cd %s; fasta2DB preads preads_norm.fasta" % pread_dir)
        os.system("cd %s; DBsplit -x500 -s400 preads" % pread_dir)
        os.system("cd %s; HPCdaligner -v -dal24 -t32 -h60 -e.96 -l500 -s1000  preads > run_jobs.sh" % pread_dir)
        os.system("cd %s; touch rdb_build_done" % pread_dir)

    rdb_build_done = makePypeLocalFile( os.path.join( pread_dir, "rdb_build_done") ) 

    db_file = makePypeLocalFile(os.path.join( pread_dir, "%s.db" % "preads" ))
    #### run daligner
    daligner_tasks = create_daligner_tasks( pread_dir, "preads", db_file, rdb_build_done, config, pread_aln= True) 
    wf.addTasks(daligner_tasks)
    #wf.refreshTargets(updateFreq = 45) # larger number better for more jobs
    wf.refreshTargets(updateFreq = 1) # larger number better for more jobs
    merge_tasks, consensus_tasks = create_merge_tasks( pread_dir, "preads", db_file, config )
    wf.addTasks( merge_tasks )
    wf.refreshTargets(updateFreq = 15) #all            



    p_las_merge_done = makePypeLocalFile( os.path.join( pread_dir, "p_las_merge_done") )
    ovlp_filter_done = makePypeLocalFile( os.path.join( pread_dir, "las_files", "ovlp_filter_done") )
    if not os.path.exists(fn(p_las_merge_done)):
        os.system("touch %s" % fn(p_las_merge_done))


    @PypeTask( inputs = {"p_las_merge_done":p_las_merge_done}, 
               outputs =  {"ovlp_filter_done":ovlp_filter_done},
               parameters = {"wd": os.path.join( pread_dir, "las_files" ),
                             "config": config},
               TaskType = PypeThreadTaskBase,
               URL = "task://localhost/ovlf" )
    def run_ovlp_filter_task(self):

        wd = self.parameters["wd"]
        config = self.parameters["config"]
        install_prefix = config["install_prefix"]
        script_dir = os.path.join( wd )
        script_fn =  os.path.join( script_dir , "ovlp_filter.sh" )
        
        script = []
        script.append( "source {install_prefix}/bin/activate".format(install_prefix = install_prefix) )
        script.append( "cd %s" % pread_dir )
        script.append( "DB2Falcon preads")
        script.append( "cd %s/las_files" % pread_dir )
        script.append( """parallel -j 24 "LA4Falcon -mo -H4000 {}  | overlap_filter_step1.py > {}.rc" ::: *.las""" )
        script.append( """cat *.rc > rc_out_all""" )
        script.append( """rm *.rc""" )
        script.append( """parallel -j 24 "LA4Falcon -mo -H4000 {}  | overlap_filter_step2.py > {}.ovl" ::: *.las""")
        script.append( """touch %s\n""" % fn(self.ovlp_filter_done))

        with open(script_fn, "w") as script_file:
            script_file.write("\n".join(script))

        job_name = self.URL.split("/")[-1]
        job_name += "-"+str(uuid.uuid1())[:8]
        job_data = {"job_name": job_name,
                    "cwd": wd,
                    "sge_option": " -pe smp 24 -q huasm",
                    "script_fn": script_fn }
        run_script(job_data, job_type = "SGE")
        wait_for_file( fn(self.ovlp_filter_done), task=self, job_name=job_name )
    
    wf.addTask( run_ovlp_filter_task )
    wf.refreshTargets(updateFreq = 1) #all            

    
    falcon_asm_done = makePypeLocalFile( os.path.join( falcon_asm_dir, "falcon_asm_done") )
    @PypeTask( inputs = {"ovlp_filter_done":ovlp_filter_done}, 
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
        script_dir = os.path.join( wd )
        script_fn =  os.path.join( script_dir ,"run_falcon_asm.sh" )
        
        script = []
        script.append( "source {install_prefix}/bin/activate".format(install_prefix = install_prefix) )
        script.append( "cd %s" % wd )
        script.append( "ln -sf %s/preads4falcon.fasta preads.fa" % self.parameters["pread_dir"])
        script.append( "ln -sf %s/las_files/rc_out_all ." % self.parameters["pread_dir"])
        script.append( "cat %s/las_files/*.ovl > preads.ovl" % self.parameters["pread_dir"])
        script.append( """falcon_asm_s.py preads.ovl preads.fa""" )
        script.append( """falcon_fixasm.py""" )
        script.append( """falcon_dedup.py""" )
        script.append( """falcon_ucns_data.py | falcon_utgcns.py > pa_cns.fa""")
        script.append( """touch %s\n""" % fn(self.falcon_asm_done))

        with open(script_fn, "w") as script_file:
            script_file.write("\n".join(script))

        job_name = self.URL.split("/")[-1]
        job_name += "-"+str(uuid.uuid1())[:8]
        job_data = {"job_name": job_name,
                    "cwd": wd,
                    "sge_option": " -pe smp 16 -q huasm",
                    "script_fn": script_fn }
        run_script(job_data, job_type = "SGE")
        wait_for_file( fn(self.ovlp_filter_done), task=self, job_name=job_name )
    
    wf.addTask( run_falcon_asm_task )
    wf.refreshTargets(updateFreq = 1) #all            
