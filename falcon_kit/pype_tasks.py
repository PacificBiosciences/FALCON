# PypeTask functions now need to be module-level.
from . import run_support as support
from pypeflow.simple_pwatcher_bridge import fn # not really needed
import collections
import logging
import os.path
LOG = logging.getLogger(__name__)

def system(call, check=False):
    LOG.debug('$(%s)' %repr(call))
    rc = os.system(call)
    msg = "Call %r returned %d." % (call, rc)
    if rc:
        LOG.warning(msg)
        if check:
            raise Exception(msg)
    else:
        LOG.debug(msg)
    return rc

def remove(*fns):
    for fn in fns:
        if os.path.exists(fn):
            os.remove(fn)
        assert not os.path.exists(fn)

# Someday, we can drop mkdir() in these.
def mkdir(d):
    if not os.path.isdir(d):
        os.makedirs(d)

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
    LOG.info('Report inputs: {}'.format(repr(kwds)))
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
    #LOG.warning('path={!r}, result={}'.format(
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
        LOG.debug('symlink {!r} -> {!r}'.format(src, tgt))
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

def task_daligner_gather(self):
    """Find all .las leaves so far.
    """
    out_dict = self.inputs
    gathered_fn = fn(self.gathered)
    nblock = self.parameters['nblock']
    LOG.debug('nblock=%d, out_dir:\n%s'%(nblock, out_dict))
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
