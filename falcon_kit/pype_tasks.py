# PypeTask functions now need to be module-level.
from . import run_support as support
from . import bash # for scattering
#from pypeflow.simple_pwatcher_bridge import fn # not really needed
import collections
import json
import logging
import os.path
LOG = logging.getLogger(__name__)

def fn(p): return p

def system(call, check=False):
    LOG.debug('$(%s)' %repr(call))
    rc = os.system(call)
    msg = 'Call %r returned %d.' % (call, rc)
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
    support.make_fofn_abs(fn(self.i_fofn), fn(self.o_fofn))

def task_make_fofn_abs_preads(self):
    script_fn = 'noop.sh'
    open(script_fn, 'w').write('echo NOOP preads')
    self.generated_script_fn = script_fn
    support.make_fofn_abs(fn(self.i_fofn), fn(self.o_fofn))

def task_build_rdb(self):
    input_fofn_fn = fn(self.input_fofn)
    job_done = fn(self.rdb_build_done)
    db = fn(self.raw_reads_db)
    run_jobs = fn(self.run_jobs)
    remove(job_done, db, run_jobs)
    work_dir = self.parameters['work_dir']
    config = self.parameters['config']

    script_fn = os.path.join( work_dir, 'prepare_rdb.sh' )
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
    work_dir = self.parameters['work_dir']
    config = self.parameters['config']

    script_fn = os.path.join( work_dir, 'prepare_pdb.sh' )
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
    wd = self.parameters['wd']
    #self.las_fofn # TODO: Are there any implicit dependencies, or can we drop this?
    job_done = fn(self.db2falcon_done)
    preads4falcon_fn = fn(self.preads4falcon)
    preads_db = fn(self.preads_db)
    config = self.parameters['config']
    script_dir = os.path.join(wd)
    script_fn = os.path.join(script_dir ,'run_db2falcon.sh')
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
    wd = self.parameters['wd']
    #self.db2falcon_done
    db_file = fn(self.db_file)
    job_done = fn(self.falcon_asm_done)
    config = self.parameters['config']
    pread_dir = self.parameters['pread_dir']
    preads4falcon_fn = fn(self.preads4falcon)
    las_fofn_fn = fn(self.las_fofn)
    script_dir = os.path.join(wd)
    script_fn =  os.path.join(script_dir ,'run_falcon_asm.sh')
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
    daligner_script = self.parameters['daligner_script']
    job_uid = self.parameters['job_uid']
    cwd = self.parameters['cwd']
    db_prefix = self.parameters['db_prefix']
    config = self.parameters['config']
    script_dir = os.path.join(cwd)
    script_fn = os.path.join(script_dir , 'rj_%s.sh' % (job_uid))
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
    job_done = fn(self.job_done)
    gathered_las_fn = fn(self.gathered_las)
    script = self.parameters['merge_script']
    job_id = self.parameters['job_id'] # aka 'block'
    cwd = self.parameters['cwd']

    gathered_dict = read_gathered_las(gathered_las_fn)
    las_paths = gathered_dict[job_id]
    for las_path in las_paths:
        src = os.path.relpath(las_path, cwd)
        tgt = os.path.join(cwd, os.path.basename(las_path))
        LOG.debug('symlink {!r} -> {!r}'.format(src, tgt))
        os.symlink(src, tgt)

    config = self.parameters['config']

    script_dir = os.path.join( cwd )
    script_fn =  os.path.join( script_dir , 'rp_%05d.sh' % (job_id))
    args = {
        'script': script,
        'config': config,
        'job_done': job_done,
        'script_fn': script_fn,
    }
    support.run_las_merge(**args)
    self.generated_script_fn = script_fn

def task_run_consensus(self):
    las_fn = fn(self.las)
    db_fn = fn(self.db)
    out_file_fn = fn(self.out_file)
    out_done = 'out.done' #fn(self.out_done)
    job_id = self.parameters['job_id']
    cwd = self.parameters['cwd']
    config = self.parameters['config']
    prefix = self.parameters['prefix']
    p_id = int(job_id)
    script_dir = os.path.join(cwd)
    script_fn = os.path.join(script_dir, 'c_%05d.sh' % (p_id))
    #merge_job_dir = os.path.dirname(merged_las_fn)
    #las_fn = os.path.abspath('{merge_job_dir}/{prefix}.{job_id}.las'.format(**locals()))
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

def task_daligner_scatter(self):
    run_jobs_fn = self.run_jobs_fn
    db_build_done = self.db_build_done
    scatter_fn = self.scatter_fn
    par = self.parameters
    db_prefix = par['db_prefix']
    nblock = par['nblock']
    config = par['config']
    pread_aln = par['pread_aln'] # False  for raw_reads
    skip_checks = config.get('skip_checks')
    tasks = []
    LOG.info('Skip LAcheck after daligner? {}'.format(skip_checks))
    func = task_run_daligner
    func_name = '{}.{}'.format(func.__module__, func.__name__)
    for job_uid, script in bash.scripts_daligner(run_jobs_fn, db_prefix, db_build_done, nblock, pread_aln, skip_checks):
        job_done_fn = 'job_%s_done' %job_uid
        parameters =  {'daligner_script': script,
                       'job_uid': job_uid,
                       'config': config,
                       'sge_option': config['sge_option_da'],
                       'db_prefix': db_prefix}
        inputs = {'db_build_done': db_build_done}
        outputs = {'job_done': job_done_fn}
        python_function = func_name,
        URL = 'task://localhost/d_%s_%s' %(job_uid, db_prefix)
        daligner_task = {
                'inputs': inputs,
                'outputs': outputs,
                'parameters': parameters,
                'python_function': python_function,
                'URL': URL,
        }
        tasks.append(daligner_task)
    content = json.dumps(tasks, sort_keys=True, indent=4, separators=(',', ': '))
    open(scatter_fn, 'w').write(content)

def task_merge_scatter(self):
    run_jobs_fn = self.run_jobs
    gathered_las_fn = self.gathered_las
    scatter_fn = self.scattered
    par = self.parameters
    db_prefix = par['db_prefix']
    config = par['config']
    func = task_run_las_merge
    func_name = '{}.{}'.format(func.__module__, func.__name__)

    merge_scripts = bash.scripts_merge(config, db_prefix, run_jobs_fn)
    tasks = []
    for p_id, merge_script, merged_las_fn in merge_scripts:
        parameters =  {'merge_script': merge_script,
                       'job_id': p_id,
                       'config': config,
                       'sge_option': config['sge_option_la'],
                      }
        job_done_fn = 'm_%05d_done' % p_id
        inputs = {'gathered_las': gathered_las_fn}
        outputs = {'job_done': job_done_fn, # probably not needed anymore
                   'merged_las': merged_las_fn,
        }
        python_function = func_name,
        URL = 'task://localhost/m_%05d_%s' %(p_id, db_prefix)
        task_desc = {
                'inputs': inputs,
                'outputs': outputs,
                'parameters': parameters,
                'python_function': python_function,
                'URL': URL,
        }
        tasks.append(task_desc)

    content = json.dumps(tasks, sort_keys=True, indent=4, separators=(',', ': '))
    open(scatter_fn, 'w').write(content)
def task_consensus_scatter(self):
    scattered_fn = self.scattered
    gathered_fn = self.gathered
    db_fn = self.db
    wd = os.path.dirname(scattered_fn)
    par = self.parameters
    db_prefix = par['db_prefix']
    config = par['config']

    func = task_run_consensus
    func_name = '{}.{}'.format(func.__module__, func.__name__)
    basedir = os.path.dirname(wd) # by convention, since we want to preseve some old paths for now

    p_ids_merge_las = read_gathered_las(gathered_fn)
    tasks = []
    for p_id, las_fns in p_ids_merge_las.iteritems():
        assert len(las_fns) == 1, repr(las_fns)
        las_fn = las_fns[0] # since we know each merge-task is for a single block
        cns_label = 'cns_%05d' %int(p_id)
        #out_done_fn = '%s_done' % cns_label
        out_file_fn = '%s.fasta' % cns_label

        parameters =  {#'cwd': rdir,
                       'job_id': p_id,
                       'prefix': db_prefix,
                       'config': config,
                       'sge_option': config['sge_option_cns'],
        }
        inputs =  {'las': las_fn,
                   'db': db_fn,
        }
        outputs = {'out_file': out_file_fn,
                   #'out_done': out_done_fn,
        }
        python_function = func_name,
        URL = 'task://localhost/%s' %cns_label
        task_desc = {
                'inputs': inputs,
                'outputs': outputs,
                'parameters': parameters,
                'python_function': python_function,
                'URL': URL,
        }
        tasks.append(task_desc)
    content = json.dumps(tasks, sort_keys=True, indent=4, separators=(',', ': '))
    open(scattered_fn, 'w').write(content)

def task_daligner_gather(self):
    """Find all .las leaves so far.
    """
    out_dict = self.inputs
    gathered_fn = fn(self.gathered)
    nblock = self.parameters['nblock']
    LOG.debug('nblock=%d, out_dir:\n%s'%(nblock, out_dict))
    job_rundirs = [os.path.dirname(fn(dal_done)) for dal_done in out_dict.values()]
    with open(gathered_fn, 'w') as ofs:
        for block, las_path in support.daligner_gather_las(job_rundirs):
            ofs.write('{} {}\n'.format(block, las_path))

def task_cns_gather(self):
    fofn_fn = fn(self.preads_fofn)
    with open(fofn_fn,  'w') as f:
        for filename in sorted(fn(plf) for plf in self.inputs.itervalues()):
            print >>f, filename

def task_merge_gather(self):
    fofn_fn = fn(self.las_fofn)
    with open(fofn_fn,  'w') as f:
        # The keys are p_ids.
        for filename in sorted(fn(plf) for plf in self.inputs.itervalues()):
            print >>f, filename
    fopfn_fn = fn(self.las_fopfn)
    with open(fopfn_fn,  'w') as f:
        # The keys are p_ids.
        for filename, p_id in sorted((fn(plf), p_id) for (p_id, plf) in self.inputs.iteritems()):
            print >>f, p_id, filename
    #wdir = os.path.dirname(las_fofn_fn)
    #pread_dir = os.path.dirname(wdir) # by convention, for now
    # Generate las.fofn in run-dir. # No longer needed!
    #system('find {}/m_*/ -name "preads.*.las" >| {}'.format(pread_dir, las_fofn_fn))


def task_dump_rawread_ids(self):
    rawread_db = fn(self.rawread_db)
    rawread_id_file = fn(self.rawread_id_file)
    system("DBshow -n %s | tr -d '>' | LD_LIBRARY_PATH= awk '{print $1}' > %s" % (rawread_db, rawread_id_file))

def task_dump_pread_ids(self):
    pread_db = fn(self.pread_db)
    pread_id_file = fn(self.pread_id_file)
    system("DBshow -n %s | tr -d '>' | LD_LIBRARY_PATH= awk '{print $1}' > %s" % (pread_db, pread_id_file))

def task_generate_read_to_ctg_map(self):
    from .fc_asm_graph import AsmGraph
    rawread_id_file = fn(self.rawread_id_file)
    pread_id_file = fn(self.pread_id_file)
    read_to_contig_map = fn(self.read_to_contig_map)

    pread_did_to_rid = open(pread_id_file).read().split('\n')
    rid_to_oid = open(rawread_id_file).read().split('\n')

    asm_G = AsmGraph(fn(self.sg_edges_list),
                        fn(self.utg_data),
                        fn(self.ctg_paths))

    pread_to_contigs = {}

    with open(read_to_contig_map, 'w') as f:
        for ctg in asm_G.ctg_data:
            if ctg[-1] == 'R':
                continue
            ctg_g = asm_G.get_sg_for_ctg(ctg)
            for n in ctg_g.nodes():
                pid = int(n.split(':')[0])

                rid = pread_did_to_rid[pid].split('/')[1]
                rid = int(int(rid)/10)
                oid = rid_to_oid[rid]
                k = (pid, rid, oid)
                pread_to_contigs.setdefault(k, set())
                pread_to_contigs[k].add(ctg)


        for k in pread_to_contigs:
            pid, rid, oid = k
            for ctg in list(pread_to_contigs[ k ]):
                print >>f, '%09d %09d %s %s' % (pid, rid, oid, ctg)
