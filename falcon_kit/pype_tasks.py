from __future__ import absolute_import
from __future__ import print_function


from future.utils import viewitems
from future.utils import itervalues
# PypeTask functions now need to be module-level.
from . import run_support as support
from . import bash  # for scattering
# from pypeflow.simple_pwatcher_bridge import fn # not really needed
import collections
import json
import logging
import os.path
LOG = logging.getLogger(__name__)


#TASK_LAS_MERGE_SCATTER_SCRIPT = """\
#python -m falcon_kit.mains.las_merge_scatter --db-prefix={params.db_prefix} --stage={params.stage} --run-jobs-fn={input.run_jobs} --gathered-las-fn={input.gathered_las} --wildcards={params.wildcards} --scattered-fn={output.scattered}
#"""
TASK_LAS_MERGE_SPLIT_SCRIPT = """\
python -m falcon_kit.mains.las_merge_split --wildcards={params.wildcards} --db-prefix={params.db_prefix} --run-jobs-fn={input.run_jobs} --gathered-las-fn={input.gathered_las} --split-fn={output.split} --bash-template-fn={output.bash_template}
"""
TASK_LAS_MERGE_SCRIPT = """\
# Note: HPC.daligner chooses a merged filename in its generated script, so we will symlink to it.
python -m falcon_kit.mains.las_merge --las-paths-fn={input.las_paths} --merge-script-fn={input.merge_script} --las-merged-fn-fn={input.merged_las_json} --las-merged-symlink-fn={output.merged_las} --job-done-fn={output.job_done} --p-id-fn={output.p_id} --p-id-num={params.p_id_num}
"""
TASK_LAS_MERGE_GATHER_SCRIPT = """\
python -m falcon_kit.mains.las_merge_gather --gathered-fn={input.gathered} --p-id2las-fn={output.p_id2las} --las-fn={output.las}
"""
#TASK_CONSENSUS_SCATTER_SCRIPT = """\
#python -m falcon_kit.mains.consensus_scatter --las-fopfn-fn={input.las_fopfn} --db-fn={input.raw_reads_db} --length-cutoff-fn={input.length_cutoff} --config-fn={input.config} --wildcards={params.wildcards} --scattered-fn={output.scattered}
#"""
TASK_CONSENSUS_SPLIT_SCRIPT = """\
python -m falcon_kit.mains.consensus_split --wildcards={params.wildcards} --p-id2las-fn={input.p_id2las} --db-fn={input.raw_reads_db} --length-cutoff-fn={input.length_cutoff} --config-fn={input.config} --split-fn={output.split} --bash-template-fn={output.bash_template}
"""
TASK_CONSENSUS_TASK_SCRIPT = """\
python -m falcon_kit.mains.consensus_task --nproc={params.pypeflow_nproc} --las-fn={input.las} --db-fn={input.db} --length-cutoff-fn={input.length_cutoff} --config-fn={input.config} --fasta-fn={output.fasta}
"""
TASK_CONSENSUS_GATHER_SCRIPT = """\
python -m falcon_kit.mains.consensus_gather_fasta_fofn --gathered-fn={input.gathered} --preads-fofn-fn={output.preads_fofn}
"""
TASK_REPORT_PRE_ASSEMBLY_SCRIPT = """\
python -m falcon_kit.mains.task_report_pre_assembly --config-fn={input.config} --length-cutoff-fn={input.length_cutoff} --raw-reads-db-fn={input.raw_reads_db} --preads-fofn-fn={input.preads_fofn} --pre-assembly-report-fn={output.pre_assembly_report}
"""

# Old
TASK_BUILD_RDB_SCRIPT = """\
python -m falcon_kit.mains.build_rdb --input-fofn-fn={input.raw_reads_fofn} --config-fn={input.config} --run-jobs-fn={output.run_jobs} --job-done-fn={output.db_build_done}
touch {output.db_build_done}
"""
# Old
TASK_BUILD_PDB_SCRIPT = """\
python -m falcon_kit.mains.build_pdb --input-fofn-fn={input.preads_fofn} --config-fn={input.config} --run-jobs-fn={output.run_jobs} --job-done-fn={output.db_build_done}
# TODO: Verify that input.preads_db exists.
touch {output.db_build_done}
"""

TASK_DB_BUILD_SCRIPT = """\
python -m falcon_kit.mains.dazzler --config-fn={input.config} --db-fn={output.db}  build --input-fofn-fn={input.input_fofn} --length-cutoff-fn={output.length_cutoff}
# TODO: Verify that db exists.
#ln -sf {output.length_cutoff} length_cutoff
"""
TASK_DB_TAN_SPLIT_SCRIPT = """\
python -m falcon_kit.mains.dazzler --config={input.config} --db={input.db}  tan-split --split={output.split} --bash-template={output.bash_template}
"""
TASK_DB_TAN_APPLY_SCRIPT = """\
python -m falcon_kit.mains.dazzler --config={input.config} --db={input.db}  tan-apply --script={input.script} --job-done={output.job_done}
"""
TASK_DB_TAN_COMBINE_SCRIPT = """\
python -m falcon_kit.mains.dazzler --config={input.config} --db={input.db}  tan-combine --gathered={input.gathered} --new-db={output.new_db}
"""
TASK_DB_DALIGNER_SPLIT_SCRIPT = """\
python -m falcon_kit.mains.dazzler --config={input.config} --db={input.db} --nproc={params.pypeflow_nproc}  daligner-split --wildcards={params.wildcards} --length-cutoff-fn={input.length_cutoff} --split-fn={output.split} --bash-template-fn={output.bash_template}
"""
TASK_DB_DALIGNER_APPLY_SCRIPT = """\
python -m falcon_kit.mains.dazzler --config={input.config} --db={input.db}  daligner-apply --script={input.script} --job-done={output.job_done}
"""
TASK_DB_DALIGNER_COMBINE_SCRIPT = """\
python -m falcon_kit.mains.dazzler --config={input.config} --db={input.db}  daligner-combine --gathered={input.gathered} --las-paths-fn={output.las_paths}
"""
TASK_DB_LAMERGE_SPLIT_SCRIPT = """\
python -m falcon_kit.mains.dazzler --config={input.config}                  merge-split --db-prefix={params.db_prefix} --las-paths={input.las_paths} --split-fn={output.split} --bash-template-fn={output.bash_template}
"""
TASK_DB_LAMERGE_APPLY_SCRIPT = """\
python -m falcon_kit.mains.dazzler --config={input.config}                  merge-apply --las-paths={input.las_paths} --las-fn={output.las_fn}
"""
TASK_DB_LAMERGE_COMBINE_SCRIPT = """\
python -m falcon_kit.mains.dazzler --config={input.config}                  merge-combine --gathered={input.gathered} --las-paths-fn={output.las_paths} --block2las-fn={output.block2las}
"""
#TASK_DALIGNER_SCATTER_SCRIPT = """\
#python -m falcon_kit.mains.daligner_scatter --run-jobs-fn={input.run_jobs} --db-prefix={params.db_prefix} --db-fn={input.db} --skip-checks={params.skip_checks} --pread-aln={params.pread_aln} --stage={params.stage} --wildcards={params.wildcards} --scattered-fn={output.scattered}
#"""
TASK_DALIGNER_SPLIT_SCRIPT = """\
python -m falcon_kit.mains.daligner_split --nproc={params.pypeflow_nproc} --wildcards={params.wildcards} --db-prefix={params.db_prefix} --skip-checks={params.skip_checks} --pread-aln={params.pread_aln} --run-jobs-fn={input.run_jobs} --db-fn={input.db} --split-fn={output.split} --bash-template-fn={output.bash_template}
"""
TASK_DALIGNER_SCRIPT = """\
# Note: HPC.daligner chooses a merged filename in its generated script, so we will symlink to it.
python -m falcon_kit.mains.daligner --daligner-settings-fn={input.daligner_settings} --daligner-script-fn={input.daligner_script} --job-done-fn={output.job_done}
"""

TASK_DALIGNER_FIND_LAS_SCRIPT = """\
python -m falcon_kit.mains.daligner_gather_las_list --gathered-fn={input.gathered} --las-paths-fn={output.las_paths}
"""
TASK_DUMP_RAWREAD_IDS_SCRIPT = """\
DBshow -n {input.rawread_db} | tr -d '>' | LD_LIBRARY_PATH= awk '{{print $1}}' > {output.rawread_id_file}
"""
TASK_DUMP_PREAD_IDS_SCRIPT = """\
DBshow -n {input.pread_db} | tr -d '>' | LD_LIBRARY_PATH= awk '{{print $1}}' > {output.pread_id_file}
"""
TASK_GENERATE_READ_TO_CTG_MAP_SCRIPT = """\
python -m falcon_kit.mains.generate_read_to_ctg_map --rawread-id={input.rawread_id_file} --pread-id={input.pread_id_file} --sg-edges-list={input.sg_edges_list} --utg-data={input.utg_data} --ctg-paths={input.ctg_paths} --output={output.read_to_contig_map}
"""
TASK_RUN_DB_TO_FALCON_SCRIPT = """\
# Given preads.db,
# write preads4falcon.fasta (implicitly) in CWD.
time DB2Falcon -U {input.preads_db}
[ -f {output.preads4falcon} ] || exit 1
touch {output.job_done}
"""
TASK_RUN_FALCON_ASM_SCRIPT = """\
# Given, las_fofn.json,
# write preads.ovl:

# mobs uses binwrappers, so it does not see our "entry-points".
# So, after dropping "src/py_scripts/*.py", we can call these via python -m:

time python -m falcon_kit.mains.ovlp_filter --db {input.db_file} --las-fofn {input.las_fofn} {params.overlap_filtering_setting} --min-len {params.length_cutoff_pr} --out-fn preads.ovl

ln -sf {input.preads4falcon_fasta} ./preads4falcon.fasta

# Given preads.ovl,
# write sg_edges_list, c_path, utg_data, ctg_paths.
time python -m falcon_kit.mains.ovlp_to_graph {params.fc_ovlp_to_graph_option} --overlap-file preads.ovl >| fc_ovlp_to_graph.log

# Given sg_edges_list, utg_data, ctg_paths, preads4falcon.fasta,
# write p_ctg.fa and a_ctg_all.fa,
# plus a_ctg_base.fa, p_ctg_tiling_path, a_ctg_tiling_path, a_ctg_base_tiling_path:
time python -m falcon_kit.mains.graph_to_contig

# Given a_ctg_all.fa, write a_ctg.fa:
time python -m falcon_kit.mains.dedup_a_tigs

# Collect all info needed to format the GFA-1 and GFA-2 representations of
# the assembly graphs.
time python -m falcon_kit.mains.collect_pread_gfa >| asm.gfa.json
time python -m falcon_kit.mains.collect_pread_gfa --add-string-graph >| sg.gfa.json
time python -m falcon_kit.mains.collect_contig_gfa >| contig.gfa.json

# Output the assembly pread graph.
time python -m falcon_kit.mains.gen_gfa_v1 asm.gfa.json >| asm.gfa
time python -m falcon_kit.mains.gen_gfa_v2 asm.gfa.json >| asm.gfa2

# Output the string graph.
time python -m falcon_kit.mains.gen_gfa_v1 sg.gfa.json >| sg.gfa
time python -m falcon_kit.mains.gen_gfa_v2 sg.gfa.json >| sg.gfa2

# Output the contig graph with associate contigs attached to each primary contig.
time python -m falcon_kit.mains.gen_gfa_v2 contig.gfa.json >| contig.gfa2

#rm -f ./preads4falcon.fasta

touch {output.falcon_asm_done}
"""


def fn(p): return p


def system(call, check=False):
    LOG.debug('$(%s)' % repr(call))
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

    script_fn = os.path.join(work_dir, 'prepare_rdb.sh')
    args = {
        'input_fofn_fn': input_fofn_fn,
        'config': config,
        'job_done': job_done,
        'script_fn': script_fn,
        'run_jobs_fn': run_jobs,
    }
    support.build_rdb(**args)
    self.generated_script_fn = script_fn


# essential the same as build_rdb() but the subtle differences are tricky to consolidate to one function
def task_build_pdb(self):
    input_fofn_fn = fn(self.preads_fofn)
    job_done = fn(self.pdb_build_done)
    db = fn(self.preads_db)
    run_jobs = fn(self.run_jobs)
    remove(job_done, db, run_jobs)
    work_dir = self.parameters['work_dir']
    config = self.parameters['config']

    script_fn = os.path.join(work_dir, 'prepare_pdb.sh')
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
    # self.las_fofn # TODO: Are there any implicit dependencies, or can we drop this?
    job_done = fn(self.db2falcon_done)
    preads4falcon_fn = fn(self.preads4falcon)
    preads_db = fn(self.preads_db)
    config = self.parameters['config']
    script_dir = os.path.join(wd)
    script_fn = os.path.join(script_dir, 'run_db2falcon.sh')
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
    # self.db2falcon_done
    db_file = fn(self.db_file)
    job_done = fn(self.falcon_asm_done)
    config = self.parameters['config']
    pread_dir = self.parameters['pread_dir']
    preads4falcon_fn = fn(self.preads4falcon)
    las_fofn_fn = fn(self.las_fofn)
    script_dir = os.path.join(wd)
    script_fn = os.path.join(script_dir, 'run_falcon_asm.sh')
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
    genome_length = int(cfg.get('genome_size', 0))  # different name in falcon
    length_cutoff = int(cfg['length_cutoff'])
    # Update length_cutoff if auto-calc (when length_cutoff is negative).
    # i_length_cutoff_fn was created long ago, so no filesystem issues.
    length_cutoff = support.get_length_cutoff(
        length_cutoff, i_length_cutoff_fn)
    cwd = self.parameters['cwd']
    script_fn = os.path.join(cwd, 'run_report_pre_assembly.sh')
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
    script_fn = os.path.join(script_dir, 'rj_%s.sh' % (job_uid))
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
    # LOG.warning('path={!r}, result={}'.format(
    #    path, pprint.pformat(result)))
    return result


def task_run_las_merge(self):
    job_done = fn(self.job_done)
    gathered_las_fn = fn(self.gathered_las)
    script = self.parameters['merge_script']
    job_id = self.parameters['job_id']  # aka 'block'
    cwd = self.parameters['cwd']

    gathered_dict = read_gathered_las(gathered_las_fn)
    las_paths = gathered_dict[job_id]
    for las_path in las_paths:
        assert os.path.isabs(las_path)
        if os.path.commonprefix([las_path, cwd]) == '/':
            src = las_path
        else:
            src = os.path.relpath(las_path, cwd)
        tgt = os.path.join(cwd, os.path.basename(las_path))
        LOG.debug('symlink {!r} <- {!r}'.format(src, tgt))
        if os.path.lexists(tgt):
            os.unlink(tgt)
        os.symlink(src, tgt)

    config = self.parameters['config']

    script_dir = os.path.join(cwd)
    script_fn = os.path.join(script_dir, 'rp_%05d.sh' % (job_id))
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
    out_done = 'out.done'  # fn(self.out_done)
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
    pread_aln = par['pread_aln']  # False  for raw_reads
    skip_checks = config.get('skip_checks')
    tasks = []
    LOG.info('Skip LAcheck after daligner? {}'.format(skip_checks))
    func = task_run_daligner
    func_name = '{}.{}'.format(func.__module__, func.__name__)
    for job_uid, script in bash.scripts_daligner(run_jobs_fn, db_prefix, db_build_done, nblock, pread_aln, skip_check=skip_checks):
        job_done_fn = 'job_%s_done' % job_uid
        parameters = {'daligner_script': script,
                      'job_uid': job_uid,
                      'config': config,
                      'sge_option': config['sge_option_da'],
                      'db_prefix': db_prefix}
        inputs = {'db_build_done': db_build_done}
        outputs = {'job_done': job_done_fn}
        python_function = func_name,
        URL = 'task://localhost/d_%s_%s' % (job_uid, db_prefix)
        daligner_task = {
            'inputs': inputs,
            'outputs': outputs,
            'parameters': parameters,
            'python_function': python_function,
            'URL': URL,
        }
        tasks.append(daligner_task)
    content = json.dumps(tasks, sort_keys=True, indent=4,
                         separators=(',', ': '))
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
        parameters = {'merge_script': merge_script,
                      'job_id': p_id,
                      'config': config,
                      'sge_option': config['sge_option_la'],
                      }
        job_done_fn = 'm_%05d_done' % p_id
        inputs = {'gathered_las': gathered_las_fn}
        outputs = {'job_done': job_done_fn,  # probably not needed anymore
                   'merged_las': merged_las_fn,
                   }
        python_function = func_name,
        URL = 'task://localhost/m_%05d_%s' % (p_id, db_prefix)
        task_desc = {
            'inputs': inputs,
            'outputs': outputs,
            'parameters': parameters,
            'python_function': python_function,
            'URL': URL,
        }
        tasks.append(task_desc)

    content = json.dumps(tasks, sort_keys=True, indent=4,
                         separators=(',', ': '))
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
    # by convention, since we want to preseve some old paths for now
    basedir = os.path.dirname(wd)

    p_ids_merge_las = read_gathered_las(gathered_fn)
    tasks = []
    for (p_id, las_fns) in viewitems(p_ids_merge_las):
        assert len(las_fns) == 1, repr(las_fns)
        # since we know each merge-task is for a single block
        las_fn = las_fns[0]
        cns_label = 'cns_%05d' % int(p_id)
        #out_done_fn = '%s_done' % cns_label
        out_file_fn = '%s.fasta' % cns_label

        parameters = {  # 'cwd': rdir,
            'job_id': p_id,
            'prefix': db_prefix,
            'config': config,
            'sge_option': config['sge_option_cns'],
        }
        inputs = {'las': las_fn,
                  'db': db_fn,
                  }
        outputs = {'out_file': out_file_fn,
                   #'out_done': out_done_fn,
                   }
        python_function = func_name,
        URL = 'task://localhost/%s' % cns_label
        task_desc = {
            'inputs': inputs,
            'outputs': outputs,
            'parameters': parameters,
            'python_function': python_function,
            'URL': URL,
        }
        tasks.append(task_desc)
    content = json.dumps(tasks, sort_keys=True, indent=4,
                         separators=(',', ': '))
    open(scattered_fn, 'w').write(content)


def task_daligner_gather(self):
    """Find all .las leaves so far.
    """
    out_dict = self.inputs
    gathered_fn = fn(self.gathered)
    nblock = self.parameters['nblock']
    LOG.debug('nblock=%d, out_dir:\n%s' % (nblock, out_dict))
    job_rundirs = [os.path.dirname(fn(dal_done))
                   for dal_done in itervalues(out_dict)]
    with open(gathered_fn, 'w') as ofs:
        for block, las_path in support.daligner_gather_las(job_rundirs):
            ofs.write('{} {}\n'.format(block, las_path))


def task_cns_gather(self):
    fofn_fn = fn(self.preads_fofn)
    with open(fofn_fn,  'w') as f:
        for filename in sorted(fn(plf) for plf in itervalues(self.inputs)):
            print(filename, file=f)


def task_merge_gather(self):
    fofn_fn = fn(self.las_fofn)
    with open(fofn_fn,  'w') as f:
        # The keys are p_ids.
        for filename in sorted(fn(plf) for plf in itervalues(self.inputs)):
            print(filename, file=f)
    fopfn_fn = fn(self.las_fopfn)
    with open(fopfn_fn,  'w') as f:
        # The keys are p_ids.
        for (filename, p_id) in sorted((fn(plf), p_id) for (p_id, plf) in viewitems(self.inputs)):
            print(p_id, filename, file=f)
    #wdir = os.path.dirname(las_fofn_fn)
    # pread_dir = os.path.dirname(wdir) # by convention, for now
    # Generate las.fofn in run-dir. # No longer needed!
    #system('find {}/m_*/ -name "preads.*.las" >| {}'.format(pread_dir, las_fofn_fn))


def task_dump_rawread_ids(self):
    rawread_db = fn(self.rawread_db)
    rawread_id_file = fn(self.rawread_id_file)
    input = object()
    input.rawread_db = rawread_db
    output = object()
    output.rawread_id_file = rawread_id_file
    system(TASK_DUMP_RAWREAD_IDS_SCRIPT.format(**locals()))


def task_dump_pread_ids(self):
    pread_db = fn(self.pread_db)
    pread_id_file = fn(self.pread_id_file)
    input = object()
    input.pread_db = pread_db
    output = object()
    output.pread_id_file = pread_id_file
    system(TASK_DUMP_PREAD_IDS_SCRIPT.format(**locals()))


def task_generate_read_to_ctg_map(self):
    input = object()
    input.rawread_id_file = fn(self.rawread_id_file)
    input.pread_id_file = fn(self.pread_id_file)
    input.sg_edges_list = fn(self.sg_edges_list)
    input.utg_data = fn(self.utg_data)
    input.ctg_paths = fn(self.ctg_paths)
    output = object()
    output.read_to_contig_map = fn(self.read_to_contig_map)
    system(TASK_GENERATE_READ_TO_CTG_MAP_SCRIPT.format(**locals()))
