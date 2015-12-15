"""Most bash-scripting is generated here.
"""
from . import functional
import os

BASH='/bin/bash'


def mkdir(d):
    if not os.path.isdir(d):
        os.makedirs(d)

def write_script_and_wrapper(script, wrapper_fn, job_done, job_exit):
    """
    job_done/_exit should be either abspath or relative to dir of wrapper_fn.
    """
    wdir = os.path.dirname(wrapper_fn)
    #mkdir(wdir) # To avoid races, callers must do this.
    root, ext = os.path.splitext(os.path.basename(wrapper_fn))
    sub_script_bfn = root + '.sub' + ext

    # Try to avoid 'text file busy' in open():
    #os.listdir(wdir)
    #os.system('touch {}'.format(os.path.join(wdir, sub_script_bfn)))
    os.system('rm -f {}'.format(os.path.join(wdir, sub_script_bfn)))

    with open(os.path.join(wdir, sub_script_bfn), 'w') as ofs:
        # We use shebang + chmod so we can see the sub-script in 'top'.
        # In order to avoid '/bin/bash: bad interpreter: Text file busy',
        # we 'touch' the sub-script after chmod.
        #   http://superuser.com/questions/934300/bin-bash-bad-interpreter-text-file-busy-even-though-the-file-editor-closed
        ofs.write('#!{}\n'.format(BASH))
        ofs.write('set -vex\n')
        ofs.write(script)
    wrapper = """
set -vex
cd {wdir}
trap 'touch {job_exit}' EXIT
ls -il {sub_script_bfn}
hostname
chmod +x {sub_script_bfn}
touch {sub_script_bfn}
ls -il {sub_script_bfn}
time ./{sub_script_bfn}
touch {job_done}
""".format(**locals())
    with open(wrapper_fn, 'w') as ofs:
        ofs.write(wrapper)
    return job_done, job_exit

def script_build_rdb(config, input_fofn_fn, run_jobs_fn):
    """
    raw_reads.db will be output into CWD, and might already exist.
    run_jobs_bfn will be output into CWD.
    """
    last_block = 1
    new_db = True
    if os.path.exists("raw_reads.db"):
        with open("raw_reads.db") as f:
            for l in f:
                l = l.strip().split()
                if l[0] == "blocks" and l[1] == "=":
                    last_block = int(l[2])
                    new_db = False
                    break

    DBsplit = 'DBsplit {pa_DBsplit_option} raw_reads'.format(**config) if new_db else ''
    openending = config["openending"]
    if openending == True:
        count = """$(cat raw_reads.db | awk '$1 == "blocks" {print $3-1}')"""
    else:
        count = """$(cat raw_reads.db | awk '$1 == "blocks" {print $3}')"""
    params = dict(config)
    params.update(locals())
    script = """\
fasta2DB -v raw_reads -f{input_fofn_fn}
{DBsplit}
LB={count}
HPCdaligner {pa_HPCdaligner_option} -H{length_cutoff} raw_reads {last_block}-$LB >| {run_jobs_fn}
""".format(**params)
    return script

def script_build_pdb(config, input_fofn_fn, run_jobs_fn):
    params = dict(config)
    params.update(locals())
    script = """\
fasta2DB -v preads -f{input_fofn_fn}
DBsplit {ovlp_DBsplit_option} preads
HPCdaligner {ovlp_HPCdaligner_option} -H{length_cutoff_pr} preads >| {run_jobs_fn}
""".format(**params)
    return script

def script_run_DB2Falcon(config):
    """Run in pread_dir.
    """
    params = dict(config)
    params.update(locals())
    script = """\
# Given preads.db,
# write preads4falcon.fasta, in 1-preads_ovl:
time DB2Falcon -U preads
""".format(**params)
    return script

def scripts_daligner(run_jobs_fn, db_prefix, rdb_build_done, pread_aln=False):
    """Yield job_uid, bash
    """
    scripts = {}
    xform_script = functional.get_script_xformer(pread_aln)
    db_dir = os.path.dirname(run_jobs_fn)
    job_descs = functional.get_daligner_job_descriptions(open(run_jobs_fn), db_prefix)
    for i, (desc, bash) in enumerate(job_descs.iteritems()):
        job_uid = '%04x' %i
        daligner_cmd = xform_script(bash)
        bash = """
db_dir={db_dir}
ln -sf ${{db_dir}}/.{db_prefix}.bps .
ln -sf ${{db_dir}}/.{db_prefix}.idx .
ln -sf ${{db_dir}}/{db_prefix}.db .
{daligner_cmd}
#rm -f *.C?.las
#rm -f *.N?.las
""".format(db_dir=db_dir, db_prefix=db_prefix, daligner_cmd=daligner_cmd)
        yield job_uid, bash

def scripts_merge(config, db_prefix, run_jobs_fn):
    """Yield p_id, bash
    """
    mjob_data = {}
    with open(run_jobs_fn) as f:
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
        bash_lines = mjob_data[p_id]

        #support.make_dirs("%s/m_%05d" % (wd, p_id))
        #support.make_dirs("%s/preads" % (wd))
        #support.make_dirs("%s/las_files" % (wd))
        #merge_script_file = os.path.abspath("%s/m_%05d/m_%05d.sh" % (wd, p_id, p_id))

        script = []
        for line in bash_lines:
            script.append(line.replace('&&', ';'))
        script.append("mkdir -p ../las_files")
        script.append("ln -sf ../m_%05d/%s.%d.las ../las_files" % (p_id, db_prefix, p_id))
        script.append("ln -sf ./m_%05d/%s.%d.las .. " % (p_id, db_prefix, p_id))
        yield p_id, '\n'.join(script + [''])

def script_run_consensus(config, db_fn, las_fn, out_file_bfn):
    """config: falcon_sense_option, length_cutoff
    """
    params = dict(config)
    params.update(locals())
    if config["falcon_sense_skip_contained"]:
        pipe = """LA4Falcon -H{length_cutoff} -fso {db_fn} {las_fn} | """
    else:
        pipe = """LA4Falcon -H{length_cutoff}  -fo {db_fn} {las_fn} | """
    pipe += """fc_consensus {falcon_sense_option} >| {out_file_bfn}"""

    script = """
set -o pipefail
%s
""" %pipe
    return script.format(**params)

def script_run_falcon_asm(config, pread_dir, db_file):
    params = dict(config)
    params.update(locals())
    script = """\
# Generate las.fofn:
find {pread_dir}/las_files -name "*.las" >| las.fofn

# Given, las.fofn,
# write preads.ovl:
time fc_ovlp_filter --db {db_file} --fofn las.fofn {overlap_filtering_setting} --min_len {length_cutoff_pr} >| preads.ovl

ln -sf {pread_dir}/preads4falcon.fasta .
# TODO: Figure out which steps need preads4falcon.fasta.

# Given preads.ovl,
# write sg_edges_list, c_path, utg_data, ctg_paths.
time fc_ovlp_to_graph preads.ovl --min_len {length_cutoff_pr} >| fc_ovlp_to_graph.log

# Given sg_edges_list, utg_data, ctg_paths,
# Write p_ctg.fa and a_ctg_all.fa,
# plus a_ctg_base.fa, p_ctg_tiling_path, a_ctg_tiling_path, a_ctg_base_tiling_path:
time fc_graph_to_contig

# Given a_ctg_all.fa, write a_ctg.fa:
time fc_dedup_a_tigs
"""
    return script.format(**params)

