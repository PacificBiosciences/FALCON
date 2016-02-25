"""Most bash-scripting is generated here.
"""
from . import functional
import os
import traceback

BASH='/bin/bash'
BUG_avoid_Text_file_busy=True
# http://stackoverflow.com/questions/1384398/usr-bin-perl-bad-interpreter-text-file-busy/

def mkdir(d):
    if not os.path.isdir(d):
        os.makedirs(d)

def make_executable(path):
    """http://stackoverflow.com/questions/12791997/how-do-you-do-a-simple-chmod-x-from-within-python
    """
    mode = os.stat(path).st_mode
    mode |= (mode & 0444) >> 2    # copy R bits to X
    os.chmod(path, mode)

def write_script_and_wrapper(script, wrapper_fn, job_done):
    """
    Write script to a filename based on wrapper_fn, in same directory.
    Write wrapper to call it,
     and to write job_done on success and job_done.exit on any exit.
    'job_done' should be either abspath or relative to dir of wrapper_fn.
    """
    job_exit = job_done + '.exit'
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
    make_executable(os.path.join(wdir, sub_script_bfn))
    if BUG_avoid_Text_file_busy:
        exe = BASH
    else:
        # We prefer to run via shebang b/c we want the script-name to appear to 'top',
        # but some users have a problem with that, e.g.
        #   https://github.com/PacificBiosciences/FALCON/issues/269
        # Another idea never worked reliably:
        #chmod +x {sub_script_bfn}
        #touch {sub_script_bfn}
        # We are trying to avoid this problem:
        #   /bin/bash: bad interpreter: Text file busy
        exe = ''

    wrapper = """
set -vex
cd {wdir}
trap 'touch {job_exit}' EXIT
ls -il {sub_script_bfn}
hostname
ls -il {sub_script_bfn}
time {exe} ./{sub_script_bfn}
touch {job_done}
"""
    wrapper = wrapper.format(**locals())
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
fasta2DB -pfakemoviename -v raw_reads -f{input_fofn_fn}
{DBsplit}
LB={count}
rm -f {run_jobs_fn}
HPCdaligner {pa_HPCdaligner_option} -H{length_cutoff} raw_reads {last_block}-$LB >| {run_jobs_fn}
""".format(**params)
    return script

def script_build_pdb(config, input_fofn_fn, run_jobs_fn):
    params = dict(config)
    params.update(locals())
    script = """\
fasta2DB -pfakemoviename -v preads -f{input_fofn_fn}
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
    try:
        job_descs = functional.get_daligner_job_descriptions(open(run_jobs_fn), db_prefix)
    except Exception:
        raise Exception('Could not parse job descriptions from file "{}":\n{}'.format(
            run_jobs_fn, traceback.format_exc()))
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

def script_run_falcon_asm(config, las_fofn_fn, preads4falcon_fasta_fn, db_file_fn):
    params = dict(config)
    params.update(locals())
    script = """\
# Given, las.fofn,
# write preads.ovl:
time fc_ovlp_filter --db {db_file_fn} --fofn {las_fofn_fn} {overlap_filtering_setting} --min_len {length_cutoff_pr} >| preads.ovl

ln -sf {preads4falcon_fasta_fn} ./preads4falcon.fasta

# Given preads.ovl,
# write sg_edges_list, c_path, utg_data, ctg_paths.
time fc_ovlp_to_graph {fc_ovlp_to_graph_option} preads.ovl >| fc_ovlp_to_graph.log

# Given sg_edges_list, utg_data, ctg_paths, preads4falcon.fasta,
# write p_ctg.fa and a_ctg_all.fa,
# plus a_ctg_base.fa, p_ctg_tiling_path, a_ctg_tiling_path, a_ctg_base_tiling_path:
time fc_graph_to_contig

rm -f ./preads4falcon.fasta

# Given a_ctg_all.fa, write a_ctg.fa:
time fc_dedup_a_tigs
"""
    return script.format(**params)

