"""Most bash-scripting is generated here.
"""
from __future__ import absolute_import


from future.utils import viewitems
from . import functional
import functools
import getpass
import logging
import md5
import os
import tempfile
import traceback

LOG = logging.getLogger(__name__)
BASH = '/bin/bash'
BUG_avoid_Text_file_busy = True
# http://stackoverflow.com/questions/1384398/usr-bin-perl-bad-interpreter-text-file-busy/


def mkdir(d):
    if not os.path.isdir(d):
        os.makedirs(d)


def make_executable(path):
    """http://stackoverflow.com/questions/12791997/how-do-you-do-a-simple-chmod-x-from-within-python
    """
    mode = os.stat(path).st_mode
    mode |= (mode & 0o444) >> 2    # copy R bits to X
    os.chmod(path, mode)


def write_sub_script(ofs, script):
    # We use shebang + chmod so we can see the sub-script in 'top'.
    # In order to avoid '/bin/bash: bad interpreter: Text file busy',
    # we 'touch' the sub-script after chmod.
    #   http://superuser.com/questions/934300/bin-bash-bad-interpreter-text-file-busy-even-though-the-file-editor-closed
    ofs.write('#!{}\n'.format(BASH))
    ofs.write('set -vex\n')
    ofs.write(script)

    if BUG_avoid_Text_file_busy:
        exe = BASH
    else:
        # We prefer to run via shebang b/c we want the script-name to appear to 'top',
        # but some users have a problem with that, e.g.
        #   https://github.com/PacificBiosciences/FALCON/issues/269
        # Another idea never worked reliably:
        # chmod +x {sub_script_bfn}
        # touch {sub_script_bfn}
        # We are trying to avoid this problem:
        #   /bin/bash: bad interpreter: Text file busy
        exe = ''
    return exe


def write_script(script, script_fn, job_done_fn=None):
    if job_done_fn:
        script += '\ntouch {}\n'.format(job_done_fn)
    with open(script_fn, 'w') as ofs:
        exe = write_sub_script(ofs, script)


def write_script_and_wrapper_top(script, wrapper_fn, job_done):
    """NOT USED.
    Write script to a filename based on wrapper_fn, in same directory.
    Write wrapper to call it,
     and to write job_done on success and job_done.exit on any exit.
    'job_done' should be either abspath or relative to dir of wrapper_fn.
    """
    job_exit = job_done + '.exit'
    wdir = os.path.abspath(os.path.dirname(wrapper_fn))
    # mkdir(wdir) # To avoid races, callers must do this.

    root, ext = os.path.splitext(os.path.basename(wrapper_fn))
    sub_script_bfn = root + '.sub' + ext
    with open(os.path.join(wdir, sub_script_bfn), 'w') as ofs:
        exe = write_sub_script(ofs, script)
    make_executable(os.path.join(wdir, sub_script_bfn))

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


def select_rundir(tmpdir, wdir, script):
    user = getpass.getuser()
    digest = md5.md5(script).hexdigest()
    return '{}/{}/falcontmp/{}/{}'.format(tmpdir, user, wdir, digest)


def write_script_and_wrapper_for_tmp(tmpdir, script, wrapper_fn, job_done):
    """NOT USED. This will be done in pypeFLOW.
    """
    wdir = os.path.dirname(os.path.abspath(wrapper_fn))
    root, ext = os.path.splitext(os.path.basename(wrapper_fn))
    sub_script_bfn = root + '.xsub' + ext
    with open(os.path.join(wdir, sub_script_bfn), 'w') as ofs:
        exe = write_sub_script(ofs, script)
    make_executable(os.path.join(wdir, sub_script_bfn))

    rdir = select_rundir(tmpdir, wdir, script)
    tmp_wrapper_script = """
shopt -s dotglob
rm -rf {rdir}
mkdir -p {rdir}
cd {rdir}
for target in {wdir}/* ; do ln -s $target . ; done
time {exe} ./{sub_script_bfn}
shopt -s dotglob
for x in ./* ; do if [ -f $x -a ! -f {wdir}/$x ] ; then mv -f $x {wdir} ; fi ; done
cd {wdir}
rm -rf {rdir}
# Presumably, the parent directories might be used by other jobs.
"""
    tmp_wrapper_script = tmp_wrapper_script.format(**locals())
    return write_script_and_wrapper_top(tmp_wrapper_script, wrapper_fn, job_done)


def filter_DBsplit_option(opt):
    """We want -a by default, but if we see --no-a[ll], we will not add -a.
    """
    flags = opt.split()
    if '-x' not in opt:
        flags.append('-x70')  # daligner belches on any read < kmer length
    return ' '.join(flags)


def update_dict_entry(d, key, func):
    d[key] = func(d[key])


def get_last_block(fn):
    """From existing db, e.g. 'raw_reads.db',
    return (bool, int)
    Related to 'openending'. Not currently used.
    """
    last_block = 1
    new_db = True
    if os.path.exists(fn):
        with open(fn) as f:
            for l in f:
                l = l.strip().split()
                if l[0] == "blocks" and l[1] == "=":
                    last_block = int(l[2])
                    new_db = False
                    break
    return (new_db, last_block)


def script_build_rdb(config, input_fofn_fn, run_jobs_bfn, length_cutoff_fn='length_cutoff'):
    """
    raw_reads.db will be output into CWD, should not already exist.
    run_jobs_bfn will be output into CWD.
    """
    last_block = 1
    config = dict(config)  # copy
    update_dict_entry(config, 'pa_DBsplit_option', filter_DBsplit_option)
    DBsplit = 'DBsplit {pa_DBsplit_option} raw_reads'.format(**config)
    # if openending == True:
    #    count = """$(cat raw_reads.db | LD_LIBRARY_PATH= awk '$1 == "blocks" {print $3-1}')"""
    count = """$(cat raw_reads.db | LD_LIBRARY_PATH= awk '$1 == "blocks" {print $3}')"""
    params = dict(config)
    length_cutoff = params.get('length_cutoff')
    if int(length_cutoff) < 0:
        bash_cutoff = '$(python2.7 -m falcon_kit.mains.calc_cutoff --coverage {} {} <(DBstats -b1 {}))'.format(
            params['seed_coverage'], params['genome_size'], 'raw_reads')
    else:
        bash_cutoff = '{}'.format(length_cutoff)
    try:
        cat_fasta = functional.choose_cat_fasta(open(input_fofn_fn).read())
    except Exception:
        LOG.exception('Using "cat" by default.')
        cat_fasta = 'cat '
    DBdust = 'DBdust {} raw_reads'.format(params.get('pa_DBdust_option', ''))
    params.update(locals())
    script = """\
echo "PBFALCON_ERRFILE=$PBFALCON_ERRFILE"
set -o pipefail
rm -f raw_reads.db .raw_reads.* # in case of re-run
#fc_fasta2fasta < {input_fofn_fn} >| fc.fofn
while read fn; do  {cat_fasta} $fn | fasta2DB -v raw_reads -i${{fn##*/}}; done < {input_fofn_fn}
#cat fc.fofn | xargs rm -f
{DBsplit}
{DBdust}
LB={count}
rm -f {run_jobs_bfn}
CUTOFF={bash_cutoff}
echo -n $CUTOFF >| {length_cutoff_fn}
echo "SMRT_PYTHON_PATH_PREPEND=$SMRT_PYTHON_PATH_PREPEND"
echo "PATH=$PATH"
which HPC.daligner
HPC.daligner -P. {pa_HPCdaligner_option} -mdust -H$CUTOFF raw_reads {last_block}-$LB >| {run_jobs_bfn}
""".format(**params)
    return script
    # Note: We dump the 'length_cutoff' file for later reference within the preassembly report
    # of pbsmrtpipe HGAP.
    # However, it is not a true dependency because we still have a workflow that starts
    # from 'corrected reads' (preads), in which case build_rdb is not run.


def script_build_pdb(config, input_fofn_bfn, run_jobs_bfn):
    # "bfn" used to mean base-file-name, but I do not remember why. ~cd
    last_block = 1
    count = """$(cat preads.db | LD_LIBRARY_PATH= awk '$1 == "blocks" {print $3}')"""
    params = dict(config)
    update_dict_entry(params, 'ovlp_DBsplit_option', filter_DBsplit_option)

    # We need to run DBdust for preads too, for consistency.
    # We will always use defaults for preads, since they should be small and clean.
    DBdust = 'DBdust preads'

    params.update(locals())
    script = """\
python -m falcon_kit.mains.copy_fofn --in={input_fofn_bfn} --out=preads.fofn --abs
rm -f preads.db .preads.* # in case of re-run
while read fn; do fasta2DB -v preads $fn; done < preads.fofn
DBsplit {ovlp_DBsplit_option} preads
{DBdust}
LB={count}
which HPC.daligner
HPC.daligner -P. {ovlp_HPCdaligner_option} -mdust -H{length_cutoff_pr} preads {last_block}-$LB >| {run_jobs_bfn}
""".format(**params)
    return script


def script_run_DB2Falcon(config, preads4falcon_fn, preads_db):
    """Run in pread_dir.
    """
    params = dict(config)
    params.update(locals())
    script = """\
# Given preads.db,
# write preads4falcon.fasta (implicitly) in CWD.
time DB2Falcon -U {preads_db}
[ -f preads4falcon.fasta ] || exit 1
""".format(**params)
    return script


def scripts_daligner(run_jobs_fn, db_prefix, rdb_build_done, nblock, pread_aln=False, skip_check=False):
    """Yield job_uid, bash
    """
    scripts = {}
    xform_script = functional.get_script_xformer(pread_aln)
    db_dir = os.path.dirname(run_jobs_fn)
    get_daligner_job_descriptions = (
        functional.get_daligner_job_descriptions_sans_LAcheck if skip_check else
        functional.get_daligner_job_descriptions)
    single = (nblock == 1)
    try:
        job_descs = get_daligner_job_descriptions(
            open(run_jobs_fn), db_prefix, single)
    except Exception:
        raise Exception('Could not parse job descriptions from file "{}":\n{}'.format(
            run_jobs_fn, traceback.format_exc()))
    for i, (desc, bash) in enumerate(viewitems(job_descs)):
        job_uid = '%04x' % i
        daligner_cmd = xform_script(bash)
        bash = """
db_dir={db_dir}
ln -sf ${{db_dir}}/.{db_prefix}.bps .
ln -sf ${{db_dir}}/.{db_prefix}.idx .
ln -sf ${{db_dir}}/{db_prefix}.db .
ln -sf ${{db_dir}}/.{db_prefix}.dust.anno .
ln -sf ${{db_dir}}/.{db_prefix}.dust.data .
rm -f *.las
{daligner_cmd}
rm -f *.C?.las *.C?.S.las *.C??.las *.C??.S.las *.C???.las *.C???.S.las
rm -f *.N?.las *.N?.S.las *.N??.las *.N??.S.las *.N???.las *.N???.S.las
""".format(db_dir=db_dir, db_prefix=db_prefix, daligner_cmd=daligner_cmd)
        yield job_uid, bash


def scripts_merge(config, db_prefix, run_jobs_fn):
    """Yield p_id, bash
    p_id can be '.123' or just ''
    """
    with open(run_jobs_fn) as f:
        mjob_data = functional.get_mjob_data(f)
    bash_funcs = """
# Get the real path of a link, following all links (since 'readlink -f' is not on OSX).
# http://stackoverflow.com/questions/7665/how-to-resolve-symbolic-links-in-a-shell-script/697552
realpath() {
    local r=$1; local t=$(readlink $r)
    while [ $t ]; do
        r=$(cd $(dirname $r) && cd $(dirname $t) && pwd -P)/$(basename $t)
        t=$(readlink $r)
    done
    echo $r
}
rmfollow() {
    local path=$(realpath $1)
    rm -f "${path}"
}
"""
    # las_fns = functional.get_las_filenames(mjob_data, db_prefix) # ok, but tricky
    for p_id in mjob_data:
        #las_fn = las_fns[p_id]
        # We already know the output .las filename by convention.
        las_fn = '%s.%s.las' % (db_prefix, p_id)

        bash_lines = mjob_data[p_id]

        script = []
        for line in bash_lines:
            script.append(line)  # Assume we no longer have "&& rm".
        for line in bash_lines:
            if not line.startswith('LAmerge'):
                continue
            las_files = [
                word + '.las' for word in functional.yield_args_from_line(line)]
            #las_fn = os.path.basename(las_files[0])
            #assert las_fn == os.path.basename(las_files[0])
            script.extend('# rmfollow {}'.format(fn) for fn in las_files[1:])
            break

        content = bash_funcs + '\n'.join(script + [''])
        yield p_id, content, las_fn


def script_run_consensus(config, db_fn, las_fn, out_file_bfn):
    """config: falcon_sense_option, length_cutoff
    TODO: Delete this version. Should not be used anywhere.
    """
    params = dict(config)
    # We calculate length_cutoff again! This is because we do not want
    # to create yet another task in pbsmrtpipe.
    length_cutoff = params.get('length_cutoff')
    if int(length_cutoff) < 0:
        bash_cutoff = '$(python2.7 -m falcon_kit.mains.calc_cutoff --coverage {} {} <(DBstats -b1 {}))'.format(
            params['seed_coverage'], params['genome_size'], db_fn)
    else:
        bash_cutoff = '{}'.format(length_cutoff)
    params.update(locals())
    LA4Falcon_flags = 'P' if params.get('LA4Falcon_preload') else ''
    if config["falcon_sense_skip_contained"]:
        LA4Falcon_flags += 'fso'
    elif config["falcon_sense_greedy"]:
        LA4Falcon_flags += 'fog'
    else:
        LA4Falcon_flags += 'fo'
    if LA4Falcon_flags:
        LA4Falcon_flags = '-' + ''.join(set(LA4Falcon_flags))
    run_consensus = "LA4Falcon -H$CUTOFF %s {db_fn} {las_fn} | python -m falcon_kit.mains.consensus {falcon_sense_option} >| {out_file_bfn}" % LA4Falcon_flags

    if config.get('dazcon', False):
        run_consensus = """
which dazcon
dazcon {pa_dazcon_option} -s {db_fn} -a {las_fn} >| {out_file_bfn}
"""

    script = """
set -o pipefail
CUTOFF=%(bash_cutoff)s
%(run_consensus)s
""" % (locals())
    return script.format(**params)


def script_run_falcon_asm(config, las_fofn_fn, preads4falcon_fasta_fn, db_file_fn):
    params = dict(config)
    params.update(locals())
    script = """\
# Given, las.fofn,
# write preads.ovl:

# mobs uses binwrappers, so it does not see our "entry-points".
# So, after dropping "src/py_scripts/*.py", we can call these via python -m:

time python -m falcon_kit.mains.ovlp_filter --db {db_file_fn} --las-fofn {las_fofn_fn} {overlap_filtering_setting} --min_len {length_cutoff_pr} --out-fn preads.ovl

ln -sf {preads4falcon_fasta_fn} ./preads4falcon.fasta

# Given preads.ovl,
# write sg_edges_list, c_path, utg_data, ctg_paths.
time python -m falcon_kit.mains.ovlp_to_graph {fc_ovlp_to_graph_option} --overlap-file preads.ovl >| fc_ovlp_to_graph.log

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
"""
    return script.format(**params)


def script_run_report_pre_assembly(i_raw_reads_db_fn, i_preads_fofn_fn, genome_length, length_cutoff, o_json_fn):
    params = dict()
    params.update(locals())
    script = """\
python2.7 -m falcon_kit.mains.report_pre_assembly --genome-length {genome_length} --length-cutoff {length_cutoff} --db {i_raw_reads_db_fn} --preads-fofn {i_preads_fofn_fn} --out {o_json_fn}
"""
    return script.format(**params)
