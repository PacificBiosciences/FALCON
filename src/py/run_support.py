import ConfigParser
import logging
import logging.config
import os
import StringIO
import sys
import tempfile
import time
import uuid

job_type = None
logger = None

def _prepend_env_paths(content, names):
    """
    E.g.
      names = ['PATH', 'PYTYHONPATH']
      content =
        echo hi
      =>
        export PATH=current:path:${PATH}
        export PYTHON=current:path:${PYTHONPATH}
        echo hi
    """
    export_env_vars = ['export %(k)s=%(v)s:${%(k)s}' %dict(
        k=name, v=os.environ.get(name, '')) for name in names]
    return '\n'.join(export_env_vars + [content])

def update_env_in_script(fn, names):
    """Modify fn using on prepend_env_paths().
    """
    with open(fn) as ifs:
        content = ifs.read()
    content = _prepend_env_paths(content, names)
    with open(fn, 'w') as ofs:
        ofs.write(content)

def use_tmpdir_for_files(basenames, src_dir, link_dir):
    """Generate script to copy db files to tmpdir (for speed).
    - Choose tmp_dir, based on src_dir name.
    - rsync basenames into tmp_dir  # after 'flock', per file
    - symlink from link_dir into tmp_dir.
    Return list of script lines, sans linefeed.
    """
    script = list()
    unique = os.path.abspath(src_dir).replace('/', '_')
    root = tempfile.gettempdir()
    tmp_dir = os.path.join(root, 'falcon', unique)
    script.append('mkdir -p %s' %tmp_dir)
    for basename in basenames:
        src = os.path.join(src_dir, basename)
        dst = os.path.join(tmp_dir, basename)
        rm_cmd = 'rm -f %s' %basename
        # Wait on lock for up to 10 minutes, in case of very large files.
        rsync_cmd = "flock -w 600 %s.lock -c 'rsync -av %s %s'" %(dst, src, dst)
        ln_cmd = 'ln -sf %s %s' %(dst, basename)
        script.extend([rm_cmd, rsync_cmd, ln_cmd])
    return script

def make_job_data(url, script_fn):
    """Choose defaults.
    Run in same directory as script_fn.
    Base job_name on script_fn.
    """
    wd = os.path.dirname(script_fn)
    job_name = '{0}-{1}-{1}'.format(
            os.path.basename(script_fn),
            url.split("/")[-1],
            str(uuid.uuid4())[:8],
            )
    job_data = {"job_name": job_name,
                "cwd": wd,
                "script_fn": script_fn }
    return job_data

def validate_config_dict(cd):
    pass

def get_config(config):
    """Temporary version for pbsmrtpipe.
    This will add missing (but curently required) options and use
    get_dict_from_old_falcon_cfg() below.
    The plan is to pass a simpler config from pbsmrtpipe,
    but that will be in a different commit.
    Side-effect: Update 'config'.
    """
    section = 'General'
    def add(name, val):
        if not config.has_option(section, name):
            config.set(section, name, val)
    add('input_fofn', 'NA')
    add('target', 'assembly')
    add('sge_option_da', 'NA')
    add('sge_option_la', 'NA')
    add('sge_option_pda', 'NA')
    add('sge_option_pla', 'NA')
    add('sge_option_fc', 'NA')
    add('sge_option_cns', 'NA')
    return get_dict_from_old_falcon_cfg(config)

def parse_config(config_fn):
    config = ConfigParser.ConfigParser()
    config.read(config_fn)
    return config

def get_dict_from_old_falcon_cfg(config):
    global job_type  # TODO: Stop using global for wait_for_file().
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
            msg = """ Target has to be "overlapping", "pre-assembly" or "assembly" in this verison. You have an unknown target %s in the configuration file.  """ % target
            raise Exception(msg)
    else:
        logger.info(""" No target specified, assuming "assembly" as target """)
        target = "assembly"

    if config.has_option('General', 'use_tmpdir'):
        use_tmpdir = config.getboolean('General','use_tmpdir')
    else:
        use_tmpdir = False

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
                   "falcon_sense_skip_contained": falcon_sense_skip_contained,
                   "use_tmpdir": use_tmpdir,
                   }

    hgap_config["install_prefix"] = sys.prefix

    return hgap_config

default_logging_config = """
[loggers]
keys=root,pypeflow,fc_run

[handlers]
keys=stream,file_pypeflow,file_fc

[formatters]
keys=form01,form02

[logger_root]
level=NOTSET
handlers=stream

[logger_pypeflow]
level=DEBUG
handlers=file_pypeflow
qualname=pypeflow
propagate=1

[logger_fc_run]
level=NOTSET
handlers=file_fc
qualname=.
propagate=1

[handler_stream]
class=StreamHandler
level=INFO
formatter=form02
args=(sys.stderr,)

[handler_file_pypeflow]
class=FileHandler
level=DEBUG
formatter=form01
args=('pypeflow.log',)

[handler_file_fc]
class=FileHandler
level=DEBUG
formatter=form01
args=('fc.log',)

[formatter_form01]
format=%(asctime)s - %(name)s - %(levelname)s - %(message)s

[formatter_form02]
format=[%(levelname)s]%(message)s
"""

def setup_logger(logging_config_fn):
    """See https://docs.python.org/2/library/logging.config.html
    """
    logging.Formatter.converter = time.gmtime # cannot be done in .ini

    if logging_config_fn:
        logger_fileobj = open(logging_config_fn)
    else:
        logger_fileobj = StringIO.StringIO(default_logging_config)
    defaults = {
    }
    logging.config.fileConfig(logger_fileobj, defaults=defaults, disable_existing_loggers=False)

    global logger
    logger = logging.getLogger("fc_run")
    return logger

def make_fofn_abs(i_fofn_fn, o_fofn_fn):
    """Copy i_fofn to o_fofn, but with relative filenames expanded for CWD.
    """
    assert os.path.abspath(o_fofn_fn) != os.path.abspath(i_fofn_fn)
    with open(i_fofn_fn) as ifs, open(o_fofn_fn, 'w') as ofs:
        for line in ifs:
            ifn = line.strip()
            if not ifn: continue
            abs_ifn = os.path.abspath(ifn)
            ofs.write('%s\n' %abs_ifn)
    #return o_fofn_fn

def make_dirs(d):
    if not os.path.isdir(d):
        os.makedirs(d)

def build_rdb(input_fofn_fn, work_dir, config, job_done, script_fn, run_jobs_fn):
    length_cutoff = config["length_cutoff"]
    pa_HPCdaligner_option = config["pa_HPCdaligner_option"]
    pa_DBsplit_option = config["pa_DBsplit_option"]
    openending = config["openending"]

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
        script_file.write("set -vex\n")
        script_file.write("trap 'touch {job_done}.exit' EXIT\n".format(job_done = job_done))
        script_file.write("cd {work_dir}\n".format(work_dir = work_dir))
        script_file.write("hostname\n")
        script_file.write("date\n")
        #script_file.write("for f in `cat {input_fofn_fn}`; do fasta2DB raw_reads $f; done\n".format(input_fofn_fn = input_fofn_fn))
        script_file.write("fasta2DB -v raw_reads -f{input_fofn_fn}\n".format(input_fofn_fn = input_fofn_fn))
        if new_db  == True:
            script_file.write("DBsplit %s raw_reads\n" % pa_DBsplit_option)
            script_file.write("date\n")
        if openending == True:
            script_file.write("""LB=$(cat raw_reads.db | awk '$1 == "blocks" {print $3-1}')\n""")
        else:
            script_file.write("""LB=$(cat raw_reads.db | awk '$1 == "blocks" {print $3}')\n""")
        script_file.write("HPCdaligner %s -H%d raw_reads %d-$LB > %s\n" %(
            pa_HPCdaligner_option, length_cutoff, last_block, run_jobs_fn))
        script_file.write("date\n")
        script_file.write("touch {job_done}\n".format(job_done = job_done))

def build_pdb(input_fofn_fn, work_dir, config, job_done, script_fn, run_jobs_fn):
    length_cutoff = config["length_cutoff_pr"]
    ovlp_HPCdaligner_option = config["ovlp_HPCdaligner_option"]
    ovlp_DBsplit_option = config["ovlp_DBsplit_option"]

    with open(script_fn,"w") as script_file:
        script_file.write("set -vex\n")
        script_file.write("trap 'touch {job_done}.exit' EXIT\n".format(job_done = job_done))
        script_file.write("cd {work_dir}\n".format(work_dir = work_dir))
        script_file.write("hostname\n")
        script_file.write("date\n")
        script_file.write("fasta2DB -v preads -f{input_fofn_fn}\n".format(input_fofn_fn = input_fofn_fn))
        script_file.write("DBsplit -x%d %s preads\n" % (length_cutoff, ovlp_DBsplit_option))
        script_file.write("HPCdaligner %s -H%d preads > %s\n" %(
            ovlp_HPCdaligner_option, length_cutoff, run_jobs_fn))
        script_file.write("touch {job_done}\n".format(job_done = job_done))

def run_falcon_asm(pread_dir, db_file, config, job_done, script_fn):
    wd = os.path.dirname(script_fn)
    overlap_filtering_setting = config["overlap_filtering_setting"]
    length_cutoff_pr = config["length_cutoff_pr"]

    script = []
    script.append( "set -vex" )
    script.append( "trap 'touch %s.exit' EXIT" % job_done )
    script.append( "cd %s" % pread_dir )
    script.append("date")
    # Given preads.db,
    # write preads4falcon.fasta, in 1-preads_ovl:
    script.append( "DB2Falcon -U preads")
    script.append("date")
    script.append( "cd %s" % wd )
    # Generate las.fofn:
    script.append( """find %s/las_files -name "*.las" > las.fofn """ % pread_dir )
    # Given, las.fofn,
    # write preads.ovl:
    script.append( """fc_ovlp_filter --db %s --fofn las.fofn %s --min_len %d > preads.ovl""" %\
            (db_file, overlap_filtering_setting, length_cutoff_pr) )
    script.append("date")
    script.append( "ln -sf %s/preads4falcon.fasta ." % pread_dir)
    # TODO: Figure out which steps need preads4falcon.fasta.

    # Given preads.ovl,
    # write sg_edges_list, c_path, utg_data, ctg_paths.
    script.append( """fc_ovlp_to_graph preads.ovl --min_len %d > fc_ovlp_to_graph.log""" % length_cutoff_pr) # TODO: drop this logfile
    script.append("date")
    # Given sg_edges_list, utg_data, ctg_paths,
    # Write p_ctg.fa and a_ctg_all.fa,
    # plus a_ctg_base.fa, p_ctg_tiling_path, a_ctg_tiling_path, a_ctg_base_tiling_path:
    script.append( """fc_graph_to_contig""" )
    script.append("date")
    # Given a_ctg_all.fa, write a_ctg.fa:
    script.append( """fc_dedup_a_tigs""" )
    script.append("date")
    script.append( """touch %s""" % job_done)

    with open(script_fn, "w") as script_file:
        script_file.write("\n".join(script) + '\n')

def run_daligner(daligner_cmd, db_prefix, nblock, config, job_done, script_fn):
    cwd = os.path.dirname(script_fn)

    script = []
    script.append( "set -vex" )
    script.append( "trap 'touch {job_done}.exit' EXIT".format(job_done = job_done) )
    script.append( "cd %s" % cwd )
    script.append( "hostname" )
    script.append( "date" )
    if config['use_tmpdir']:
        basenames = [pattern.format(db_prefix) for pattern in ('.{}.idx', '.{}.bps', '{}.db')]
        dst_dir = os.path.abspath(cwd)
        src_dir = os.path.abspath(os.path.dirname(cwd)) # by convention
        script.extend(use_tmpdir_for_files(basenames, src_dir, dst_dir))
    script.append( "time "+ daligner_cmd )
    script.append( "date" )
    script.append( "touch {job_done}".format(job_done = job_done) )

    with open(script_fn,"w") as script_file:
        script_file.write("\n".join(script) + '\n')

def run_las_merge(p_script_fn, job_done, config, script_fn):
    cwd = os.path.dirname(script_fn)
    script = []
    script.append( "set -vex" )
    script.append( "trap 'touch {job_done}.exit' EXIT".format(job_done = job_done) )
    script.append( "cd %s" % cwd )
    script.append( "hostname" )
    script.append( "time bash %s" % p_script_fn )
    script.append( "touch {job_done}".format(job_done = job_done) )

    with open(script_fn,"w") as script_file:
        script_file.write("date\n")
        script_file.write("\n".join(script) + '\n')
        script_file.write("date\n")

def run_consensus(job_id, out_file_fn, prefix, config, job_done, script_fn):
    cwd = os.path.dirname(script_fn)
    falcon_sense_option = config["falcon_sense_option"]
    length_cutoff = config["length_cutoff"]

    script = []
    script.append("set -vex")
    script.append("set -o pipefail")
    script.append("trap 'touch {job_done}.exit' EXIT".format(job_done = job_done))
    script.append("cd ..")
    script.append( "date" )
    pipe = ''
    if config["falcon_sense_skip_contained"]:
        pipe += """LA4Falcon -H%d -fso %s las_files/%s.%d.las | """ % (length_cutoff, prefix, prefix, job_id)
    else:
        pipe += """LA4Falcon -H%d -fo %s las_files/%s.%d.las | """ % (length_cutoff, prefix, prefix, job_id)
    pipe += """fc_consensus %s > %s""" % (falcon_sense_option, out_file_fn)
    script.append(pipe)
    script.append("date")
    script.append("touch {job_done}".format(job_done = job_done))

    c_script_fn = os.path.join(cwd, "cp_%05d.sh" % job_id)
    with open(c_script_fn, "w") as f:
        f.write('\n'.join(script + ['']))

    script = []
    script.append( "set -vex" )
    script.append( "cd %s" % cwd )
    script.append( "hostname" )
    script.append( "date" )
    script.append( "time bash %s" %os.path.basename(c_script_fn) )
    script.append( "date" )

    with open(script_fn,"w") as f:
        f.write("\n".join(script + ['']))
