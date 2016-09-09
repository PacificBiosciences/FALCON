from . import bash
import ConfigParser
import json
import logging
import logging.config
import os
import re
import StringIO
import sys
import tempfile
import time
import uuid

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

def update_HPCdaligner_option(option):
    if '-dal' in option:
        logger.warning('HPC.daligner option "-dal" has changed to "-B". Correcting this for you.')
        option = option.replace('-dal', '-B')
    if '-deg' in option:
        logger.warning('HPC.daligner option "-deg" has changed to "-D". Correcting this for you.')
        option = option.replace('-deg', '-D')
    return option

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

def dict2config(jdict, section):
    config = ConfigParser.ConfigParser()
    if not config.has_section(section):
        config.add_section(section)
    for k,v in jdict.iteritems():
        config.set(section, k, str(v))
    return config

def parse_config(config_fn):
    ext = os.path.splitext(config_fn)[1]
    if ext in ('.json', '.js'):
        jdict = json.loads(open(config_fn).read())
        config = dict2config(jdict, "General")
    else:
        config = ConfigParser.ConfigParser()
        config.readfp(open(config_fn))
    return config

def get_dict_from_old_falcon_cfg(config):
    job_type = "SGE"
    section = 'General'
    if config.has_option(section, 'job_type'):
        job_type = config.get(section, 'job_type')

    watcher_type = 'fs_based'
    if config.has_option(section, 'watcher_type'):
        watcher_type = config.get(section, 'watcher_type')

    watcher_directory = 'mypwatcher'
    if config.has_option(section, 'watcher_directory'):
        watcher_directory = config.get(section, 'watcher_directory')

    pa_concurrent_jobs = 8
    if config.has_option(section, 'pa_concurrent_jobs'):
        pa_concurrent_jobs = config.getint(section, 'pa_concurrent_jobs')

    cns_concurrent_jobs = 8
    if config.has_option(section, 'cns_concurrent_jobs'):
        cns_concurrent_jobs = config.getint(section, 'cns_concurrent_jobs')

    ovlp_concurrent_jobs = 8
    if config.has_option(section, 'ovlp_concurrent_jobs'):
        ovlp_concurrent_jobs = config.getint(section, 'ovlp_concurrent_jobs')

    #appending = False
    #if config.has_option(section, 'appending'):
    #    appending = config.get(section, 'appending')
    #    if appending == "True":
    #        appending = True

    #openending = False
    #if config.has_option(section, 'openending'):
    #    openending = config.get(section, 'openending')
    #    if openending == "True":
    #        openending = True

    input_type = "raw"
    if config.has_option(section, 'input_type'):
        input_type = config.get(section, 'input_type')

    overlap_filtering_setting =  """--max_diff 1000 --max_cov 1000 --min_cov 2"""
    if config.has_option(section, 'overlap_filtering_setting'):
        overlap_filtering_setting = config.get(section, 'overlap_filtering_setting')

    pa_HPCdaligner_option = """-v -D24 -t16 -e.70 -l1000 -s100"""
    if config.has_option(section, 'pa_HPCdaligner_option'):
        pa_HPCdaligner_option = config.get(section, 'pa_HPCdaligner_option')

    ovlp_HPCdaligner_option = """ -v -D24 -t32 -h60 -e.96 -l500 -s1000"""
    if config.has_option(section, 'ovlp_HPCdaligner_option'):
        ovlp_HPCdaligner_option = config.get(section, 'ovlp_HPCdaligner_option')

    pa_HPCdaligner_option = update_HPCdaligner_option(pa_HPCdaligner_option)
    ovlp_HPCdaligner_option = update_HPCdaligner_option(ovlp_HPCdaligner_option)

    pa_DBsplit_option = """ -x500 -s200"""
    if config.has_option(section, 'pa_DBsplit_option'):
        pa_DBsplit_option = config.get(section, 'pa_DBsplit_option')

    skip_checks = True
    if config.has_option(section, 'skip_checks'):
        skip_checks = config.getboolean(section, 'skip_checks')

    dust = False
    if config.has_option(section, 'dust'):
        dust = config.getboolean(section, 'dust')

    pa_DBdust_option = "-w128 -t2.5 -m20"
    if config.has_option(section, 'pa_DBdust_option'):
        pa_DBdust_option = config.get(section, 'pa_DBdust_option')

    dazcon = False
    if config.has_option(section, 'dazcon'):
        dazcon = config.getboolean(section, 'dazcon')

    pa_dazcon_option = "-j 4 -x -l 500"
    if config.has_option(section, 'pa_dazcon_option'):
        pa_dazcon_option = config.get(section, 'pa_dazcon_option')

    ovlp_DBsplit_option = """ -x500 -s200"""
    if config.has_option(section, 'ovlp_DBsplit_option'):
        ovlp_DBsplit_option = config.get(section, 'ovlp_DBsplit_option')

    falcon_sense_option = """ --output_multi --min_idt 0.70 --min_cov 2 --max_n_read 1800 --n_core 6"""
    if config.has_option(section, 'falcon_sense_option'):
        falcon_sense_option = config.get(section, 'falcon_sense_option')
    if 'local_match_count' in falcon_sense_option or 'output_dformat' in falcon_sense_option:
        raise Exception('Please remove obsolete "--local_match_count_*" or "--output_dformat"' +
                        ' from "falcon_sense_option" in your cfg: %s' %repr(falcon_sense_option))

    falcon_sense_skip_contained = False
    if config.has_option(section, 'falcon_sense_skip_contained'):
        falcon_sense_skip_contained = config.get(section, 'falcon_sense_skip_contained')
        if falcon_sense_skip_contained in ["True", "true", "1"]:
            falcon_sense_skip_contained = True
        else:
            falcon_sense_skip_contained = False

    genome_size = 0
    if config.has_option(section, 'genome_size'):
        genome_size = config.getint(section, 'genome_size')

    seed_coverage = 20
    if config.has_option(section, 'seed_coverage'):
        seed_coverage = config.getfloat(section, 'seed_coverage')

    length_cutoff = -1
    if config.has_option(section, 'length_cutoff'):
        length_cutoff = config.getint(section, 'length_cutoff')
    if length_cutoff < 0:
        if genome_size < 1:
            raise Exception('Must specify either length_cutoff>0 or genome_size>0')

    length_cutoff_pr = config.getint(section, 'length_cutoff_pr')
    input_fofn_fn = config.get(section, 'input_fofn')

    # This one depends on length_cutoff_pr for its default.
    fc_ovlp_to_graph_option = ''
    if config.has_option(section, 'fc_ovlp_to_graph_option'):
        fc_ovlp_to_graph_option = config.get(section, 'fc_ovlp_to_graph_option')
    if '--min_len' not in fc_ovlp_to_graph_option:
        fc_ovlp_to_graph_option += ' --min_len %d' %length_cutoff_pr

    bestn = 12
    if config.has_option(section, 'bestn'):
        bestn = config.getint(section, 'bestn')

    if config.has_option(section, 'target'):
        target = config.get(section, 'target')
        if target not in ["overlapping", "pre-assembly", "assembly"]:
            msg = """ Target has to be "overlapping", "pre-assembly" or "assembly" in this verison. You have an unknown target %s in the configuration file.  """ % target
            raise Exception(msg)
    else:
        logger.info(""" No target specified, assuming "assembly" as target """)
        target = "assembly"

    if config.has_option(section, 'stop_all_jobs_on_failure'):
        stop_all_jobs_on_failure = config.getboolean(section, 'stop_all_jobs_on_failure')
    else:
        # Good default. Rarely needed, since we already stop early if *all* tasks fail
        # in a given refresh.
        stop_all_jobs_on_failure = False
    if config.has_option(section, 'use_tmpdir'):
        tmpdir = config.get(section, 'use_tmpdir')
        if '/' in tmpdir:
            tempfile.tempdir = tmpdir
            use_tmpdir = True
        else:
            use_tmpdir = config.getboolean(section, 'use_tmpdir')
    else:
        use_tmpdir = False

    TEXT_FILE_BUSY = 'avoid_text_file_busy'
    if config.has_option(section, TEXT_FILE_BUSY):
        bash.BUG_avoid_Text_file_busy = config.getboolean(section, TEXT_FILE_BUSY)

    hgap_config = {#"input_fofn_fn" : input_fofn_fn, # deprecated
                   "input_fofn" : input_fofn_fn,
                   "target" : target,
                   "job_type" : job_type,
                   "input_type": input_type,
                   #"openending": openending,
                   "pa_concurrent_jobs" : pa_concurrent_jobs,
                   "ovlp_concurrent_jobs" : ovlp_concurrent_jobs,
                   "cns_concurrent_jobs" : cns_concurrent_jobs,
                   "overlap_filtering_setting": overlap_filtering_setting,
                   "genome_size" : genome_size,
                   "seed_coverage" : seed_coverage,
                   "length_cutoff" : length_cutoff,
                   "length_cutoff_pr" : length_cutoff_pr,
                   "sge_option_da": config.get(section, 'sge_option_da'),
                   "sge_option_la": config.get(section, 'sge_option_la'),
                   "sge_option_pda": config.get(section, 'sge_option_pda'),
                   "sge_option_pla": config.get(section, 'sge_option_pla'),
                   "sge_option_fc": config.get(section, 'sge_option_fc'),
                   "sge_option_cns": config.get(section, 'sge_option_cns'),
                   "pa_HPCdaligner_option": pa_HPCdaligner_option,
                   "ovlp_HPCdaligner_option": ovlp_HPCdaligner_option,
                   "pa_DBsplit_option": pa_DBsplit_option,
                   "skip_checks": skip_checks,
                   "dust": dust,
                   "pa_DBdust_option": pa_DBdust_option,
                   "dazcon": dazcon,
                   "pa_dazcon_option": pa_dazcon_option,
                   "ovlp_DBsplit_option": ovlp_DBsplit_option,
                   "fc_ovlp_to_graph_option": fc_ovlp_to_graph_option,
                   "falcon_sense_option": falcon_sense_option,
                   "falcon_sense_skip_contained": falcon_sense_skip_contained,
                   "stop_all_jobs_on_failure": stop_all_jobs_on_failure,
                   "use_tmpdir": use_tmpdir,
                   "watcher_type": watcher_type,
                   "watcher_directory": watcher_directory,
                   TEXT_FILE_BUSY: bash.BUG_avoid_Text_file_busy,
                   }
    provided = dict(config.items(section))
    unused = set(provided) - set(k.lower() for k in hgap_config)
    if unused:
        import warnings
        warnings.warn("Unexpected keys in input config: %s" %repr(unused))

    hgap_config["install_prefix"] = sys.prefix

    return hgap_config

default_logging_config = """
[loggers]
keys=root

[handlers]
keys=stream,file_all

[formatters]
keys=form01,form02

[logger_root]
level=NOTSET
handlers=stream,file_all

[handler_stream]
class=StreamHandler
level=INFO
formatter=form02
args=(sys.stderr,)

[handler_file_all]
class=FileHandler
level=DEBUG
formatter=form01
args=('all.log', 'w')

[formatter_form01]
format=%(asctime)s - %(name)s - %(levelname)s - %(message)s

[formatter_form02]
format=[%(levelname)s]%(message)s
"""

def _setup_logging(logging_config_fn):
    """See https://docs.python.org/2/library/logging.config.html
    """
    logging.Formatter.converter = time.gmtime # cannot be done in .ini

    if logging_config_fn:
        if logging_config_fn.endswith('.json'):
            logging.config.dictConfig(json.loads(open(logging_config_fn).read()))
            #print repr(logging.Logger.manager.loggerDict) # to debug
            return
        logger_fileobj = open(logging_config_fn)
    else:
        logger_fileobj = StringIO.StringIO(default_logging_config)
    defaults = {
    }
    logging.config.fileConfig(logger_fileobj, defaults=defaults, disable_existing_loggers=False)

def setup_logger(logging_config_fn):
    global logger
    try:
        _setup_logging(logging_config_fn)
        logger = logging.getLogger("fc_run")
        logger.info('Setup logging from file "{}".'.format(logging_config_fn))
    except Exception:
        logging.basicConfig()
        logger = logging.getLogger()
        logger.exception('Failed to setup logging from file "{}". Using basicConfig().'.format(logging_config_fn))
    try:
        import logging_tree
        logger.info(logging_tree.format.build_description())
    except ImportError:
        pass

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

def get_nblock(db_file):
    """Return #blocks in dazzler-db.
    """
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

def daligner_gather_las(job_rundirs):
    """Return list of (block, las_fn).
    """
    # Could be L1.* or preads.*
    re_las = re.compile(r'\.(\d*)(\.\d*)?\.las$')
    for job_rundir in job_rundirs:
        # 'out' sub-dir by convention. See run_daligner() above. (Improve that someday.)
        for las_fn in os.listdir(job_rundir):
            mo = re_las.search(las_fn)
            if not mo:
                continue
            block = int(mo.group(1)) # We will merge in the m_* dir of the left block.
            yield block, os.path.join(job_rundir, las_fn)

def get_length_cutoff(length_cutoff, fn):
    if length_cutoff < 0:
        try:
            length_cutoff = int(open(fn).read().strip())
            logger.info('length_cutoff=%d from %r' %(length_cutoff, fn))
        except Exception:
            logger.exception('Unable to read length_cutoff from "%s".' %fn)
    return length_cutoff # possibly updated

def build_rdb(input_fofn_fn, config, job_done, script_fn, run_jobs_fn):
    run_jobs_fn = os.path.basename(run_jobs_fn)
    script = bash.script_build_rdb(config, input_fofn_fn, run_jobs_fn)
    bash.get_write_script_and_wrapper(config)(script, script_fn, job_done)

def build_pdb(input_fofn_fn, config, job_done, script_fn, run_jobs_fn):
    run_jobs_fn = os.path.basename(run_jobs_fn)
    script = bash.script_build_pdb(config, input_fofn_fn, run_jobs_fn)
    bash.get_write_script_and_wrapper(config)(script, script_fn, job_done)

def run_db2falcon(config, job_done, script_fn):
    script = bash.script_run_DB2Falcon(config)
    bash.get_write_script_and_wrapper(config)(script, script_fn, job_done)

def run_falcon_asm(config, las_fofn_fn, preads4falcon_fasta_fn, db_file_fn, job_done, script_fn):
    script = bash.script_run_falcon_asm(config, las_fofn_fn, preads4falcon_fasta_fn, db_file_fn)
    bash.get_write_script_and_wrapper(config)(script, script_fn, job_done)

def run_daligner(daligner_script, db_prefix, config, job_done, script_fn):
    bash.get_write_script_and_wrapper(config)(daligner_script, script_fn, job_done)

def run_las_merge(script, job_done, config, script_fn):
    bash.get_write_script_and_wrapper(config)(script, script_fn, job_done)

def run_consensus(db_fn, las_fn, out_file_fn, config, job_done, script_fn):
    script = bash.script_run_consensus(config, db_fn, las_fn, os.path.basename(out_file_fn))
    bash.get_write_script_and_wrapper(config)(script, script_fn, job_done)
