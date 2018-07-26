from __future__ import absolute_import


from future.utils import viewitems

from . import bash, functional
from .functional import cfg_tobool
from .io import NativeIO
from .util.system import (make_fofn_abs, make_dirs, cd)
import json
import logging
import logging.config
import os
import re
import io
import sys
import tempfile
import time
import uuid

logger = logging.getLogger(__name__)

from ConfigParser import SafeConfigParser as ConfigParser


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
    export_env_vars = ['export %(k)s=%(v)s:${%(k)s}' % dict(
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
    """NOT USED. Kept only for reference. This will be done in pypeFLOW.

    Generate script to copy db files to tmpdir (for speed).
    - Choose tmp_dir, based on src_dir name.
    - rsync basenames into tmp_dir  # after 'flock', per file
    - symlink from link_dir into tmp_dir.
    Return list of script lines, sans linefeed.
    """
    script = list()
    unique = os.path.abspath(src_dir).replace('/', '_')
    root = tempfile.gettempdir()
    tmp_dir = os.path.join(root, 'falcon', unique)
    script.append('mkdir -p %s' % tmp_dir)
    for basename in basenames:
        src = os.path.join(src_dir, basename)
        dst = os.path.join(tmp_dir, basename)
        rm_cmd = 'rm -f %s' % basename
        # Wait on lock for up to 10 minutes, in case of very large files.
        rsync_cmd = "flock -w 600 %s.lock -c 'rsync -av %s %s'" % (
            dst, src, dst)
        ln_cmd = 'ln -sf %s %s' % (dst, basename)
        script.extend([rm_cmd, rsync_cmd, ln_cmd])
    return script


def make_job_data(url, script_fn):
    """Choose defaults.
    Run in same directory as script_fn.
    Base job_name on script_fn.
    """
    wd = os.path.dirname(script_fn)
    job_name = '{0}-{1}-{2}'.format(
        os.path.basename(script_fn),
        url.split("/")[-1],
        str(uuid.uuid4())[:8],
    )
    job_data = {"job_name": job_name,
                "cwd": wd,
                "script_fn": script_fn}
    return job_data


def update_HPCdaligner_option(option):
    if '-dal' in option:
        logger.warning(
            'HPC.daligner option "-dal" has changed to "-B". Correcting this for you.')
        option = option.replace('-dal', '-B')
    if '-deg' in option:
        logger.warning(
            'HPC.daligner option "-deg" has changed to "-D". Correcting this for you.')
        option = option.replace('-deg', '-D')
    return option


def clean_falcon_options(fc):
    """Update some values in fc.
    Replace _ with - in a couple places.
    """
    keys = ('falcon_sense_option', 'overlap_filtering_setting', 'fc_ovlp_to_graph_option',
    )
    for key in keys:
        update_dash_flags(fc, key)
    for dk in ('pa_HPCdaligner_option', 'ovlp_HPCdaligner_option'):
        if dk in fc:
            fc[dk] = update_HPCdaligner_option(fc[dk])


def get_config(config):
    """
    This is only for the call from pbsmrtpipe:
       upport.get_config(support.parse_config(fn))
    We have changed parse_config() to return a dict.
    So this is a no-op.
    """
    cfg = dict(config) # already a dict now
    return cfg


def dict2config(jdict, section):
    config = ConfigParser()
    if not config.has_section(section):
        config.add_section(section)
    for (k, v) in viewitems(jdict):
        config.set(section, k, str(v))
    return config


def parse_config(config_fn):
    """Deprecated.
    Called from pbsmrtpipe, for now.
    """
    return parse_cfg_file(config_fn)

def parse_cfg_file(config_fn):
    """Return as dict.
    """
    with open(config_fn) as stream:
        ext = os.path.splitext(config_fn)[1]
        if ext in ('.json', '.js'):
            config = json.loads(stream.read())
        else:
            # Parse sections (and case-sensitively), into sub-dicts.
            config = parse_cfg_with_sections(stream)
    update_defaults(config['General'])
    # Copy General section to top, for now.
    #for key, val in config['General'].items():
    #    config[key] = val
    ##cfg.update(config.get('General', {}))
    check_config_sections(config) # Ensure that the right sections exist.
    update_job_sections(config)
    return config

def process_job_defaults(job_defaults):
    key = 'use_tmpdir'
    use_tmpdir = job_defaults.get(key, '')
    if '/' in use_tmpdir:
        tempfile.tempdir = use_tmpdir
    else:
        if use_tmpdir.lower().startswith('t'):
            use_tmpdir = tempfile.gettempdir()
        else:
            use_tmpdir = False
        job_defaults[key] = use_tmpdir

def update_job_defaults_section(config):
    """For backwards compatibility with stuff from 'General' section.
    """
    General = config['General']
    job_defaults = config['job.defaults']

    if 'njobs' in General:
        logger.warning('"njobs" belongs in the [job.defaults] section.')
    if 'pwatcher_type' in General:
        logger.warning('Please specify "pwatcher_type" only in the [job.defaults] section, not in [General].')
    if 'job_type' in General:
        logger.warning('Please specify "job_type" only in the [job.defaults] section, not in [General].')
    if 'stop_all_jobs_on_failure' in General:
        logger.warning('Please specify "stop_all_jobs_on_failure" only in the [job.defaults] section, not in [General].')
    if 'use_tmpdir' in General:
        logger.warning('Please specify "use_tmpdir" only in the [job.defaults] section, not in [General].')
    if 'job_name_style' in General:
        logger.warning('Please specify "job_name_style" only in the [job.defaults] section, not in [General].')
    if 'job_queue' in General:
        logger.warning('Please specify "JOB_QUEUE" only in the [job.defaults] section, not as "job_queue" in [General].')
    if 'sge_option' in General:
        logger.warning('Please specify "JOB_OPTS" in the [job.defaults] section, not as "sge_option" in [General].')

    pwatcher_type = General.get('pwatcher_type', 'fs_based') #, config.get('pwatcher_type')))
    job_type = job_defaults.get('job_type', General.get('job_type', '')).lower()
    job_queue = General.get('job_queue', '')
    sge_option = General.get('sge_option', '')

    if 'pwatcher_type' not in job_defaults:
        job_defaults['pwatcher_type'] = pwatcher_type
    else:
        pwatcher_type = job_defaults['pwatcher_type']
    if 'submit' not in config['job.defaults']:
        if 'blocking' == pwatcher_type:
            if not job_queue or ' ' not in job_queue:
                raise Exception('pwatcher_type=blocking, but "submit" is not in [job.defaults] section.')
            config['job.defaults']['submit'] = job_queue
            logger.warning('Please set "submit" in [job.defaults] section. (For now, we will use "job_queue" from [General], which was a hack.)')
        elif 'fs_based' == pwatcher_type or 'network_based' == pwatcher_type:
            if not job_type:
                logger.error('job.defaults.submit is not set; pwatcher_type={}; but job_type is not set. Maybe try "job_type=local" first.'.format(pwatcher_type))
                job_type = 'local'
                job_defaults['job_type'] = job_type
            allowed_job_types = ['sge', 'pbs', 'torque', 'slurm', 'lsf', 'local']
            assert job_type in allowed_job_types, 'job_type={} not in {}'.format(
                    job_type, allowed_job_types)
            if job_queue and 'JOB_QUEUE' not in config['job.defaults']:
                job_defaults['JOB_QUEUE'] = job_queue
        else:
            raise Exception('Unknown pwatcher_type={}'.format(pwatcher_type))
    #assert 'submit' in config['job.defaults'], repr(config)
    if sge_option and 'JOB_OPTS' not in config['job.defaults']:
        job_defaults['JOB_OPTS'] = sge_option
    if 'njobs' not in job_defaults:
        config['job.defaults']['njobs'] = int(General.get('default_concurrent_jobs', 8)) # GLOBAL DEFAULT CONCURRENCY
        msg = 'Please supply a default for "njobs" (aka concurrency) in section [job.defaults]. For now, we will use {}'.format(
                config['job.defaults']['njobs'])
        logger.warning(msg)
    def update_if_if(key):
        if key not in job_defaults:
            if key in General:
                job_defaults[key] = General[key]
                logger.warning('Found "{}" from [General] section; should be in [job.defaults] instead.'.format(key))
    update_if_if('job_name_style')
    update_if_if('stop_all_jobs_on_failure')
    update_if_if('use_tmpdir')

    legacy_names = [
            'pwatcher_type', 'pwatcher_directory',
            'job_type', 'job_queue', 'job_name_style',
            'use_tmpdir',
    ]
    def update_if_missing(name, sub_dict):
        if General.get(name) and name not in sub_dict:
            sub_dict[name] = General[name]
    for name in legacy_names:
        update_if_missing(name, config['job.defaults'])
    process_job_defaults(job_defaults)

def update_job_sections(config):
    """More for backwards compatibility with stuff from 'General' section.
    """
    update_job_defaults_section(config)
    General = config['General']

    # Update a few where the names change and the section is non-default.
    def update_step_job_opts(name):
        if General.get('sge_option_'+name) and 'JOB_OPTS' not in config['job.step.'+name]:
            config['job.step.'+name]['JOB_OPTS'] = General['sge_option_'+name]
    def update_step_njobs(name):
        if General.get(name+'_concurrent_jobs') and 'njobs' not in config['job.step.'+name]:
            config['job.step.'+name]['njobs'] = int(General[name+'_concurrent_jobs'])
    for name in ['da', 'la', 'pda', 'pla', 'cns', 'fc', 'asm']:
        update_step_job_opts(name)
        update_step_njobs(name)
    # Prefer 'asm' to 'fc'.
    asm = dict(config['job.step.asm'])
    config['job.step.asm'] = config['job.step.fc']
    del config['job.step.fc']
    config['job.step.asm'].update(asm)

def parse_cfg_with_sections(stream):
    """Return as dict of dict of ...
    """
    #Experimental:
    """
    ConfigParser sections become sub-sub sections when separated by dots.

        [foo.bar]
        baz = 42

    is equivalent to JSON

        {"foo": {"bar": {"baz": 42}}}
    """
    content = stream.read()
    result = dict()
    try:
        jdict = json.loads(NativeIO(content).read())
        return jdict
    except ValueError:
        pass #logger.exception('Could not parse stream as JSON.')
    try:
        config = ConfigParser() #strict=False?
        config.optionxform = str
        config.readfp(NativeIO(content))
        sections = config.sections()
        for sec in sections:
            result[sec] = dict(config.items(sec))
        return result
    except:
        raise


def check_config_sections(cfg):
    """And ensure these all exist.
    """
    allowed_sections = set(['General',
            'job.step.da', 'job.step.pda',
            'job.step.la', 'job.step.pla',
            'job.step.cns', 'job.step.fc',
            'job.step.asm',
            'job.defaults',
    ])
    all_sections = set(k for k,v in cfg.items() if isinstance(v, dict))
    unexpected = all_sections - allowed_sections
    if unexpected:
        msg = 'You have {} unexpected cfg sections: {}'.format(
            len(unexpected), unexpected)
        raise Exception(msg)
    # Guarantee they all exist.
    for sec in allowed_sections:
        if sec not in cfg:
            cfg[sec] = dict()

def update_dash_flags(cfg, key):
    if key not in cfg:
        return
    val = cfg[key]
    cfg[key] = new_val = functional.dash_flags(cfg[key])
    if val != new_val:
        msg = '''\
Option contains flags with "_":
 "{key}={val}". Those should be "-", as in
 "{key}={new_val}". Auto-replaced.'''.format(**locals())
        logger.warning(msg)


TEXT_FILE_BUSY = 'avoid_text_file_busy'

def update_defaults(cfg):
    """cfg is probably the General sub-dict.
    """
    def set_default(key, val):
        if key not in cfg:
            cfg[key] = val
    set_default('input_type', 'raw')
    set_default('overlap_filtering_setting', '--max-diff 1000 --max-cov 1000 --min-cov 2')
    set_default('pa_HPCdaligner_option', '-v -D24 -t16 -e.70 -l1000 -s100')
    set_default('ovlp_HPCdaligner_option', '-v -D24 -t32 -h60 -e.96 -l500 -s1000')
    set_default('pa_DBsplit_option', '-x500 -s200')
    set_default('skip_checks', False)
    set_default('pa_DBdust_option', '') # Gene recommends the defaults. I have tried -w128 -t2.5 -m20
    set_default('dazcon', False)
    set_default('pa_dazcon_option', '-j 4 -x -l 500')
    set_default('ovlp_DBsplit_option', '-x500 -s200')
    set_default('falcon_sense_option', '--output-multi --min-idt 0.70 --min-cov 2 --max-n-read 1800')
    set_default('falcon_sense_skip_contained', False)
    set_default('falcon_sense_greedy', False)
    set_default('LA4Falcon_preload', '')
    set_default('fc_ovlp_to_graph_option', '')
    set_default('genome_size', 0)
    set_default('seed_coverage', 20)
    set_default('length_cutoff', -1)
    set_default('length_cutoff_pr', 0)
    set_default('bestn', 12)
    set_default('target', 'assembly')
    set_default(TEXT_FILE_BUSY, bash.BUG_avoid_Text_file_busy)

    for bool_key in ('skip_checks', 'dazcon', 'falcon_sense_skip_contained', 'falcon_sense_greedy', 'LA4Falcon_preload', TEXT_FILE_BUSY):
        cfg[bool_key] = functional.cfg_tobool(cfg.get(bool_key, False))

    if 'dust' in cfg:
        logger.warning(
            "The 'dust' option is deprecated and ignored. We always run DBdust now. Use pa_DBdust_option to override its default arguments.")

    bash.BUG_avoid_Text_file_busy = cfg[TEXT_FILE_BUSY]

    clean_falcon_options(cfg)

    falcon_sense_option = cfg['falcon_sense_option']
    if 'local_match_count' in falcon_sense_option or 'output_dformat' in falcon_sense_option:
        raise Exception('Please remove obsolete "--local_match_count_*" or "--output_dformat"' +
                        ' from "falcon_sense_option" in your cfg: %s' % repr(falcon_sense_option))
    length_cutoff = int(cfg['length_cutoff'])
    if length_cutoff < 0:
        genome_size = int(cfg['genome_size'])
        if genome_size < 1:
            raise Exception(
                'Must specify either length_cutoff>0 or genome_size>0')

    # This one depends on length_cutoff_pr for its default.
    fc_ovlp_to_graph_option = cfg['fc_ovlp_to_graph_option']
    if '--min_len' not in fc_ovlp_to_graph_option and '--min-len' not in fc_ovlp_to_graph_option:
        length_cutoff_pr = cfg['length_cutoff_pr']
        fc_ovlp_to_graph_option += ' --min_len {}'.format(length_cutoff_pr)
        cfg['fc_ovlp_to_graph_option'] = fc_ovlp_to_graph_option

    target = cfg['target']
    if target not in ["overlapping", "pre-assembly", "assembly"]:
        msg = """ Target has to be "overlapping", "pre-assembly" or "assembly" in this verison. You have an unknown target {!r} in the configuration file.  """.format(target)
        raise Exception(msg)

    possible_extra_keys = [
            'sge_option', 'default_concurrent_jobs',
            'pwatcher_type', 'pwatcher_directory',
            'job_type', 'job_queue', 'job_name_style',
            'use_tmpdir',
    ]
    for step in ['da', 'la', 'pda', 'pla', 'fc', 'cns', 'asm']:
        sge_option_key = 'sge_option_' + step
        possible_extra_keys.append(sge_option_key)
        concurrent_jobs_key = step + '_concurrent_jobs'
        possible_extra_keys.append(concurrent_jobs_key)
    extra = list()
    for key in possible_extra_keys:
        if key in cfg:
            extra.append(key)
    if extra:
        extra.sort()
        msg = 'You have several old-style options. These should be provided in the `[job.defaults]` or `[job.step.*]` sections, and possibly renamed. See https://github.com/PacificBiosciences/FALCON/wiki/Configuration\n {}'.format(extra)
        logger.warning(msg)

    check_unexpected_keys(cfg)

def check_unexpected_keys(cfg):
    # Warn on unused variables.
    expected = (TEXT_FILE_BUSY,
        'input_fofn',
        'input_type',
        'overlap_filtering_setting',
        'pa_HPCdaligner_option',
        'ovlp_HPCdaligner_option',
        'pa_DBsplit_option',
        'skip_checks',
        'pa_DBdust_option',
        'dazcon',
        'pa_dazcon_option',
        'ovlp_DBsplit_option',
        'falcon_sense_option',
        'falcon_sense_skip_contained',
        'falcon_sense_greedy',
        'LA4Falcon_preload',
        'fc_ovlp_to_graph_option',
        'genome_size',
        'seed_coverage',
        'length_cutoff',
        'length_cutoff_pr',
        'bestn',
        'target',
    )
    unused = set(cfg.keys()) - set(expected)
    if unused:
        logger.warning("Unexpected keys in input config: {}".format(unused))


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
format=%(asctime)s - %(name)s:%(lineno)d - %(levelname)s - %(message)s

[formatter_form02]
format=[%(levelname)s]%(message)s
"""


def _setup_logging(logging_config_fn):
    """See https://docs.python.org/2/library/logging.config.html
    """
    logging.Formatter.converter = time.gmtime  # cannot be done in .ini

    if logging_config_fn:
        if logging_config_fn.endswith('.json'):
            logging.config.dictConfig(
                json.loads(open(logging_config_fn).read()))
            # print repr(logging.Logger.manager.loggerDict) # to debug
            return
        logger_fileobj = open(logging_config_fn)
    else:
        logger_fileobj = NativeIO(default_logging_config)
    defaults = {
    }
    logging.config.fileConfig(
        logger_fileobj, defaults=defaults, disable_existing_loggers=False)


def setup_logger(logging_config_fn):
    global logger
    try:
        _setup_logging(logging_config_fn)
        logger = logging.getLogger("fc_run")
        logger.info('Setup logging from file "{}".'.format(logging_config_fn))
    except Exception:
        logging.basicConfig()
        logger = logging.getLogger()
        logger.exception(
            'Failed to setup logging from file "{}". Using basicConfig().'.format(logging_config_fn))
    try:
        import logging_tree
        logger.info(logging_tree.format.build_description())
    except ImportError:
        pass

    return logger


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
            # We will merge in the m_* dir of the left block.
            block = int(mo.group(1))
            yield block, os.path.join(job_rundir, las_fn)


def get_length_cutoff(length_cutoff, fn):
    if length_cutoff < 0:
        length_cutoff = int(open(fn).read().strip())
        logger.info('length_cutoff=%d from %r' % (length_cutoff, fn))
    return length_cutoff  # possibly updated


def build_rdb(input_fofn_fn, config, job_done, script_fn, run_jobs_fn):
    run_jobs_fn = os.path.basename(run_jobs_fn)
    script = bash.script_build_rdb(config, input_fofn_fn, run_jobs_fn)
    bash.write_script(script, script_fn, job_done)


def build_pdb(input_fofn_fn, config, job_done, script_fn, run_jobs_fn):
    run_jobs_fn = os.path.basename(run_jobs_fn)
    script = bash.script_build_pdb(config, input_fofn_fn, run_jobs_fn)
    bash.write_script(script, script_fn, job_done)


def run_db2falcon(config, preads4falcon_fn, preads_db, job_done, script_fn):
    script = bash.script_run_DB2Falcon(config, preads4falcon_fn, preads_db)
    bash.write_script(script, script_fn, job_done)


def run_falcon_asm(config, las_fofn_fn, preads4falcon_fasta_fn, db_file_fn, job_done, script_fn):
    script = bash.script_run_falcon_asm(
        config, las_fofn_fn, preads4falcon_fasta_fn, db_file_fn)
    bash.write_script(script, script_fn, job_done)


def run_report_pre_assembly(i_raw_reads_db_fn, i_preads_fofn_fn, genome_length, length_cutoff, o_json_fn, job_done, script_fn):
    script = bash.script_run_report_pre_assembly(
        i_raw_reads_db_fn, i_preads_fofn_fn, genome_length, length_cutoff, o_json_fn)
    bash.write_script(script, script_fn, job_done)


def run_daligner(daligner_script, db_prefix, config, job_done, script_fn):
    bash.write_script(daligner_script, script_fn, job_done)


def run_las_merge(script, job_done, config, script_fn):
    bash.write_script(script, script_fn, job_done)


def run_consensus(db_fn, las_fn, out_file_fn, config, job_done, script_fn):
    script = bash.script_run_consensus(
        config, db_fn, las_fn, os.path.basename(out_file_fn))
    bash.write_script(script, script_fn, job_done)
