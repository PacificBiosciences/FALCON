'''
    # DAMASKER options
    """
    Example config usage:
    pa_use_tanmask = true
    pa_use_repmask = true
    pa_HPCtanmask_option =
    pa_repmask_levels = 2
    pa_HPCrepmask_1_option = -g1 -c20 -mtan
    pa_HPCrepmask_2_option = -g10 -c15 -mtan -mrep1
    pa_damasker_HPCdaligner_option = -mtan -mrep1 -mrep10
    """
    pa_use_tanmask = False
    if config.has_option(section, 'pa_use_tanmask'):
        pa_use_tanmask = config.getboolean(section, 'pa_use_tanmask')

    pa_HPCtanmask_option = ""
    if config.has_option(section, 'pa_HPCtanmask_option'):
        pa_HPCtanmask_option = config.get(section, 'pa_HPCtanmask_option')

    pa_use_repmask = False
    if config.has_option(section, 'pa_use_repmask'):
        pa_use_repmask = config.getboolean(section, 'pa_use_repmask')

    pa_repmask_levels = 0   # REPmask tool can be used multiple times.
    if config.has_option(section, 'pa_repmask_levels'):
        pa_repmask_levels = config.getint(section, 'pa_repmask_levels')

    pa_HPCrepmask_1_option = """ -g1 -c20 -mtan"""
    if config.has_option(section, 'pa_HPCrepmask_1_option'):
        pa_HPCrepmask_1_option = config.get(section, 'pa_HPCrepmask_1_option')

    pa_HPCrepmask_2_option = """ -g10 -c15 -mtan -mrep1"""
    if config.has_option(section, 'pa_HPCrepmask_2_option'):
        pa_HPCrepmask_2_option = config.get(section, 'pa_HPCrepmask_2_option')

    pa_damasker_HPCdaligner_option = """ -mtan -mrep1 -mrep10"""    # Repeat masks need to be passed to Daligner.
    if config.has_option(section, 'pa_damasker_HPCdaligner_option'):
        pa_damasker_HPCdaligner_option = config.get(
            section, 'pa_damasker_HPCdaligner_option')
    # End of DAMASKER options.

    # Note: We dump the 'length_cutoff' file for later reference within the preassembly report
    # of pbsmrtpipe HGAP.
    # However, it is not a true dependency because we still have a workflow that starts
    # from 'corrected reads' (preads), in which case build_db is not run.

    Catrack -dfv raw_reads.db tan
'''
from __future__ import absolute_import

import argparse
import collections
import glob
import logging
import os
import re
import sys
from ..util.io import yield_validated_fns
from .. import io, functional
from .. import(
        bash,  # for write_sub_script
        pype_tasks,  # for TASKS
)

LOG = logging.getLogger()
WAIT = 20 # seconds to wait for file to exist

def filter_DBsplit_option(opt):
    """We want -a by default, but if we see --no-a[ll], we will not add -a.
    """
    flags = opt.split()
    if '-x' not in opt:
        flags.append('-x70')  # daligner belches on any read < kmer length
    return ' '.join(flags)
def script_build_db(config, input_fofn_fn, db):
    """
    db (e.g. 'raw_reads.db') will be output into CWD, should not already exist.
    'dust' track will also be generated.
    """
    params = dict(config)
    try:
        cat_fasta = functional.choose_cat_fasta(open(input_fofn_fn).read())
    except Exception:
        LOG.exception('Using "cat" by default.')
        cat_fasta = 'cat '
    DBdust = 'DBdust {} {}'.format(config.get('pa_DBdust_option', ''), db)
    params.update(locals())
    script = """\
echo "PBFALCON_ERRFILE=$PBFALCON_ERRFILE"
set -o pipefail
rm -f {db}.db .{db}.* # in case of re-run
#fc_fasta2fasta < {input_fofn_fn} >| fc.fofn
while read fn; do  {cat_fasta} $fn | fasta2DB -v {db} -i${{fn##*/}}; done < {input_fofn_fn}
#cat fc.fofn | xargs rm -f
{DBdust}
""".format(**params)
    return script
def script_length_cutoff(config, db, length_cutoff_fn='length_cutoff'):
    params = dict(config)
    length_cutoff = config['user_length_cutoff']
    if int(length_cutoff) < 0:
        bash_cutoff = '$(python2.7 -m falcon_kit.mains.calc_cutoff --coverage {} {} <(DBstats -b1 {}))'.format(
            params['seed_coverage'], params['genome_size'], db)
    else:
        bash_cutoff = '{}'.format(length_cutoff)
    params.update(locals())
    return """
CUTOFF={bash_cutoff}
echo -n $CUTOFF >| {length_cutoff_fn}
""".format(**params)
def script_DBsplit(config, db):
    params = dict(config)
    params.update(locals())
    DBsplit_opt = filter_DBsplit_option(config['DBsplit_opt'])
    params.update(locals())
    return """
DBsplit -f {DBsplit_opt} {db}
#LB=$(cat {db}.db | LD_LIBRARY_PATH= awk '$1 == "blocks" {{print $3}}')
#echo -n $LB >| db_block_count
""".format(**params)
def build_db(config, input_fofn_fn, db_fn, length_cutoff_fn):
    LOG.info('Building rdb from {!r}, to write {!r}'.format(
        input_fofn_fn, db_fn))
    db = os.path.splitext(db_fn)[0]

    # First, make FOFN relative to thisdir.
    my_input_fofn_fn = 'my.' + os.path.basename(input_fofn_fn)
    with open(my_input_fofn_fn, 'w') as stream:
        for fn in yield_validated_fns(input_fofn_fn):
            stream.write(fn)
            stream.write('\n')
    script = ''.join([
        script_build_db(config, my_input_fofn_fn, db),
        script_DBsplit(config, db),
        script_length_cutoff(config, db, length_cutoff_fn),
    ])
    script_fn = 'build_db.sh'
    with open(script_fn, 'w') as ofs:
        exe = bash.write_sub_script(ofs, script)
    io.syscall('bash -vex {}'.format(script_fn))

def script_HPC_daligner(config, db, length_cutoff_fn, tracks, prefix='daligner-jobs'):
    params = dict(config)
    masks = ' '.join('-m{}'.format(track) for track in tracks)
    params.update(locals())
    symlink(length_cutoff_fn, 'CUTOFF')
    return """
#LB=$(cat db_block_count)
CUTOFF=$(cat CUTOFF)
rm -f daligner-jobs.*
echo "SMRT_PYTHON_PATH_PREPEND=$SMRT_PYTHON_PATH_PREPEND"
echo "PATH=$PATH"
which HPC.daligner
HPC.daligner -P. {daligner_opt} {masks} -H$CUTOFF -f{prefix} {db} >| run_jobs.sh
    """.format(**params)
def script_HPC_TANmask(config, db, prefix):
    assert prefix and '/' not in prefix
    params = dict(config)
    params.update(locals())
    return """
rm -f {prefix}.*
rm -f .{db}.*.tan.anno
rm -f .{db}.*.tan.data
echo "SMRT_PYTHON_PATH_PREPEND=$SMRT_PYTHON_PATH_PREPEND"
echo "PATH=$PATH"
which HPC.TANmask
HPC.TANmask -P. {TANmask_opt} -v -f{prefix} {db}
    """.format(**params)

def symlink(actual, symbolic=None, force=True):
    """Symlink into cwd, relatively.
    symbolic name is basename(actual) if not provided.
    If not force, raise when already exists and does not match.
    But ignore symlink to self.
    """
    symbolic = os.path.basename(actual) if not symbolic else symbolic
    if os.path.abspath(actual) == os.path.abspath(symbolic):
        LOG.warning('Cannot symlink {!r} as {!r}, itself.'.format(actual, symbolic))
        return
    rel = os.path.relpath(actual)
    if force:
        LOG.info('ln -sf {} {}'.format(rel, symbolic))
        if os.path.lexists(symbolic):
            if os.readlink(symbolic) == rel:
                return
            else:
                os.unlink(symbolic)
    else:
        LOG.info('ln -s {} {}'.format(rel, symbolic))
        if os.path.lexists(symbolic):
            if os.readlink(symbolic) != rel:
                msg = '{!r} already exists as {!r}, not {!r}'.format(
                        symbolic, os.readlink(symbolic), rel)
                raise Exception(msg)
            else:
                LOG.info('{!r} already points to {!r}'.format(symbolic, rel))
                return
    os.symlink(rel, symbolic)

def symlink_db(db_fn):
    """Symlink everything that could be related to this Dazzler DB.
    """
    db_dirname, db_basename = os.path.split(db_fn)
    dbname = os.path.splitext(db_basename)[0]

    fn = os.path.join(db_dirname, dbname + '.db')
    symlink(fn)
    possible_suffixes = (
            '.idx', '.bps',
            '.dust.data', '.dust.anno',
            '.tan.data', '.tan.anno',
    )
    for suffix in possible_suffixes:
        fn = os.path.join(db_dirname, '.' + dbname + suffix)
        if os.path.exists(fn):
            symlink(fn)
    return dbname

def tan_split(config, config_fn, db_fn, uows_fn, bash_template_fn):
    with open(bash_template_fn, 'w') as stream:
        stream.write(pype_tasks.TASK_DB_TAN_APPLY_SCRIPT)
    # TANmask would put track-files in the DB-directory, not '.',
    # so we need to symlink everything first.
    db = symlink_db(db_fn)

    script = ''.join([
        script_HPC_TANmask(config, db, prefix='tan-jobs'),
    ])
    script_fn = 'split_db.sh'
    with open(script_fn, 'w') as ofs:
        exe = bash.write_sub_script(ofs, script)
    io.syscall('bash -vex {}'.format(script_fn))

    # We now have files like tan-jobs.01.OVL
    # We need to parse that one. (We ignore the others.)
    lines = open('tan-jobs.01.OVL').readlines()

    re_block = re.compile(r'{}(\.\d+|)'.format(db))

    def get_blocks(line):
        """Return ['.1', '.2', ...]
        """
        return [mo.group(1) for mo in re_block.finditer(line)]

    scripts = list()
    for line in lines:
        if line.startswith('#'):
            continue
        if not line.strip():
            continue
        blocks = get_blocks(line)
        assert blocks, 'No blocks found in {!r} from {!r}'.format(line, 'tan-jobs.01.OVL')
        las_files = ' '.join('TAN.{db}{block}.las'.format(db=db, block=block) for block in blocks)
        script_lines = [
            line,
            'LAcheck {} {}\n'.format(db, las_files),
            'TANmask {} {}\n'.format(db, las_files),
            'rm -f {}\n'.format(las_files),
        ]
        if len(blocks) == 1:
            # special case
            assert blocks == [''], "{!r} != ['']".format(blocks)
            script_lines.append('mv .{db}.tan.data .{db}.1.tan.data\n'.format(db=db))
            script_lines.append('mv .{db}.tan.anno .{db}.1.tan.anno\n'.format(db=db))
        scripts.append(''.join(script_lines))

    jobs = list()
    for i, script in enumerate(scripts):
        job_id = 'tan_{:03d}'.format(i)
        script_dir = os.path.join('tan-scripts', job_id)
        script_fn = os.path.join(script_dir, 'run_datander.sh')
        io.mkdirs(script_dir)
        with open(script_fn, 'w') as stream:
            stream.write('{}\n'.format(script))
        # Record in a job-dict.
        job = dict()
        job['input'] = dict(
                config=config_fn,
                db=db_fn,
                script=os.path.abspath(script_fn),
        )
        job['output'] = dict(
                job_done = 'job.done'
        )
        job['params'] = dict(
        )
        job['wildcards'] = dict(
                tan0_id=job_id,
        )
        jobs.append(job)
    io.serialize(uows_fn, jobs)

def tan_apply(db_fn, script_fn, job_done_fn):
    # datander would put track-files in the DB-directory, not '.',
    # so we need to symlink everything first.
    db = symlink_db(db_fn)

    symlink(script_fn)
    io.syscall('bash -vex {}'.format(os.path.basename(script_fn)))
    io.touch(job_done_fn)

def tan_combine(db_fn, gathered_fn, new_db_fn):
    new_db = os.path.splitext(new_db_fn)[0]
    db = symlink_db(db_fn)
    assert db == new_db, 'basename({!r})!={!r}, but old and new DB names must match.'.format(db_fn, new_db_fn)

    # Remove old, in case of resume.
    io.syscall('rm -f .{db}.*.tan.anno .{db}.*.tan.data'.format(db=db))

    gathered = io.deserialize(gathered_fn)
    gathered_dn = os.path.dirname(gathered_fn)

    # Create symlinks for all track-files.
    for job in gathered:
        done_fn = job['job_done']
        done_dn = os.path.dirname(done_fn)
        if not os.path.isabs(done_dn):
            LOG.info('Found relative done-file: {!r}'.format(done_fn))
            done_dn = os.path.join(gathered_dn, done_dn)
        annos = glob.glob('{}/.{}.*.tan.anno'.format(done_dn, db))
        datas = glob.glob('{}/.{}.*.tan.data'.format(done_dn, db))
        assert len(annos) == len(datas), 'Mismatched globs:\n{!r}\n{!r}'.format(annos, datas)
        for fn in annos + datas:
            symlink(fn, force=False)
    cmd = 'Catrack -vdf {} tan'.format(db)
    io.syscall(cmd)

def split_db(config_fn, daligner_uows_fn, db_fn):
    raise NotImplementedError('Who dey?')
def get_tracks(db_fn):
    db_dirname, db_basename = os.path.split(db_fn)
    dbname = os.path.splitext(db_basename)[0]
    fns = glob.glob('{}/.{}.*.anno'.format(db_dirname, dbname))
    re_anno = re.compile(r'\.{}\.([^\.]+)\.anno'.format(dbname))
    tracks = [re_anno.search(fn).group(1) for fn in fns]
    return tracks

def daligner_split(config, config_fn, db_fn, nproc, wildcards, length_cutoff_fn, split_fn, bash_template_fn):
    with open(bash_template_fn, 'w') as stream:
        stream.write(pype_tasks.TASK_DB_DALIGNER_APPLY_SCRIPT)
    db = os.path.splitext(db_fn)[0]
    dbname = os.path.basename(db)

    tracks = get_tracks(db_fn)

    script = ''.join([
        script_HPC_daligner(config, db, length_cutoff_fn, tracks, 'daligner-jobs'),
    ])
    script_fn = 'split_db.sh'
    with open(script_fn, 'w') as ofs:
        exe = bash.write_sub_script(ofs, script)
    io.syscall('bash -vex {}'.format(script_fn))

    # We now have files like tan-jobs.01.OVL
    # We need to parse that one. (We ignore the others.)
    lines = open('daligner-jobs.01.OVL').readlines()

    preads_aln = True if dbname == 'preads' else False
    xformer = functional.get_script_xformer(preads_aln)
    LOG.debug('preads_aln={!r} (True => use daligner_p)'.format(preads_aln))

    scripts = list()
    for line in lines:
        if line.startswith('#'):
            continue
        if not line.strip():
            continue
        line = xformer(line) # Use daligner_p for preads.
        scripts.append(line)
    """
    Special case:
        # Daligner jobs (1)
        daligner raw_reads raw_reads && mv raw_reads.raw_reads.las raw_reads.las
    In that case, the "block" name is empty. (See functional.py)
    We will rename the file. (LAmerge on a single input is a no-op, which is fine.)
    """
    if len(scripts) == 1:
        script = scripts[0]
        re_script = re.compile(r'(mv\b.*\S+\s+)(\S+)$') # no trailing newline, for now
        mo = re_script.search(script)
        if not mo:
            msg = 'Only 1 line in daligner-jobs.01.OVL, but\n {!r} did not match\n {!r}.'.format(
                re_script.pattern, script)
            LOG.warning(msg)
        else:
            new_script = re_script.sub(r'\1{dbname}.1.{dbname}.1.las'.format(dbname=dbname), script, 1)
            msg = 'Only 1 line in daligner-jobs.01.OVL:\n {!r} matches\n {!r}. Replacing with\n {!r}.'.format(
                re_script.pattern, script, new_script)
            LOG.warning(msg)
            scripts = [new_script]

    for i, script in enumerate(scripts):
        LAcheck = 'LAcheck -vS {} *.las'.format(db)
        script += '\n' + LAcheck + '\n'
        scripts[i] = script

    jobs = list()
    for i, script in enumerate(scripts):
        job_id = 'j_{:04d}'.format(i)
        script_dir = os.path.join('daligner-scripts', job_id)
        script_fn = os.path.join(script_dir, 'run_daligner.sh')
        io.mkdirs(script_dir)
        with open(script_fn, 'w') as stream:
            stream.write(script)
        # Record in a job-dict.
        job = dict()
        job['input'] = dict(
                config=config_fn,
                db=db_fn,
                script=os.path.abspath(script_fn),
        )
        job['output'] = dict(
                job_done = 'daligner.done'
        )
        job['params'] = dict(
        )
        job['wildcards'] = {wildcards: job_id}
        jobs.append(job)
    io.serialize(split_fn, jobs)

def daligner_apply(db_fn, script_fn, job_done_fn):
    symlink(script_fn)
    io.syscall('bash -vex {}'.format(os.path.basename(script_fn)))
    io.touch(job_done_fn)

class MissingLas(Exception):
    pass

def is_perfect_square(n):
    import math
    root = round(math.sqrt(n))
    return n == root*root

def daligner_combine(db_fn, gathered_fn, las_paths_fn):
    """Merge all .las pair-files from gathered daligner runs.
    Write simple las_paths_fn and archaic p_id2las_fn.
    """
    gathered = io.deserialize(gathered_fn)
    d = os.path.abspath(os.path.realpath(os.path.dirname(gathered_fn)))
    def abspath(fn):
        if os.path.isabs(fn):
            return fn # I expect this never to happen though.
        return os.path.join(d, fn)
    job_done_fns = list()
    for job_output in gathered:
        for fn in job_output.values():
            abs_fn = abspath(fn)
            job_done_fns.append(abs_fn)
    #import pprint
    #LOG.info('job_done_fns: {}'.format(pprint.pformat(job_done_fns)))
    job_rundirs = sorted(os.path.dirname(fn) for fn in job_done_fns)
    import time
    def find_las_paths():
        las_paths = list()
        for uow_dir in job_rundirs:
            # We could assert the existence of a job_done file here.
            d = os.path.abspath(uow_dir)
            las_path_glob = glob.glob(os.path.join(d, '*.las'))
            LOG.debug('dir={!r}, glob={!r}'.format(d, las_path_glob))
            if not las_path_glob:
                LOG.info('No .las found in {!r}. Sleeping {} seconds before retrying.'.format(d, WAIT))
                time.sleep(WAIT)
                las_path_glob = glob.glob(os.path.join(d, '*.las'))
                LOG.debug('dir={!r}, glob={!r}'.format(d, las_path_glob))
            if not las_path_glob:
                msg = 'No .las files found in daligner dir {!r}.'.format(d)
                raise MissingLas(msg)
            las_paths.extend(las_path_glob)
        n = len(las_paths)
        if not is_perfect_square(n):
            msg = '{} is not a perfect square. We must be missing some .las files.'.format(n)
            raise MissingLas(msg)
        return las_paths
    try:
        las_paths = sorted(find_las_paths(), key=lambda fn: os.path.basename(fn))
    except MissingLas:
        LOG.info('Not enough .las found from {!r}. Sleeping {} seconds before retrying.'.format(gathered_fn, WAIT))
        time.sleep(WAIT)
        las_paths = find_las_paths()
    io.serialize(las_paths_fn, las_paths)


def merge_split(config_fn, dbname, las_paths_fn, split_fn, bash_template_fn):
    with open(bash_template_fn, 'w') as stream:
        stream.write(pype_tasks.TASK_DB_LAMERGE_APPLY_SCRIPT)

    las_paths = io.deserialize(las_paths_fn)

    re_las_pair = re.compile(r'{db}\.(\d+)\.{db}\.(\d+)\.las$'.format(db=dbname))
    las_map = collections.defaultdict(list)
    for path in las_paths:
        mo = re_las_pair.search(path)
        if not mo:
            msg = '{!r} does not match regex {!r}'.format(
                path, re_las_pair.pattern)
            raise Exception(msg)
        a, b = int(mo.group(1)), int(mo.group(2))
        las_map[a].append(path)

    jobs = list()
    for i, block in enumerate(las_map):
        job_id = 'm_{:05d}'.format(i)

        # Write the las files for this job.
        input_dir = os.path.join('merge-scripts', job_id)
        las_paths_fn = os.path.join('.', input_dir, 'las-paths.json')
        io.mkdirs(input_dir)
        las_paths = las_map[block]
        io.serialize(las_paths_fn, las_paths)

        # Record in a job-dict.
        las_fn = '{}.{}.las'.format(dbname, block)
        job = dict()
        job['input'] = dict(
                config=config_fn,
                #db=db_fn,
                #script=os.path.abspath(script_fn),
                las_paths=las_paths_fn,
        )
        job['output'] = dict(
                #job_done = 'daligner.done'
                las_fn=las_fn
        )
        job['params'] = dict(
        )
        job['wildcards'] = dict(
                mer0_id=job_id,
        )
        jobs.append(job)
    io.serialize(split_fn, jobs)

def ichunked(seq, chunksize):
    """Yields items from an iterator in iterable chunks.
    https://stackoverflow.com/a/1335572
    """
    from itertools import chain, islice
    it = iter(seq)
    while True:
        yield chain([it.next()], islice(it, chunksize-1))

def merge_apply(las_paths_fn, las_fn):
    """Merge the las files into one, a few at a time.
    This replaces the logic of HPC.daligner.
    """
    io.rm_force(las_fn)

    #all_las_paths = rel_to(io.deserialize(las_paths_fn), os.path.dirname(las_paths_fn))
    all_las_paths = io.deserialize(las_paths_fn)

    # Create symlinks, so system calls will be shorter.
    all_syms = list()
    for fn in all_las_paths:
        symlink(fn)
        all_syms.append(os.path.basename(fn))
    curr_paths = sorted(all_syms)

    # Merge a few at-a-time.
    at_a_time = 250 # max is 252 for LAmerge
    level = 1
    while len(curr_paths) > 1:
        level += 1
        next_paths = list()
        for i, paths in enumerate(ichunked(curr_paths, at_a_time)):
            tmp_las = 'L{}.{}.las'.format(level, i+1)
            paths_arg = ' '.join(paths)
            cmd = 'LAmerge -v {} {}'.format(tmp_las, paths_arg)
            io.syscall(cmd)
            next_paths.append(tmp_las)
        curr_paths = next_paths

    io.syscall('mv -f {} {}'.format(curr_paths[0], 'keep-this'))
    io.syscall('#rm -f *.las')
    io.syscall('mv -f {} {}'.format('keep-this', las_fn))

def merge_combine(gathered_fn, las_paths_fn, block2las_fn):
    gathered = io.deserialize(gathered_fn)
    d = os.path.abspath(os.path.realpath(os.path.dirname(gathered_fn)))
    def abspath(fn):
        if os.path.isabs(fn):
            return fn # I expect this never to happen though.
        return os.path.join(d, fn)
    las_fns = list()
    for job_output in gathered:
        assert len(job_output) == 1, 'len(job_output) == {} != 1'.format(len(job_output))
        for fn in job_output.values():
            abs_fn = abspath(fn)
            las_fns.append(abs_fn)
    #import pprint
    #LOG.info('job_done_fns: {}'.format(pprint.pformat(job_done_fns)))

    import time
    #job_rundirs = sorted(os.path.dirname(fn) for fn in job_done_fns)
    #las_paths = list()
    #for uow_dir in job_rundirs:
    #    # We could assert the existence of a job_done file here.
    #    d = os.path.abspath(uow_dir)
    #    las_path_glob = glob.glob(os.path.join(d, '*.las'))
    #    LOG.debug('dir={!r}, glob={!r}'.format(d, las_path_glob))
    #    if not las_path_glob:
    #        time.sleep(WAIT)
    #        las_path_glob = glob.glob(os.path.join(d, '*.las'))
    #        LOG.debug('dir={!r}, glob={!r}'.format(d, las_path_glob))
    #    if not las_path_glob:
    #        msg = 'Missing .las file. Skipping block for dir {}'.format(d)
    #        LOG.error(msg)
    #    if len(las_path_glob) > 1:
    #        msg = 'Expected exactly 1 .las in {!r}, but found {}:\n {!r}'.format(
    #            d, len(las_path_glob), las_path_glob)
    #        LOG.warning(mgs)
    #    las_paths.extend(las_path_glob)

    las_paths = list()
    for las_fn in sorted(las_fns):
        if not os.path.exists(las_fn):
            msg = 'Did not find las-file {!r}. Waiting {} seconds.'.format(las_fn, WAIT)
            LOG.info(las_fn)
            time.sleep(WAIT)
            if not os.path.exists(las_fn):
                msg = 'Did not find las-file {!r}, even after waiting {} seconds. Maybe retry later?'.format(las_fn, WAIT)
                raise Exception(msg)
                #LOG.warning(las_fn)
                #continue
        las_paths.append(las_fn)

    # Map block nums to .las files.
    re_block = re.compile(r'\.(\d+)\.las$')
    block2las = dict()
    for fn in las_paths:
        mo = re_block.search(fn)
        if not mo:
            msg = 'Las file {!r} did not match regex {!r}.'.format(
                fn, re_block.pattern)
            raise Exception(msg)
        block = int(mo.group(1))
        block2las[block] = fn

    # Verify sequential block nums.
    blocks = sorted(block2las.keys())
    expected = list(range(1, len(blocks)+1))
    if blocks != expected:
        msg = '{!r} has {} .las files, but their block-numbers are not sequential: {!r}'.format(
            gathered_fn, len(blocks), blocks)
        raise Exception(msg)
    # Serialize result, plus an archaric file.
    io.serialize(las_paths_fn, sorted(block2las.values()))
    io.serialize(block2las_fn, block2las)


def setup_logging(log_level):
    hdlr = logging.StreamHandler(sys.stderr)
    hdlr.setLevel(log_level)
    hdlr.setFormatter(logging.Formatter('[%(levelname)s]%(message)s'))
    LOG.addHandler(hdlr)
    LOG.setLevel(logging.NOTSET)
    LOG.info('Log-level: {}'.format(log_level))

def cmd_build(args):
    ours = get_ours(args.config_fn, args.db_fn)
    build_db(ours, args.input_fofn_fn, args.db_fn, args.length_cutoff_fn)
def cmd_tan_split(args):
    ours = get_ours(args.config_fn, args.db_fn)
    tan_split(ours, args.config_fn, args.db_fn, args.split_fn, args.bash_template_fn)
def cmd_tan_apply(args):
    tan_apply(args.db_fn, args.script_fn, args.job_done_fn)
def cmd_tan_combine(args):
    tan_combine(args.db_fn, args.gathered_fn, args.new_db_fn)
def cmd_daligner_split(args):
    ours = get_ours(args.config_fn, args.db_fn)
    daligner_split(ours, args.config_fn, args.db_fn, args.nproc, args.wildcards, args.length_cutoff_fn, args.split_fn, args.bash_template_fn)
def cmd_daligner_apply(args):
    daligner_apply(args.db_fn, args.script_fn, args.job_done_fn)
def cmd_daligner_combine(args):
    daligner_combine(args.db_fn, args.gathered_fn, args.las_paths_fn)
def cmd_merge_split(args):
    merge_split(args.config_fn, args.db_prefix, args.las_paths_fn, args.split_fn, args.bash_template_fn)
def cmd_merge_apply(args):
    merge_apply(args.las_paths_fn, args.las_fn)
def cmd_merge_combine(args):
    merge_combine(args.gathered_fn, args.las_paths_fn, args.block2las_fn)

options_note = """

For raw_reads.db, we also look for the following config keys:

- pa_DBsplit_option
- pa_HPCdaligner_option
- pa_HPCTANmask_option
- length_cutoff: -1 => calculate based on "genome_size" and "seed_coverage" config.
- seed_coverage
- genome_size

For preads.db, these are named:

- ovlp_DBsplit_option
- ovlp_HPCdaligner_option
- ovlp_HPCTANmask_option
- length_cutoff_pr
"""

def get_ours(config_fn, db_fn):
    ours = dict()
    config = io.deserialize(config_fn)
    ours['genome_size'] = int(config['genome_size'])
    ours['seed_coverage'] = float(config['seed_coverage'])
    if os.path.basename(db_fn).startswith('preads'):
        ours['DBsplit_opt'] = config.get('ovlp_DBsplit_option', '')
        ours['daligner_opt'] = config.get('ovlp_HPCdaligner_option', '')
        ours['TANmask_opt'] = config.get('ovlp_HPCTANmask_option', '')
        ours['user_length_cutoff'] = int(config.get('length_cutoff_pr', '0'))
    else:
        ours['DBsplit_opt'] = config.get('pa_DBsplit_option', '')
        ours['daligner_opt'] = config.get('pa_HPCdaligner_option', '')
        ours['TANmask_opt'] = config.get('pa_HPCTANmask_option', '')
        ours['user_length_cutoff'] = int(config.get('length_cutoff', '0'))

    LOG.info('config({!r}):\n{}'.format(config_fn, config))
    LOG.info('our subset of config:\n{}'.format(ours))
    return ours

def add_build_arguments(parser):
    parser.add_argument(
        '--input-fofn-fn', required=True,
        help='input. User-provided file of fasta filename. Relative paths are relative to directory of FOFN.',
    )
    parser.add_argument(
        '--length-cutoff-fn', required=True,
        help='output. Simple file of a single integer, either calculated or specified by --user-length-cutoff.'
    )
def add_tan_split_arguments(parser):
    parser.add_argument(
        '--split-fn', default='tan-mask-uows.json',
        help='output. Units-of-work from TANmask, for datander.',
    )
    parser.add_argument(
        '--bash-template-fn', default='bash-template.sh',
        help='output. Script to apply later.',
    )
def add_tan_apply_arguments(parser):
    parser.add_argument(
        '--script-fn', required=True,
        help='input. Script to run datander.',
    )
    parser.add_argument(
        '--job-done-fn', default='job.done',
        help='output. Sentinel.',
    )
def add_tan_combine_arguments(parser):
    parser.add_argument(
        '--gathered-fn', required=True,
        help='input. List of sentinels. Produced by gen_parallel_tasks() gathering. The tan-track files are next to these.',
    )
    parser.add_argument(
        '--new-db-fn', required=True,
        help='output. This must match the input DB name. It will be symlinked, except the new track files.',
    )
def add_daligner_split_arguments(parser):
    parser.add_argument(
        '--wildcards', default='dummy_wildcard',
        help='Comma-separated string of keys to be subtituted into output paths for each job, if any. (Helps with snakemake and pypeflow; not needed in pbsmrtpipe, since outputs are pre-determined.)',
    )
    parser.add_argument(
        '--length-cutoff-fn', required=True,
        help='input. Simple file of a single integer, either calculated or specified by --user-length-cutoff.'
    )
    parser.add_argument(
        '--split-fn', default='daligner-mask-uows.json',
        help='output. Units-of-work from TANmask, for dadalignerder.',
    )
    parser.add_argument(
        '--bash-template-fn', default='bash-template.sh',
        help='output. Script to apply later.',
    )
def add_daligner_apply_arguments(parser):
    parser.add_argument(
        '--script-fn', required=True,
        help='input. Script to run dadalignerder.',
    )
    parser.add_argument(
        '--job-done-fn', default='job.done',
        help='output. Sentinel.',
    )
def add_daligner_combine_arguments(parser):
    parser.add_argument(
        '--gathered-fn', required=True,
        help='input. List of sentinels. Produced by gen_parallel_tasks() gathering. The .las files are next to these.',
    )
    parser.add_argument(
        '--las-paths-fn', required=True,
        help='output. JSON list of las paths.')
def add_merge_split_arguments(parser):
    parser.add_argument(
        '--db-prefix', required=True,
        help='DB is named "prefix.db", and the prefix is expected to match .las files',
    )
    parser.add_argument(
        '--las-paths-fn',
        help='input. foo.a.foo.b.las files from daligner.',
    )
    parser.add_argument(
        '--split-fn', default='merge-mask-uows.json',
        help='output. Units-of-work from TANmask, for damergeder.',
    )
    parser.add_argument(
        '--bash-template-fn', default='bash-template.sh',
        help='output. Script to apply later.',
    )
def add_merge_apply_arguments(parser):
    parser.add_argument(
        '--las-paths-fn', required=True,
        help='input. JSON list of las paths to merge. These must be .las "pairs" (foo.a.foo.b.las). The a-blocks must all be the same. Ultimately, we will generate a single .las from these, named after the a-block.')
    parser.add_argument(
        '--las-fn', required=True,
        help='output. The merged las-file.',
    )
def add_merge_combine_arguments(parser):
    parser.add_argument(
        '--gathered-fn', required=True,
        help='input. List of sentinels. Produced by gen_parallel_tasks() gathering. The .las files are next to these.',
    )
    parser.add_argument(
        '--las-paths-fn', required=True,
        help='output. JSON list of las paths.')
    parser.add_argument(
        '--block2las-fn', required=True,
        help='output. JSON dict of block (int) to las.')

def add_split_arguments(parser):
    parser.add_argument(
        '--run-jobs-fn',
        help='Output. Produced by HPC.daligner.',
    )
    parser.add_argument(
        '--rep-mask-uows-fn', default='rep-mask-uows.json',
        help='output. Units-of-work for REPmask.',
    )
    parser.add_argument(
        '--daligner-uows-fn', default='daligner-uows.json',
        help='output. Units-of-work for daligner.',
    )
    parser.add_argument(
        '--las-merge-uows-fn', default='las-merge-uows.json',
        help='output. Units-of-work for LAmerge.',
    )

class HelpF(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass

def parse_args(argv):
    description = 'Basic daligner steps: build; split into units-of-work; combine results and prepare for next step.'
    epilog = 'These tasks perform the split/apply/combine strategy (of which map/reduce is a special case).' + options_note
    parser = argparse.ArgumentParser(
        description=description,
        epilog=epilog,
        formatter_class=HelpF,
    )
    parser.add_argument(
        '--log-level', default='INFO',
        help='Python logging level.',
    )
    parser.add_argument(
        '--nproc', type=int, default=0,
        help='ignored for now, but non-zero will mean "No more than this."',
    )
    parser.add_argument(
        '--config-fn', required=True,
        help='Input. JSON of user-configuration. (This is probably the [General] section.)',
    )
    parser.add_argument(
        '--db-fn', default='raw_reads.db',
        help='Input or Output. Dazzler DB. (Dot-files are implicit.)',
    )

    help_build = 'build Dazzler DB for raw_reads; calculate length-cutoff for HGAP seed reads; split Dazzler DB into blocks; run DBdust to mask low-complexity regions'
    
    help_tan_split = 'generate units-of-work for datander, via HPC.TANmask'
    help_tan_apply = 'run datander and TANmask as a unit-of-work (according to output of HPC.TANmask), and remove .las files'
    help_tan_combine = 'run Catrack on partial track-files, to produce a "tan" mask'
    help_daligner_split = 'generate units-of-work for daligner, via HPC.daligner'
    help_daligner_apply = 'run daligner as a unit-of-work (according to output of HPC.TANmask)'
    help_daligner_combine = 'generate a file of .las files, plus units-of-work for LAmerge (alias for merge-split)'
    help_merge_split = 'generate a file of .las files, plus units-of-work for LAmerge (alias for daligner-combine)'
    help_merge_apply = 'run merge as a unit-of-work, and (possibly) remove .las files'
    help_merge_combine = 'generate a file of .las files'

    subparsers = parser.add_subparsers(help='sub-command help')

    parser_build = subparsers.add_parser('build',
            formatter_class=HelpF,
            description=help_build,
            help=help_build)
    add_build_arguments(parser_build)
    parser_build.set_defaults(func=cmd_build)

    parser_tan_split = subparsers.add_parser('tan-split',
            formatter_class=HelpF,
            description=help_tan_split,
            epilog='',
            help=help_tan_split)
    add_tan_split_arguments(parser_tan_split)
    parser_tan_split.set_defaults(func=cmd_tan_split)

    parser_tan_apply = subparsers.add_parser('tan-apply',
            formatter_class=HelpF,
            description=help_tan_apply,
            epilog='',
            help=help_tan_apply)
    add_tan_apply_arguments(parser_tan_apply)
    parser_tan_apply.set_defaults(func=cmd_tan_apply)

    parser_tan_combine = subparsers.add_parser('tan-combine',
            formatter_class=HelpF,
            description=help_tan_combine,
            epilog='The result will be mostly symlinks, plus new tan-track files. To use these as a mask, subsequent steps will need to add "-mtan".',
            help=help_tan_combine)
    add_tan_combine_arguments(parser_tan_combine)
    parser_tan_combine.set_defaults(func=cmd_tan_combine)

    parser_daligner_split = subparsers.add_parser('daligner-split',
            formatter_class=HelpF,
            description=help_daligner_split,
            epilog='HPC.daligner will be passed mask flags for any mask tracks which we glob. Note: if db is named "preads.db", we run daligner_p instead of daligner.',
            help=help_daligner_split)
    add_daligner_split_arguments(parser_daligner_split)
    parser_daligner_split.set_defaults(func=cmd_daligner_split)

    parser_daligner_apply = subparsers.add_parser('daligner-apply',
            formatter_class=HelpF,
            description=help_daligner_apply,
            epilog='',
            help=help_daligner_apply)
    add_daligner_apply_arguments(parser_daligner_apply)
    parser_daligner_apply.set_defaults(func=cmd_daligner_apply)

    parser_daligner_combine = subparsers.add_parser('daligner-combine',
            formatter_class=HelpF,
            description=help_daligner_combine,
            epilog='',
            help=help_daligner_combine)
    add_daligner_combine_arguments(parser_daligner_combine)
    parser_daligner_combine.set_defaults(func=cmd_daligner_combine)

    parser_merge_split = subparsers.add_parser('merge-split',
            formatter_class=HelpF,
            description=help_merge_split,
            epilog='HPC.merge will be passed mask flags for any mask tracks which we glob.',
            help=help_merge_split)
    add_merge_split_arguments(parser_merge_split)
    parser_merge_split.set_defaults(func=cmd_merge_split)

    parser_merge_apply = subparsers.add_parser('merge-apply',
            formatter_class=HelpF,
            description=help_merge_apply,
            epilog='Ultimately, there will be 1 .las per block.',
            help=help_merge_apply)
    add_merge_apply_arguments(parser_merge_apply)
    parser_merge_apply.set_defaults(func=cmd_merge_apply)

    parser_merge_combine = subparsers.add_parser('merge-combine',
            formatter_class=HelpF,
            description=help_merge_combine,
            epilog='',
            help=help_merge_combine)
    add_merge_combine_arguments(parser_merge_combine)
    parser_merge_combine.set_defaults(func=cmd_merge_combine)

    args = parser.parse_args(argv[1:])
    return args


def main(argv=sys.argv):
    args = parse_args(argv)
    setup_logging(args.log_level)
    args.func(args)


if __name__ == '__main__':  # pragma: no cover
    main()
