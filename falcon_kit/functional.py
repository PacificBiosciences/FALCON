"""Purely functional code.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function


from future.utils import viewitems
from .io import NativeIO as StringIO
import collections
import logging
import re

#LOG = logging.getLogger(__name__)


def _verify_pairs(pairs1, pairs2):
    if pairs1 != pairs2:  # pragma: no cover
        print('pair2dali:', pairs1)
        print('pair2sort:', pairs2)
        print('dali-sort:', set(pairs1) - set(pairs2))
        print('sort-dali:', set(pairs2) - set(pairs1))
        print('pair2dali:', len(pairs1))
        print('pair2sort:', len(pairs2))
        assert pairs1 == pairs2


def skip_LAcheck(bash):
    def lines():
        for line in StringIO(bash):
            if 'LAcheck' in line:
                yield 'set +e\n'
                yield line
                yield 'set -e\n'
            else:
                yield line
    return ''.join(lines())


def get_daligner_job_descriptions_sans_LAcheck(run_jobs_stream, db_prefix, single=False):
    """Strip LAcheck (somehow) from each bash script.
    (For now, we will run it but not fail on error.)
    """
    descs = get_daligner_job_descriptions(run_jobs_stream, db_prefix, single)
    result = {}
    for (k, v) in viewitems(descs):
        bash = skip_LAcheck(v)
        bash = bash.replace(
            'LAsort', 'python2.7 -m falcon_kit.mains.LAsort {}'.format(db_prefix))
        bash = bash.replace(
            'LAmerge', 'python2.7 -m falcon_kit.mains.LAmerge {}'.format(db_prefix))
        result[k] = bash
    return result


def get_daligner_job_descriptions(run_jobs_stream, db_prefix, single=False):
    """Return a dict of job-desc-tuple -> HPCdaligner bash-job.

    Comments are ignored.

    E.g., each item will look like:
      ('.2', '.1', '.2', '.3'): 'daligner

    Rationale
    ---------
    For i/o efficiency, we want to daligner calls with LAsort/LAmerge lines. But
    Gene has done this himself too. So now, we only want the daligner calls here.

    Later, we will do the extra LAmerge lines, grouped by A-read.
    """
    re_block_dali = re.compile(r'%s(\.\d+|)' % db_prefix)

    def blocks_dali(line):
        """Return ['.1', '.2', ...]
        Can return [''] if only 1 block.
        """
        return [mo.group(1) for mo in re_block_dali.finditer(line)]
    # X == blocks[0]; A/B/C = blocks[...]

    lines = [line.strip() for line in run_jobs_stream]
    # in case caller passed filename, not stream
    assert any(len(l) > 1 for l in lines), repr('\n'.join(lines))

    lines_dali = [l for l in lines if l.startswith(
        'daligner')]  # could be daligner_p
    result = {}
    for dali in lines_dali:
        id = tuple(blocks_dali(dali))
        early_checks = [
            "LAcheck -v {db_prefix} *.las".format(db_prefix=db_prefix)]
        script = '\n'.join([dali] + early_checks) + '\n'
        result[id] = script
    return result


re_first_block_las = re.compile(r'^(?:\S+)(?:\s+-\S+)*\s+[^\.]+\.(\d+|)')


def first_block_las(line):
    """
    >>> first_block_las('LAsort -v -a foo.1.foo.1.C0')
    1
    """
    mo = re_first_block_las.search(line)
    try:
        return int(mo.group(1))
    except Exception as e:
        raise Exception('Pattern {!r} does not match line {!r}: {}'.format(
            re_first_block_las.pattern, line, e))


def get_las_filenames(mjob_data, db_prefix):
    """Given result of get_mjob_data(),
    return {int: final_las_filename}
    """
    # This is our best guess.
    # (We might not even need this, since we know the output filename of each merge-task by convention.)
    # Eventually, we need to re-write HPC.daligner.
    result = {}
    re_LAmerge = re.compile(r'^LAmerge\s+(?:\-\S+\s+)(\S+)')
    re_LAcheck = re.compile(r'^LAcheck\s+(?:\-\S+\s+)\S+\s+(\S+)')
    for (p_id, bash_lines) in viewitems(mjob_data):
        if not bash_lines:
            # The daligner+LAsort+LAmerge job produced the final .las
            # for this block. We will symlink it later.
            las_fn = '{}.{}.las'.format(db_prefix, p_id)
            result[p_id] = las_fn
            continue
        # Find the last line which can tell us the final .las name.
        i = len(bash_lines) - 1
        while bash_lines[i].split()[0] not in ('LAmerge', 'LAcheck'):
            i -= 1
        # Now we will raise an exception if there were none. But in theory, there
        # must be at least an LAsort.
        first_word = bash_lines[i].split()[0]
        if first_word == 'LAmerge':
            regex = re_LAmerge
        elif first_word == 'LAcheck':
            regex = re_LAcheck
        else:
            raise Exception('first_word={!r} in line#{} of {!r}'.format(
                first_word, i, bash_lines))
        mo = regex.search(bash_lines[i])
        if not mo:
            raise Exception('Regex {!r} failed on {!r}'.format(
                regex.pattern, bash_lines[i]))
        las_fn = mo.group(1) + '.las'
        result[p_id] = las_fn
    return result


def get_mjob_data(run_jobs_stream):
    """Given output of HPC.daligner,
    return {int: [bash-lines]}
    """
    f = run_jobs_stream

    # Strip either '&& rm ...' or '; rm ...' ?
    #re_strip_rm = re.compile(r'^(.*) ((\&\&)|;) .*$')

    # Copied from scripts_merge()
    mjob_data = {}
    for l in f:
        l = l.strip()
        first_word = l.split()[0]
        if first_word not in ("LAsort", "LAmerge", "rm"):
            continue
        if first_word in ["LAsort"]:
            # We now run this part w/ daligner, but we still need
            # a small script for some book-keeping.
            p_id = first_block_las(l)
            mjob_data.setdefault(p_id, [])
            # mjob_data[p_id].append(  " ".join(l) ) # Already done w/ daligner!
            raise Exception('We do not expect to see LAsort at all anymore.')
        elif first_word in ["LAmerge"]:
            p_id = first_block_las(l)
            mjob_data.setdefault(p_id, [])
            # l = re_strip_rm.sub(r'\1', l) # (LAmerge && rm) rm is very safe if we run in /tmp
            mjob_data[p_id].append(l)
            #LOG.info('{}: {}'.format(p_id, l))
        elif first_word in ["rm"]:
            p_id = first_block_las(l)
            mjob_data.setdefault(p_id, [])
            mjob_data[p_id].append(l)
            #LOG.info('rm{}: {}'.format(p_id, l))
    #for key, data in mjob_data.items():
    #    mjob_data[key] = '\n'.join(data)
    return mjob_data


def yield_args_from_line(bash_line):
    """Given a line of LAmerge, etc.,
    return [output_las_fn, input_las_fn0, input_las_fn1, ...]
    """
    for word in bash_line.split():
        if word.startswith('-') or word in ('LAcheck', 'LAmerge', 'LAsort'):
            continue
        yield word


_re_sub_daligner = re.compile(r'^daligner\b', re.MULTILINE)


def xform_script_for_preads(script):
    daligner_exe = 'daligner_p'
    # , flags=re.MULTILINE) # flags in py2.7
    return _re_sub_daligner.sub(daligner_exe, script)


def xform_script_for_raw_reads(script):
    return script


def get_script_xformer(pread_aln):
    if pread_aln:
        return xform_script_for_preads
    else:
        return xform_script_for_raw_reads


class GenomeCoverageError(Exception):
    pass


def calc_cutoff_from_reverse_sorted_readlength_counts(rl_counts, target):
    """Return first read_len which gives at least 'target' bases.
    """
    total = sum(pair[0] * pair[1] for pair in rl_counts)
    subtotal = 0
    if target > total:
        msg = 'Not enough reads available for desired genome coverage (bases needed={} > actual={})'.format(
            target, total)
        raise GenomeCoverageError(msg)
    cutoff = 0
    for (rl, count) in rl_counts:
        subtotal += rl * count
        if subtotal >= target:
            cutoff = rl
            break
    else:  # pragma: no cover
        msg = 'Impossible target (probably a bug): target={target}, subtotal={subtotal}, total={total}'.format(
            locals())
        raise Exception(msg)
    return cutoff


def num2int(num):
    """
    >>> num2int('1,000,000')
    1000000
    """
    return int(num.replace(',', ''))


def get_reverse_sorted_readlength_counts_from_DBstats(DBstats_output):
    """Return pairs of (readlength, count).
        Bin:      Count  % Reads  % Bases     Average
    169,514:          1      0.0      0.0      169514
    ...
    ->
    [(169514, 1), ...]
    """
    rl_counts = list()
    lines = DBstats_output.splitlines()
    re_stat = re.compile(
        r'^\s*(?P<bin>\S+):\s+(?P<count>\S+)\s+\S+\s+\S+\s+\S+\s*$')
    for line in lines:
        match = re_stat.search(line)
        if not match:
            continue
        rl = num2int(match.group('bin'))
        count = num2int(match.group('count'))
        rl_counts.append((rl, count))
    return rl_counts


def calc_cutoff(target, DBstats_output):
    """Calculate the length_cutoff needed for at least 'target' bases.
    DBstats_output: ASCII output of 'DBstats -b1 DB',
    """
    rl_counts = get_reverse_sorted_readlength_counts_from_DBstats(
        DBstats_output)
    return calc_cutoff_from_reverse_sorted_readlength_counts(rl_counts, target)


def parse_2columns_of_ints(data):
    r"""Given 2 columns of integers,
    space- and line-delimited,
    yield tuples.

    >>> tuple(parse_2columns_of_ints("1 2\n3 4"))
    ((1, 2), (3, 4))
    """
    for line in data.splitlines():
        line = line.strip()
        if not line:
            continue
        yield tuple(int(x) for x in line.split())


def weighted_average(cols):
    """Given tuples of (weight, value),
    return weighted average.

    >>> weighted_average(((100, 1), (200, 2), (100, 5)))
    2.5
    """
    return sum(w * v for (w, v) in cols) / sum(w for (w, v) in cols)


def parsed_readlengths_from_dbdump_output(output):
    """Given output text from the DBump command,
    yield all read-lengths.
    """
    re_length = re.compile('^L\s+\d+\s+(\d+)\s+(\d+)$')
    for line in output.splitlines():
        mo = re_length.search(line)
        if mo:
            beg, end = mo.group(1, 2)
            beg = int(beg)
            end = int(end)
            yield end - beg


def mapped_readlengths_from_dbdump_output(output):
    """Given output text from the DBump command,
    return dict of (id => read-length).
    There will be alternate lines like these:
      R #
      L # # #
    https://dazzlerblog.wordpress.com/command-guides/dazz_db-command-guide/
    """
    lengths = dict()
    re_rid = re.compile('^R\s+(\d+)$')
    re_length = re.compile('^L\s+(\d+)\s+(\d+)\s+(\d+)$')
    for line in output.splitlines():
        mo = re_rid.search(line)
        if mo:
            idx = int(mo.group(1))
            continue
        mo = re_length.search(line)
        if mo:
            well, beg, end = mo.group(1, 2, 3)
            well = int(idx)
            beg = int(beg)
            end = int(end)
            lengths[idx] = (end - beg)
            continue
    return lengths


def average_difference(dictA, dictB):
    """Return the average difference of
    values in dictA minus dictB, only
    using values in dictA.
    If a value is missing from dictB, raise Exception.
    """
    total_diff = 0.0
    for (k, va) in viewitems(dictA):
        vb = dictB[k]
        total_diff += (va - vb)
    return total_diff / len(dictA)


def calc_metric_fragmentation(perl_counts_output):
    # """perl -e 'while (<>) { if ( m{>[^/]+/(\d+)\d/} ) { $id{$1}++; } }; while (my ($k, $v) = each %%id) { $counts{$v}++; }; while (my ($k, $v) = each %%counts) { print "$v $k\n"; };' %s""" %(fastas)
    cols = tuple(parse_2columns_of_ints(perl_counts_output))
    avg = weighted_average(cols)
    return avg


def calc_metric_truncation(dbdump_output, length_pairs_output):
    # """perl -e 'while (<>) { if ( m{>[^/]+/0*(\d+)\d/(\d+)_(\d+)} ) { $lengths{$1} += ($3 - $2); } }; while (my ($k, $v) = each %%lengths) { print "$k $v\n"; };' %s""" %(fastas)
    cols = tuple(parse_2columns_of_ints(length_pairs_output))
    pread_lengths = dict((k, v) for (k, v) in cols)
    orig_lengths = mapped_readlengths_from_dbdump_output(dbdump_output)
    avg = -average_difference(pread_lengths, orig_lengths)
    return avg


def choose_cat_fasta(fofn):
    """Given the contents of a fasta FOFN,
    return a command to write the contents of a fasta to stdout,
    keeping the original file.
    Raise Exception on error.

    >>> choose_cat_fasta('abc.gz')
    'zcat '
    >>> choose_cat_fasta('abc.dexta')
    'undexta -vkU -w60 -i < '
    >>> choose_cat_fasta('abc')
    'cat '
    """
    first_line = fofn.splitlines()[0]
    if first_line.endswith('.gz'):
        return 'zcat '
    elif first_line.endswith('.dexta'):
        return 'undexta -vkU -w60 -i < '
    else:
        return 'cat '


re_underscore_flag = re.compile(r'(--[\w-]+)(_)')
def dash_flags(val):
    """
    >>> dash_flags('--foo_bar --one_two_three')
    '--foo-bar --one-two-three'
    >>> dash_flags('')
    ''
    """
    while True:
        # Repeat until settled, as there might be multiple _ in the same flag.
        new_val = re_underscore_flag.sub(r'\1-', val)
        if new_val == val:
            return new_val
        val = new_val


def cfg_tobool(v):
    """
    >>> cfg_tobool('yes')
    True
    >>> cfg_tobool('true')
    True
    >>> cfg_tobool('T')
    True
    >>> cfg_tobool('1')
    True
    >>> cfg_tobool('no')
    False
    >>> cfg_tobool('false')
    False
    >>> cfg_tobool('F')
    False
    >>> cfg_tobool('0')
    False
    >>> cfg_tobool('')
    False
    """
    if v in (True, False, None):
        return v
    if not v:
        return False
    if v.upper()[0] in ('T', 'Y'):
        return True
    if v.upper()[0] in ('F', 'N'):
        return False
    return bool(int(v))


# https://stackoverflow.com/questions/3387691/how-to-perfectly-override-a-dict
# We derived from dict instead of from MutableMapping to json.dumps() works.

_RaiseKeyError = object() # singleton for no-default behavior

class LowerDict(dict):  # dicts take a mapping or iterable as their optional first argument
    __slots__ = () # no __dict__ - that would be redundant
    def __init__(self):
        # No args allowed, to keep it simple.
        super(LowerDict, self).__init__(self)
    def __getitem__(self, k):
        return super(LowerDict, self).__getitem__(k.lower())
    def __setitem__(self, k, v):
        return super(LowerDict, self).__setitem__(k.lower(), v)
    def __delitem__(self, k):
        return super(LowerDict, self).__delitem__(k.lower())
    def get(self, k, default=None):
        return super(LowerDict, self).get(k.lower(), default)
    def setdefault(self, k, default=None):
        return super(LowerDict, self).setdefault(k.lower(), default)
    def pop(self, k, v=_RaiseKeyError):
        if v is _RaiseKeyError:
            return super(LowerDict, self).pop(k.lower())
        return super(LowerDict, self).pop(k.lower(), v)
    #def update(self, mapping=(), **kwargs):
    #    super(LowerDict, self).update(self._process_args(mapping, **kwargs))
    def __contains__(self, k):
        return super(LowerDict, self).__contains__(k.lower())
    #def copy(self): # don't delegate w/ super - dict.copy() -> dict :(
    #    return type(self)(self)
    @classmethod
    def fromkeys(cls, keys, v=None):
        return super(LowerDict, cls).fromkeys((k.lower() for k in keys), v)
    def __repr__(self):
        return '{0}({1})'.format(type(self).__name__, super(LowerDict, self).__repr__())


__loop_set = set()

def toLowerDict(cfg):
    """Change key-names to be lower-case, at all levels of dict cfg.
    Then, return the case-insensitive LowerDict, substituted recursively.
    """
    if isinstance(cfg, LowerDict):
        return cfg
    if id(cfg) in __loop_set:
        # Prevent infinite loop.
        raise Exception('Already ran update_lowercase({}) (len(set)=={}):\n  {}'.format(
            id(cfg), len(__loop_set), cfg))
    __loop_set.add(id(cfg))

    low = LowerDict()

    for k,v in list(cfg.items()):
        if isinstance(v, dict):
            v = toLowerDict(v) # RECURSION
        if k in low:
            msg = 'Collision for "{}" in dict:\n{}'.format(k, cfg)
            if v != low[k]:
                raise Exception(msg)
        low[k] = v
    return low
