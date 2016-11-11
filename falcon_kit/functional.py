"""Purely functional code.
"""
from __future__ import division
import collections
import re
import StringIO

def _verify_pairs(pairs1, pairs2):
    if pairs1 != pairs2:
        print('pair2dali:', pairs1)
        print('pair2sort:', pairs2)
        print('dali-sort:', set(pairs1) - set(pairs2))
        print('sort-dali:', set(pairs2) - set(pairs1))
        print('pair2dali:', len(pairs1))
        print('pair2sort:', len(pairs2))
        assert pairs1 == pairs2

def skip_LAcheck(bash):
    def lines():
        for line in StringIO.StringIO(bash):
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
    for k,v in descs.iteritems():
        result[k] = skip_LAcheck(v)
    return result

def get_daligner_job_descriptions(run_jobs_stream, db_prefix, single=False):
    """Return a dict of job-desc-tuple -> HPCdaligner bash-job.

    Comments and lines starting with LAmerge are ignored.

    E.g., each item will look like:
      ('.2', '.1', '.2', '.3'): 'daligner ...; LAsort ...; LAmerge ...; LAcheck ...; rm ...'

    Rationale
    ---------
    For i/o efficiency, we will combine daligner calls with LAsort lines, which include 0-level merge.
    Example:
      daligner -v -t16 -H12000 -e0.7 -s1000 raw_reads.2 raw_reads.1 raw_reads.2
    That would be combined with two LAsort lines:
      LAsort -v raw_reads.2.raw_reads.1.C0 ...
      LAsort -v raw_reads.2.raw_reads.2.C0 ...
    For each returned job, the result of
      daligner X A B C; LAsort*
    will then be
      L1.X.A, L1.X.B, and L1.X.C
    where A, B, or C could be X.
    (In the example, X=2 A=1 B=2.)

    Oddly, b/c of a recent change by GM, if there is only 1 block,
    then the merged file is named differently. To achieve that, use single=True.
    """
    re_block_dali = re.compile(r'%s(\.\d+|)' %db_prefix)
    def blocks_dali(line):
        """Return ['.1', '.2', ...]
        Can return [''] if only 1 block.
        """
        return [mo.group(1) for mo in re_block_dali.finditer(line)]
    # X == blocks[0]; A/B/C = blocks[...]

    re_pair_sort = re.compile(r'%s(\.\d+|)\.%s(\.\d+|)' %(db_prefix, db_prefix))
    def LAsort_pair(line):
        """Return [('.1', '.1'), ('.1', '.2'), ('.2', '.1'), ...]
        Can return [('', '')] if only 1 block.
        """
        mo = re_pair_sort.search(line)
        if not mo:
            raise Exception('Pattern {!r} does not match line {!r}'.format(
                re_pair_sort.pattern, line))
        return mo.group(1, 2)

    lines = [line.strip() for line in run_jobs_stream]
    assert any(len(l) > 1 for l in lines), repr(lines) # in case caller passed filename, not stream
    lines_dali = [l for l in lines if l.startswith('daligner')] # could be daligner_p
    lines_sort = [l for l in lines if l.startswith('LAsort')]
    pair2dali = {}
    for line in lines_dali:
        blocks = blocks_dali(line)
        block0 = blocks[0]
        for block in blocks[1:]:
            pair = (block0, block)
            pair2dali[pair] = line
            if block != block0:
                # Then we have a reverse comparison too.
                # https://dazzlerblog.wordpress.com/2014/07/10/dalign-fast-and-sensitive-detection-of-all-pairwise-local-alignments/
                rpair = (block, block0)
                pair2dali[rpair] = line
    pair2sort = {}
    for line in lines_sort:
        pair = LAsort_pair(line)
        pair2sort[pair] = line
    _verify_pairs(sorted(pair2dali.keys()), sorted(pair2sort.keys()))
    dali2pairs = collections.defaultdict(set)
    total_pairs = 0
    for pair, dali in pair2dali.items():
        dali2pairs[dali].add(pair)
        total_pairs += 1
    if single:
        assert total_pairs == 1, 'In single-mode, but total_pair={}'.format(total_pairs)
        assert list(pair2dali.keys())[0] == ('.1', '.1'), repr(pair2dali)
    result = {}
    for dali, pairs in dali2pairs.items():
        pairs = list(pairs)
        pairs.sort( key=lambda k: ( (int(k[0][1:]) if k[0].startswith('.') else 0), (int(k[1][1:]) if k[1].startswith('.') else 0) ) )
        sorts = [pair2sort[pair] for pair in pairs]
        id = tuple(blocks_dali(dali))
        early_checks = [ "LAcheck -v {db_prefix} *.las".format( db_prefix = db_prefix )  ]
        if single:
            checks = [ "LAcheck -vS {db_prefix} {db_prefix}.1".format( db_prefix = db_prefix )  ]
        else:
            checks = [ "LAcheck -vS {db_prefix} L1{p1}{p2}".format( db_prefix = db_prefix, p1=pair[0], p2=pair[1]) for pair in pairs ]

        script = '\n'.join([dali] + early_checks + sorts + checks) + '\n'
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
    for p_id, bash_lines in mjob_data.iteritems():
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
                re_las_name.pattern, bash_lines[i]))
        las_fn = mo.group(1) + '.las'
        result[p_id] = las_fn
    return result

def get_mjob_data(run_jobs_stream):
    """Given output of HPC.daligner,
    return {int: [bash-lines]}
    """
    f = run_jobs_stream

    ## Strip either '&& rm ...' or '; rm ...' ?
    #re_strip_rm = re.compile(r'^(.*) ((\&\&)|;) .*$')

    # Copied from scripts_merge()
    mjob_data = {}
    for l in f:
        l = l.strip()
        first_word = l.split()[0]
        if first_word not in ("LAsort", "LAmerge"):
            continue
        if first_word in ["LAsort"]:
            # We now run this part w/ daligner, but we still need
            # a small script for some book-keeping.
            p_id = first_block_las(l)
            mjob_data.setdefault( p_id, [] )
            #mjob_data[p_id].append(  " ".join(l) ) # Already done w/ daligner!
        elif first_word in ["LAmerge"]:
            p_id = first_block_las(l)
            mjob_data.setdefault( p_id, [] )
            #l = re_strip_rm.sub(r'\1', l) # rm is very safe if we run in /tmp
            mjob_data[p_id].append(l)
    return mjob_data

_re_sub_daligner = re.compile(r'^daligner\b', re.MULTILINE)
def xform_script_for_preads(script):
    daligner_exe = 'daligner_p'
    return _re_sub_daligner.sub(daligner_exe, script) #, flags=re.MULTILINE) # flags in py2.7

def xform_script_for_raw_reads(script):
    return script

def get_script_xformer(pread_aln):
    if pread_aln:
        return xform_script_for_preads
    else:
        return xform_script_for_raw_reads

def calc_cutoff_from_reverse_sorted_readlength_counts(rl_counts, target):
    """Return first read_len which gives at least 'target' bases.
    """
    total = sum(pair[0]*pair[1] for pair in rl_counts)
    subtotal = 0
    assert target <= total, 'Not enough genome coverage (target={} < actual={})'.format(target, total)
    cutoff = 0
    for (rl, count) in rl_counts:
        subtotal += rl*count
        if subtotal >= target:
            cutoff = rl
            break
    else:
        raise Exception('Impossible target: target={target}, subtotal={subtotal}, total={total}'.format(locals()))
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
    re_stat = re.compile(r'^\s*(?P<bin>\S+):\s+(?P<count>\S+)\s+\S+\s+\S+\s+\S+\s*$')
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
    rl_counts = get_reverse_sorted_readlength_counts_from_DBstats(DBstats_output)
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
        if not line: continue
        yield tuple(int(x) for x in line.split())

def weighted_average(cols):
    """Given tuples of (weight, value),
    return weighted average.

    >>> weighted_average(((100, 1), (200, 2), (100, 5)))
    2.5
    """
    return sum(w*v for (w,v) in cols) / sum(w for (w,v) in cols)

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
    for k,va in dictA.iteritems():
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
    pread_lengths = dict((k,v) for (k,v) in cols)
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
