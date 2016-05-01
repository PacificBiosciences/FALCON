"""Purely functional code.
"""
import collections
import re

def _verify_pairs(pairs1, pairs2):
    if pairs1 != pairs2:
        print('pair2dali:', pairs1)
        print('pair2sort:', pairs2)
        print('pair2dali:', len(pairs1))
        print('pair2sort:', len(pairs2))
        assert pairs1 == pairs2

def get_daligner_job_descriptions(run_jobs_stream, db_prefix):
    """Return a dict of job-desc-tuple -> HPCdaligner bash-job.

    Comments and lines starting with LAmerge are ignored.

    E.g., each item will look like:
      ('.2', '.1', '.2', '.3'): 'daligner ...; LAsort ...; LAmerge ...; rm ...'

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
      L.1.X.A, L.1.X.B, and L.1.X.C
    where A, B, or C could be X.
    (In the example, X=2 A=1 B=2.)

    Oddly, b/c of a recent change by GM, if there is only 1 block, then the suffix string is empty.
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
        for block in blocks:
            pair = (blocks[0], block)
            pair2dali[pair] = line
            if block != blocks[0]:
                # Then we have a reverse comparison too.
                # https://dazzlerblog.wordpress.com/2014/07/10/dalign-fast-and-sensitive-detection-of-all-pairwise-local-alignments/
                rpair = (block, blocks[0])
                pair2dali[rpair] = line
    pair2sort = {}
    for line in lines_sort:
        pair = LAsort_pair(line)
        pair2sort[pair] = line
    _verify_pairs(sorted(pair2dali.keys()), sorted(pair2sort.keys()))
    dali2pairs = collections.defaultdict(set)
    for pair, dali in pair2dali.items():
        dali2pairs[dali].add(pair)
    result = {}
    for dali, pairs in dali2pairs.items():
        sorts = [pair2sort[pair] for pair in sorted(pairs, key=lambda k: (
            (int(k[0][1:]) if k[0].startswith('.') else 0),
            (int(k[1][1:]) if k[1].startswith('.') else 0)
        ))]
        id = tuple(blocks_dali(dali))
        script = '\n'.join([dali] + sorts) + '\n'
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

def get_mjob_data(run_jobs_stream):
    """Given output of HPC.daligner,
    return {int: [bash-lines]}
    """
    f = run_jobs_stream

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
    assert target <= total, (target, total)
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
