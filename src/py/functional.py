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

    E.g., each item will look like:
      (2, 1, 2, 3): 'daligner ...; LAsort ...; LAmerge ...; rm ...'

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

    Comments and lines starting with LAmerge are ignored.
    """
    re_block_dali = re.compile(r'%s\.(\d+)' %db_prefix)
    def blocks_dali(line):
        return [mo.group(1) for mo in re_block_dali.finditer(line)]
    # X == blocks[0]; A/B/C = blocks[...]

    re_pair_sort = re.compile(r'%s\.(\d+)\.%s\.(\d+)' %(db_prefix, db_prefix))
    def LAsort_pair(line):
        return re_pair_sort.search(line).group(1, 2)

    lines = [line.strip() for line in run_jobs_stream]
    assert any(len(l) > 1 for l in lines) # in case caller passed filename, not stream
    lines_dali = [l for l in lines if l.startswith('daligner')] # could be daligner_p
    lines_sort = [l for l in lines if l.startswith('LAsort')]
    pair2dali = {}
    for line in lines_dali:
        blocks = blocks_dali(line)
        for block in blocks[1:]:
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
        sorts = [pair2sort[pair] for pair in sorted(pairs, key=lambda k: (int(k[0]), int(k[1])))]
        id = tuple(map(int, blocks_dali(dali)))
        script = '\n'.join([dali] + sorts) + '\n'
        result[id] = script
    return result

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
