import helpers
from nose.tools import assert_equal, assert_raises, eq_
import falcon_kit.functional as f
import StringIO
import collections
import os

thisdir = os.path.dirname(os.path.abspath(__file__))
example_HPCdaligner_fn = os.path.join(thisdir, 'HPCdaligner_synth0.sh')
example_HPCdaligner_small_fn = os.path.join(thisdir, 'HPCdaligner_synth0_preads.sh')

def test_get_daligner_job_descriptions():
    example_HPCdaligner = open(example_HPCdaligner_fn)
    result = f.get_daligner_job_descriptions(
            example_HPCdaligner, 'raw_reads')
    assert result
    helpers.equal_multiline(result[('.1', '.1')], "daligner -v -h1 -t16 -H1 -e0.7 -l1 -s1000 raw_reads.1 raw_reads.1\nLAcheck -v raw_reads *.las\nLAsort -v raw_reads.1.raw_reads.1.C0 raw_reads.1.raw_reads.1.N0 && LAmerge -v L1.1.1 raw_reads.1.raw_reads.1.C0.S raw_reads.1.raw_reads.1.N0.S && rm raw_reads.1.raw_reads.1.C0.S.las raw_reads.1.raw_reads.1.N0.S.las\nLAcheck -vS raw_reads L1.1.1\n")
    helpers.equal_multiline(result[('.2', '.1', '.2')], "daligner -v -h1 -t16 -H1 -e0.7 -l1 -s1000 raw_reads.2 raw_reads.1 raw_reads.2\nLAcheck -v raw_reads *.las\nLAsort -v raw_reads.1.raw_reads.2.C0 raw_reads.1.raw_reads.2.N0 && LAmerge -v L1.1.2 raw_reads.1.raw_reads.2.C0.S raw_reads.1.raw_reads.2.N0.S && rm raw_reads.1.raw_reads.2.C0.S.las raw_reads.1.raw_reads.2.N0.S.las\nLAsort -v raw_reads.2.raw_reads.1.C0 raw_reads.2.raw_reads.1.N0 && LAmerge -v L1.2.1 raw_reads.2.raw_reads.1.C0.S raw_reads.2.raw_reads.1.N0.S && rm raw_reads.2.raw_reads.1.C0.S.las raw_reads.2.raw_reads.1.N0.S.las\nLAsort -v raw_reads.2.raw_reads.2.C0 raw_reads.2.raw_reads.2.N0 && LAmerge -v L1.2.2 raw_reads.2.raw_reads.2.C0.S raw_reads.2.raw_reads.2.N0.S && rm raw_reads.2.raw_reads.2.C0.S.las raw_reads.2.raw_reads.2.N0.S.las\nLAcheck -vS raw_reads L1.1.2\nLAcheck -vS raw_reads L1.2.1\nLAcheck -vS raw_reads L1.2.2\n")
    eq_(len(result), 2)

def test_get_daligner_job_descriptions_small():
    # when there is only 1 block, a special case
    example_HPCdaligner = open(example_HPCdaligner_small_fn)
    result = f.get_daligner_job_descriptions(
            example_HPCdaligner, 'preads', single=True)
    assert result
    helpers.equal_multiline(result[('.1', '.1')], "daligner -v -h1 -t50 -H1 -e0.99 -l1 -s1000 preads.1 preads.1\nLAcheck -v preads *.las\nLAsort -v preads.1.preads.1.C0 preads.1.preads.1.N0 preads.1.preads.1.C1 preads.1.preads.1.N1 preads.1.preads.1.C2 preads.1.preads.1.N2 preads.1.preads.1.C3 preads.1.preads.1.N3 && LAmerge -v preads.1 preads.1.preads.1.C0.S preads.1.preads.1.N0.S preads.1.preads.1.C1.S preads.1.preads.1.N1.S preads.1.preads.1.C2.S preads.1.preads.1.N2.S preads.1.preads.1.C3.S preads.1.preads.1.N3.S\nLAcheck -vS preads preads.1\n")
    eq_(len(result), 1)

example_se161 = os.path.join(thisdir, 'se161.sh')

def test_get_daligner_job_descriptions_se161():
    example_HPCdaligner = open(example_se161)
    result = f.get_daligner_job_descriptions(
            example_HPCdaligner, 'raw_reads', single=False)
    assert result

def test_get_mjob_data():
    example_HPCdaligner = open(example_HPCdaligner_fn)
    result = f.get_mjob_data(
            example_HPCdaligner)
    assert result
    eq_(result[1], ['LAmerge -v raw_reads.1 L1.1.1 L1.1.2 && rm L1.1.1.las L1.1.2.las'])
    eq_(result[2], ['LAmerge -v raw_reads.2 L1.2.1 L1.2.2 ; rm L1.2.1.las L1.2.2.las'])

def test_skip_LAcheck():
    orig = """set -e
hello there
LAcheck foo bar
middle
LAcheck -vS foo bar
goodbye
"""
    got = f.skip_LAcheck(orig)
    expected = """set -e
hello there
set +e
LAcheck foo bar
set -e
middle
set +e
LAcheck -vS foo bar
set -e
goodbye
"""
    eq_(got, expected)

def test_first_block_las():
    line = 'LAsort -v -a -q foo.1.foo.1.C0'
    result = f.first_block_las(line)
    eq_(result, 1)
    line = 'LAsort foo.1.foo.1.C0'
    result = f.first_block_las(line)
    eq_(result, 1)
    # We expect HPC.daligner to use block numbers always.
    line = 'LAsort -v -a -q foo.foo.C0'
    assert_raises(Exception, f.first_block_las, line)
    line = 'LAsort foo.foo.C0'
    assert_raises(Exception, f.first_block_las, line)

def test_xform_script_for_preads():
    # Technically, we never have more than one daligner in a script, but that
    # could change in pbsmrtpipe, since it limits the number of chunks.
    script = 'daligner x y\nLAsort a b\ndaligner x1 y1\n'
    expected = 'daligner_p x y\nLAsort a b\ndaligner_p x1 y1\n'
    result = f.xform_script_for_preads(script)
    eq_(result, expected)

    script = 'daligner x y\nLAsort a b\ndaligner x1 y1\n'
    expected = script # no-op
    result = f.xform_script_for_raw_reads(script)
    eq_(result, expected)

    eq_(f.get_script_xformer(True), f.xform_script_for_preads)
    eq_(f.get_script_xformer(False), f.xform_script_for_raw_reads)

def test_calc_cutoff_from_reverse_sorted_readlength_counts():
    read_lens = [1,2,2,3]
    rs_rl_counts = [(3, 1), (2, 2), (1, 1)]
    assert collections.Counter(read_lens) == dict(rs_rl_counts)
    def check(target, expected):
        got = f.calc_cutoff_from_reverse_sorted_readlength_counts(rs_rl_counts, target)
        assert_equal(expected, got)
    for n, expected in (
            (8, 1),
            (7, 2),
            (4, 2),
            (3, 3),
            (1, 3),
            (0, 3),
        ):
        yield check, n, expected
    assert_raises(Exception, check, 9, None)

capture = """
Statistics for all reads of length 500 bases or more

        240,548 reads        out of         309,410  ( 77.7%)
  1,117,649,525 base pairs   out of   1,320,493,403  ( 84.6%)

          4,646 average read length
          6,418 standard deviation

  Base composition: 0.264(A) 0.247(C) 0.279(G) 0.211(T)

  Distribution of Read Lengths (Bin size = 1)

        Bin:      Count  % Reads  % Bases     Average
    169,514:          1      0.0      0.0      169514
    169,513:          0      0.0      0.0      169514

"""
def test_get_reverse_sorted_readlength_counts_from_DBstats():
    got = f.get_reverse_sorted_readlength_counts_from_DBstats(capture)
    expected = [(169514, 1), (169513, 0)]
    assert_equal(expected, got)

def test_num2int():
    assert_equal(5, f.num2int('5'))
    assert_equal(1000, f.num2int('1,000'))
    assert_equal(1000000, f.num2int('1,000,000'))

partial_capture = """
        Bin:      Count  % Reads  % Bases     Average
          4:          2      0.0      0.0      xxx
          3:          0      0.0      0.0      xxx
          2:          3      0.0      0.0      xxx
          1:          8      0.0      0.0      xxx
"""
def test_calc_cutoff():
    target = 14
    expected = 2
    got = f.calc_cutoff(target, partial_capture)
    eq_(expected, got)

sample_DBdump_output = """+ R 2
+ M 0
+ H 400
@ H 8
R 0
H 8 m000_000
L 10 1899 3899
R 1
H 8 m000_000
L 11 2080 4080
R 2
H 8 m000_000
L 12 0 2500
"""
def test_parsed_readlengths_from_dbdump_output():
    lengths = list(f.parsed_readlengths_from_dbdump_output(sample_DBdump_output))
    helpers.equal_list(lengths, [2000, 2000, 2500])
def test_mapped_readlengths_from_dbdump_output():
    mapped = f.mapped_readlengths_from_dbdump_output(sample_DBdump_output)
    helpers.equal_dict(mapped, {0: 2000, 1: 2000, 2: 2500})

def test_average_difference():
    dictA = {1: 50, 2:60}
    dictB = {1: 55, 2:65, 3: 70}
    avg = f.average_difference(dictA, dictB)
    eq_(avg, -5.0)

    dictB = {1: 55}
    assert_raises(Exception, f.average_difference, dictA, dictB)

sample_perl_output = """
0 1900
1 1950
"""
def test_calc_metric_truncation():
    # and prove that dbdump can have extra reads, ignored.
    trunc = f.calc_metric_truncation(sample_DBdump_output, sample_perl_output)
    eq_(trunc, 75.0)

sample_perl_counts_output = """
100 1
200 2
100 5
"""
def test_calc_metric_fragmentation():
    frag = f.calc_metric_fragmentation(sample_perl_counts_output)
    eq_(frag, 2.5)
