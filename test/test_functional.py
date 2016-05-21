from nose.tools import assert_equal, assert_raises, eq_
import falcon_kit.functional as f
import StringIO
import collections
import os

thisdir = os.path.dirname(os.path.abspath(__file__))
example_HPCdaligner_fn = os.path.join(thisdir, 'HPCdaligner_synth0.sh')

def test_get_daligner_job_descriptions():
    example_HPCdaligner = open(example_HPCdaligner_fn)
    result = f.get_daligner_job_descriptions(
            example_HPCdaligner, 'raw_reads')
    assert result
    eq_(result[('.1', '.1')], "daligner -v -h1 -t16 -H1 -e0.7 -l1 -s1000 raw_reads.1 raw_reads.1\nLAsort -v raw_reads.1.raw_reads.1.C0 raw_reads.1.raw_reads.1.N0 && LAmerge -v L1.1.1 raw_reads.1.raw_reads.1.C0.S raw_reads.1.raw_reads.1.N0.S && rm raw_reads.1.raw_reads.1.C0.S.las raw_reads.1.raw_reads.1.N0.S.las\nLAcheck raw_reads.db L1.1.1.las || exit 1\n")
    eq_(result[('.2', '.1', '.2')], "daligner -v -h1 -t16 -H1 -e0.7 -l1 -s1000 raw_reads.2 raw_reads.1 raw_reads.2\nLAsort -v raw_reads.1.raw_reads.2.C0 raw_reads.1.raw_reads.2.N0 && LAmerge -v L1.1.2 raw_reads.1.raw_reads.2.C0.S raw_reads.1.raw_reads.2.N0.S && rm raw_reads.1.raw_reads.2.C0.S.las raw_reads.1.raw_reads.2.N0.S.las\nLAsort -v raw_reads.2.raw_reads.1.C0 raw_reads.2.raw_reads.1.N0 && LAmerge -v L1.2.1 raw_reads.2.raw_reads.1.C0.S raw_reads.2.raw_reads.1.N0.S && rm raw_reads.2.raw_reads.1.C0.S.las raw_reads.2.raw_reads.1.N0.S.las\nLAsort -v raw_reads.2.raw_reads.2.C0 raw_reads.2.raw_reads.2.N0 && LAmerge -v L1.2.2 raw_reads.2.raw_reads.2.C0.S raw_reads.2.raw_reads.2.N0.S && rm raw_reads.2.raw_reads.2.C0.S.las raw_reads.2.raw_reads.2.N0.S.las\nLAcheck raw_reads.db L1.1.2.las || exit 1\nLAcheck raw_reads.db L1.2.1.las || exit 1\nLAcheck raw_reads.db L1.2.2.las || exit 1\n")
    eq_(len(result), 2)

def test_get_mjob_data():
    example_HPCdaligner = open(example_HPCdaligner_fn)
    result = f.get_mjob_data(
            example_HPCdaligner)
    assert result
    eq_(result[1], ['LAmerge -v raw_reads.1 L1.1.1 L1.1.2 && rm L1.1.1.las L1.1.2.las'])
    eq_(result[2], ['LAmerge -v raw_reads.2 L1.2.1 L1.2.2 ; rm L1.2.1.las L1.2.2.las'])

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
