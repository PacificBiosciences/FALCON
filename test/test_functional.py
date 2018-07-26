import helpers
import pytest
import falcon_kit.functional as f
import collections
import os
import re

thisdir = os.path.dirname(os.path.abspath(__file__))
example_HPCdaligner_fn = os.path.join(thisdir, 'HPCdaligner_synth0_new.sh')
example_HPCdaligner_small_fn = os.path.join(
    thisdir, 'HPCdaligner_synth0_preads.sh')


def eq_(a, b):
    assert a == b


def test_get_daligner_job_descriptions():
    example_HPCdaligner = open(example_HPCdaligner_fn)
    result = f.get_daligner_job_descriptions(
        example_HPCdaligner, 'raw_reads')
    assert result
    import sys, pprint
    sys.stderr.write(pprint.pformat(result))
    helpers.equal_multiline(result[('.1', '.1')], 'daligner -v -w1 -h1 -t50 -H2000 -e0.99 -l1 -s1000 -P=. -mdust raw_reads.1 raw_reads.1\nLAcheck -v raw_reads *.las\n')
    helpers.equal_multiline(result[('.2', '.1', '.2')], 'daligner -v -w1 -h1 -t50 -H2000 -e0.99 -l1 -s1000 -P=. -mdust raw_reads.2 raw_reads.1 raw_reads.2\nLAcheck -v raw_reads *.las\n')
    eq_(len(result), 2)


def test_get_daligner_job_descriptions_with_bad_arg():
    with pytest.raises(AssertionError) as excinfo:
        f.get_daligner_job_descriptions(
            'fake_filename.txt', 'raw_reads')
    assert r"f\na\nk\ne" in str(excinfo.value)


def test_get_daligner_job_descriptions_small():
    # when there is only 1 block, a special case
    example_HPCdaligner = open(example_HPCdaligner_small_fn)
    result = f.get_daligner_job_descriptions(
        example_HPCdaligner, 'preads', single=True)
    assert result
    helpers.equal_multiline(result[('.1', '.1')], 'daligner -v -h1 -t50 -H1 -e0.99 -l1 -s1000 preads.1 preads.1\nLAcheck -v preads *.las\n')
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
    eq_(result[1], [
        'LAmerge -v raw_reads.1 raw_reads.1.raw_reads.1 raw_reads.1.raw_reads.2',
        'rm raw_reads.1.raw_reads.1.las raw_reads.1.raw_reads.2.las'])
    eq_(result[2], [
        'LAmerge -v raw_reads.2 raw_reads.2.raw_reads.1 raw_reads.2.raw_reads.2',
        'rm raw_reads.2.raw_reads.1.las raw_reads.2.raw_reads.2.las'])


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
    with pytest.raises(Exception):
        f.first_block_las(line)
    line = 'LAsort foo.foo.C0'
    with pytest.raises(Exception):
        f.first_block_las(line)


def test_xform_script_for_preads():
    # Technically, we never have more than one daligner in a script, but that
    # could change in pbsmrtpipe, since it limits the number of chunks.
    script = 'daligner x y\nLAsort a b\ndaligner x1 y1\n'
    expected = 'daligner_p x y\nLAsort a b\ndaligner_p x1 y1\n'
    result = f.xform_script_for_preads(script)
    eq_(result, expected)

    script = 'daligner x y\nLAsort a b\ndaligner x1 y1\n'
    expected = script  # no-op
    result = f.xform_script_for_raw_reads(script)
    eq_(result, expected)

    eq_(f.get_script_xformer(True), f.xform_script_for_preads)
    eq_(f.get_script_xformer(False), f.xform_script_for_raw_reads)


read_lens = [1, 2, 2, 3]
rs_rl_counts = [(3, 1), (2, 2), (1, 1)]
assert collections.Counter(read_lens) == dict(rs_rl_counts)


def check(target, expected):
    got = f.calc_cutoff_from_reverse_sorted_readlength_counts(
        rs_rl_counts, target)
    eq_(expected, got)


@pytest.mark.parametrize("n_target, n_expected", [
    (8, 1),
    (7, 2),
    (4, 2),
    (3, 3),
    (1, 3),
    (0, 3),
])
def test_calc_cutoff_from_reverse_sorted_readlength_counts(n_target, n_expected):
    check(n_target, n_expected)


def test_calc_cutoff_from_reverse_sorted_readlength_counts_raises():
    with pytest.raises(Exception):
        check(9, None)


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
    eq_(expected, got)


def test_num2int():
    eq_(5, f.num2int('5'))
    eq_(1000, f.num2int('1,000'))
    eq_(1000000, f.num2int('1,000,000'))


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


def test_calc_cutoff_bad_coverage():
    target = 23  # > 22 available
    expected_message = 'Not enough reads available for desired genome coverage (bases needed=23 > actual=22)'
    with pytest.raises(f.GenomeCoverageError) as ctx:
        f.calc_cutoff(target, partial_capture)
    assert expected_message == str(ctx.value)


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
    lengths = list(f.parsed_readlengths_from_dbdump_output(
        sample_DBdump_output))
    helpers.equal_list(lengths, [2000, 2000, 2500])


def test_mapped_readlengths_from_dbdump_output():
    mapped = f.mapped_readlengths_from_dbdump_output(sample_DBdump_output)
    helpers.equal_dict(mapped, {0: 2000, 1: 2000, 2: 2500})


def test_average_difference():
    dictA = {1: 50, 2: 60}
    dictB = {1: 55, 2: 65, 3: 70}
    avg = f.average_difference(dictA, dictB)
    eq_(avg, -5.0)

    dictB = {1: 55}
    with pytest.raises(Exception):
        f.average_difference(dictA, dictB)


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


def test_args_from_line():
    line = 'LAmerge -v preads.1 preads.1.preads.1.C0.S preads.1.preads.1.N0.S'
    expected = ['preads.1', 'preads.1.preads.1.C0.S', 'preads.1.preads.1.N0.S']
    result = list(f.yield_args_from_line(line))
    helpers.equal_list(result, expected)
    bash_lines = [line]
    las_files = [word + '.las' for word in f.yield_args_from_line(
        line) for line in bash_lines if line.startswith('LAmerge')][1:]


from falcon_kit.util import io


def test_splitlines_iter():
    for text in ['', 'a', 'a\n', 'a\nb', 'a\nb\nc', 'a\nb\nc\n', '\n', '\na', '\na\n']:
        assert list(io.splitlines_iter(text)) == list(text.splitlines())


def test_Readers():
    reader = io.CapturedProcessReaderContext('echo "hi\nthere"')
    with reader:
        assert ['hi', 'there'] == list(reader.readlines())
    reader = io.StreamedProcessReaderContext('echo "hi\nthere"')
    with reader:
        assert ['hi', 'there'] == list(reader.readlines())


from falcon_kit.mains import consensus_task
# These are now redudant with the doctests, but I guess they don't hurt anything.
# And I'm not sure whether SonarQube sees doctest results. (I think so.) ~cd

def test_get_falcon_sense_option():
    assert consensus_task.get_falcon_sense_option('', 11) == ' --n-core=11'
    assert consensus_task.get_falcon_sense_option('--n-core=24', 10) == ' --n-core=10'

def test_get_pa_dazcon_option():
    assert consensus_task.get_pa_dazcon_option('', 12) == ' -j 12'
    assert consensus_task.get_pa_dazcon_option('-j  48', 13) == ' -j 13'

def test_get_option_with_proper_nproc():
    regexp = re.compile(r'-j[^\d]*(\d+)')
    assert consensus_task.get_option_with_proper_nproc(regexp, 'foo -j 5', 'baz', nproc=7, cpu_count=6) == ('foo ', 5)
    assert consensus_task.get_option_with_proper_nproc(regexp, 'foo -j 5', 'baz', nproc=3, cpu_count=4) == ('foo ', 3)
    assert consensus_task.get_option_with_proper_nproc(regexp, 'foo -j 5', 'baz', nproc=3, cpu_count=2) == ('foo ', 2)
    assert consensus_task.get_option_with_proper_nproc(regexp, 'foo', 'baz', nproc=3, cpu_count=3) == ('foo', 3)
