from nose.tools import eq_
import falcon_kit.functional as f
import StringIO
import os

thisdir = os.path.dirname(os.path.abspath(__file__))
example_HPCdaligner = open(os.path.join(thisdir, 'HPCdaligner_synth0.sh'))

def test_get_daligner_job_descriptions():
    result = f.get_daligner_job_descriptions(
            example_HPCdaligner, 'raw_reads')
    assert result
    eq_(result[(1, 1)], "daligner -v -h1 -t16 -H1 -e0.7 -l1 -s1000 raw_reads.1 raw_reads.1\nLAsort -v raw_reads.1.raw_reads.1.C0 raw_reads.1.raw_reads.1.N0 && LAmerge -v L1.1.1 raw_reads.1.raw_reads.1.C0.S raw_reads.1.raw_reads.1.N0.S && rm raw_reads.1.raw_reads.1.C0.S.las raw_reads.1.raw_reads.1.N0.S.las\n")
    eq_(result[(2, 1, 2)], "daligner -v -h1 -t16 -H1 -e0.7 -l1 -s1000 raw_reads.2 raw_reads.1 raw_reads.2\nLAsort -v raw_reads.1.raw_reads.2.C0 raw_reads.1.raw_reads.2.N0 && LAmerge -v L1.1.2 raw_reads.1.raw_reads.2.C0.S raw_reads.1.raw_reads.2.N0.S && rm raw_reads.1.raw_reads.2.C0.S.las raw_reads.1.raw_reads.2.N0.S.las\nLAsort -v raw_reads.2.raw_reads.1.C0 raw_reads.2.raw_reads.1.N0 && LAmerge -v L1.2.1 raw_reads.2.raw_reads.1.C0.S raw_reads.2.raw_reads.1.N0.S && rm raw_reads.2.raw_reads.1.C0.S.las raw_reads.2.raw_reads.1.N0.S.las\nLAsort -v raw_reads.2.raw_reads.2.C0 raw_reads.2.raw_reads.2.N0 && LAmerge -v L1.2.2 raw_reads.2.raw_reads.2.C0.S raw_reads.2.raw_reads.2.N0.S && rm raw_reads.2.raw_reads.2.C0.S.las raw_reads.2.raw_reads.2.N0.S.las\n")
    eq_(len(result), 2)
