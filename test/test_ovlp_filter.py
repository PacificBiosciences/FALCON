from __future__ import unicode_literals
import falcon_kit.mains.ovlp_filter as mod


def assert_equal(expected, got):
    assert expected == got


def test_help():
    """Can be called 'pytest' or something, but reports
    proper help message otherwise.
    """
    try:
        mod.main(['prog', '--help'])
    except SystemExit:
        pass


def test_several():
    expected = ['000000001', '000000002', '000000017', '000000028']
    data = """\
000000000 000000001 -1807 100.00 0 181 1988 1988 0 0 1807 1989 overlap
000000000 000000002 -823 99.88 0 0 823 1988 0 1166 1989 1989 overlap
000000000 000000003 -50 99.94 0 0 50 1988 0 0 50 50 overlap
000000000 000000017 -61 98.36 0 0 61 1988 0 1928 1989 1989 overlap
000000000 000000028 -1952 79.95 0 0 1952 1988 0 37 1989 1989 overlap
000000001 000000000 -1807 100.00 0 0 1807 1989 0 181 1988 1988 overlap
000000001 000000002 -642 99.84 0 0 642 1989 0 1347 1989 1989 overlap
000000002 000000000 -823 99.88 0 1166 1989 1989 0 0 823 1988 overlap
000000002 000000001 -642 99.84 0 1347 1989 1989 0 0 642 1989 overlap
000000003 000000000 -50 99.94 0 0 50 50 0 0 50 1988 overlap
000000017 000000000 -61 98.36 0 1928 1989 1989 0 0 61 1988 overlap
000000028 000000000 -1952 79.95 0 37 1989 1989 0 0 1952 1988 overlap
"""
    readlines = data.strip().splitlines
    max_diff, max_ovlp, min_ovlp, min_len = 1000, 1000, 1, 1
    got = mod.filter_stage1(readlines, max_diff, max_ovlp, min_ovlp, min_len)
    assert_equal(expected, got)


def test_one_not_ignored():
    """This is the same as a line dropped in the earlier test.
    """
    expected = []
    data = """\
000000003 000000000 -50 99.94 0 0 50 50 0 0 50 1988 overlap
"""
    readlines = data.strip().splitlines
    max_diff, max_ovlp, min_ovlp, min_len = 1000, 1000, 1, 1
    got = mod.filter_stage1(readlines, max_diff, max_ovlp, min_ovlp, min_len)
    assert_equal(expected, got)


def test_one_line_ignored():
    """This is the same as a line kept in the earlier test.
    """
    expected = ['000000017']
    data = """\
000000017 000000000 -61 98.36 0 1928 1989 1989 0 0 61 1988 overlap
"""
    readlines = data.strip().splitlines
    max_diff, max_ovlp, min_ovlp, min_len = 1000, 1000, 1, 1
    got = mod.filter_stage1(readlines, max_diff, max_ovlp, min_ovlp, min_len)
    assert_equal(expected, got)


def test_empty():
    expected = []
    data = """\
"""
    readlines = data.strip().splitlines
    max_diff, max_ovlp, min_ovlp, min_len = 1000, 1000, 1, 1
    got = mod.filter_stage1(readlines, max_diff, max_ovlp, min_ovlp, min_len)
    assert_equal(expected, got)
