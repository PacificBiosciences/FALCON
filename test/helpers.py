"""
The equal_*() funcs are not really needed with pytest,
but they do not hurt.
"""
from builtins import bytes
from falcon_kit.io import NativeIO as StringIO
import difflib
import os.path


def equal_list(a, b):
    assert set(a) ^ set(b) == set()


def equal_dict(a, b):
    equal_list(sorted(a.keys()), sorted(b.keys()))
    for k in list(a.keys()):
        assert \
            a[k] == b[k], 'Inequal at k={!r} ({!r}!={!r})'.format(k, a[k], b[k])


def equal_multiline(a, b):
    alines = a.splitlines()
    blines = b.splitlines()
    equal_list(alines, blines)


def get_test_data_dir():
    return os.path.join(os.path.dirname(__file__), '..', 'test_data')


def assert_filecmp(got, expected_path):
    #result = [line.strip() for line in StringIO(got)] # fails on unicode
    #result = [line.strip() for line in io.TextIOWrapper(io.StringIO(got))] # for unicode
    result = [line.strip() for line in StringIO(bytes(got, 'utf-8'))]
    expected = [line.strip() for line in open(expected_path)]
    diffs = list(difflib.context_diff(result, expected, fromfile='got', tofile='expected'))
    if diffs:
        assert False, 'context_diff:\n' + '\n'.join(diffs[:30])
