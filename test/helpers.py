"""
The equal_*() funcs are not really needed with pytest,
but they do not hurt.
"""
import os.path


def equal_list(a, b):
    assert set(a) ^ set(b) == set()


def equal_dict(a, b):
    equal_list(sorted(a.keys()), sorted(b.keys()))
    for k in a.keys():
        assert \
            a[k] == b[k], 'Inequal at k={!r} ({!r}!={!r})'.format(k, a[k], b[k])


def equal_multiline(a, b):
    alines = a.splitlines()
    blines = b.splitlines()
    equal_list(alines, blines)


def get_test_data_dir():
    return os.path.join(os.path.dirname(__file__), '..', 'test_data')
