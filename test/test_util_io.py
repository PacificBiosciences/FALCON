from __future__ import unicode_literals
import time
import pytest
from builtins import str
import falcon_kit.util.io as M


def test_io_se1331():
    """Regression test for unicode conversion slow-down.
    """
    x = ''
    cmd = 'seq 20000'
    beg = time.clock()
    reader = M.CapturedProcessReaderContext(cmd)
    with reader:
        for line in reader.readlines():
            x += line
    end = time.clock()
    assert end-beg < 1


@pytest.mark.parametrize('Context',
        [M.CapturedProcessReaderContext, M.StreamedProcessReaderContext])
def test_str_type(Context):
    cmd = 'seq 2'
    reader = Context(cmd)
    with reader:
        lines = list(reader.readlines())
        assert isinstance(lines[0], str)
        assert lines == ['1', '2']
