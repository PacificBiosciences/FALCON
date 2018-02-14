from __future__ import unicode_literals
import time
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
