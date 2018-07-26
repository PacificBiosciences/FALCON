"""Job pools for multiprocessing.
"""
from __future__ import absolute_import

from builtins import map
from builtins import object
import multiprocessing


class FakePool(object):
    """Fake version of multiprocessing.Pool
    """

    def map(self, func, iterable, chunksize=None):
        return list(map(func, iterable))

    def imap(self, func, iterable, chunksize=None):
        return map(func, iterable)

    def terminate(self):
        pass

    def __init__(self, initializer=None, initargs=[], *args, **kwds):
        if initializer:
            initializer(*initargs)


def Pool(processes, *args, **kwds):
    """Pool factory.
    If 'not processes', return our FakePool;
    otherwise, a multiprocessing.Pool.
    """
    if processes:
        return multiprocessing.Pool(processes, *args, **kwds)
    else:
        return FakePool(*args, **kwds)
