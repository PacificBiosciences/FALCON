import redis
import sys
from pbcore.io import FastaReader


r = redis.StrictRedis(host='localhost', port=6379, db=0)

class RedisList(object):

    def __init__(self, rs):
        self._rs = rs
        self.id_ = "pid:" + str( id(self) )

    def append(self, value):
        self._rs.rpush( self.id_, value)

    def __len__(self):
        return self._rs.llen( self.id_ )

    def __getitem__(self, i):
        return self._rs.lrange( self.id_, i, i)

    def pylist(self):
        return self._rs.lrange( self.id_, 0, -1)

    def __del__(self):
        self._rs.delete(self.id_)

class RedisDict(object):

    def __init__(self, rs):
        self._rs = rs
        self.id_ = "pid:" + str( id(self) )

    def __setitem__(self, key, value):
        self._rs.hset( self.id_, key, value )

    def __getitem__(self, key):
        return self._rs.hget( self.id_, key )

    def __delitem__(self, key):
        return self._rs.hdel( self.id_, key)


    def __len__(self):
        return self._rs.hlen( self.id_ )
    
    def keys(self):
        return self._rs.hgetall( self.id_ ).keys()

    def values(self):
        return self._rs.hgetall( self.id_ ).values()

    def pydict(self):
        return self._rs.hgetall( self.id_ )

    def __del__(self):
        self._rs.delete(self.id_)

def test_list():
    x = RedisList(r)
    x.append( "1" )
    x.append( "2" )
    print len(x)
    print x.pylist()
    del x

    y = RedisDict(r)
    y["a"] = "b"
    y["b"] = 1
    print y["a"]
    del y["a"]
    print y.values()
    print y.keys()
    print y.pydict()
    del y

if __name__ == "__main__":
    test_list()
