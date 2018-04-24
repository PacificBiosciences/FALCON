from __future__ import absolute_import
from __future__ import division

from pypeflow.io import (
        syscall, capture, cd,
        mkdirs, symlink, rm, touch, filesize, exists_and_not_empty) # needed?
import io
import logging
import os
import pprint
import sys

if sys.version_info >= (3, 0):
    NativeIO = io.StringIO
else:
    NativeIO = io.BytesIO

LOG = logging.getLogger()


def log(*msgs):
    LOG.debug(' '.join(repr(m) for m in msgs))


def validate_config(config):
    # This simple and quick check catches common problems early.
    # This code might go somewhere else someday.
    smrt_bin_cmds = [
        'blasr', 'samtools', 'pbalign', 'variantCaller',
    ]
    path_cmds = [
        'nucmer',
        'show-coords',
        'fc_rr_hctg_track2.exe',
    ]
    LOG.info('PATH={}'.format(os.environ['PATH']))
    try:
        capture('which which')
    except Exception:
        LOG.warning('Could not find "which" command. Skipping checks for "blasr", etc.')
        return
    for cmd in smrt_bin_cmds + path_cmds:
        syscall('which ' + cmd)
    syscall('show-coords -h')
    syscall('nucmer --version')


def update_env_from_config(config, fn):
    LOG.info('From {!r}, config={}'.format(fn, pprint.pformat(config)))
    smrt_bin = config.get('smrt_bin')
    if smrt_bin:
        PATH = '{}:{}'.format(os.environ['PATH'], smrt_bin)
        os.environ['PATH'] = PATH
    validate_config(config)


def eng(number):
    return '{:.1f}MB'.format(number / 2**20)


def read_as_msgpack(stream):
    import msgpack
    content = stream.read()
    log('  Read {} as msgpack'.format(eng(len(content))))
    return msgpack.loads(content)


def read_as_json(stream):
    import json
    content = stream.read()
    log('  Read {} as json'.format(eng(len(content))))
    return json.loads(content)


def write_as_msgpack(stream, val):
    import msgpack
    content = msgpack.dumps(val)
    log('  Serialized to {} as msgpack'.format(eng(len(content))))
    stream.write(content)


def write_as_json(stream, val):
    import json
    content = json.dumps(val, indent=2, separators=(',', ': '))
    log('  Serialized to {} as json'.format(eng(len(content))))
    stream.write(content)


def deserialize(fn):
    log('Deserializing from {!r}'.format(fn))
    with open(fn) as ifs:
        log('  Opened for read: {!r}'.format(fn))
        if fn.endswith('.msgpack'):
            val = read_as_msgpack(ifs)
        elif fn.endswith('.json'):
            val = read_as_json(ifs)
        else:
            raise Exception('Unknown extension for {!r}'.format(fn))
    log('  Deserialized {} records'.format(len(val)))
    return val


def serialize(fn, val):
    """Assume dirname exists.
    """
    log('Serializing {} records'.format(len(val)))
    mkdirs(os.path.dirname(fn))
    with open(fn, 'w') as ofs:
        log('  Opened for write: {!r}'.format(fn))
        if fn.endswith('.msgpack'):
            write_as_msgpack(ofs, val)
        elif fn.endswith('.json'):
            write_as_json(ofs, val)
            ofs.write('\n') # for vim
        else:
            raise Exception('Unknown extension for {!r}'.format(fn))


def yield_bam_fn(input_bam_fofn_fn):
    LOG.info('Reading BAM names from FOFN {!r}'.format(input_bam_fofn_fn))
    fofn_basedir = os.path.normpath(os.path.dirname(input_bam_fofn_fn))

    def abs_fn(maybe_rel_fn):
        if os.path.isabs(maybe_rel_fn):
            return maybe_rel_fn
        else:
            return os.path.join(fofn_basedir, maybe_rel_fn)
    for row in open(input_bam_fofn_fn):
        yield abs_fn(row.strip())


def yield_abspath_from_fofn(fofn_fn):
    """Yield each filename.
    Relative paths are resolved from the FOFN directory.
    """
    try:
        basedir = os.path.dirname(fofn_fn)
        for line in open(fofn_fn):
            fn = line.strip()
            if not os.path.isabs(fn):
                fn = os.path.abspath(os.path.join(basedir, fn))
            yield fn
    except Exception:
        LOG.error('Problem resolving paths in FOFN {!r}'.format(fofn_fn))
        raise


def rmdirs(*dirnames):
    for d in dirnames:
        assert os.path.normpath(d.strip()) not in ['.', '', '/']
    syscall('rm -rf {}'.format(' '.join(dirnames)))

def rmdir(d):
    rmdirs(d)

def rm_force(*fns):
    for fn in fns:
        if os.path.exists(fn):
            os.unlink(fn)
