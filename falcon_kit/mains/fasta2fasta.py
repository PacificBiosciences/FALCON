#!/usr/bin/env python2.7
"""A pre-processor for DAZZ_DB/fasta2DB.

Since fasta2DB has several constraints
(a single movie per fasta, limited line-width, actual filenames),
we write intermediate fasta
files to disk. To reduce disk I/O, we can also compress them.

Currently, we ignore zmw numbers and instead use a global counter.

Inputs may be compressed, and may be either fasta or fastq.
(For now, we ignore QVs.)
"""
from __future__ import absolute_import


from future.utils import itervalues
from builtins import object
from ..util.system import abs_fns
import argparse
import glob
import gzip
import logging
import os
import re
import sys

log = logging.getLogger()

DNA_BASES = ['A', 'C', 'G', 'T']
COMPLEMENT = {
    'A': 'T',
    'C': 'G',
    'G': 'C',
    'T': 'A',
}


def complement(x): return (COMPLEMENT[base] for base in x)


zmw_counter = None


def WriteSplit(write, seq, split=8000):
    i = 0
    while i < len(seq):
        slice = seq[i:i + split]
        write(slice)
        write('\n')
        i += split


def parse_header(header, zmw_counter=None):
    """
    >>> parse_header('>mine foo bar', 1)
    ('mine', 1, 'foo bar', 2)
    >>> parse_header('>mine/123/5_75 foo bar')
    ('mine', 123, '5_75 foo bar', None)

    For now, ignore the zmw and instead use a global counter.
    """
    if '/' in header:
        parts = header[1:].split('/')
    else:
        parts = header[1:].split(None, 1)
    movie = parts[0]
    if zmw_counter is None:
        zmw = int(parts[1])
    else:
        zmw = zmw_counter
        zmw_counter += 1
    if len(parts) > 1:
        extra = parts[-1]
    else:
        extra = ''
    return movie, zmw, extra, zmw_counter


re_range = re.compile('^(\d+)_(\d+)\s*(.*)$')


def process_fasta(ifs, movie2write):
    header = ifs.readline().strip()
    if header[0] != '>':
        raise Exception('{!r} is not a fasta file.'.format(ifs.name))
    while header:
        global zmw_counter
        movie, zmw, extra, zmw_counter = parse_header(header, zmw_counter)
        write = movie2write[movie]
        # log.info('header={!r}'.format(header))
        seq = ''
        line = ifs.readline().strip()
        while line and not line.startswith('>'):
            seq += line.strip()
            line = ifs.readline().strip()
        length = len(seq)
        #log.info('seq:{!r}...({})'.format(seq[:5], length))
        beg, end = 0, length
        mo = re_range.search(extra)
        if mo:
            beg, end, extra = mo.groups()
            beg = int(beg)
            end = int(end)
            if (end - beg) != length:
                end = beg + length
                # Probably never happens tho.
        if extra:
            extra = ' ' + extra
        new_header = '>{movie}/{zmw}/{beg}_{end}{extra}\n'.format(**locals())
        write(new_header)
        WriteSplit(write, seq)
        header = line


def process_fastq(ifs, movie2write):
    header = ifs.readline().strip()
    if header[0] != '@':
        raise Exception('{!r} is not a fastq file.'.format(ifs.name))
    while header:
        global zmw_counter
        movie, zmw, extra, zmw_counter = parse_header(header, zmw_counter)
        write = movie2write[movie]
        # log.info('header={!r}'.format(header))
        seq = ifs.readline().strip()
        header2 = ifs.readline().strip()
        quals = ifs.readline().strip()
        length = len(seq)
        #log.info('seq:{!r}...({})'.format(seq[:5], length))
        new_header = '>{movie}/{zmw}/0_{length} {extra}\n'.format(**locals())
        write(new_header)
        WriteSplit(write, seq)
        header = ifs.readline().strip()


def process_try_both(ifs, movie2write):
    try:
        process_fasta(ifs, movie2write)
    except Exception:
        log.exception('bad fasta: {!r}; trying as fastq...'.format(ifs.name))
        process_fastq(ifs, movie2write)


def process(ifn, movie2write):
    root, ext = os.path.splitext(ifn)
    if ifn.endswith('.gz'):
        Open = gzip.GzipFile
        ext = os.path.splitext(root)[1]
    elif ifn.endswith('.bz2'):
        import bz2
        Open = bz2.BZ2File
        ext = os.path.splitext(root)[1]
    else:
        Open = open

    log.info('ext={!r}'.format(ext))
    if ext in ('.fasta', '.fa'):
        func = process_fasta
    elif ext in ('.fastq', '.fq'):
        func = process_fastq
    else:
        func = process_try_both

    with Open(ifn) as ifs:
        func(ifs, movie2write)


class WriterMap(object):
    def basenames(self):
        return list(self.__obn2movie.keys())

    def close(self):
        for ofs in itervalues(self.__movie2ofs):
            ofs.close()

    def __getitem__(self, movie):
        """Get or create a 'write' function.
        """
        ofs = self.__movie2ofs.get(movie)
        if ofs is None:
            obn = self.__basename(movie)
            self.__obn2movie[obn] = movie
            if os.path.exists(obn):
                log.info('Over-writing {!r}'.format(obn))
            else:
                log.info('Creating {!r}'.format(obn))
            ofs = self.__open(obn, mode='w')
            self.__movie2ofs[movie] = ofs
        return ofs.write

    def __init__(self, Basename, Open):
        self.__obn2movie = dict()
        self.__movie2ofs = dict()
        self.__basename = Basename
        self.__open = Open


def get_writer(Gzip=False):
    if Gzip:
        def Basename(movie): return movie + '.fasta.gz'
        import functools
        Open = functools.partial(gzip.GzipFile, compresslevel=1)
        # A little better, a little slower:
        #import bz2
        #Open = bz2.BZ2File
        #Basename = lambda movie: movie + '.fasta.bz2'
    else:
        def Basename(movie): return movie + '.fasta'
        Open = open
    movie2write = WriterMap(Basename, Open)
    return movie2write


def fixall(ifns, Gzip=False):
    """Given an iterator of input absolute filenames (fasta or fastq),
    return a list of output basenames of resulting .fasta(.gz) files, relative to CWD.
    """
    if Gzip:
        open = gzip.GzipFile
    movie2write = get_writer(Gzip)
    for ifn in ifns:
        process(ifn, movie2write)
    movie2write.close()
    return movie2write.basenames()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--gzip',
                        action='store_true',
                        help='Compress intermediate fasta with gzip. (Not currently implemented.)')
    parser.add_argument('--zmw-start',
                        type=int,
                        help='Ignore the zmw number in the fasta header. Instead, use a global counter, starting at this numer.')
    # parser.add_argument('--clean',
    #    action='store_true',
    #    help='Remove intermediate fasta when done.')
    # parser.add_argument('--fofn',
    #    help='Dump intermediate FOFN. This can be used directly by "fasta2DB foo -ffofn" if fasta are uncompressed.')
    # parser.add_argument('--fasta2DB',
    #    help='Pass these arguments along to fasta2DB. These should exclude fasta inputs.')
    #global ARGS
    ARGS = parser.parse_args()
    global zmw_counter
    zmw_counter = ARGS.zmw_start
    for obn in fixall(abs_fns(sys.stdin, os.getcwd()), Gzip=ARGS.gzip):
        sys.stdout.write('{}\n'.format(os.path.abspath(obn)))


if __name__ == "__main__":
    logging.basicConfig()
    log.setLevel(logging.DEBUG)
    main()
