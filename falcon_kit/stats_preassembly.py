#!/usr/bin/env python
""" PreAssembly Report.

See FALCON-pbsmrtpipe/pbfalcon/report_preassembly.py for XML version.
"""
# Copied from
#   http://swarm/files/depot/branches/springfield/S2.3/software/smrtanalysis/bioinformatics/tools/pbreports/pbreports/report/preassembly.py
from __future__ import absolute_import
from __future__ import division
from .FastaReader import FastaReader
import collections
import itertools
import logging
import os
import pprint

log = logging.getLogger(__name__)
__version__ = '0.1'


# Copied from pbreports/util.py
# We want to avoid a dependency on pbreports b/c it needs matplotlib.
def get_fasta_readlengths(fasta_file):
    """
    Get a sorted list of contig lengths
    :return: (tuple)
    """
    lens = []
    with FastaReader(fasta_file) as f:
        for record in f:
            lens.append(len(record.sequence))
    lens.sort()
    return lens


class FastaContainer(object):

    def __init__(self, nreads, total, file_name):
        self.nreads = nreads
        self.total = total
        self.file_name = file_name

    @staticmethod
    def from_file(file_name):
#        nreads, total = _compute_values(file_name)
        read_lens = get_fasta_readlengths(file_name)
        nreads = len(read_lens)
        total = sum(read_lens)
        return FastaContainer(nreads, total, file_name)

    def __str__(self):
        return "N {n} Total {t} File: {f}".format(n=self.nreads, t=self.total, f=self.file_name)

def _validate_file(file_name):
    if os.path.isfile(file_name):
        return os.path.abspath(file_name)
    else:
        msg = "Unable to find {f}".format(f=file_name)
        log.error(msg)
        raise IOError(msg)

def cutoff_reads(read_lens, min_read_len):
    return [rl for rl in read_lens if rl >= min_read_len]

def read_len_above(read_lens, threshold):
    subtotal = 0
    # Reverse-order calculation is faster.
    for irev, rl in enumerate(reversed(read_lens)):
        subtotal += rl
        if subtotal >= threshold:
            return rl

def percentile(read_lens, p):
    # TODO: Fix this when p=1.0
    return read_lens[int(len(read_lens)*p)]

def stats_from_sorted_readlengths(read_lens):
    nreads = len(read_lens)
    total = sum(read_lens)
    n50 = read_len_above(read_lens, int(total * 0.50))
    p95 = percentile(read_lens, 0.95)
    #alt_n50 = pbreports.util.compute_n50(read_lens)
    #log.info('our n50=%s, pbreports=%s' %(n50, alt_n50)) # Ours is more correct when median is between 2 reads.
    stats = collections.namedtuple('FastaStats', ['nreads', 'total', 'n50', 'p95'])
    return stats(nreads=nreads, total=total, n50=n50, p95=p95)

def read_lens_from_fofn(fofn_fn):
    fns = [fn.strip() for fn in open(fofn_fn) if fn.strip()]
    # get_fasta_readlengths() returns sorted, so sorting the chain is roughly linear.
    return list(sorted(itertools.chain.from_iterable(get_fasta_readlengths(fn) for fn in fns)))

def stats_dict(stats_raw_reads, stats_seed_reads, stats_corrected_reads, genome_length, length_cutoff):
    """All inputs are paths to fasta files.
    genome_length and length_cutoff can be None.
    """
    log.info('stats for raw reads:       %s' %repr(stats_raw_reads))
    log.info('stats for seed reads:      %s' %repr(stats_seed_reads))
    log.info('stats for corrected reads: %s' %repr(stats_corrected_reads))

    kwds = {}
    kwds['genome_length'] = 0 if genome_length is None else genome_length
    kwds['length_cutoff'] = 0 if length_cutoff is None else length_cutoff
    kwds['raw_reads'] = stats_raw_reads.nreads
    kwds['raw_bases'] = stats_raw_reads.total
    kwds['raw_mean'] = stats_raw_reads.total / stats_raw_reads.nreads
    kwds['raw_n50'] = stats_raw_reads.n50
    kwds['raw_p95'] = stats_raw_reads.p95
    kwds['raw_coverage'] = stats_raw_reads.total / genome_length
    kwds['seed_reads'] = stats_seed_reads.nreads
    kwds['seed_bases'] = stats_seed_reads.total
    kwds['seed_mean'] = stats_seed_reads.total / stats_seed_reads.nreads
    kwds['seed_n50'] = stats_seed_reads.n50
    kwds['seed_p95'] = stats_seed_reads.p95
    kwds['seed_coverage'] = stats_seed_reads.total / genome_length
    kwds['preassembled_reads'] = stats_corrected_reads.nreads
    kwds['preassembled_bases'] = stats_corrected_reads.total
    kwds['preassembled_mean'] = stats_corrected_reads.total / stats_corrected_reads.nreads
    kwds['preassembled_n50'] = stats_corrected_reads.n50
    kwds['preassembled_p95'] = stats_corrected_reads.p95
    kwds['preassembled_coverage'] = stats_corrected_reads.total / genome_length
    kwds['preassembled_yield'] = stats_corrected_reads.total / stats_seed_reads.total
    return kwds

def make_dict(
        i_preads_fofn_fn,
        i_raw_reads_fofn_fn,
        genome_length,
        length_cutoff,
    ):
    raw_reads = read_lens_from_fofn(i_raw_reads_fofn_fn)
    stats_raw_reads = stats_from_sorted_readlengths(raw_reads)

    seed_reads = cutoff_reads(raw_reads, length_cutoff)
    stats_seed_reads = stats_from_sorted_readlengths(seed_reads)

    preads = read_lens_from_fofn(i_preads_fofn_fn)
    stats_preads = stats_from_sorted_readlengths(preads)
    report_dict = stats_dict(
            stats_raw_reads=stats_raw_reads,
            stats_seed_reads=stats_seed_reads,
            stats_corrected_reads=stats_preads,
            genome_length=genome_length,
            length_cutoff=length_cutoff,
    )
    return report_dict
