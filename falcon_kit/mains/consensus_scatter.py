from __future__ import unicode_literals
import argparse
import collections
import logging
import os
import sys
from .. import io
from .. import bash

LOG = logging.getLogger()

def corrected_relpath(p, was_rel_to):
    if os.path.isabs(p):
        return p
    #LOG.warning('{},{},{}'.format(p, was_rel_to, os.path.relpath(os.path.join(was_rel_to, p))))
    return os.path.normpath(os.path.relpath(os.path.join(was_rel_to, p)))

def read_gathered_las(path):
    """Return dict of block->[las_paths].
    For now, these are ws separated on each line of input.
    """
    result = collections.defaultdict(list)
    dn = os.path.normpath(os.path.dirname(path))
    with open(path) as ifs:
        for line in ifs:
            block, las_path = line.split()
            result[int(block)].append(corrected_relpath(las_path, dn))
    #import pprint
    #LOG.warning('path={!r}, result={}'.format(
    #   path, pprint.pformat(result)))
    return result


def run(las_fopfn_fn, db_fn, length_cutoff_fn, config_fn, wildcards, scattered_fn):
    LOG.info('Scattering las from {!r} (based on {!r}) into {!r}.'.format(
        las_fopfn_fn, db_fn, scattered_fn))

    wildcards = wildcards.split(',')
    basedir = os.path.dirname(os.path.abspath(scattered_fn))
    rootdir = os.path.dirname(os.path.dirname(basedir)) # for now
    jobs = list()
    p_ids_merge_las = read_gathered_las(las_fopfn_fn)
    tasks = []
    for p_id, las_fns in p_ids_merge_las.iteritems():
        assert len(las_fns) == 1, repr(las_fns)
        # since we know each merge-task is for a single block
        las_fn = las_fns[0]
        cns_id = 'cns_%05d' % int(p_id)
        cns_id2 = cns_id
        ##out_done_fn = '%s_done' % cns_label
        #out_file_fn = '%s.fasta' % cns_label
        symlinked_las_fn = '{rootdir}/0-rawreads/cns-scatter/{cns_id}/merged.{cns_id2}.las'.format(**locals())
        io.mkdirs(os.path.normpath(os.path.dirname(symlinked_las_fn)))
        src = os.path.relpath(las_fn,
            os.path.normpath(os.path.dirname(symlinked_las_fn)))
        io.symlink(src, symlinked_las_fn)

        # Record in a job-dict.
        job = dict()
        job['input'] = dict(
                las = symlinked_las_fn,
                db = db_fn,
                length_cutoff = length_cutoff_fn,
                config = config_fn,
        )
        job['output'] = dict(
                fasta =
                '{rootdir}/0-rawreads/consensus/{cns_id}/consensus.{cns_id2}.fasta'.format(**locals()),
        )
        job['params'] = dict(
        )
        job['wildcards'] = {wildcards[0]: cns_id, wildcards[1]: cns_id}
        jobs.append(job)

    io.serialize(scattered_fn, jobs)


class HelpF(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass


def parse_args(argv):
    description = 'Prepare for parallel consensus jobs.'
    epilog = ''
    parser = argparse.ArgumentParser(
        description=description,
        epilog=epilog,
        formatter_class=HelpF,
    )
    parser.add_argument(
        '--las-fopfn-fn',
        help='Input. FOFN for las files.)',
    )
    parser.add_argument(
        '--db-fn',
        help='Input. Dazzler DB of raw_reads.',
    )
    parser.add_argument(
        '--length-cutoff-fn',
        help='Input. Contains a single integer, the length-cutoff.',
    )
    parser.add_argument(
        '--config-fn',
        help='Input. JSON of relevant configuration (currently from General section of full-prog config).',
    )
    parser.add_argument(
        '--wildcards',
        help='To be used in substitutions',
    )
    parser.add_argument(
        '--scattered-fn',
        help='Output. JSON list of jobs, where each is a dict of input/output/params/wildcards.',
    )
    args = parser.parse_args(argv[1:])
    return args


def main(argv=sys.argv):
    args = parse_args(argv)
    logging.basicConfig(level=logging.INFO)
    run(**vars(args))


if __name__ == '__main__':  # pragma: no cover
    main()
