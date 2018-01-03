import argparse
import collections
import logging
import os
import sys
from .. import io
from .. import bash

LOG = logging.getLogger()

def read_gathered_las(path):
    """Return dict of block->[las_paths].
    For now, these are ws separated on each line of input.
    """
    result = collections.defaultdict(list)
    with open(path) as ifs:
        for line in ifs:
            block, las_path = line.split()
            result[int(block)].append(las_path)
    # LOG.warning('path={!r}, result={}'.format(
    #    path, pprint.pformat(result)))
    return result

def run(run_jobs_fn, gathered_las_fn, db_prefix, stage, wildcards, scattered_fn):
    LOG.info('Scattering las from {!r} (based on {!r}) into {!r}.'.format(
        gathered_las_fn, run_jobs_fn, scattered_fn))
    config = dict() # not used anyway
    merge_scripts = list(bash.scripts_merge(config, db_prefix, run_jobs_fn))
    gathered_dict = read_gathered_las(gathered_las_fn)
    gathered_dict_dir = os.path.normpath(os.path.dirname(gathered_las_fn))

    basedir = os.path.dirname(os.path.abspath(scattered_fn))
    rootdir = os.path.dirname(os.path.dirname(basedir)) # for now
    jobs = list()
    for p_id, merge_script, merged_las_fn in merge_scripts:
        job_id = 'm_%05d' %p_id
        # Write the scattered inputs.
        merge_script_fn = '{rootdir}/{stage}/merge-scripts/{job_id}/merge-script.sh'.format(**locals())
        las_paths_fn = '{rootdir}/{stage}/merge-scripts/{job_id}/las_paths.json'.format(**locals())
        merged_las_json_fn = '{rootdir}/{stage}/merge-scripts/{job_id}/merged_las.json'.format(**locals())
        io.mkdirs(os.path.dirname(merge_script_fn), os.path.dirname(las_paths_fn), os.path.dirname(merged_las_json_fn))
        with open(merge_script_fn, 'w') as stream:
            stream.write(merge_script)
        # las_paths must be updated for new relative paths
        las_paths_dir = os.path.normpath(os.path.dirname(las_paths_fn))
        las_paths = list()
        for p in gathered_dict[p_id]:
            updated_relative_path = os.path.relpath(os.path.join(gathered_dict_dir, p), las_paths_dir)
            las_paths.append(updated_relative_path)
        io.serialize(las_paths_fn, las_paths)
        io.serialize(merged_las_json_fn, merged_las_fn)

        # Record in a job-dict.
        job = dict()
        job['input'] = dict(
                las_paths = las_paths_fn,
                merge_script = merge_script_fn,
                merged_las_json = merged_las_json_fn,
        )
        job['output'] = dict(
                merged_las = '{rootdir}/{stage}/{job_id}/merged.las'.format(**locals()),
                job_done = '{rootdir}/{stage}/{job_id}/merge.done'.format(**locals()),
        )
        job['params'] = dict(
                actual_merged_las = merged_las_fn,
        )
        job['wildcards'] = {wildcards: job_id}
        jobs.append(job)

    io.serialize(scattered_fn, jobs)


class HelpF(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass


def parse_args(argv):
    description = 'Prepare for parallel las-merge jobs.'
    epilog = ''
    parser = argparse.ArgumentParser(
        description=description,
        epilog=epilog,
        formatter_class=HelpF,
    )
    parser.add_argument(
        '--run-jobs-fn',
        help='Input. Result of HPC.daligner.',
    )
    parser.add_argument(
        '--gathered-las-fn',
        help='Input. (Not sure of content yet.)',
    )
    parser.add_argument(
        '--db-prefix', default='raw_reads',
        help='Either preads or raw_reads.',
    )
    parser.add_argument(
        '--stage', default='0-rawreads',
        help='Either 0-rawreads or 1-preads_ovl, for now.',
    )
    parser.add_argument(
        '--wildcards',
        help='To be used in substitutions',
    )
    # Do we need config too?
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
