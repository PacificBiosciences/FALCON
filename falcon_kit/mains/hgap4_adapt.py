"""Given a full HGAP4 run,

generate directories and symlinks to make it look like
a pypeflow run.

Then, fc_run/fc_unzip/fc_quiver can operate on it. In fact, fc_run
should be already satisfied.

One caveat: At the moment, parts of falcon-unzip actually write into the
falcon dirs. We should fix that. But for now, we create writable run-dirs.
We do *not* write into the HGAP4 run-dir.
"""
from __future__ import absolute_import


from future.utils import viewitems
from ..util.system import (cd, touch, make_dirs)
import argparse
import contextlib
import glob
import json
import logging
import os
import sys

LOG = logging.getLogger(__name__)

"""
Note that, even though this program is designed to let unzip run,
it has nothing to do with unzip. It merely mimics falcon, so that
the falcon jobs appear as fully satisfied to pypeflow. That is why
this program is in this repo rather than in FALCON_unzip.

However, if HGAP4/pbsmrtpipe-tasks change, then this would need to
be updated.
"""

"""Post-FALCON steps:

pbcoretools.tasks.fasta2referenceset-0
pbalign.tasks.pbalign-0                 *******
pbreports.tasks.summarize_coverage-0
genomic_consensus.tasks.variantcaller-0 *******
pbreports.tasks.coverage_report_hgap-0
genomic_consensus.tasks.gff2bed-0
pbcoretools.tasks.contigset2fasta-0
pbreports.tasks.polished_assembly-0
pbreports.tasks.mapping_stats_hgap-0
falcon_ns.tasks.task_report_preassembly_yield-0

Or with chunking:

pbcoretools.tasks.fasta2referenceset-0
pbcoretools.tasks.subreadset_align_scatter-1
pbalign.tasks.pbalign-2
pbalign.tasks.pbalign-1
.pbcoretools.tasks.subreadset_align_scatter-b473df0f-c3d5-46ab-8c45-be5054ea0dbd-gathered-pipeline.chunks.json
pbcoretools.tasks.gather_alignmentset-1
pbreports.tasks.summarize_coverage-0
pbcoretools.tasks.alignment_contig_scatter-1
pbreports.tasks.coverage_report_hgap-0
genomic_consensus.tasks.variantcaller-2
genomic_consensus.tasks.variantcaller-1
.pbcoretools.tasks.alignment_contig_scatter-7eda161b-3ed9-4891-97f3-300dd975407a-gathered-pipeline.chunks.json
pbcoretools.tasks.gather_gff-1
pbreports.tasks.mapping_stats_hgap-0
pbcoretools.tasks.gather_fastq-1
pbcoretools.tasks.gather_contigset-1
pbcoretools.tasks.gather_vcf-1
genomic_consensus.tasks.gff2bed-0
pbreports.tasks.polished_assembly-0
pbcoretools.tasks.contigset2fasta-0
"""

@contextlib.contextmanager
def mkcd(newdir):
    make_dirs(newdir)
    with cd(newdir):
        yield


def symlink(jo):
    """Caller should first cd into link-dir.
    """
    def assert_exists(path):
        assert os.path.exists(path), 'File does not exist: {!r}'.format(path)

    def assert_dir(path):
        assert os.path.isdir(path), 'Not a directory: {!r}'.format(path)
    assert_dir(jo)
    taskdir = os.path.join(jo, 'tasks')
    assert_dir(taskdir)

    def touch_done():
        """Standard pypeflow convention.
        """
        touch('run.sh.done')

    def abstdir(basetdir):
        return os.path.abspath(os.path.join(jo, 'tasks', basetdir))

    def link(targetdir, basename, linkname=None):
        if not linkname:
            linkname = basename
        reldir = os.path.relpath(targetdir)
        target = os.path.join(reldir, basename)
        assert_exists(os.path.abspath(target))
        if os.path.lexists(linkname):
            if os.readlink(linkname) == target:
                return
            os.unlink(linkname)
        LOG.info('link {!r} to {}/{}'.format(linkname, reldir, basename))
        os.symlink(target, linkname)

    # Define task symlinkers

    def task_make_fofn_abs_raw():
        """deprecated
        "input_fofn" from cfg
        """
        rdir = abstdir('falcon_ns2.tasks.task_falcon_make_fofn_abs-0')
        with mkcd('0-rawreads/raw-fofn-abs/'):
            link(rdir, 'file.fofn', 'input.fofn')
            # touch('input.fofn')
            touch_done()

    def task_build_rdb():
        rdir = abstdir('falcon_ns2.tasks.task_falcon0_dazzler_build_raw-0')
        with mkcd('0-rawreads/build/'):
            #touch('length_cutoff', 'rdb_build_done', 'run_jobs.sh', 'raw_reads.db')
            link(rdir, 'raw_reads.db')
            link(rdir, '.raw_reads.bps')
            link(rdir, '.raw_reads.idx')
            link(rdir, '.raw_reads.dust.data')
            link(rdir, '.raw_reads.dust.anno')
            touch_done()

    def task_tan_split():
        rdir = abstdir('falcon_ns2.tasks.task_falcon0_dazzler_tan_split-0')
        with mkcd('0-rawreads/tan-split/'):
            #link(rdir, 'split.json', 'tan-uows.json')
            with open('tan-uows.json', 'w') as stream:
                data = dict()
                stream.write(json.dumps(data))
            link(rdir, 'bash_template.txt', 'bash_template.sh')
            touch_done()

    def task_tan_gathered():
        with mkcd('0-rawreads/tan-gathered/'):
            touch_done()

    def task_tan_combine():
        rdir = abstdir('falcon_ns2.tasks.task_falcon0_dazzler_tan_combine-0')
        with mkcd('0-rawreads/tan-combine/'):
            link(rdir, 'raw_reads.db')
            link(rdir, '.raw_reads.bps')
            link(rdir, '.raw_reads.idx')
            link(rdir, '.raw_reads.dust.data')
            link(rdir, '.raw_reads.dust.anno')
            link(rdir, '.raw_reads.tan.data')
            link(rdir, '.raw_reads.tan.anno')
            touch_done()

    def task_daligner_split():
        rdir = abstdir('falcon_ns2.tasks.task_falcon0_dazzler_daligner_split-0')
        with mkcd('0-rawreads/daligner-split/'):
            #link(rdir, 'split.json', 'all-units-of-work.json')
            with open('all-units-of-work.json', 'w') as stream:
                data = dict()
                stream.write(json.dumps(data))
            link(rdir, 'bash_template.txt', 'bash_template.sh')
            touch_done()

    def task_daligner_gathered():
        with mkcd('0-rawreads/daligner-gathered/'):
            touch_done()

    def task_daligner_combine():
        rdir = abstdir('falcon_ns2.tasks.task_falcon0_dazzler_daligner_combine-0')
        with mkcd('0-rawreads/daligner-combine/'):
            link(rdir, 'las_paths.json', 'gathered-las.json')
            touch_done()

    def task_lamerge_split():
        rdir = abstdir('falcon_ns2.tasks.task_falcon0_dazzler_lamerge_split-0')
        with mkcd('0-rawreads/las-merge-split/'):
            #link(rdir, 'split.json', 'all-units-of-work.json')
            with open('all-units-of-work.json', 'w') as stream:
                data = dict()
                stream.write(json.dumps(data))
            link(rdir, 'bash_template.txt', 'las-merge-bash-template.sh')
            touch_done()

    def task_lamerge_gathered():
        with mkcd('0-rawreads/las-merge-gathered/'):
            touch_done()

    def task_lamerge_combine():
        # falcon_unzip/rr_hctg_track.py looks at las-merge-combine/las_paths.json, with abspaths
        rdir = abstdir(
            'falcon_ns2.tasks.task_falcon0_dazzler_lamerge_combine-0')
        with mkcd('0-rawreads/las-merge-combine/'):
            link(rdir, 'las_paths.json', 'las_fofn.json') # unzip/quiver, for now
            link(rdir, 'las_paths.json')
            link(rdir, 'block2las.json')
            touch_done()

    def task_cns_split():
        rdir = abstdir('falcon_ns2.tasks.task_falcon0_dazzler_cns_split-0')
        with mkcd('0-rawreads/cns-split/'):
            #link(rdir, 'split.json', 'all-units-of-work.json')
            with open('split.json', 'w') as stream:
                data = dict()
                stream.write(json.dumps(data))
            #link(rdir, 'bash_template.txt', 'bash_template.sh')
            touch_done()

    def task_cns_gather():
        #rdir = abstdir('falcon_ns2.tasks.task_falcon0_dazzler_cns_split-0')
        #rdir = abstdir('falcon_ns2.tasks.task_falcon0_run_cns_post_gather-0')
        with mkcd('0-rawreads/cns-gather/'):
            touch_done()

    #def task_cns_combine():
    #    rdir = abstdir('falcon_ns2.tasks.task_falcon0_dazzler_cns_combine-0')
    #    with mkcd('0-rawreads/cns-combine/'):
    #        touch_done()

    def task_preads():
        rdir = abstdir('falcon_ns2.tasks.task_falcon0_run_cns_post_gather-0')
        with mkcd('0-rawreads/preads/'):
            link(rdir, 'input-preads.fofn', 'input_preads.fofn')
            touch_done()

    def task_report_pre_assembly():
        """
        """
        with mkcd('0-rawreads/report/'):
            # touch('pre_assembly_stats.json')
            touch_done()

    def task_build_pdb():
        """
        """
        rdir = abstdir('falcon_ns2.tasks.task_falcon1_build_pdb-0')
        with mkcd('1-preads_ovl/build/'):
            #touch('pdb_build_done', 'run_jobs.sh', 'preads.db')
            link(rdir, 'preads.db')
            link(rdir, '.preads.bps')
            link(rdir, '.preads.idx')
            link(rdir, '.preads.dust.data')
            link(rdir, '.preads.dust.anno')
            touch_done()

    def task_daligner_split1():
        #rdir = abstdir('falcon_ns2.tasks.task_falcon1_dazzler_daligner_split-0')
        with mkcd('1-preads_ovl/daligner-split/'):
            #link(rdir, 'split.json', 'all-units-of-work.json')
            with open('all-units-of-work.json', 'w') as stream:
                data = dict()
                stream.write(json.dumps(data))
            #link(rdir, 'bash_template.txt', 'bash_template.sh')
            touch_done()

    def task_daligner_gathered1():
        with mkcd('1-preads_ovl/daligner-gathered/'):
            touch_done()

    def task_daligner_combine1():
        rdir = abstdir('falcon_ns2.tasks.task_falcon1_run_daligner_find_las-0')
        with mkcd('1-preads_ovl/daligner-combine/'):
            #link(rdir, 'las_paths.json', 'gathered-las.json')
            link(rdir, 'gathered-las.json', 'gathered-las.json')
            touch_done()

    def task_lamerge_split1():
        #rdir = abstdir('falcon_ns2.tasks.task_falcon1_dazzler_lamerge_split-0')
        with mkcd('1-preads_ovl/las-merge-split/'):
            #link(rdir, 'split.json', 'all-units-of-work.json')
            with open('all-units-of-work.json', 'w') as stream:
                data = dict()
                stream.write(json.dumps(data))
            #link(rdir, 'bash_template.txt', 'las-merge-bash-template.sh')
            touch_done()

    def task_lamerge_gathered1():
        with mkcd('1-preads_ovl/las-merge-gathered/'):
            touch_done()

    def task_lamerge_combine1():
        rdir = abstdir(
            'falcon_ns2.tasks.task_falcon1_run_las_merge_post_gather-0')
        #    'falcon_ns2.tasks.task_falcon1_dazzler_lamerge_combine-0')
        with mkcd('1-preads_ovl/las-merge-combine/'):
            link(rdir, 'las-fofn.json', 'las_fofn.json') # unzip/quiver, for now
            link(rdir, 'las-fofn.json', 'las_paths.json')
            link(rdir, 'p_id2las.json', 'block2las.json')
            touch_done()

    def task_run_db2falcon():
        rdir = abstdir('falcon_ns2.tasks.task_falcon1_run_db2falcon-0')
        with mkcd('1-preads_ovl/db2falcon/'):
            link(rdir, 'preads4falcon.fasta')
            touch_done()

    def task_run_falcon_asm():
        rdir = abstdir('falcon_ns2.tasks.task_falcon2_run_falcon_asm-0')
        with mkcd('2-asm-falcon/'):
            # workflow depends on:
            touch('falcon_asm_done')
            # get_read_ctg_map needs:
            link(rdir, 'sg_edges_list')
            link(rdir, 'utg_data')
            link(rdir, 'ctg_paths')
            # fetch_reads needs:
            link(rdir, 'p_ctg.fa')
            link(rdir, 'a_ctg.fa')
            link(rdir, 'p_ctg_tiling_path')
            link(rdir, 'a_ctg_tiling_path')

            touch_done()

    #task_make_fofn_abs_raw()
    task_build_rdb()
    task_tan_split()
    task_tan_gathered()
    task_tan_combine()
    task_daligner_split()
    task_daligner_gathered()
    task_daligner_combine()
    task_lamerge_split()
    task_lamerge_gathered()
    task_lamerge_combine()
    task_cns_split()
    task_cns_gather()
    #task_cns_combine()
    task_preads()
    task_report_pre_assembly()
    task_build_pdb()
    task_daligner_split1()
    task_daligner_gathered1()
    task_daligner_combine1()
    task_lamerge_split1()
    task_lamerge_gathered1()
    task_lamerge_combine1()
    task_run_db2falcon()
    task_run_falcon_asm()

    def dump_fc_run(fn):
        input_fofn = os.path.join(abstdir('falcon_ns2.tasks.task_falcon_make_fofn_abs-0'), 'file.fofn')
        length_cutoff = int(open(os.path.join(abstdir('falcon_ns2.tasks.task_falcon0_dazzler_build_raw-0'), 'length_cutoff.txt')).read())
        with open(fn, 'w') as stream:
            p = lambda x: stream.write(x + '\n')
            p('[General]')
            p('input_fofn = {}'.format(input_fofn))
            p('length_cutoff = {}'.format(length_cutoff))
            p('[Unzip]')
            p('input_fofn = {}'.format(input_fofn))
            p('input_bam_fofn = {} # You need to find this!'.format('input_bam.fofn'))
            p('[job.defaults]')
            p('pwatcher_type = blocking')
            #p('submit = /bin/bash -c "${JOB_SCRIPT}"')
            p('submit = /bin/bash -c "${JOB_SCRIPT}" > "${JOB_STDOUT}" 2> "${JOB_STDERR}"')

    dump_fc_run('fc_run.generated.cfg')


def get_parser():
    class HelpF(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
        pass
    description = 'Given a full HGAP4 run, generate directories and symlinks to make it look like a pypeflow run.'
    epilog = """
Typically:
    mkdir mydir/
    cd mydir/
    python -m falcon_kit.mains.hgap4_adapt --job-output-dir=../job_output/

    fc_run fc_run.cfg          -- (A)
    fc_unzip.py fc_unzip.cfg   -- (B)
    fc_quiver.py fc_unzip.cfg  -- (C)

You need to create/modify the .cfg files.

(A) This should be a no-op, and you do not need to run this at all. Just a sanity check.
It will tell you that everything is already satisfied. But it
cannot work unless you provide `input.fofn` (which can be empty) and set it to `input_fofn`
in your .cfg.

(B)/(C) These will need both `input_fofn` and `input_bam_fofn`. The latter
should name actual BAM files to use for Quiver (also for partitioning for pbalign).

For more help on .cfg files, see
* https://github.com/PacificBiosciences/FALCON/wiki
* https://github.com/PacificBiosciences/FALCON_unzip/wiki
* https://github.com/PacificBiosciences/FALCON-integrate/wiki
"""

    parser = argparse.ArgumentParser(
        description=description,
        epilog=epilog,
        formatter_class=HelpF)
    parser.add_argument('--job-output-dir', default='job_output',
                        help='Directory of HGAP4 job_output. (A symlink or relative path is fine.) Task-dirs are under here in "tasks/"')
    return parser


def main(argv=sys.argv):
    args = get_parser().parse_args(argv[1:])
    symlink(args.job_output_dir)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    main()
