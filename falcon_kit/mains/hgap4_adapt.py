"""Given a full HGAP4 run,

generate directories and symlinks to make it look like
a pypeflow run.

Then, fc_run/fc_unzip/fc_quiver can operate on it. In fact, fc_run
should be already satisfied.

One caveat: At the moment, parts of falcon-unzip actually write into the
falcon dirs. We should fix that. But for now, we create writable run-dirs.
We do *not* write into the HGAP4 run-dir.
"""
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

"""HGAP4 with --force-chunk-mode

pbcoretools.tasks.filterdataset-0
pbcoretools.tasks.subreadset_zmw_scatter-1
pbcoretools.tasks.bam2fasta-2
pbcoretools.tasks.bam2fasta-1
.pbcoretools.tasks.subreadset_zmw_scatter-f87695df-095d-44f0-958c-380106fc60
pbcoretools.tasks.gather_fasta-1
pbcoretools.tasks.fasta2fofn-0
falcon_ns.tasks.task_falcon_gen_config-0
falcon_ns.tasks.task_falcon_config-0
falcon_ns.tasks.task_falcon_make_fofn_abs-0
falcon_ns.tasks.task_falcon0_build_rdb-0
pbfalcon.tasks.scatter0_run_daligner_jobs-1
.pbfalcon.tasks.scatter0_run_daligner_jobs-c71da9cc-ee0e-41ca-93b3-63b5d3f14857-gathered-pipeline.chunks.json
falcon_ns.tasks.task_falcon0_run_daligner_jobs-1
pbfalcon.tasks.gather0_run_daligner_jobs-1
falcon_ns.tasks.task_falcon0_rm_las-0
falcon_ns.tasks.task_falcon0_run_merge_consensus_jobs-0
pbfalcon.tasks.scatter_run_scripts_in_json-1
falcon_ns.tasks.task_falcon0_merge-1
.pbfalcon.tasks.scatter_run_scripts_in_json-94c652c2-e483-4cc3-81e8-264e5de74da3-gathered-pipeline.chunks.json
pbcoretools.tasks.gather_txt-1
pbfalcon.tasks.scatter_run_scripts_in_json_2-1
.pbfalcon.tasks.scatter_run_scripts_in_json_2-23ee9e9f-415d-4312-90ba-5744320f4f7b-gathered-pipeline.chunks.json
falcon_ns.tasks.task_falcon0_cons-1
pbcoretools.tasks.gather_txt-2
falcon_ns.tasks.task_report_preassembly_yield-0
falcon_ns.tasks.task_falcon1_rm_las-0
falcon_ns.tasks.task_falcon1_build_pdb-0
pbfalcon.tasks.scatter1_run_daligner_jobs-1
.pbfalcon.tasks.scatter1_run_daligner_jobs-53468780-ee9d-44e9-a9cb-60676f3ae7b4-gathered-pipeline.chunks.json
pbfalcon.tasks.gather1_run_daligner_jobs-1
falcon_ns.tasks.task_falcon1_merge-0
falcon_ns.tasks.task_falcon1_run_merge_consensus_jobs-0
falcon_ns.tasks.task_falcon1_db2falcon-0
falcon_ns.tasks.task_falcon2_rm_las-0
falcon_ns.tasks.task_falcon2_run_asm-0
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


"""In job_output/tasks/ dir, without chunking:

pbcoretools.tasks.filterdataset-0
pbcoretools.tasks.bam2fasta-0                    *******
pbcoretools.tasks.fasta2fofn-0
falcon_ns.tasks.task_falcon_gen_config-0
falcon_ns.tasks.task_falcon_config-0
falcon_ns.tasks.task_falcon_make_fofn_abs-0
falcon_ns.tasks.task_falcon0_build_rdb-0
falcon_ns.tasks.task_falcon0_run_daligner_jobs-0 *******
falcon_ns.tasks.task_falcon0_run_merge_consensus_jobs-0
falcon_ns.tasks.task_falcon0_rm_las-0
falcon_ns.tasks.task_falcon0_merge-0
falcon_ns.tasks.task_falcon0_cons-0
falcon_ns.tasks.task_falcon1_rm_las-0
falcon_ns.tasks.task_falcon1_build_pdb-0
falcon_ns.tasks.task_falcon1_run_daligner_jobs-0 *******
falcon_ns.tasks.task_falcon1_merge-0
falcon_ns.tasks.task_falcon1_run_merge_consensus_jobs-0
falcon_ns.tasks.task_falcon1_db2falcon-0
falcon_ns.tasks.task_falcon2_rm_las-0
falcon_ns.tasks.task_falcon2_run_asm-0


0-rawreads/
0-rawreads/raw-fofn-abs
0-rawreads/daligner-scatter
0-rawreads/job_0000
0-rawreads/raw-gather
0-rawreads/merge-scatter
0-rawreads/m_00001
0-rawreads/merge-gather
0-rawreads/cns-scatter
0-rawreads/preads
0-rawreads/preads/cns_00001
0-rawreads/report

1-preads_ovl/
1-preads_ovl/daligner-scatter
1-preads_ovl/job_0000
1-preads_ovl/gathered-las
1-preads_ovl/merge-scatter
1-preads_ovl/m_00001
1-preads_ovl/merge-gather
1-preads_ovl/db2falcon

2-asm-falcon/
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
        """
        "input_fofn" from cfg
        """
        with mkcd('0-rawreads/raw-fofn-abs/'):
            # touch('input.fofn')
            touch_done()

    def task_build_rdb():
        rdir = abstdir('falcon_ns.tasks.task_falcon0_build_rdb-0')
        with mkcd('0-rawreads/'):
            #touch('length_cutoff', 'rdb_build_done', 'run_jobs.sh', 'raw_reads.db')
            link(rdir, 'raw_reads.db')
            link(rdir, '.raw_reads.bps')
            link(rdir, '.raw_reads.idx')
            link(rdir, '.raw_reads.dust.data')
            link(rdir, '.raw_reads.dust.anno')
            touch_done()

    def task_daligner_scatter():
        """
        """
        with mkcd('0-rawreads/daligner-scatter/'):
            data = dict()
            with open('scattered.json', 'w') as stream:
                stream.write(json.dumps(data))
            touch_done()

    def create_daligner_tasks():
        """
        0-rawreads/job_0000 ***
        """
        #create_daligner_tasks('0-rawreads', '0-rawreads/daligner-scatter/scattered.json')
    def task_daligner_gather():
        """
        """
        with mkcd('0-rawreads/raw-gather/'):
            # touch('gathered_las.txt')
            touch_done()

    def task_merge_scatter():
        """
        """
        with mkcd('0-rawreads/merge-scatter/'):
            data = dict()
            with open('scattered.json', 'w') as stream:
                stream.write(json.dumps(data))
            touch_done()

    def create_merge_tasks():
        """
        0-rawreads/m_00001 ***
        """
        #create_merge_tasks('0-rawreads', '0-rawreads/merge-scatter/scattered.json')
    def task_merge_gather():
        # falcon_unzip/rr_hctg_track.py expects files under 0-rawreads/m*/raw_reads.*.las
        dn2fn = dict()
        rdir = abstdir(
            'falcon_ns.tasks.task_falcon0_run_merge_consensus_jobs-0/')
        with cd(rdir):
            for path in glob.glob(os.path.join('m*/raw_reads.*.las')):
                dn, fn = os.path.split(path)
                dn2fn[dn] = fn
        with mkcd('0-rawreads/'):
            for dn, fn in dn2fn.items():
                with mkcd(dn):
                    link(os.path.join(rdir, dn), fn)
        with mkcd('0-rawreads/merge-gather/'):
            #touch('las.fofn', 'las.fopfn')
            touch_done()

    def task_consensus_scatter():
        """
        """
        with mkcd('0-rawreads/cns-scatter/'):
            data = dict()
            with open('scattered.json', 'w') as stream:
                stream.write(json.dumps(data))
            touch_done()

    def create_consensus_tasks():
        """
        0-rawreads/preads/cns_00001 ***
        """
        #create_consensus_tasks('0-rawreads', '0-rawreads/cns-scatter/scattered.json')
    def task_consensus_gather():
        """
        """
        with mkcd('0-rawreads/preads'):
            # touch('input_preads.fofn')
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
        rdir = abstdir('falcon_ns.tasks.task_falcon1_build_pdb-0')
        with mkcd('1-preads_ovl/'):
            #touch('pdb_build_done', 'run_jobs.sh', 'preads.db')
            link(rdir, 'preads.db')
            link(rdir, '.preads.bps')
            link(rdir, '.preads.idx')
            link(rdir, '.preads.dust.data')
            link(rdir, '.preads.dust.anno')
            touch_done()

    def task_daligner_scatter1():
        """
        """
        with mkcd('1-preads_ovl/daligner-scatter/'):
            data = dict()
            with open('scattered.json', 'w') as stream:
                stream.write(json.dumps(data))
            touch_done()

    def create_daligner_tasks1():
        """
        1-preads_ovl/job_0000 ***
        """
        #create_daligner_tasks('1-preads', '1-preads/daligner-scatter/scattered.json')
    def task_daligner_gather1():
        """
        """
        with mkcd('1-preads_ovl/gathered-las/'):
            # touch('gathered_las.txt')
            touch_done()

    def task_merge_scatter1():
        """
        """
        with mkcd('1-preads_ovl/merge-scatter/'):
            data = dict()
            with open('scattered.json', 'w') as stream:
                stream.write(json.dumps(data))
            touch_done()

    def create_merge_tasks1():
        """
        1-preads_ovl/m_00001 ***
        """
        #create_merge_tasks('1-preads_ovl', '1-preads_ovl/merge-scatter/scattered.json')
    def task_merge_gather1():
        # falcon/pr_ctg_track.py expects files under 1-preads_ovl/m*/preads.*.las
        dn2fn = dict()
        rdir = abstdir(
            'falcon_ns.tasks.task_falcon1_run_merge_consensus_jobs-0/')
        with cd(rdir):
            for path in glob.glob(os.path.join('m*/preads.*.las')):
                dn, fn = os.path.split(path)
                dn2fn[dn] = fn
        with mkcd('1-preads_ovl/'):
            for dn, fn in dn2fn.items():
                with mkcd(dn):
                    link(os.path.join(rdir, dn), fn)
        #rdir = abstdir('falcon_ns.tasks.task_falcon1_merge-0')
        rdir = abstdir(
            'falcon_ns.tasks.task_falcon1_run_merge_consensus_jobs-0')
        with mkcd('1-preads_ovl/merge-gather/'):
            #touch('las.fofn', 'las.fopfn')
            link(rdir, 'file.fofn', 'las.fofn')
            touch_done()

    def task_run_db2falcon():
        rdir = abstdir('falcon_ns.tasks.task_falcon1_db2falcon-0')
        rdirx = abstdir(
            'falcon_ns.tasks.task_falcon1_run_merge_consensus_jobs-0')
        with mkcd('1-preads_ovl/db2falcon/'):
            #touch('db2falcon', 'preads4falcon.fasta')
            link(rdirx, 'preads4falcon.fasta')
            touch_done()

    def task_run_falcon_asm():
        """falcon_ns.tasks.task_falcon2_run_asm-0
          (falcon_ns.tasks.task_falcon2_rm_las-0)
        """
        rdir = abstdir('falcon_ns.tasks.task_falcon2_run_asm-0')
        with mkcd('2-asm-falcon/'):
            # workflow depends on:
            touch('falcon_asm_done')
            # get_read_ctg_map needs:
            link(rdir, 'sg_edges_list')
            link(rdir, 'utg_data')
            link(rdir, 'ctg_paths')
            # fetch_reads needs:
            link(rdir, 'p_ctg.fa')

            touch_done()

    task_make_fofn_abs_raw()
    task_build_rdb()
    task_daligner_scatter()
    create_daligner_tasks()
    task_daligner_gather()
    task_merge_scatter()
    create_merge_tasks()
    task_merge_gather()
    task_consensus_scatter()
    create_consensus_tasks()
    task_consensus_gather()
    task_report_pre_assembly()
    task_build_pdb()
    task_daligner_scatter1()
    create_daligner_tasks1()
    task_daligner_gather1()
    task_merge_scatter1()
    create_merge_tasks1()
    task_merge_gather1()
    task_run_db2falcon()
    task_run_falcon_asm()


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
