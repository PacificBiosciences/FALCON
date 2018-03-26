import pytest
import StringIO
import falcon_kit.run_support as mod


def parse_config(content):
    """Used by tests.
    (Clean this code up later.)
    """
    config = dict()
    stream = StringIO.StringIO
    cfg2 = mod.parse_cfg_with_sections(stream(content))
    mod.update_config_from_sections(config, cfg2)
    mod.update_job_sections(config)
    return config


FC_RUN_CFG = """\
[General]
job_name_style = 1
#use_tmpdir = true
#job_type = local
#job_type = sge
stop_all_jobs_on_failure = true

#skip_checks = true
#use_tmpdir = /scratch
pwatcher_type = blocking
job_type = string
job_queue = bash -C ${CMD} >| ${STDOUT_FILE} 2>| ${STDERR_FILE}
job_queue = bash -C ${CMD}
# By dropping STD*_FILE, we see all output on the console.
# That helps debugging in TravisCI/Bamboo.


# list of files of the initial bas.h5 files
input_fofn = input.fofn
#input_fofn = preads.fofn

input_type = raw
#input_type = preads

# The length cutoff used for seed reads used for initial mapping
#length_cutoff = 1
genome_size = 5000
seed_coverage = 20

# The length cutoff used for seed reads usef for pre-assembly
length_cutoff_pr = 1


#job_queue = production
sge_option_da = -pe smp 8 -q %(job_queue)s
sge_option_la = -pe smp 2 -q %(job_queue)s
sge_option_pda = -pe smp 8 -q %(job_queue)s
sge_option_pla = -pe smp 2 -q %(job_queue)s
sge_option_fc = -pe smp 24 -q %(job_queue)s
sge_option_cns = -pe smp 8 -q %(job_queue)s

default_concurrent_jobs = 64
da_concurrent_jobs = 32
la_concurrent_jobs = 32
cns_concurrent_jobs = 32
pda_concurrent_jobs = 32
pla_concurrent_jobs = 32

pa_HPCdaligner_option =   -v -B4 -t50 -h1 -e.99 -w1 -l1 -s1000
ovlp_HPCdaligner_option = -v -B4 -t50 -h1 -e.99 -l1 -s1000

#pa_DBsplit_option =   -a -x5 -s.00065536
pa_DBsplit_option =   -a -x5 -s.065536
#pa_DBsplit_option =   -a -x5 -s1
ovlp_DBsplit_option = -a -x5 -s50

LA4Falcon_preload = true
falcon_sense_option = --output_multi --min_idt 0.70 --min_cov 1 --max_n_read 20000 --n_core 0
#--min_cov_aln 1 --min_len_aln 40

overlap_filtering_setting = --max_diff 10000 --max_cov 100000 --min_cov 1 --min_len 1 --bestn 1000 --n_core 0
#dazcon = 1


[job.defaults]
submit = bash -C ${CMD} 2> ${STDERR_FILE}.${NPROC}.txt

[job.step.da]
NPROC = 256

[job.step.pda]
njobs = 42
"""

def test_update_config_from_sections():
    config = parse_config(FC_RUN_CFG)

    assert int(config['General']['default_concurrent_jobs']) == 64
    assert int(config['General']['da_concurrent_jobs']) == 32
    assert int(config['job.defaults']['njobs']) == 64
    assert int(config['job.step.da']['njobs']) == 32
    assert int(config['job.step.da']['NPROC']) == 256
    assert int(config['job.step.pda']['njobs']) == 42

def test_update_config_from_sections_foo():
    config = dict()
    cfg2 = mod.parse_cfg_with_sections(StringIO.StringIO(FC_RUN_CFG))

    import collections
    cfg2['foo'] = collections.defaultdict(list)
    with pytest.raises(Exception) as excinfo:
        mod.update_config_from_sections(config, cfg2)
    assert 'You have 1 unexpected cfg sections' in str(excinfo.value)


CFG_SANS_SUBMIT_LOCAL = """\
[General]
pwatcher_type = fs_based
job_type = local
[job.defaults]
# empty
"""
def test_update_job_sections0():
    config = parse_config(CFG_SANS_SUBMIT_LOCAL)
    assert 'submit' not in config['job.defaults']
    #assert config['job.defaults']['submit'] == 'bash -C ${CMD}'

CFG_SANS_SUBMIT_SANS_JOB_TYPE = """\
[General]
pwatcher_type = fs_based
#job_type = sge
[job.defaults]
# empty
"""
def test_update_job_sections1():
    with pytest.raises(Exception) as excinfo:
        parse_config(CFG_SANS_SUBMIT_SANS_JOB_TYPE)
    assert 'but General.job_type is not set' in str(excinfo.value)

CFG_SANS_SUBMIT_UNKNOWN_JOB_TYPE = """\
[General]
pwatcher_type = fs_based
job_type = foo
[job.defaults]
# empty
"""
def test_update_job_sections2():
    with pytest.raises(Exception) as excinfo:
        config = parse_config(CFG_SANS_SUBMIT_UNKNOWN_JOB_TYPE)
    assert 'job_type=foo not in ' in str(excinfo.value)

CFG_SANS_SUBMIT_UNKNOWN_PWATCHER_TYPE = """\
[General]
pwatcher_type = foo
job_type = sge
[job.defaults]
# empty
"""
def test_update_job_sections3():
    with pytest.raises(Exception) as excinfo:
        config = parse_config(CFG_SANS_SUBMIT_UNKNOWN_PWATCHER_TYPE)
    assert 'Unknown pwatcher_type' in str(excinfo.value)
