
![Falcon Logo](./falcon_icon2.png =128x128)

# Falcon Genome Assembly Tool Kit Manual for v0.2 Branch

Author: Jason Chin
Date: Dec 4, 2014

## Overview of the Hierarchical Genome Assembly Process

A "Hierarchical Genome Assembly Process" is constituted of the following steps
for generating a genome assembly from a set of sequencing reads:

* Raw sub-reads overlapping for error correction

* Pre-assembly and error correction

* Overlapping detection of the error corrected reads

* Overlap filtering

* Constructing graph from overlaps

* Constructing contig from graph

Each of the step is accomplished with different command line tools implementing
different set algorithms to accomplish the work.  Also, the computation
requirement are also quite different for each steps.  The manual assumes the
user has a reasonable amount computation resources. For example, to assemble
100M size genome with a reasonable amount of time, one might need at least 32
core cpus and 128Gb RAM. The code is written with the assumption of a cluster
computation environment. One need a job queue for long last scripting job and
cpu-rich computational job queues

With a file that contains a set of sub-reads, the script `fc_run.py` can drive
the workflow managing checking the data dependency and submitting the jobs for
each step and generating a draft assembly from the giving data. 

`fc_run.py` is the workflow driving script needs to be run on a machine which
allow long last time through the period of time of the whole assembly process.
It takes a configuration file as single input.  The input files of the raw
sequence data is included in the configuration files.

The configuration file can be used for controlling the computation resource used
and important algorithm parameters for reaching optimized assembly according to
the input data set. Unfortunately, there is no magic way to guess what the right
options are as the available computational resource from place to place and the
scope of a sequencing project varies from case to case.  The best way to tune
the parameter is to understand some assembly theory, a little bit of the
implementation so one can guess the impact of changing certain option correctly.
It is also very important to do quick look at the read length distribution and
overall coverage and adjust the options accordingly.

In this manual, we will go over the hierarchical genome assembly process and the
`fc_run.py` option choice side-by-side.

## Installation

Well, this is a TODO item for now.

## Running `fc_run.py` 

It is easy to start the assembly process like this,

```
    fc_run.py fc_run.cfg

```


Here is an example of `fc_run.cfg` for a small E. coli assembly,

```

    [General]
    # list of files of the initial subread fasta files
    input_fofn = input.fofn
    #input_fofn = preads.fofn

    input_type = raw
    #input_type = preads

    # The length cutoff used for seed reads used for initial mapping
    length_cutoff = 12000

    # The length cutoff used for seed reads usef for pre-assembly
    length_cutoff_pr = 12000

    # Cluster queue setting
    sge_option_da = -pe smp 8 -q bigmem
    sge_option_la = -pe smp 2 -q bigmem
    sge_option_pda = -pe smp 8 -q bigmem 
    sge_option_pla = -pe smp 2 -q bigmem
    sge_option_fc = -pe smp 24 -q bigmem
    sge_option_cns = -pe smp 8 -q bigmem

    # concurrency settgin
    pa_concurrent_jobs = 32
    ovlp_concurrent_jobs = 32

    # overlapping options for Daligner
    pa_HPCdaligner_option =  -v -dal4 -t16 -e.70 -l1000 -s1000 
    ovlp_HPCdaligner_option = -v -dal4 -t32 -h60 -e.96 -l500 -s1000 

    pa_DBsplit_option = -x500 -s50
    ovlp_DBsplit_option = -x500 -s50

    # error correction consensus optione
    falcon_sense_option = --output_multi --min_idt 0.70 --min_cov 4 --local_match_count_threshold 2 --max_n_read 200 --n_core 6

    # overlap filtering options
    overlap_filtering_setting = --max_diff 100 --max_cov 100 --min_cov 20 --bestn 10 
```

## Raw sub-reads overlapping for error correction

In this version of the Falcon kit, the overlapping is done with a modified
version of Gene Myers' Daligner (http://dazzlerblog.wordpress.com). The forked
version can be found at https://github.com/cschin/DALIGNER . Most changes from
the original Gene Myers' code is on adapting some I/O for the downstream
bioinformatics process. You can just do a simple `git diff` to see the
difference. (Isn't open source great!?)  

The option `input_fofn` points to the file that contains all input data.
`fasta2DB` from `Daligner` is called within `fc_run.py`. (This is I/O intensive
and it will be run from the computer node where you execute `fc_run.py`. If this
is an issue in your cluster, you will have to modify the code to wrap the
related operation into a script that can be submitted in your job management
system.)

This version of `fc_run.py` supports running assembly from error corrected
reads.  If you set the option `input_type = preads` rather than `input_type =
raw`, `fc_run.py` will assume the fasta files in `input_fofn` are all
error-corrected reads and it will ignore any error correction step and go
directly into the final assembly overlapping step.

You will need to decide the length cutoff. Typically, it will be nice to chose
the threshold at the point you can get longest 15x to 20x for genome assembly.
However, if the computational resource is abundant and you might find other
applications of error corrected reads, you can set lower length cutoff to get
more error corrected reads for your applications. 

The option `length_cutoff` controls the cutoff used during the error correction
process and `length_cutoff_pr` controls the cutoff used for the later assembly
overlapping step.  In the final assembly, more reads may not lead to a better
assembly due to some of the reads can be noisy and create false links in the
assembly graph. Sometimes, it might make sense to try different
`length_cutoff_pr` as it is relative cheap for computation than the first
overlapping step for error correction. One strategy is to chose smaller
`length_cutoff` and do the computation once. Later, we can use different
`length_cutoff_pr` for getting better assembly. 

The option `pa_concurrent_jobs` controls the number of concurrent jobs that can
be submitted by `fc_run.py`.  `sge_option_da` and `sge_option_la` controls the
job queue and the number of slots of the `daligner` jobs.  The default number of
thread used by `daligner` is 4. However, depending on the cluster configuration
and the amount of memoery of the computational nodes, you might want to use more
than 4 slots.  The best to chose the right number is to consult your local HPC
gurus and do some small experiments firest.

The total numebr of jobs that is run is determined how one "splits" the sequence
database. You should read Gene Myers's blog ( http://dazzlerblog.wordpress.com )
carefully to know how to tune the option `pa_DBsplit_option` and
`pa_HPCdaligner_option`.  Generally, for large genome, you should use `-s400`
(400Mb sequence per block) in `pa_DBsplit_option`. This will make smaller number
of jobs but each job runs longer. However, if you have job queue system which
has limit of how long a job can run, it might be desirable to have smaller
number for the `-s` option.

Another parameter affects the total number of jobs is the `-dal` option in
`pa_HPCdaligner_option`. The number for the `-dal` option determines how many
blocks are compared to each in single jobs. Larger number gives larger jobs but
smaller amount of total jobs. Smaller number gives smaller jobs but you have to
submit more jobs to your cluster.

In this workflow, the trace point generated by `daligner` is not used. ( Well,
to be efficient, one should use the trace points but one have to know how to
pull them out correctly first. )  The `-s1000` in `pa_HPCdaligner_option` makes
the trace points sparse to save some disk space (not much though).  We also
ignore all reads less than 1kb by specifying `-l1000`.

## Pre-assembly and error correction

The output of `daligner` is a set of `.las` files that contains information of
the alignments between the reads.  Such information is dumped as sequences for
error correction by a binary executable `LA4Falcon` to `fc_consensus.py`. The
`fc_consensus.py` dose the work to generate consensus. (The alignments for
generating consensus are done with back-end code written in C for speed.)

The `fc_consensus.py` has many options. You can use the `falcon_sense_option` to
control it.  In most of case, the `--min_cov` and `--max_n_read` are the most
important options.  `--min_cov` controls when a seed read getting trimmed or
broken due to low coverage.  `--max_n_read` put a cap on the number of reads
used for error correction. In high repetitive genome, you will need to put
smaller `--max_n_read` to make sure the consensus code does not waste time
aligning repeats. The longest proper overlaps are used for correction to reduce the
probability of collapsed repeats.

## Overlapping detection of the error corrected reads

This part is pretty much the same as the first overlapping stage, although some
"hacks" are necessary as `daligner` only take native raw reads as default.  
`fc_run.py` generates a fasta file of error corrects where the fasta header is
parse-able by `daligner`.  The following parameters control the computation
process for this step:

```
    sge_option_pda = -pe smp 8 -q bigmem 
    sge_option_pla = -pe smp 2 -q bigmem
    ovlp_concurrent_jobs = 32
    ovlp_DBsplit_option = -x500 -s50
    ovlp_HPCdaligner_option = -v -dal4 -t32 -h60 -e.96 -l500 -s1000 
```

The setting is mostly parallel to the first overlapping step. The major
difference is the `-e` option in `ovlp_HPCdaligner_option`. The error rate is
much lower now so we expect much higher correlation between the p-reads.

Glossary
---------

sub-read : 

pre-assembly:

error correction:

overlap:

proper overlap:

string graph:

contig:

primary contig:

associate contig:

p-read: Pre-assembled Reads, error corrected reads through the pre-assembly process. 


