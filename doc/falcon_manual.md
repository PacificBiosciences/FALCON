
![Falcon Logo](https://github.com/PacificBiosciences/FALCON/blob/master/doc/falcon_icon2.png)

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

Here is the environment this example is used to construct a user space
python `virtualenv` for running Falcon.

```
$ uname -a
Linux login14-biofx 3.13.0-24-generic #47-Ubuntu SMP Fri May 2 23:30:00 UTC 2014 x86_64 GNU/Linux

$ gcc --version
gcc (Ubuntu 4.8.2-19ubuntu1) 4.8.2
Copyright (C) 2013 Free Software Foundation, Inc.
This is free software; see the source for copying conditions.  There is NO
warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

$ python --version
Python 2.7.6 

# I do plan to move the Python 3 sometime.
```

Running `virtualenv` in a shared directory across the computational nodes:

```
$ cd installation_dir
$ virtualenv --no-site-packages  --always-copy   $PWD/fc_env
New python executable in /home/UNIXHOME/jchin/task2014/falcon_pb_github/fc_env/bin/python
Installing setuptools, pip...done.
```

Activate the virtual environment,
```
$ . $PWD/fc_env/bin/activate
(fc_env) $ which python
installation_dir/fc_env/bin/python
```

We will need two python packges `pypeFLOW` and `FALCON` and we need compile the
DALINGER code and put the binary executables into the virtual environment.

Here is a simple shell script that creates the environment and installs the
packages:
``` 
virtualenv --no-site-packages  --always-copy   $PWD/fc_env
. $PWD/fc_env/bin/activate
git clone https://github.com/cschin/pypeFLOW 
cd pypeFLOW

cd ..
git clone https://github.com/PacificBiosciences/FALCON.git
cd FALCON
python setup.py install

cd ..
git clone https://github.com/cschin/DAZZ_DB.git
cd DAZZ_DB/
make
cp DBrm DBshow DBsplit DBstats fasta2DB ../fc_env/bin/ 

cd ..
git clone https://github.com/cschin/DALIGNER.git
cd DALIGNER
make
cp daligner daligner_p DB2Falcon HPCdaligner LA4Falcon LAmerge LAsort  ../fc_env/bin
cd ..
```

Here is a quick testing to try it out for a simple E. coli assembly

```
. $PWD/fc_env/bin/activate
cd ecoli_test/
mkdir data
cd data
wget https://www.dropbox.com/s/tb78i5i3nrvm6rg/m140913_050931_42139_c100713652400000001823152404301535_s1_p0.1.subreads.fasta
wget https://www.dropbox.com/s/v6wwpn40gedj470/m140913_050931_42139_c100713652400000001823152404301535_s1_p0.2.subreads.fasta
wget https://www.dropbox.com/s/j61j2cvdxn4dx4g/m140913_050931_42139_c100713652400000001823152404301535_s1_p0.3.subreads.fasta
cd ..
find $PWD/data -name "*.fasta" > input.fofn
cp ../FALCON/examples/fc_run_ecoli.cfg .
cp ../../ecoli_test/input.fofn .
fc_run.py fc_run_ecoli.cfg 
```


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


## Overlap filtering

Not all overlap are "independent" so it is possible to impose some filtering
step to reduce computation and assembly graph complexity. For example, if a read
is fully contained in another read, the overlap information between these two
read does not provide extra information for re-constructing the genome. Also,
due to the transitive property of the overlapping relations, many overlap
information can be simply inferred.  In fact, the first stage for constructing
contigs are to remove the transitive reducible edges. It means that we might
just needs the "best `n` overlaps" in the 5' or 3' ends.  The `--bestn`
parameter in `overlap_filtering_setting` option can be used to control the
maximum overlap reported for each read.

Another useful heuristics is to only keep reads that has average 5' and 3'
coverage.  If a read is ended in a repeat regions, it might have higher than
normal coverages of the end which is a repeat.  Such reads does not provide much
value for uniquely resolving the related repeats.  We can filter them out and
hopefully there are reads which span through the repeats and have "normal"
unique anchors on both ends.  Also, if the coverage is too low on one end of a
read, it could be just too many errors or sequencing artifacts over there.  Such
reads create "spurs" in the assembly graph which are typically filtered out
anyway. The `--max_cov` and `--min_cov` are used for filtering reads that have
too high or too low overlaps. 

The filtering scripts also allows filtering out some "split" reads.  If a read
have very unequal coverages between the 5' and 3' ends, it can be also the signal that
one end is a repeat.  The `--max_diff` parameter can be used to filter out the
reads where one ends has much more coverage than the other end.

What is the right numbers used for these parameters?  These parameters may the
most tricky ones to be set right.  If the overall coverage of the error
corrected reads longer than the length cut off is known and reasonable high
(e.g. greater than 20x), it might be safe to set `min_cov` to be 5, `max_cov` to
be three times of the average coverage and the `max_diff` to be twice of the
average coverage.  However, in low coverage case, it might better to set
`min_cov` to be one or two.  A helper script called `fc_ovlp_stats.py` can help
to dump the number of the 3' and 5' overlap of a given length cutoff, you can
plot the distribution of the number of overlaps to make a better decision. 

One can also set the `max_diff` and `max_cov` to be really high to avoid any
filtering if that is preferred in some cases.  

This filtering process will certainly filter out information about high copy
repeats.  Namely, those repeats will likely to be filtered out totally and do
not appear in the final assembly.  If you are interested in those repeats even
though they may not be able to placed within some longer contig, you will
probably want to avoid filtering them out or process them differently.  In
general, it might be more efficient and useful to process those repeats
separately. Including them in the assembly process typically does not help much
for getting better contiguity and maybe messy for post-processing with current
algorithms. I think it is a very interesting but also very challenging
bioinformatics topic on how to process these repeats better for improving
assembly beside understand the nature of these repeats. 


## Constructing graph from overlaps

Given the overlapping data, the string graph is constructed by
`fc_ovlp_to_graph.py` using the default parameters. `fc_ovlp_to_graph.py`
generated several files representing the final string graph of the assembly. The
final `ctg_path` contain the information of the graph of each contig. A contig
is a linear of path of simple paths and compound paths. "Compound paths" are
those subgraph that is not simple but have unique inlet and outlet after graph
edge reduction. They can be induced by genome polymorphism or sequence errors.
By explicitly encoding such information in the graph output, we can examine the
sequences again to classify them later.

(TODO: write more details about the layout rule and how it is useful for
polyploid assembly.) 


## Constructing contig from graph

The final step to create draft contigs is to find a single path of each contig
graph and to generate sequences accordingly.  In the case that a contig graph is
not a simple path, we find the end-to-end path that has the most overlapped
bases. This is called as the `primary contig`.  For each compound path within 
the graph, if an alternative path different from primary one is possible, we
will construct the associated contig.  In the case which the associated contigs
are induced by sequencing error, the identity of the alternative contig and the
primary contig will be high ( > 99% identity most of time). In the case where
there are true structure polymorphism, there are typically bigger difference
between the associated contigs and the primary contigs.   

The script "fc_graph_to_contig.py" takes the sequence data and graph output to
construct contigs. It generated all associated contigs at this moment. Some
post-processing procedure to de-duplicate some of the associated contigs induced
by errors will be developed in the future. ( You are encourage to find some
creative way to solve this problem for sure. )

# Working directory structure, job recovery and trouble shooting

The code is designed to work in single directory.  The typical layout of a
working directory looks like this:

```
$ ls -l
total 56
drwxr-xr-x 84 jchin Domain Users  8192 Nov 30 12:30 0-rawreads
drwxr-xr-x 18 jchin Domain Users  4096 Nov 30 12:33 1-preads_ovl
drwxr-xr-x  2 jchin Domain Users  4096 Nov 30 12:44 2-asm-falcon
-rwxr-xr-x  1 jchin Domain Users  1041 Nov 30 12:13 fc_run.cfg
-rw-r--r--  1 jchin Domain Users   378 Nov 29 23:20 input.fofn
drwxr-xr-x  2 jchin Domain Users  4096 Nov 30 12:13 scripts
drwxr-xr-x  2 jchin Domain Users 24576 Nov 30 12:33 sge_log
```

A typical input.fofn looks like this:

```
/mnt/data/Analysis_Results/m140913_050931_42139_c100713652400000001823152404301535_s1_p0.1.subreads.fasta
/mnt/data/Analysis_Results/m140913_050931_42139_c100713652400000001823152404301535_s1_p0.2.subreads.fasta
/mnt/data/Analysis_Results/m140913_050931_42139_c100713652400000001823152404301535_s1_p0.3.subreads.fasta
```

## Inside `0-rawreads` directory
The directory `0-rawreads` includes all the scripts and data for overlapping the
raw sequences. It contains various `job_*` and `m_*` directories:

For example, if we divide the E. coli data into 20 chunks, the directory look
like this,

```
cns_done   job_00011  job_00024  job_00037  job_00050  m_00003  m_00016
da_done    job_00012  job_00025  job_00038  job_00051  m_00004  m_00017
job_00000  job_00013  job_00026  job_00039  job_00052  m_00005  m_00018
job_00001  job_00014  job_00027  job_00040  job_00053  m_00006  m_00019
job_00002  job_00015  job_00028  job_00041  job_00054  m_00007  m_00020
job_00003  job_00016  job_00029  job_00042  job_00055  m_00008  preads
job_00004  job_00017  job_00030  job_00043  job_00056  m_00009  prepare_db.sh
job_00005  job_00018  job_00031  job_00044  job_00057  m_00010  raw_reads.db
job_00006  job_00019  job_00032  job_00045  job_00058  m_00011  rdb_build_done
job_00007  job_00020  job_00033  job_00046  job_00059  m_00012  run_jobs.sh
job_00008  job_00021  job_00034  job_00047  las_files  m_00013
job_00009  job_00022  job_00035  job_00048  m_00001    m_00014
job_00010  job_00023  job_00036  job_00049  m_00002    m_00015
```

The `job_*` directories store the output for each `daligner` job. The `m_*`
directories are the working directory for merging jobs. There are some empty
files which the name is ended with `done`. The time stamp of those files are
used to track the stage of the workflow.  You can modify the time stamps and
re-satrt the `fc_run.py` to trigger doing the computation for certain part of
the workflow.  However, it is not recommended unless you have some time to read
the source code of `fc_run.py` to know how the dependences inside the workflow
are tracked. (For example, if you `touch rdb_build_done` after a successfully
assembly and re-run `fc_run.py`, since all intermediate processes depends on the
file and the `rdb_build_done` is newer than any of the intermediate files, it
will trigger the `fc_run.py` to repeat the whole assembly process again.)

The `las_files` stores the final alignment information. If you do not plan to
re-run the jobs but like to know how the alignments look like, you can delete
all `job_*` and `m_*` directory but keep the `las_files` and `preads`
directories.

The major output of this step is stored in `0-rawreads/preads`. The `out.%04d.fa`
inside `0-rawreads/preads` are the fasta files of the output reads. The sentinel
file `cns_done` will be created if this step is successfully finished.

## Inside `1-preads_ovl` directory

This directory store the data of the p-read vs. p-read overlapping. It is
overall similar to the `0-rawreads` directory, but without the consensus step.
The major output are the `*.las` files inside `1-preads_ovl/` directory.

## The `2-asm-falcon` directory

This is final output directory. It contains the information the assembly graph
and the draft contigs.  The detail will be describe in the `Graph output format`
section. 

# Options for `fc_ovlp_to_graph.py` and assembly graph output format

Here is the usage information for running `fc_ovlp_to_graph.py`:

```
usage: fc_ovlp_to_graph.py [-h] [--min_len MIN_LEN] [--min_idt MIN_IDT]
                           [--lfc]
                           overlap_file

a example string graph assembler that is desinged for handling diploid genomes

positional arguments:
  overlap_file       a file that contains the overlap information.

optional arguments:
  -h, --help         show this help message and exit
  --min_len MIN_LEN  minimum length of the reads to be considered for
                     assembling
  --min_idt MIN_IDT  minimum alignment identity of the reads to be considered
                     for assembling
  --lfc              use local flow constraint method rather than best overlap
                     method to resolve knots in string graph
```

In some case, you might want to lower the `min_idt` to keep more overlap or
increase `min_len` to reduce the number of overlap used for constructing the
contig after the overlap filtering step.  The `--lfc` toggles the rule for
resolving local knots in the graph. If `--lfc` is not specified, "the best
overlapped edge" will be kept when there are multiple in- or out- edges from a
node while the others will be removed.  


The first stage of the assembly is to construct the initial string graph and
classify each edges in the string graph.  `sg_edges_list` contained the
information of the information of the edges in the full string graph and the
classification. For example, 5 edges are shown in the five lines of the file
below

```
$ head -5 sg_edges_list 
000017363:B 000007817:E 000007817 10841 28901 10841 99.52 TR
000015379:E 000004331:B 000004331 6891 0 18178 99.35 TR
000006813:B 000000681:E 000000681 7609 23795 7616 99.72 TR
000002258:E 000002505:B 000002505 5850 0 17215 99.62 TR
000013449:B 000012565:B 000012565 3317 0 20570 99.72 G
```

The first two columns indicates the in and out node of the edge.  The node
notation contain two files operated by `:`.  The first field is the read
identifier. The second field is either `B` or `E`. `B` is the 5' end of the read
and `E` is the 3' end of the reads. The next three field indicates the
corresponding sequences of the edges. In this example, the edge in the first
line contains the sequence from read `000007817` base `[10841, 28901)`. If the
second coordinate is smaller than the first one, it means the corresponded
sequence is reverse complimented.  The next two column are the number of
overlapped base and the overlap identity. The final column is the
classification.  Currently, there are 4 different types `G`, `TR`, `R`, and `S`.
An edge with type "`G`" is used for the final string graph. A "`TR`" means the edge is
transitive reducible.  "`R`" means the edge is removed during the local repeat
resolution and "`S`" means the edge is likely to be a "spur" which only one ends
is connected.


The initial string graph is further to be simplified into a set of "unitig" edges.
The `utg_data` file contains the details of each unitig. Each line in the file
represents a unitig. The first three files are "start node", "via node", and "end
node".  Two unitgs might have the same "start node" and "end node", so we need
another "via node" to uniquely identify the unitigs. Here is an example of the
`utg_data` files:

```
$ head -10 utg_data
000015696:B 000009028:B 000016941:B contained 16438 134865 000015696:B~000006612:B~000002456:B~000014643:B~000007407:B~000015939:E~000009028:B~000016941:B
000010623:B 000015633:B 000014991:B contained 30158 18666 000010623:B~000015633:B~000014991:B
000015636:B 000002245:B 000010757:E contained 15402 40356 000015636:B~000002245:B~000010757:E
000014184:E NA 000012028:E compound 14895 56928 000014184:E~000012765:E~000012028:E|000014184:E~000007953:B~000012028:E
000010757:B NA 000015636:E compound 15402 40356 000010757:B~000002245:E~000015636:E|000010757:B~000014783:E~000015636:E
000014184:E 000007953:B 000012028:E contained 14792 32932 000014184:E~000007953:B~000012028:E
000010623:B NA 000014991:B compound 30148 163627 000010623:B~000015633:B~000014991:B|000010623:B~000001407:B~000014991:B
000012028:B 000012765:B 000014184:B contained 19137 56928 000012028:B~000000382:E~000012765:B~000014184:B
000016941:B 000003353:B 000008783:B simple 88381 615439 000016941:B~000003353:B~000010261:B~000011789:E~000017006:B~000016307:B~...
000014991:B 000013790:E 000002926:B simple 392373 2274104 000014991:B~000013790:E~000004614:B~000003329:B~000004898:B~000000461:B~000017105:E~... 
```

The forth field indicates the type of the unitigs, the fifth field is the
estimate length of the untig and the six field is the total number of overlapped
bases in the unitig.  There are three kinds of unitigs: "simple", "contained",
and "compound".  "Simple" unitigs are those unitigs which are just a simple path
(every node has one in- and one out-edge except the begining and ending nodes of
the path.) It is represented by a list of nodes which each node is separated by
`~` characters in the 7th column. The "contained" contigs are simple path but
those unitigs are also part of other "compound" paths. The "compound" unitigs
represents bubble-like subgraph in the graph. While it is not "simple", it has
well defined in- and out- nodes and they are treated as a single unit when the
contigs are constructed. The structure inside a "compound" unitig can be from
biological nature or sequencing/alignment errors.  Each edge in the "compound"
unitig subgraph are encoded explicitly as a collection of simple contained
unitigs in the 7th column. The contained unitigs within a compound unitig 
are separated by the `|` character.

The file `ctg_paths` encodes the graph for each contig after the unitigs are
analyzed and put into contigs. Each line has 7 columns.  The first column is the
contig ID. The contig ID are just the serial numbers followed by `R` or `F`. Two
contigs with same serial number but different endings are "dual" to each other.
Namely, they are constructed from "dual" edges and they are mostly reverse
complemented to each other except near the ends of the contigs.  The second
column is the type of contig. If a unitig is circular (the beginning node and
the ending node are the same), then it will be marked as "`ctg_circular`".
Everything else will be "`ctg_linear`".  In some case, even a contig is marked
as "`ctg_linear`", it can be still a circular contig if the beginning node and
the ending node are the same but it is not a "simple" path.  One can detect that
by checking the beginning and ending nodes if necessary. 

The third field indicates the first unitig in the contig in the form of
`begin_node~via_node~end_node`. The fourth field is the ending node of the
contig. The 5th and 6th fields are the estimated length and the overlapped based
of the contig respectively.  The final column are the unitigs in the contig. The
three node format unitig IDs are separated by `|`.

# De-dup example and strategy

TODO

# Low coverage assembly

TODO

# Assembly layout rule
TODO


# CPU usages

This command helps to calculate the user cpu time.
```
 find . -name "*log" | xargs cat | grep sys | awk -F "u" '{print $1}' | awk '{s+=$1/3600};{print s}'
```

# Applications of the assembly graph

TODO

# Reproducibility and replicability 


# C code for sequence alignment and consensus
---------------------------------------------

Several C code files for implementing sequence matching, alignment and consensus:

```
    kmer_lookup.c  # kmer match code for quickly identify potential hits
    DW_banded.c    # function for detailed sequence alignment
                   # It is based on Eugene Myers' Paper 
                   # "AnO(ND) difference algorithm and its variations", 1986, 
                   # http://dx.doi.org/10.1007/BF01840446
    falcon.c       # functions for generating consensus sequences for a set of multiple sequence alginment
    common.h       # header file for common declaration
```

A python wrapper library using Python's ctypes to call the C functions: falcon_kit.py


TDOD

# Other TODOs

* Incremental overlapping
* Pre-processing repeat for overlapping

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

associated contig:

p-read: Pre-assembled Reads, error corrected reads through the pre-assembly process. 

compound path:

simple path:
