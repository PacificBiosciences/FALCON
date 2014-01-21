Running an Amazon EC2 instance that has HBAR-DTK + Falcon pre-installed
=======================================================================

1. Install the latest verison of StarCluster
```
    git clone https://github.com/jtriley/StarCluster.git
    cd StarCluster
    python setup.py install #better in virtualenv
```
The stable version of StarCluster does not support the `c3` instance.  For
assembly, using one node of `c3.8xlarge` instance is more convenient. In my
test, I can finish single E. coli genome within almost one hour. Namely, one can
assembly a bacteria genome in less then 5 bucks.

2. Use the `StarCluster.cfg` as the configuration file for `StarCluster` to
setup a `falcon` cluster

3. Start the cluster 
```
    starcluster start falcon
```

4. login to the cluster
```
    starcluster sshmaster falcon
```

5. set up the SGE
```
    cd /home/sge_setup
    bash sge_setup.sh
```

6. There is alreay an existing assembly results in `/home/Ecoli_ASM/`. Here I
show how to reproduce it. First, create a new assembly working directory in
`/mnt`, set it up and run HBAR_WF3.py to get preassembled reads
```
    cd /mnt
    mkdir test_asm
    cd test_asm
    cp /home/Ecoli_ASM/HBAR.cfg .
    cp /home/Ecoli_ASM/input.fofn .
    source /home/HBAR_ENV/bin/activate
    HBAR_WF3.py HBAR.cfg
```

7. The next part of the assembly does not start automatically yet. The detail
steps are in the `run_asm.sh` script and one can use to get contigs and
consensus. 
```
    cp /home/Ecoli_ASM/run_asm.sh .
    bash run_asm.sh
```
The consensus result is in `/mnt/consensus.fasta`. Since we did not do any
consensus after the unitig step. One more run of quiver consensus may further
improve the final assembly accuracy.

8. A yeast (S. cerevisiae W303) data set is also included in the AMI. One can try
to assemble it with a larger cluster setting.


9. Here is the result of a timing test:
```
    (HBAR_ENV)root@master:/mnt/test_asm# time HBAR_WF3.py HBAR.cfg
    
    Your job 1 ("mapping_task_q00002_t000011416727c") has been submitted
    Your job 2 ("qf_task_q00002a3e75f4c") has been submitted
    Your job 3 ("mapping_task_q00003_t00001b667b504") has been submitted
    Your job 4 ("qf_task_q000036974ef22") has been submitted
    Your job 5 ("mapping_task_q00001_t000017bf52d9c") has been submitted
    Your job 6 ("qf_task_q000010b31d960") has been submitted
    Your job 7 ("pa_task_000001ee38aee") has been submitted
    
    
    
    real    26m51.030s
    user    1m10.152s
    sys     0m11.993s
    
    (HBAR_ENV)root@master:/mnt/test_asm# time bash run_asm.sh
    [WARNING] This .cmp.h5 file lacks some of the QV data tracks that are required for optimal performance of the Quiver algorithm.  For optimal results use the ResequencingQVs workflow in SMRTPortal with bas.h5 files from an instrument using software version 1.3.1 or later.

    real    13m2.945s
    user    244m44.322s
    sys     2m7.032s
```
For a better results, one might run quiver twice. It is possible to do the whole assembly within one hour (~ 26 + 13 * 2 = 52 minutes). With the overhead on setting up, file transfer, etc., one can assembly a bacteria genome in EC2 less than 5 bucks in principle.


--
Jason Chin, 01/18/2014

