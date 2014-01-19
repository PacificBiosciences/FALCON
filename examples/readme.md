Running an AWS instance that has HBAR-DTK + Falcon pre-installed
================================================================

1. Install the latest verison of StarCluster
```
    git clone https://github.com/jtriley/StarCluster.git
    cd StarCluster
    python setup.py install #better in virtualenv
```
The stable versio nof StarCluster does not support the `c3` instance.  For
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

6. There is alreay an existing assembly results iVn `/home/Ecoli_ASM/`. Here I
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

7. The next part of the assembly is not started automatically yet. The detail
steps are in the `run_asm.sh` script and one can use to get contigs and
consensus. 
```
    cp /home/Ecoli_ASM/run_asm.sh .
    bash run_asm.sh
```
The consensus result is in `/mnt/consensus.fasta`. Since we did not do any
consensus after the unitig step. One more run of quiver consensus may further
improve the final assembly accuracy.

8. A yeast (S. cerevisiae W303) data is also included in the AMI. One can try
to assemble it with a larger cluster setting.

--
Jason Chin, 01/18/2014

