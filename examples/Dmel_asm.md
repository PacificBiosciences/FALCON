Dmel Assembly with FALCON on Amazon EC2
=========================================

Preparation for Running StarCluster
-----------------------------------

I use a development version of StarCluster since the stable version does
not support the kind of instance that we need to use in AWS EC2.

```
    git clone https://github.com/jtriley/StarCluster.git
    git checkout 4149bbed292b0298478756d778d8fbf1dd210daf
```

First of all, I assume you already know how to create an AWS EC2 account, and
go through the tutorial for running it on EC2.

I build a public EC2 EBS snapshot. You should create a new EBS volume using
the `PacBio_Dmel_Asm` snapshop. It already contains the raw sequence fasta files and
an assembly as example already.

Here is an example of the configuration for StarCluster::

```
    [aws info]
    aws_access_key_id = your_access_key
    aws_secret_access_key = your_secret_access_key
    aws_user_id = your_user_id

    [volume DMEL]
    volume_id=your_dmel_data_EBS_id #e.g volume_id=vol-c9df3b85
    mount_path=/mnt/dmel_asm

    [cluster falcon-pre-asm]
    keyname = starcluster
    cluster_size = 1
    cluster_user = sgeadmin
    cluster_shell = bash
    master_image_id = ami-ef3c0e86
    master_instance_type = c3.8xlarge
    node_image_id = ami-ef3c0e86
    node_instance_type = c3.8xlarge
    availability_zone = us-east-1a
    volumes = DMEL

    [cluster falcon-bigmem]
    keyname = starcluster
    cluster_size = 1
    cluster_user = sgeadmin
    cluster_shell = bash
    master_image_id = ami-73d2d21a
    master_instance_type = cr1.8xlarge
    node_image_id = ami-ef3c0e86
    node_instance_type = c3.8xlarge
    availability_zone = us-east-1a
    volumes = DMEL

    [global]
    default_template = falcon-bigmem
    ENABLE_EXPERIMENTAL=True
```

I set up two cluster configurations for different part of the assembly process.
If you want to run end-to-end in one kind of instance, you can just use the 
`falcon-bigmem` for assembly. It costs a little bit more.

The AMI images (ami-ef3c0e86 and ami-73d2d21a) are pre-built with most package
necessary for the assembly work. If you will like to build your own, you can
check with this script:

```
    https://raw.github.com/PacificBiosciences/FALCON/v0.1.1/examples/install_note.sh
```

Get preassembled reads
------------------------

"Pre-assembly" is the process to error correct PacBio reads to generate
"preassembled reads" (p-reads) which have good accuracy to be assembled by
traditional Overlap-Layout-Consensus assembly algorithms directly. In this
instruction, we use an experimental code `falcon_qrm.py` to match the reads for
error correction. It is much faster than using `blasr` for the same purpose but
it may not as robust as `blasr` to generate high quality results yet as many
statistical properties of the algorithm is not fully studied.


First, let start an EC2 cluster of one node to set up a few things by running 
following `starcluster` command:

```
    starcluster start -c falcon-pre-asm falcon
```

Once the cluster is built, one can login the master node by:

```
    starcluster sshmaster falcon
```

We will need the following steps to setup the running environment::

1. update SGE environment

```
    cd /mnt/dmel_asm/sge_setup
    bash sge_setup.sh
```

2. setup HBAD-DTK environment

```
    . /home/HBAR_ENV/bin/activate
```

3. update HBAR-DTK and falcon_asm

```
    cd /mnt/dmel_asm/packages/pbtools.hbar-dtk-0.1.5
    python setup.py install
    cd /mnt/dmel_asm/packages/falcon_kit-0.1.1
    #edit falcon_asm.py to set identity threshold for overlapping at 98%, it is done in the EBS snapshot
    python setup.py install
```

If you want to do an assembly in `/mnt/dmel_asm/new_asm/`, just clone the 
configuration in `/mnt/dmel_asm/asm_template/` to `/mnt/dmel_asm/new_asm/`:

```
    cd /mnt/dmel_asm/
    cp -a asm_template/ new_asm/
    cd new_asm/
```

An example of the assembly result can be found in `/mnt/dmel_asm/asm_example`.

You can start the pre-assembly stage by running the `HBAR_WF3.py` script as following:

```
    python HBAR_WF3.py HBAR_step1.cfg
```

It will take while to preparing the fasta files for pre-assembly. Once that is
one, SGE jobs for matching reads will be submitted. Once the SGE jobs are
submmited, you can use add more node to run the jobs concurrently to speed up
the process by issuing this command on your local host to add the nodes:

    starcluster addnode -n 15 falcon # add 15 nodes 

When all nodes are up. You can try to run the load balancer so once the jobs are
done, the node can be terminated automatically to save some money.

    starcluster loadbalance -k 9 -K -m 16 -n 1 falcon

I found I have to comment out one line of code in `starcluster/plugins/sge.py`
to make it work properly to remove unused nodes:
    
    class SGEPlugin(clustersetup.DefaultClusterSetup):
        def _remove_from_sge(self, node):
            #comment out the following line in the code
            #master.ssh.execute('qconf -de %s' % node.alias)

If you use 16 nodes, it will takes about 4 hours to finish all jobs.  If all
pre-assembly jobs finish the cluster will be terminated automatically, but the
results will be kept in the EBS volume.

The generated p-reads will be in `/mnt/dmel_asm/new_asm/2-preads-falcon/pread_*.fa`.

Assembling the p-reads
------------------------

We use a different instance type which has bigger memory to assemble the genome. We
only needs one node for the assembly part.  We still use SGE as the code was written 
to run end-to-end assembly in a general SGE cluster. First, start single node cluster by
running the commands in the local host:

```
    starcluster start -c falcon-bigmem falcon
```

Repeat the setup process:

```
    cd /mnt/dmel_asm/sge_setup
    bash sge_setup.sh

    . /home/HBAR_ENV/bin/activate

    cd /mnt/dmel_asm/packages/pbtools.hbar-dtk-0.1.5
    python setup.py install
    cd /mnt/dmel_asm/packages/falcon_kit-0.1.1
    #edit falcon_asm.py to set identity threshold for overlapping at 98%, it is done in the EBS snapshot
    python setup.py install
```

You can start the assembly stage by running the `HBAR_WF3.py` script as following:

```
    cd /mnt/dmel_asm/new_asm/
    python HBAR_WF3.py HBAR_step2.cfg
```

It takes about two hours for the assembly process to finish. The results will 
be in `/mnt/dmel_asm/new_asm/3-asm-falcon`. 

Here is a list of the output files:

```
    full_string_graph.adj  # the adjacent nodes of the edges in the full string graph
    string_graph.gexf      # the gexf file of the string graph for graph visualization
    string_graph.adj       # the adjecent nodes of the edges in the string graph after transitive reduction
    edges_list             # full edge list 
    paths                  # path for the unitigs
    unit_edges.dat         # path and sequence of the untigs
    uni_graph.gexf         # unitig graph in gexf format 
    unitgs.fa              # fasta files of the unitigs
    all_tigs_paths         # paths for all final contigs (= primary contigs + associated contigs)
    all_tigs.fa            # fasta file for all contigs
    primary_tigs_paths     # paths for all primary contigs 
    primary_tigs.fa        # fasta file fot the primary contigs
    primary_tigs_paths_c   # paths for all primary contigs, detectable mis-assemblies are broken 
    primary_tigs_c.fa      # fasta file fot the primary contigs, detectable mis-assemblies are broken
    asm_graph.gexf         # the assembly graph where the edges are the contigs
```

There might be redundant contig. The following script can be used to remove
redundant contigs:

```
    export PATH=$PATH:/home/HBAR_ENV/MUMmer3.23
    nucmer -mum all_tigs.fa all_tigs.fa -p all_tigs_self >& /dev/null
    show-coords -o -H -T all_tigs_self.delta | grep CONTAINS | awk '$7>96' | awk '{print $9}' | sort -u > all_tigs_duplicated_ids
    remove_dup_ctg.py
    cat p-tigs_nodup.fa a-tigs_nodup.fa > pa-tigs_nodup.fa
```

The non-reduant set of contigs in `pa-tigs_nodup.fa` will be suitable for further correction
by the Quvier algorithm. 

-
Jason Chin, March 9, 2014

