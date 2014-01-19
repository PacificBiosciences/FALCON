mkdir 3-asm-falcon/
cd 3-asm-falcon/
cat ../2-preads-falcon/pread_00000.fa > preads.fa
falcon_overlap.py  --min_len 10000 --n_core 10 --d_core 1 preads.fa > preads.ovlp
falcon_asm.py preads.ovlp preads.fa
falcon_fixasm.py

export PATH=$PATH:/home/HBAR_ENV/MUMmer3.23
nucmer -maxmatch all_tigs.fa all_tigs.fa -p all_tigs_self >& /dev/null
show-coords -o -H -T all_tigs_self.delta | grep CONTAINS | awk '$7>96' | awk '{print $9}' | sort -u > all_tigs_duplicated_ids
remove_dup_ctg.py
cat p-tigs_nodup.fa a-tigs_nodup.fa > pa-tigs_nodup.fa
cat p-tigs_nodup.fa a-tigs_nodup.fa > /mnt/pa-tigs_nodup.fa

find /home/data/Ecoli/ -name "*.bax.h5" > /mnt/h5_input.fofn
cd /mnt
pbalign.py --forQuiver --nproc 32  --tmpDir /mnt --maxHits 1  h5_input.fofn pa-tigs_nodup.fa output.cmp.h5 
samtools faidx pa-tigs_nodup.fa
quiver -j 24 output.cmp.h5 -r pa-tigs_nodup.fa -o variants.gff -o consensus.fasta
