mkdir ecoli_test/
cd ecoli_test/
mkdir data
cd data
wget https://www.dropbox.com/s/tb78i5i3nrvm6rg/m140913_050931_42139_c100713652400000001823152404301535_s1_p0.1.subreads.fasta
wget https://www.dropbox.com/s/v6wwpn40gedj470/m140913_050931_42139_c100713652400000001823152404301535_s1_p0.2.subreads.fasta
wget https://www.dropbox.com/s/j61j2cvdxn4dx4g/m140913_050931_42139_c100713652400000001823152404301535_s1_p0.3.subreads.fasta
cd ..
find $PWD/data -name "*.fasta" > input.fofn
cp ../FALCON/examples/fc_run_ecoli_2.cfg  .
fc_run.py fc_run_ecoli_2.cfg
