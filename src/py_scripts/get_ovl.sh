/mnt/secondary/Share/HBAR_03202013/bin/parallel -j 32 "LA4Falcon -mo -H2000 {}  | python overlap_filter_step1.py > {}.ignore" ::: *.las
rm all.ignore
cat *.ignore > all.ignore
/mnt/secondary/Share/HBAR_03202013/bin/parallel -j 32 "LA4Falcon -mo -H2000 {}  | python overlap_filter_step2.py > {}.rc" ::: *.las
cat *.rc > rc_out_all
rm *.rc
/mnt/secondary/Share/HBAR_03202013/bin/parallel -j 32 "LA4Falcon -mo -H2000 {}  | python overlap_filter_step3.py > {}.ovl" ::: *.las
