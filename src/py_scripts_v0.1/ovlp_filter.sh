source /mnt/secondary/Share/HBAR_03202013/bin/activate
parallel -j 24 "LA4Falcon -mo -H10000 {}  | python overlap_filter_step1.py > {}.ignore" ::: *.las
cat *.ignore > all.ignore
parallel -j 24 "LA4Falcon -mo -H10000 {}  | python overlap_filter_step2.py > {}.rc" ::: *.las
cat *.rc > rc_out_all
parallel -j 24 "LA4Falcon -mo -H10000 {}  | python overlap_filter_step3.py > {}.ovl" ::: *.las
