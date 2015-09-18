# Daligner jobs (2)
daligner -v -h1 -t16 -H1 -e0.7 -l1 -s1000 raw_reads.1 raw_reads.1
daligner -v -h1 -t16 -H1 -e0.7 -l1 -s1000 raw_reads.2 raw_reads.1 raw_reads.2
# Initial sort jobs (4)
LAsort -v raw_reads.1.raw_reads.1.C0 raw_reads.1.raw_reads.1.N0 && LAmerge -v L1.1.1 raw_reads.1.raw_reads.1.C0.S raw_reads.1.raw_reads.1.N0.S && rm raw_reads.1.raw_reads.1.C0.S.las raw_reads.1.raw_reads.1.N0.S.las
LAsort -v raw_reads.1.raw_reads.2.C0 raw_reads.1.raw_reads.2.N0 && LAmerge -v L1.1.2 raw_reads.1.raw_reads.2.C0.S raw_reads.1.raw_reads.2.N0.S && rm raw_reads.1.raw_reads.2.C0.S.las raw_reads.1.raw_reads.2.N0.S.las
LAsort -v raw_reads.2.raw_reads.1.C0 raw_reads.2.raw_reads.1.N0 && LAmerge -v L1.2.1 raw_reads.2.raw_reads.1.C0.S raw_reads.2.raw_reads.1.N0.S && rm raw_reads.2.raw_reads.1.C0.S.las raw_reads.2.raw_reads.1.N0.S.las
LAsort -v raw_reads.2.raw_reads.2.C0 raw_reads.2.raw_reads.2.N0 && LAmerge -v L1.2.2 raw_reads.2.raw_reads.2.C0.S raw_reads.2.raw_reads.2.N0.S && rm raw_reads.2.raw_reads.2.C0.S.las raw_reads.2.raw_reads.2.N0.S.las
# Level 1 jobs (2)
LAmerge -v raw_reads.1 L1.1.1 L1.1.2 && rm L1.1.1.las L1.1.2.las
LAmerge -v raw_reads.2 L1.2.1 L1.2.2 && rm L1.2.1.las L1.2.2.las
