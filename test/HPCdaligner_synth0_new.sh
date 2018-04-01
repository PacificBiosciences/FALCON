# Daligner jobs (2)
daligner -v -w1 -h1 -t50 -H2000 -e0.99 -l1 -s1000 -P=. -mdust raw_reads.1 raw_reads.1
daligner -v -w1 -h1 -t50 -H2000 -e0.99 -l1 -s1000 -P=. -mdust raw_reads.2 raw_reads.1 raw_reads.2
# Check initial .las files jobs (2) (optional but recommended)
LAcheck -vS raw_reads raw_reads.1.raw_reads.1 raw_reads.1.raw_reads.2
LAcheck -vS raw_reads raw_reads.2.raw_reads.1 raw_reads.2.raw_reads.2
# Level 1 merge jobs (2)
LAmerge -v raw_reads.1 raw_reads.1.raw_reads.1 raw_reads.1.raw_reads.2
LAmerge -v raw_reads.2 raw_reads.2.raw_reads.1 raw_reads.2.raw_reads.2
# Check level 2 .las files jobs (2) (optional but recommended)
LAcheck -vS raw_reads raw_reads.1
LAcheck -vS raw_reads raw_reads.2
# Remove level 1 .las files (optional)
rm raw_reads.1.raw_reads.1.las raw_reads.1.raw_reads.2.las
rm raw_reads.2.raw_reads.1.las raw_reads.2.raw_reads.2.las
