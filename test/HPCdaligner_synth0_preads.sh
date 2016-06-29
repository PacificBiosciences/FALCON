# Daligner jobs (1)
daligner -v -h1 -t50 -H1 -e0.99 -l1 -s1000 preads.1 preads.1
# Initial sort jobs (1)
LAsort -v preads.1.preads.1.C0 preads.1.preads.1.N0 preads.1.preads.1.C1 preads.1.preads.1.N1 preads.1.preads.1.C2 preads.1.preads.1.N2 preads.1.preads.1.C3 preads.1.preads.1.N3 && LAmerge -v preads.1 preads.1.preads.1.C0.S preads.1.preads.1.N0.S preads.1.preads.1.C1.S preads.1.preads.1.N1.S preads.1.preads.1.C2.S preads.1.preads.1.N2.S preads.1.preads.1.C3.S preads.1.preads.1.N3.S
# Check all level 1 .las files (optional but recommended)
LAcheck -vS preads preads.1
# Remove initial .las files (optional)
rm preads.1.preads.1.C0.las preads.1.preads.1.N0.las preads.1.preads.1.C1.las preads.1.preads.1.N1.las preads.1.preads.1.C2.las preads.1.preads.1.N2.las preads.1.preads.1.C3.las preads.1.preads.1.N3.las
rm preads.1.preads.1.C0.S.las preads.1.preads.1.N0.S.las preads.1.preads.1.C1.S.las preads.1.preads.1.N1.S.las preads.1.preads.1.C2.S.las preads.1.preads.1.N2.S.las preads.1.preads.1.C3.S.las preads.1.preads.1.N3.S.las
