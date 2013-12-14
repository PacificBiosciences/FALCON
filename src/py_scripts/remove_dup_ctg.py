import pbcore.io

import sys
"""nucmer -maxmatch all_tigs.fa all_tigs.fa -p all_tigs_self >& /dev/null"""
"""show-coords -o -H -T all_tigs_self.delta | grep CONTAINS | awk '$7>96' | awk '{print $9}' | sort -u > all_tigs_duplicated_ids"""

id_to_remove = set()
with open("all_tigs_duplicated_ids") as f:
    for l in f:
        l = l.strip().split("-")
        major, minor = l[:2]
        id_to_remove.add ( (major, minor) )

f = pbcore.io.FastaReader("all_tigs.fa")
with open("a-tigs_nodup.fa", "w") as f_out:
    for r in f:
        major, minor = r.name.split()[0].split("-")[:2]
        if minor == "0000":
            continue
        if (major, minor) in id_to_remove:
            continue
        if len(r.sequence) < 500:
            continue
        print >>f_out, ">"+r.name
        print >>f_out, r.sequence

f = pbcore.io.FastaReader("primary_tigs_c.fa")
with open("p-tigs_nodup.fa", "w") as f_out:
    for r in f:
        major, minor = r.name.split()[0].split("_")[:2]
        if (major, "0000") in id_to_remove:
            continue
        if len(r.sequence) < 500:
            continue
        print >>f_out, ">"+r.name
        print >>f_out, r.sequence
