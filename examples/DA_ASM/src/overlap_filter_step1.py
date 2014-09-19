#!/usr/bin/env python2.7
import sys
from numpy import array

"""002293727 000305923 -3871 98.44 0 0 3820 4468 1 0 3871 5543 overlap"""

current_q_id = None
contained_list = set()
chimer_list = set()

for l in sys.stdin:
    l = l.strip().split()
    q_id, t_id = l[:2]

    if q_id in contained_list:
        continue
    if q_id in contained_list:
        continue

    if current_q_id == None:
        current_q_id = q_id
        overlap_data = {"5p":[], "3p":[]}
        coverage = array([0] *  int(l[7]))

    elif q_id != current_q_id:

        left = overlap_data["5p"]
        right = overlap_data["3p"]
        if len(left) != 0 and len(right) != 0:
            if min(coverage[100:-100]) < 2 and max(coverage[100:-100]) > 12: #for low coverage region, we don't call chimer
                chimer_list.add(current_q_id)

        overlap_data = {"5p":[], "3p":[]}
        current_q_id = q_id
        coverage = array([0] *  int(l[7]))

    overlap_len = -int(l[2])
    idt = float(l[3])
    q_s, q_e, q_l = int(l[5]), int(l[6]), int(l[7])
    t_s, t_e, t_l = int(l[9]), int(l[10]), int(l[11])
    if idt < 90:
        continue
    coverage[q_s+10:q_e-10] += 1

    if q_l < 500 or t_l < 500:
        continue

    if l[-1] == "contains":
        contained_list.add( t_id ) 

    if l[-1] == "contained":
        contained_list.add( q_id )
    
    if q_s == 0:
        overlap_data["5p"].append( (-overlap_len,  t_l - (t_e - t_s),  l) )
    elif q_e == q_l:
        overlap_data["3p"].append( (-overlap_len, t_l - (t_e - t_s), l) )

for r_id in contained_list:
    print r_id, "1"
for r_id in chimer_list:
    print r_id, "2"

