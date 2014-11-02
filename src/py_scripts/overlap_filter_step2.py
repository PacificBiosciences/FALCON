#!/usr/bin/env python2.7
import sys
from numpy import array, min, max, mean

"""002293727 000305923 -3871 98.44 0 0 3820 4468 1 0 3871 5543 overlap"""

ignore_ids = set()
with open("all.ignore") as f:
    for l in f:
        l = l.strip().split()
        ignore_ids.add(l[0])



current_q_id = None
contained_list = set()
chimer_list = set()
contained = False
ave_idt = 0.0
all_over_len = 0.0
overlap_data = {"5p":[], "3p":[]}
q_id = None

for l in sys.stdin:
    l = l.strip().split()
    q_id, t_id = l[:2]

    if current_q_id == None:
        current_q_id = q_id
        overlap_data = {"5p":[], "3p":[]}
        coverage = array([0] *  int(l[7]))
        contained = False
        ave_idt = 0.0
        all_over_len = 0.0

    elif q_id != current_q_id:

        left = overlap_data["5p"]
        right = overlap_data["3p"]
    
        if not contained and all_over_len > 0.01:
            print "o", current_q_id, len(left), len(right), max(coverage[:100]), max(coverage[-100:]), max(coverage[100:-100]), min(coverage[100:-100]), round(mean( coverage[50:-50] )), "%.1f" % (ave_idt/all_over_len), q_l
        else:
            print "c", current_q_id

        overlap_data = {"5p":[], "3p":[]}
        current_q_id = q_id
        coverage = array([0] *  int(l[7]))
        contained = False
        ave_idt = 0.0
        all_over_len = 0.0

    overlap_len = -int(l[2])
    idt = float(l[3])
    q_s, q_e, q_l = int(l[5]), int(l[6]), int(l[7])
    t_s, t_e, t_l = int(l[9]), int(l[10]), int(l[11])

    if idt < 90:
        continue

    if q_l < 500 or t_l < 500:
        continue

    #if l[-1] == "contains":
    #   contained_list.add( t_id ) 

    if t_id not in ignore_ids and  l[-1] == "contained":
        contained = True
        #contained_list.add( q_id )
    
    if not contained:
        if l[-1] in ("contains", "overlap"):
            ave_idt += idt * overlap_len
            all_over_len += overlap_len
            coverage[q_s+10:q_e-10] += 1
            
        if q_s == 0:
            overlap_data["5p"].append( (-overlap_len,  t_l - (t_e - t_s) ) )
        elif q_e == q_l:
            overlap_data["3p"].append( (-overlap_len, t_l - (t_e - t_s) ) )

if q_id !=  None:
    left = overlap_data["5p"]
    right = overlap_data["3p"]

    if not contained and all_over_len > 0.01:
        print "o", q_id, len(left), len(right), max(coverage[:100]), max(coverage[-100:]), max(coverage[100:-100]), min(coverage[100:-100]), round(mean( coverage[50:-50] )), "%.1f" % (ave_idt/all_over_len), q_l
    else:
        print "c", q_id

