#!/usr/bin/env python
import sys
from numpy import array


current_q_id = None
contained_list = set()
chimer_list = set()

with open("rc_out_all") as f:
    for l in f:
        l = l.strip().split()
        if l[1] == "1":
            contained_list.add(l[0])
        if l[1] == "2":
            chimer_list.add(l[0])

for l in sys.stdin:
    l = l.strip().split()
    q_id, t_id = l[:2]

    if current_q_id == None:
        current_q_id = q_id
        overlap_data = {"5p":[], "3p":[]}

    elif q_id != current_q_id:

        left = overlap_data["5p"]
        right = overlap_data["3p"]
        left.sort()
        right.sort()

        for i in xrange(len(left)):
            score, m_range, ovlp = left[i]
            print " ".join(ovlp) 
            if i >= 2 and m_range > 1000:
                break
        
        for i in xrange(len(right)):
            score, m_range, ovlp = right[i]
            print " ".join(ovlp) 
            if i >= 2 and m_range > 1000:
                break

        overlap_data = {"5p":[], "3p":[]}
        current_q_id = q_id

    overlap_len = -int(l[2])
    idt = float(l[3])
    q_s, q_e, q_l = int(l[5]), int(l[6]), int(l[7])
    t_s, t_e, t_l = int(l[9]), int(l[10]), int(l[11])

    if q_id in contained_list:
        continue
    if t_id in contained_list:
        continue
    if idt < 90:
        continue
    if q_l < 1000 or t_l < 1000:
        continue
    if q_id in chimer_list or t_id in chimer_list:
        continue
    if q_id in contained_list or t_id in contained_list:
        continue
    if q_s == 0:
        overlap_data["5p"].append( (-overlap_len,  t_l - (t_e - t_s),  l) )
    elif q_e == q_l:
        overlap_data["3p"].append( (-overlap_len, t_l - (t_e - t_s), l) )

