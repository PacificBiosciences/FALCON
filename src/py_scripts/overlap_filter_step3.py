#!/usr/bin/env python
import sys
from numpy import array


current_q_id = None
contained_list = set()
chimer_list = set()
ignore_list = set()

read_end_data = {}

with open("rc_out_all") as f:
    for l in f:
        l = l.strip().split()
        if l[0] == "c":
            contained_list.add(l[1])
        if l[0] == "o":
            left_count, right_count = int(l[2]), int(l[3])
            left_c_max, right_c_max = int(l[4]), int(l[5])
            max_c = l[6]
            min_c = l[7]
            mean_c = l[8]
            ave_idt = l[9]
            read_end_data[l[1]] = ( left_count, right_count,
                                    left_c_max, right_c_max,
                                    max_c, min_c,
                                    mean_c, ave_idt)
            if abs(left_count - right_count) > 60:
                ignore_list.add(l[1])
            if left_count > 120 or right_count > 120:
                ignore_list.add(l[1])
            if left_count < 2 or right_count < 2:
                ignore_list.add(l[1])


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
            print " ".join(ovlp), read_end_data[current_q_id] 
            if i >= 3 and m_range > 1000:
                break
        
        for i in xrange(len(right)):
            score, m_range, ovlp = right[i]
            print " ".join(ovlp), read_end_data[current_q_id]
            if i >= 3 and m_range > 1000:
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
    if q_id in ignore_list:
        continue
    if t_id in ignore_list:
        continue

    if idt < 90:
        continue
    if q_l < 1000 or t_l < 1000:
        continue
    if q_s == 0:
        overlap_data["5p"].append( (-overlap_len,  t_l - (t_e - t_s),  l) )
    elif q_e == q_l:
        overlap_data["3p"].append( (-overlap_len, t_l - (t_e - t_s), l) )

left = overlap_data["5p"]
right = overlap_data["3p"]
left.sort()
right.sort()


for i in xrange(len(left)):
    score, m_range, ovlp = left[i]
    print " ".join(ovlp), read_end_data[current_q_id] 
    if i >= 3 and m_range > 1000:
        break

for i in xrange(len(right)):
    score, m_range, ovlp = right[i]
    print " ".join(ovlp), read_end_data[current_q_id]
    if i >= 3 and m_range > 1000:
        break
