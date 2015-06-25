from falcon_kit.multiproc import Pool
import argparse
import subprocess as sp
import shlex


def filter_stats(lines, min_len):
        current_q_id = None
        contained = False
        ave_idt = 0.0
        all_over_len = 0.0
        overlap_data = {"5p":0, "3p":0}
        q_id = None
        rtn_data = []
        q_l = 0
        for l in lines:
            l = l.strip().split()
            q_id, t_id = l[:2]

            if q_id != current_q_id:
                left_count = overlap_data["5p"]
                right_count = overlap_data["3p"]
                if (current_q_id != None and
                        (left_count > 0 or right_count > 0)):
                    rtn_data.append( (current_q_id, q_l, left_count, right_count  ) )
                overlap_data = {"5p":0, "3p":0}
                current_q_id = q_id
                contained = False
                ave_idt = 0.0
                all_over_len = 0.0

            overlap_len = -int(l[2])
            idt = float(l[3])
            q_s, q_e, q_l = int(l[5]), int(l[6]), int(l[7])
            t_s, t_e, t_l = int(l[9]), int(l[10]), int(l[11])

            if q_l < min_len or t_l < min_len:
                continue

            if idt < 90:
                continue

            if not contained:
                if l[-1] in ("contains", "overlap"):
                    ave_idt += idt * overlap_len
                    all_over_len += overlap_len
                if q_s == 0:
                    overlap_data["5p"] += 1
                elif q_e == q_l:
                    overlap_data["3p"] += 1

        if q_id != None:
            left_count = overlap_data["5p"]
            right_count = overlap_data["3p"]
            if (left_count > 0 or right_count > 0):
                rtn_data.append( (q_id, q_l, left_count, right_count  ) )
            
        return rtn_data


def try_filter_stats(input_):
    fn, min_len = input_
    try:
        lines = sp.check_output(shlex.split("LA4Falcon -mo ../1-preads_ovl/preads.db %s" % fn)).splitlines()
        return fn, filter_stats(lines, min_len)
    except (KeyboardInterrupt, SystemExit):
        return


def main(*argv):
    parser = argparse.ArgumentParser(description='a simple multi-processes LAS ovelap data filter')
    parser.add_argument('--n_core', type=int, default=4,
                        help='number of processes used for generating consensus; '
                        '0 for main process only (default=%(default)s)')
    parser.add_argument('--fofn', type=str, help='file contains the path of all LAS file to be processed in parallel')
    parser.add_argument('--min_len', type=int, default=2500, help="min length of the reads")
    args = parser.parse_args(argv)
    exe_pool = Pool(args.n_core)

    file_list = open(args.fofn).read().split("\n")
    #print "all", len(contained)
    inputs = []
    for fn in file_list:
        if len(fn) != 0:
            inputs.append( (fn, args.min_len ) )
    for res in exe_pool.imap(try_filter_stats, inputs):  
        for l in res[1]:
            print " ".join([str(c) for c in l])
