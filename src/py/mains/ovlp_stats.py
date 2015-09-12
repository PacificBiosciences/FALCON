from falcon_kit.multiproc import Pool
import falcon_kit.util.io as io
import argparse
import shlex
import subprocess as sp
import sys

Reader = io.CapturedProcessReaderContext


def filter_stats(readlines, min_len):
        current_q_id = None
        contained = False
        ave_idt = 0.0
        all_over_len = 0.0
        overlap_data = {"5p":0, "3p":0}
        q_id = None
        rtn_data = []
        q_l = 0
        for l in readlines():
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


def run_filter_stats(fn, min_len):
    cmd = "LA4Falcon -mo ../1-preads_ovl/preads.db %s" % fn
    reader = Reader(cmd)
    with reader:
        return fn, filter_stats(reader.readlines, min_len)

def run_ovlp_stats(exe_pool, file_list, min_len):
    inputs = []
    for fn in file_list:
        if len(fn) != 0:
            inputs.append( (run_filter_stats, fn, min_len ) )
    for res in exe_pool.imap(io.run_func, inputs):
        for l in res[1]:
            print " ".join([str(c) for c in l])

def try_run_ovlp_stats(n_core, fofn, min_len):
    io.LOG('starting ovlp_stats')
    file_list = io.validated_fns(fofn)
    io.LOG('fofn %r: %r' %(fofn, file_list))
    n_core = min(n_core, len(file_list))
    exe_pool = Pool(n_core)
    try:
        run_ovlp_stats(exe_pool, file_list, min_len)
        io.LOG('finished ovlp_stats')
    except KeyboardInterrupt:
        io.LOG('terminating ovlp_stats workers...')
        exe_pool.terminate()

def ovlp_stats(fofn, min_len, n_core, stream, debug, silent):
    if debug:
        n_core = 0
        silent = False
    if silent:
        io.LOG = io.write_nothing
    if stream:
        global Reader
        Reader = io.StreamedProcessReaderContext
    try_run_ovlp_stats(n_core, fofn, min_len)

def parse_args(argv):
    parser = argparse.ArgumentParser(description='a simple multi-processes LAS ovelap data filter')
    parser.add_argument('--n_core', type=int, default=4,
                        help='number of processes used for generating consensus; '
                        '0 for main process only (default=%(default)s)')
    parser.add_argument('--fofn', type=str, help='file contains the path of all LAS file to be processed in parallel')
    parser.add_argument('--min_len', type=int, default=2500, help="min length of the reads")
    parser.add_argument('--stream', action='store_true', help='stream from LA4Falcon, instead of slurping all at once; can save memory for large data')
    parser.add_argument('--debug', '-g', action='store_true', help="single-threaded, plus other aids to debugging")
    parser.add_argument('--silent', action='store_true', help="suppress cmd reporting on stderr")
    return parser.parse_args(argv[1:])

def main(argv=sys.argv):
    args = parse_args(argv)
    ovlp_stats(**vars(args))
