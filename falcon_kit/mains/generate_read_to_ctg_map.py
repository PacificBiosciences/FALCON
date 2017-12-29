import argparse
import logging
import sys
from ..util import io
from ..fc_asm_graph import AsmGraph

def run(rawread_id_fn, pread_id_fn, sg_edges_list_fn, utg_data_fn, ctg_paths_fn, output_fn):
    read_to_contig_map = output_fn
    pread_did_to_rid = open(pread_id_fn).read().split('\n')
    rid_to_oid = open(rawread_id_fn).read().split('\n')

    asm_G = AsmGraph(sg_edges_list_fn,
                     utg_data_fn,
                     ctg_paths_fn)

    pread_to_contigs = {}

    with open(read_to_contig_map, 'w') as f:
        for ctg in asm_G.ctg_data:
            if ctg[-1] == 'R':
                continue
            ctg_g = asm_G.get_sg_for_ctg(ctg)
            for n in ctg_g.nodes():
                pid = int(n.split(':')[0])

                rid = pread_did_to_rid[pid].split('/')[1]
                rid = int(int(rid) / 10)
                oid = rid_to_oid[rid]
                k = (pid, rid, oid)
                pread_to_contigs.setdefault(k, set())
                pread_to_contigs[k].add(ctg)

        for k in pread_to_contigs:
            pid, rid, oid = k
            for ctg in list(pread_to_contigs[k]):
                print >>f, '%09d %09d %s %s' % (pid, rid, oid, ctg)


class HelpF(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass


def parse_args(argv):
    description = 'Generate read_to_ctg_map from rawread_id file and pread_id file'
    epilog = ''
    parser = argparse.ArgumentParser(
        description=description,
        epilog=epilog,
        formatter_class=HelpF,
    )
    parser.add_argument(
        '--rawread-id-fn',
        required=True,
        help='From TASK_DUMP_RAWREAD_IDS_SCRIPT',
    )
    parser.add_argument(
        '--pread-id-fn',
        required=True,
        help='From TASK_DUMP_PREAD_IDS_SCRIPT',
    )
    parser.add_argument(
        '--sg-edges-list-fn',
        required=True,
        help='From Falcon stage 2-asm-falcon',
    )
    parser.add_argument(
        '--utg-data-fn',
        required=True,
        help='From Falcon stage 2-asm-falcon',
    )
    parser.add_argument(
        '--ctg-paths-fn',
        required=True,
        help='From Falcon stage 2-asm-falcon',
    )
    parser.add_argument(
        '--output-fn',
        required=True,
        help='read-to-ctg-map',
    )
    args = parser.parse_args(argv[1:])
    return args


def main(argv=sys.argv):
    args = parse_args(argv)
    logging.basicConfig(level=logging.INFO)
    run(**vars(args))


if __name__ == '__main__':  # pragma: no cover
    main()
