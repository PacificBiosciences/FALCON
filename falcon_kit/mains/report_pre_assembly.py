from .. import stats_preassembly
import argparse
import json
import logging

log = logging.getLogger()

def do_report(db, preads_fofn, genome_length, length_cutoff, out):
    kwds = dict(
        i_preads_fofn_fn=preads_fofn,
        i_raw_reads_db_fn=db,
        genome_length=genome_length,
        length_cutoff=length_cutoff,
    )
    report_dict = stats_preassembly.calc_dict(**kwds)
    content = json.dumps(report_dict, sort_keys=True, indent=4, separators=(',', ': '))
    open(out, 'w').write(content)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--genome-length',
        type=int,
        required=True,
        help='Estimated number of bases in the full genome haplotype.')
    parser.add_argument('--length-cutoff',
        type=int,
        required=True,
        help='Minimum length of any seed read.')
    parser.add_argument('--db',
        required=True,
        help='Path to raw_reads.db (dazzler DB)')
    parser.add_argument('--preads-fofn',
        required=True,
        help='Path to FOFN of preads fasta files.')
    parser.add_argument('--out',
        required=True,
        help='Path to JSON output file.')
    ARGS = parser.parse_args()
    do_report(**vars(ARGS))


if __name__ == "__main__":
    logging.basicConfig()
    log.setLevel(logging.DEBUG)
    main()
