#!/usr/bin/env python

import argparse
import Bio.SeqIO
import sys

def get_args():
    parser = argparse.ArgumentParser(
        description="""Make the BoulderIO formatter config file for primer3 to direct guide sequences against a reference.

Typical usage:

./zebraguide.py sample2.fasta | primer3_core -format_output

""",
        epilog="cariaso@gmail.com Copyright 2016"
       )

    parser.add_argument('fasta', 
                        help='fasta file of references and guides')
    parser.add_argument('--size-limit', type=int, default=20,
                        help='cutoff length of references vs guides')
    parser.add_argument('--out',
                        help='output filename, defaults to stdout')
    parser.add_argument('--add', action='append', dest='primer3tags',
                        default=[],
                        help='Add these tags to each record')


    args = parser.parse_args()
    return args


def analyze(out, guide, ref, args=None):
    guideseq = guide.seq.tostring().upper()
    refseq = ref.seq.tostring().upper()

    guidestartpos = refseq.find(guideseq)

    if guidestartpos == -1:
        out.write('COMMENT=%s not found in %s\n' % (guide.name, ref.name))
        out.write('=\n')
        return

    guidendpos = guidestartpos + len(guide)
    out.write('SEQUENCE_ID=%s-%s\n' % (guide.name, ref.name))
    out.write('SEQUENCE_TEMPLATE=%s\n' % ref.seq.tostring())
    out.write('SEQUENCE_TARGET=%s,%s\n' % (guidestartpos, guidendpos))
    #out.write('PRIMER_THERMODYNAMIC_PARAMETERS_PATH=/usr/local/Cellar/primer3/2.3.6/share/primer3/primer3_config/\n')
    out.write('PRIMER_EXPLAIN_FLAG=1\n')
    if args:
        for tag in args.primer3tags:
            out.write('%s\n' % tag)

    # out.write('PRIMER_FIRST_BASE_INDEX=1\n')
    # out.write('PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT=1\n')
    # out.write('PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT=0\n')
    # out.write('PRIMER_PICK_LEFT_PRIMER=1\n')
    # out.write('PRIMER_PICK_INTERNAL_OLIGO=0\n')
    # out.write('PRIMER_PICK_RIGHT_PRIMER=1\n')
    # out.write('PRIMER_LIBERAL_BASE=1\n')
    # out.write('PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS=0\n')
    # out.write('PRIMER_LOWERCASE_MASKING=0\n')
    # out.write('PRIMER_PICK_ANYWAY=1\n')
    # out.write('PRIMER_EXPLAIN_FLAG=1\n')
    # out.write('PRIMER_TASK=generic\n')
    # out.write('PRIMER_MIN_QUALITY=0\n')
    # out.write('PRIMER_MIN_END_QUALITY=0\n')
    # out.write('PRIMER_QUALITY_RANGE_MIN=0\n')
    # out.write('PRIMER_QUALITY_RANGE_MAX=100\n')
    # out.write('PRIMER_MIN_SIZE=18\n')
    # out.write('PRIMER_OPT_SIZE=20\n')
    # out.write('PRIMER_MAX_SIZE=23\n')
    # out.write('PRIMER_MIN_TM=57.0\n')
    # out.write('PRIMER_OPT_TM=59.0\n')
    # out.write('PRIMER_MAX_TM=62.0\n')
    # out.write('PRIMER_PAIR_MAX_DIFF_TM=5.0\n')
    # out.write('PRIMER_TM_FORMULA=1\n')
    # out.write('PRIMER_PRODUCT_MIN_TM=-1000000.0\n')
    # out.write('PRIMER_PRODUCT_OPT_TM=0.0\n')
    # out.write('PRIMER_PRODUCT_MAX_TM=1000000.0\n')
    # out.write('PRIMER_MIN_GC=30.0\n')
    # out.write('PRIMER_OPT_GC_PERCENT=50.0\n')
    # out.write('PRIMER_MAX_GC=70.0\n')
    # out.write('PRIMER_PRODUCT_SIZE_RANGE=150-250 100-300 301-400 401-500 501-600 601-700 701-850 851-1000\n')
    # out.write('PRIMER_NUM_RETURN=5\n')
    # out.write('PRIMER_MAX_END_STABILITY=9.0\n')
    # out.write('PRIMER_MAX_LIBRARY_MISPRIMING=12.00\n')
    # out.write('PRIMER_PAIR_MAX_LIBRARY_MISPRIMING=20.00\n')
    # out.write('PRIMER_MAX_TEMPLATE_MISPRIMING_TH=40.00\n')
    # out.write('PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING_TH=70.00\n')
    # out.write('PRIMER_MAX_SELF_ANY_TH=45.0\n')
    # out.write('PRIMER_MAX_SELF_END_TH=35.0\n')
    # out.write('PRIMER_PAIR_MAX_COMPL_ANY_TH=45.0\n')
    # out.write('PRIMER_PAIR_MAX_COMPL_END_TH=35.0\n')
    # out.write('PRIMER_MAX_HAIRPIN_TH=24.0\n')
    # out.write('PRIMER_MAX_TEMPLATE_MISPRIMING=12.00\n')
    # out.write('PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING=24.00\n')
    # out.write('PRIMER_MAX_SELF_ANY=8.00\n')
    # out.write('PRIMER_MAX_SELF_END=3.00\n')
    # out.write('PRIMER_PAIR_MAX_COMPL_ANY=8.00\n')
    # out.write('PRIMER_PAIR_MAX_COMPL_END=3.00\n')
    # out.write('PRIMER_MAX_NS_ACCEPTED=0\n')
    # out.write('PRIMER_MAX_POLY_X=4\n')
    # out.write('PRIMER_INSIDE_PENALTY=-1.0\n')
    # out.write('PRIMER_OUTSIDE_PENALTY=0\n')
    # out.write('PRIMER_GC_CLAMP=0\n')
    # out.write('PRIMER_MAX_END_GC=5\n')
    # out.write('PRIMER_MIN_LEFT_THREE_PRIME_DISTANCE=3\n')
    # out.write('PRIMER_MIN_RIGHT_THREE_PRIME_DISTANCE=3\n')
    # out.write('PRIMER_MIN_5_PRIME_OVERLAP_OF_JUNCTION=7\n')
    # out.write('PRIMER_MIN_3_PRIME_OVERLAP_OF_JUNCTION=4\n')
    # out.write('PRIMER_SALT_MONOVALENT=50.0\n')
    # out.write('PRIMER_SALT_CORRECTIONS=1\n')
    # out.write('PRIMER_SALT_DIVALENT=1.5\n')
    # out.write('PRIMER_DNTP_CONC=0.6\n')
    # out.write('PRIMER_DNA_CONC=50.0\n')
    # out.write('PRIMER_SEQUENCING_SPACING=500\n')
    # out.write('PRIMER_SEQUENCING_INTERVAL=250\n')
    # out.write('PRIMER_SEQUENCING_LEAD=50\n')
    # out.write('PRIMER_SEQUENCING_ACCURACY=20\n')
    # out.write('PRIMER_WT_SIZE_LT=1.0\n')
    # out.write('PRIMER_WT_SIZE_GT=1.0\n')
    # out.write('PRIMER_WT_TM_LT=1.0\n')
    # out.write('PRIMER_WT_TM_GT=1.0\n')
    # out.write('PRIMER_WT_GC_PERCENT_LT=0.0\n')
    # out.write('PRIMER_WT_GC_PERCENT_GT=0.0\n')
    # out.write('PRIMER_WT_SELF_ANY_TH=0.0\n')
    # out.write('PRIMER_WT_SELF_END_TH=0.0\n')
    # out.write('PRIMER_WT_HAIRPIN_TH=0.0\n')
    # out.write('PRIMER_WT_TEMPLATE_MISPRIMING_TH=0.0\n')
    # out.write('PRIMER_WT_SELF_ANY=0.0\n')
    # out.write('PRIMER_WT_SELF_END=0.0\n')
    # out.write('PRIMER_WT_TEMPLATE_MISPRIMING=0.0\n')
    # out.write('PRIMER_WT_NUM_NS=0.0\n')
    # out.write('PRIMER_WT_LIBRARY_MISPRIMING=0.0\n')
    # out.write('PRIMER_WT_SEQ_QUAL=0.0\n')
    # out.write('PRIMER_WT_END_QUAL=0.0\n')
    # out.write('PRIMER_WT_POS_PENALTY=0.0\n')
    # out.write('PRIMER_WT_END_STABILITY=0.0\n')
    # out.write('PRIMER_PAIR_WT_PRODUCT_SIZE_LT=0.0\n')
    # out.write('PRIMER_PAIR_WT_PRODUCT_SIZE_GT=0.0\n')
    # out.write('PRIMER_PAIR_WT_PRODUCT_TM_LT=0.0\n')
    # out.write('PRIMER_PAIR_WT_PRODUCT_TM_GT=0.0\n')
    # out.write('PRIMER_PAIR_WT_COMPL_ANY_TH=0.0\n')
    # out.write('PRIMER_PAIR_WT_COMPL_END_TH=0.0\n')
    # out.write('PRIMER_PAIR_WT_TEMPLATE_MISPRIMING_TH=0.0\n')
    # out.write('PRIMER_PAIR_WT_COMPL_ANY=0.0\n')
    # out.write('PRIMER_PAIR_WT_COMPL_END=0.0\n')
    # out.write('PRIMER_PAIR_WT_TEMPLATE_MISPRIMING=0.0\n')
    # out.write('PRIMER_PAIR_WT_DIFF_TM=0.0\n')
    # out.write('PRIMER_PAIR_WT_LIBRARY_MISPRIMING=0.0\n')
    # out.write('PRIMER_PAIR_WT_PR_PENALTY=1.0\n')
    # out.write('PRIMER_PAIR_WT_IO_PENALTY=0.0\n')
    # out.write('PRIMER_INTERNAL_MIN_SIZE=18\n')
    # out.write('PRIMER_INTERNAL_OPT_SIZE=20\n')
    # out.write('PRIMER_INTERNAL_MAX_SIZE=27\n')
    # out.write('PRIMER_INTERNAL_MIN_TM=57.0\n')
    # out.write('PRIMER_INTERNAL_OPT_TM=60.0\n')
    # out.write('PRIMER_INTERNAL_MAX_TM=63.0\n')
    # out.write('PRIMER_INTERNAL_MIN_GC=20.0\n')
    # out.write('PRIMER_INTERNAL_OPT_GC_PERCENT=50.0\n')
    # out.write('PRIMER_INTERNAL_MAX_GC=80.0\n')
    # out.write('PRIMER_INTERNAL_MAX_SELF_ANY_TH=47.00\n')
    # out.write('PRIMER_INTERNAL_MAX_SELF_END_TH=47.00\n')
    # out.write('PRIMER_INTERNAL_MAX_HAIRPIN_TH=47.00\n')
    # out.write('PRIMER_INTERNAL_MAX_SELF_ANY=12.00\n')
    # out.write('PRIMER_INTERNAL_MAX_SELF_END=12.00\n')
    # out.write('PRIMER_INTERNAL_MIN_QUALITY=0\n')
    # out.write('PRIMER_INTERNAL_MAX_NS_ACCEPTED=0\n')
    # out.write('PRIMER_INTERNAL_MAX_POLY_X=5\n')
    # out.write('PRIMER_INTERNAL_MAX_LIBRARY_MISHYB=12.00\n')
    # out.write('PRIMER_INTERNAL_SALT_MONOVALENT=50.0\n')
    # out.write('PRIMER_INTERNAL_DNA_CONC=50.0\n')
    # out.write('PRIMER_INTERNAL_SALT_DIVALENT=1.5\n')
    # out.write('PRIMER_INTERNAL_DNTP_CONC=0.0\n')
    # out.write('PRIMER_INTERNAL_WT_SIZE_LT=1.0\n')
    # out.write('PRIMER_INTERNAL_WT_SIZE_GT=1.0\n')
    # out.write('PRIMER_INTERNAL_WT_TM_LT=1.0\n')
    # out.write('PRIMER_INTERNAL_WT_TM_GT=1.0\n')
    # out.write('PRIMER_INTERNAL_WT_GC_PERCENT_LT=0.0\n')
    # out.write('PRIMER_INTERNAL_WT_GC_PERCENT_GT=0.0\n')
    # out.write('PRIMER_INTERNAL_WT_SELF_ANY_TH=0.0\n')
    # out.write('PRIMER_INTERNAL_WT_SELF_END_TH=0.0\n')
    # out.write('PRIMER_INTERNAL_WT_HAIRPIN_TH=0.0\n')
    # out.write('PRIMER_INTERNAL_WT_SELF_ANY=0.0\n')
    # out.write('PRIMER_INTERNAL_WT_SELF_END=0.0\n')
    # out.write('PRIMER_INTERNAL_WT_NUM_NS=0.0\n')
    # out.write('PRIMER_INTERNAL_WT_LIBRARY_MISHYB=0.0\n')
    # out.write('PRIMER_INTERNAL_WT_SEQ_QUAL=0.0\n')
    # out.write('PRIMER_INTERNAL_WT_END_QUAL=0.0\n')
    out.write('=\n')

    
def main():
    args = get_args()

    infastafh = open(args.fasta)
    seq_stream = Bio.SeqIO.parse(infastafh, 'fasta')

    if args.out:
        outfh = file(args.out,'w')
    else:
        outfh = sys.stdout

    ref = None
    for seqobj in seq_stream:
        if len(seqobj) > args.size_limit:
            ref = seqobj
        else:
            analyze(outfh, seqobj, ref, args)

if __name__ == '__main__':
    main()
