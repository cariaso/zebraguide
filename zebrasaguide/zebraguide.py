#!/usr/bin/env python

#import biopython
import argparse
import sys
from Bio import SeqIO

def get_args():
    parser = argparse.ArgumentParser(
        description="make the config file for primer3 to direct guide sequences against a reference",
        epilog="cariaso@gmail.com Copyright 2016"
       )

    parser.add_argument('fasta', 
                        help='fasta file of references and guides')
    parser.add_argument('--size-limit', type=int, default=20,
                        help='cutoff length of references vs guides')

    args = parser.parse_args()
    return args


def analyze(guide, ref):
    print guide.name, len(guide), guide.seq.tostring()
    print ref.name, len(ref), ref.seq.tostring()
    print '----------'
    
def main():
    args = get_args()
    print(args)


    infastafh = open(args.fasta)
    fasta_sequences = SeqIO.parse(infastafh, 'fasta')
    #with open(output_file) as out_file:
    if True:
        ref = None
        guide = None
        for seqobj in fasta_sequences:
            #name, sequence = seqobj.id, seqobj.seq.tostring()
            if len(seqobj) > args.size_limit:
                ref = seqobj
            else:
                guide = seqobj
                analyze(guide, ref)

    

if __name__ == '__main__':
    main()
