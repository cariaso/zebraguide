#!/usr/bin/env python

from abifpy import Trace
import sys
import pdb
from termcolor import colored


def fastafn2seq(fn):
    infh = file(fn)
    topline = infh.readline()
    body = infh.read()
    seq = body.replace('\n','')
    return seq


ab1infn = sys.argv[1]
refinfn = sys.argv[2]
if False:
    aseq = Trace(ab1infn)
    alen = len(aseq.qual_val)
    #pdb.set_trace()
    for i in range(alen):
        viz = '|' * aseq.qual_val[i]
        print'\t'.join([str(i), aseq.seq[i],str(aseq.qual_val[i]),viz])
    for k in sorted(aseq.data):
        print k
    print'-------'
    print aseq.data['baseorder']
    print'-------'
    for i, val in enumerate(aseq.data['raw1']):
        if i > alen: break
        if val > 1:
            viz = '|' * val
        else:
            viz = '#'
        print '%d\t%d\t%s' % (i, val, viz)



import tempfile

tmpdir = tempfile.mkdtemp()
import os
import shutil
baseab1 = os.path.basename(ab1infn)
tempab1 = '%s/%s' % (tmpdir, baseab1)
shutil.copyfile(ab1infn, tempab1)
iupacfn = 'abc'
os.system('~/phred/phred -id %s -m %s' % (tmpdir, iupacfn))

seq = fastafn2seq(iupacfn)
#print seq


known = {
'A' : 'aa',
'M' : 'ac',
'R' : 'ag',
'W' : 'at',
'm' : 'ca',
'C' : 'cc',
'S' : 'cg',
'Y' : 'cT',
'r' : 'ga',
's' : 'gc',
'G' : 'gg',
'K' : 'gt',
'w' : 'ta',
'y' : 'tc',
'k' : 'tg',
'T' : 'tt',
'n' : 'NN',
}

    
#queryfn = 'query1.fa'
#queryfh = open(queryfn, 'w')
#bestfh.write('>aseq1\n')
#for i, achr in enumerate(seq):
#    alt1 = known[achr][0]
#    alt2 = known[achr][1]
    #print i, achr, alt1, alt2

#    if achr == 'n':
#        queryfh.write('[ATCG]')
#    elif alt1 == alt2:
#        queryfh.write(alt1.upper())
#    else:
#        queryfh.write('[%s]' % known[achr].upper())
#queryfh.close()





from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist

import pdb
#pdb.set_trace()

SC_perfect = 10
SC_1of2 = 8
SC_miss = -5

ref_gap_open = -10
ref_gap_extend = -5
query_gap_open = -10
query_gap_extend = -5


matrix = {

('A','A'):SC_perfect,
('A','M'):SC_1of2,
('A','R'):SC_1of2,
('A','W'):SC_1of2,
('A','C'):SC_miss,
('A','S'):SC_miss,
('A','Y'):SC_miss,
('A','G'):SC_miss,
('A','K'):SC_miss,
('A','T'):SC_miss,
('A','N'):0,


('M','A'):SC_1of2,
('M','M'):SC_perfect,
('M','R'):SC_miss,
('M','W'):SC_miss,
('M','C'):SC_1of2,
('M','S'):SC_miss,
('M','Y'):SC_miss,
('M','G'):SC_miss,
('M','K'):SC_miss,
('M','T'):SC_miss,
('M','N'):0,


('R','A'):SC_1of2,
('R','M'):SC_miss,
('R','R'):SC_perfect,
('R','W'):SC_miss,
('R','C'):SC_miss,
('R','S'):SC_miss,
('R','Y'):SC_miss,
('R','G'):SC_1of2,
('R','K'):SC_miss,
('R','T'):SC_miss,
('R','N'):0,


('W','A'):SC_1of2,
('W','M'):SC_miss,
('W','R'):SC_miss,
('W','W'):SC_perfect,
('W','C'):SC_miss,
('W','S'):SC_miss,
('W','Y'):SC_miss,
('W','G'):SC_miss,
('W','K'):SC_miss,
('W','T'):SC_1of2,
('W','N'):0,


('C','A'):SC_miss,
('C','M'):SC_1of2,
('C','R'):SC_miss,
('C','W'):SC_miss,
('C','C'):SC_perfect,
('C','S'):SC_1of2,
('C','Y'):SC_1of2,
('C','G'):SC_miss,
('C','K'):SC_miss,
('C','T'):SC_miss,
('C','N'):0,


('S','A'):SC_miss,
('S','M'):SC_miss,
('S','R'):SC_miss,
('S','W'):SC_miss,
('S','C'):SC_1of2,
('S','S'):SC_perfect,
('S','Y'):SC_miss,
('S','G'):SC_1of2,
('S','K'):SC_miss,
('S','T'):SC_miss,
('S','N'):0,


('Y','A'):SC_miss,
('Y','M'):SC_miss,
('Y','R'):SC_miss,
('Y','W'):SC_miss,
('Y','C'):SC_1of2,
('Y','S'):SC_miss,
('Y','Y'):SC_perfect,
('Y','G'):SC_miss,
('Y','K'):SC_miss,
('Y','T'):SC_1of2,
('Y','N'):0,


('G','A'):SC_miss,
('G','M'):SC_miss,
('G','R'):SC_1of2,
('G','W'):SC_miss,
('G','C'):SC_miss,
('G','S'):SC_1of2,
('G','Y'):SC_miss,
('G','G'):SC_perfect,
('G','K'):SC_1of2,
('G','T'):SC_miss,
('G','N'):0,


('K','A'):SC_miss,
('K','M'):SC_miss,
('K','R'):SC_miss,
('K','W'):SC_miss,
('K','C'):SC_miss,
('K','S'):SC_miss,
('K','Y'):SC_miss,
('K','G'):SC_1of2,
('K','K'):SC_perfect,
('K','T'):SC_1of2,
('K','N'):0,


('T','A'):SC_miss,
('T','M'):SC_miss,
('T','R'):SC_miss,
('T','W'):SC_1of2,
('T','C'):SC_miss,
('T','S'):SC_miss,
('T','Y'):SC_1of2,
('T','G'):SC_miss,
('T','K'):SC_1of2,
('T','T'):SC_perfect,
('T','N'):0,


('N','A'):0,
('N','M'):0,
('N','R'):0,
('N','W'):0,
('N','C'):0,
('N','S'):0,
('N','Y'):0,
('N','G'):0,
('N','K'):0,
('N','T'):0,
('N','N'):0,

}


refseq = fastafn2seq(refinfn)


step = 50000

print len(refseq)
smallseq = seq.upper()
steps = range(0, len(refseq),step)


for start in steps:
    bigseq = refseq.upper()[start:start+step]
    print start,'@',len(bigseq), 'x', len(smallseq)

    print
    print 'searching for'
    print smallseq
    #alns = pairwise2.align.globalds(bigseq, smallseq, matrix, gap_open, gap_extend)
    alns = pairwise2.align.localdd(bigseq, smallseq, matrix, ref_gap_open, ref_gap_extend, query_gap_open, query_gap_extend)


    from Bio.pairwise2 import format_alignment
    for a in alns[:1]:
        #print(format_alignment(*a))

        aln_top, aln_bot, aln_score, begin, end = a
        align_chars = []
        alt_chars = []
        for t,b in zip(aln_top, aln_bot):
            score = matrix.get((t,b), 0)
            if score > 1:
                compare = '='
                couldbe = known[b].upper()
                if t == couldbe[1]:
                    alts = couldbe[0]
                    alt_iupac = couldbe[0]
                else:
                    alts = couldbe
                    alt_iupac = b
                alt_chars.append(alt_iupac)
                #print 'ref=%s%s\tbest=%s:%s\talts=%s\talt_iupac=%s' % (t, compare, b, couldbe, alts, alt_iupac)
            else:
                compare = ' '
                #print 'ref=%s%s\tbest=%s' % (t, compare, b)

            align_chars.append(compare)

        print
        print 'colored alignment'
        for i, c in enumerate(aln_top):
            perfect = ''.join(align_chars[i-2:i+3]) == '====='
            good = align_chars[i] == '='
            if perfect:
                sys.stdout.write(colored(c, 'red',attrs=['bold']))
            elif good:
                sys.stdout.write(colored(c, 'grey',attrs=['bold']))
            else:
                sys.stdout.write(colored(c, 'white'))
        print
        print
        print 'reference/hits/query'
        print aln_top
        print ''.join(align_chars)
        print aln_bot
        print 'score=',aln_score, 'begin=',begin, 'end=',end
        print '\n'*4


        altseq = ''.join(alt_chars)
        print 'altseq'
        print altseq

        alt_alns = pairwise2.align.localdd(bigseq, altseq, matrix, ref_gap_open, ref_gap_extend, query_gap_open, query_gap_extend)


        for a2 in alt_alns[:1]:

            alt_aln_top, alt_aln_bot, alt_score, alt_begin, alt_end = a2
            alt_align_chars = []
            for alt_t,alt_b in zip(alt_aln_top, alt_aln_bot):
                alt_base_score = matrix.get((alt_t,alt_b), 0)
                if alt_base_score > 1:
                    alt_compare = '='
                else:
                    alt_compare = ' '
                #print '%s%s%s' % (alt_t, alt_compare, alt_b)

                alt_align_chars.append(alt_compare)
            #print alt_aln_top
            print ''.join(alt_align_chars)
            print alt_aln_bot
            print 'score=',alt_score, 'begin=',alt_begin, 'end=',alt_end
            print '\n'*4

        print altseq



sys.exit()
for i, rest in enumerate(zip(smallseq, altseq)):
    a1 = smallseq[i]
    a2 = altseq[i]
    #a1, a2 = rest
    print i, aln_top[i], align_chars[i], aln_bot[i], a1, known[a1], a2, known[a2], ';', alt_align_chars[i], alt_aln_bot[i],
    if a1 in 'ATCG':
        print '*'
    elif a1 != a2:
        print '#'
    else:
        print ' '



sys.exit()



cmd = '~/sim4.2012-10-10/sim4 %s %s' % (iupacfn, refinfn)
print cmd
os.system(cmd)
