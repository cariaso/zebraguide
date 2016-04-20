#!/usr/bin/env python

from abifpy import Trace
import sys
import pdb
ab1infn = sys.argv[1]
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
infh = file(iupacfn)
topline = infh.readline()
body = infh.read()
seq = body.replace('\n','')
print seq


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
'n' : '??',
}

    
for i, achr in enumerate(seq):
    alt1 = known[achr][0]
    alt2 = known[achr][1]
    print i, achr, alt1, alt2
