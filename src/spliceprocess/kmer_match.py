import sys

import os
import os.path
import argparse
from collections import defaultdict, OrderedDict

import pandas as pd
import pysam

import re

def zopen(fn,mode,encoding='utf-8'):
    # return gzip.GzipFile(fn,mode=mode) if fn.endswith('gz') else open(fn,mode=mode)
    return gzip.open(fn,mode=mode,encoding=encoding) if fn.endswith('gz') else open(fn,mode=mode,encoding=encoding)

def main():

    opts = argparse.ArgumentParser('very slow and simplistic, just check for the given kmers, and count each read towards the first one found')

    opts.add_argument('--libname', dest='libname')
    opts.add_argument('--in_fq', dest='in_fq')
    opts.add_argument('--in_checktbl', dest='in_checktbl')
    opts.add_argument('--out_counts', dest='out_counts')

    o = opts.parse_args()

    fq = pysam.FastxFile(o.in_fq)
    kmertbl = pd.read_table(o.in_checktbl)

    mname_re = { r['event'] : re.compile(r['kmer'].replace(' ','').upper()) 
                for _,r in kmertbl.iterrows() }

    mname_ct = { nm:0 for nm in mname_re }
    mname_ct['other'] = 0

    for rd in fq:
        match = 0
        for nm in mname_re:
            if mname_re[nm].search( rd.sequence ):
                mname_ct[nm]+=1
                match = 1
                break
        if not match: mname_ct['other']+=1            
    
    mname_ct = pd.DataFrame( {'event':mname_ct.keys(), 'count':mname_ct.values()} )
    mname_ct['libname'] = o.libname
    mname_ct.to_csv( o.out_counts, sep='\t', index=False )

if __name__=='__main__':
    main()
