import sys
import gzip
import os
import os.path
import argparse
from collections import defaultdict, OrderedDict

import pandas as pd
import pysam

import re


def classify_event(event):
    """Return the category (INCL / SKIP / OTHER) for a given event string."""
    ev = str(event).lower()
    if 'incl' in ev:
        return 'INCL'
    elif 'skip' in ev:
        return 'SKIP'
    else:
        return 'OTHER'

def main():

    opts = argparse.ArgumentParser('check for given kmers on fwd and rev read and then assign as incl, skip, or other')

    opts.add_argument('--libname', dest='libname')
    opts.add_argument('--in_fwd_fq', dest='in_fwd_fq')
    opts.add_argument('--in_rev_fq', dest='in_rev_fq')
    opts.add_argument('--in_checktbl', dest='in_checktbl')
    opts.add_argument('--out_counts', dest='out_counts')

    o = opts.parse_args()

    with pysam.FastxFile(o.in_fwd_fq) as fwd_fq, pysam.FastxFile(o.in_rev_fq) as rev_fq:

        kmertbl = pd.read_table(o.in_checktbl)
        kmertbl['category'] = kmertbl['event'].apply(lambda e: classify_event(e))
        
        fwd_mname_re = { r['category'] : re.compile(r['kmer'].replace(' ','').upper()) 
                    for _,r in kmertbl[kmertbl['whichread']=='r1'].iterrows() }

        rev_mname_re = { r['category'] : re.compile(r['kmer'].replace(' ','').upper()) 
                    for _,r in kmertbl[kmertbl['whichread']=='r2'].iterrows() }

        # caveat with other alt_acc and alt_donor is that if the exon length is <10bp then these categories will not be accurate
        # and alt donor/acceptor would just fall in the OTHER category
        mname_ct = { 'INCL':0, 'SKIP':0, 'OTHER: ALT_ACC': 0,'OTHER: ALT_DON': 0, 'OTHER: ALL':0, 'MISMATCH':0}

        for rd1, rd2 in zip(fwd_fq, rev_fq):
            fwd_bc = rd1.comment.strip().split('BC=')[1]
            rev_bc = rd2.comment.strip().split('BC=')[1]
            
            if fwd_bc != rev_bc:
                print(f'WARNING: bc mismatch for {rd1.name}...skipping')
                continue
            
            fwd_seq = rd1.sequence.upper()
            cat_fwd = None
            for nm in fwd_mname_re:
                if fwd_mname_re[nm].search( fwd_seq ):
                    cat_fwd = nm
                    break
            if not cat_fwd: 
                cat_fwd='OTHER'   

            rev_seq = rd2.sequence.upper()
            cat_rev = None
            for nm in rev_mname_re:
                if rev_mname_re[nm].search( rev_seq ):
                    cat_rev = nm
                    break
            if not cat_rev: 
                cat_rev='OTHER' 

            if cat_fwd == 'INCL' and cat_rev == 'INCL':
                mname_ct['INCL'] += 1
            elif cat_fwd =='SKIP' and cat_rev == 'SKIP':
                mname_ct['SKIP'] +=1
            elif cat_fwd=='SKIP' and cat_rev != 'SKIP':
                mname_ct['MISMATCH'] += 1
            elif cat_rev=='SKIP' and cat_fwd != 'SKIP':
                mname_ct['MISMATCH'] += 1
            elif cat_fwd=='INCL' and cat_rev=='OTHER':
                mname_ct['OTHER: ALT_DON'] += 1
                mname_ct['OTHER: ALL'] += 1
            elif cat_rev=='INCL' and cat_fwd=='OTHER':
                mname_ct['OTHER: ALT_ACC'] += 1
                mname_ct['OTHER: ALL'] += 1
            else:
                mname_ct['OTHER: ALL'] += 1
                
    mname_ct = pd.DataFrame( {'event':mname_ct.keys(), 'count':mname_ct.values()} )
    mname_ct['libname'] = o.libname
    mname_ct.to_csv( o.out_counts, sep='\t', index=False )

if __name__=='__main__':
    main()
