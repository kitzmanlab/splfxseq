import sys

import os
import os.path
import argparse
from collections import defaultdict, OrderedDict

import Bio.Seq

import gzip

import re

def zopen(fn,mode,encoding='utf-8'):
    # return gzip.GzipFile(fn,mode=mode) if fn.endswith('gz') else open(fn,mode=mode)
    return gzip.open(fn,mode=mode,encoding=encoding) if fn.endswith('gz') else open(fn,mode=mode,encoding=encoding)

def main():

    opts = argparse.ArgumentParser(description='extract bc from read names and turn into its own fsatq')

    opts.add_argument('--in_fq', dest='in_fq')
    opts.add_argument('--out_bc_fq', dest='out_bc_fq')
    
    o = opts.parse_args()

    fin_bc = zopen(o.in_fq,'rt')
    fout_bc = zopen(o.out_bc_fq,'wt')

    re_bcn = re.compile('BC=([ACGTN]+)')

    i=0

    for bcl1 in fin_bc:
        i+=1
        if i%500000==0:
            print('%d reads...'%i)

        fin_bc.readline(); fin_bc.readline() ;fin_bc.readline()
        
        m_ucb = re_bcn.search(bcl1)

        if m_ucb:
            bc = m_ucb.groups()[0]

            fout_bc.write(bcl1)
            fout_bc.write(bc)
            fout_bc.write('\n+\n%s\n'%( 'F'*len(bc) ))



if __name__=='__main__':
    main()
