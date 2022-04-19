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

    opts = argparse.ArgumentParser(description='extract upper case part (barcode) from reverse read fastq and label both the fwd and rev reads w/ it')

    opts.add_argument('--in_bc_fq', dest='in_bc_fq')
    opts.add_argument('--in_splice_fq', dest='in_splice_fq')
    opts.add_argument('--out_bc_fq', dest='out_bc_fq')
    opts.add_argument('--out_splice_fq', dest='out_splice_fq')

    opts.add_argument('--splice_umi_len',type=int,default=0,dest='splice_umi_len')

    opts.add_argument('--dont_rc_bc', default=True, action='store_false', dest='revcomp_bc')

    o = opts.parse_args()

    fin_bc = zopen(o.in_bc_fq,'rt')
    fin_splice = zopen(o.in_splice_fq,'rt')    
    fout_bc = zopen(o.out_bc_fq,'wt')
    fout_splice = zopen(o.out_splice_fq,'wt')

    re_ucb = re.compile('([ACGTN]+)')

    dorcbc = o.revcomp_bc

    doumi = o.splice_umi_len>0
    umil = o.splice_umi_len

    i=0

    for bcl1 in fin_bc:

        i+=1
        if i%500000==0:
            print('%d reads...'%i)

        bcl2,_,bcl4=fin_bc.readline(),fin_bc.readline(),fin_bc.readline()
        insl1,insl2,_,insl4=fin_splice.readline(),fin_splice.readline(),fin_splice.readline(),fin_splice.readline()

        if doumi:
            umi = insl2[:umil]
            insl2 = insl2[umil:]
            insl4 = insl4[umil:]

            if len(insl2)<20 or len(bcl2)<2:
                continue

        m_ucb = re_ucb.search(bcl2)

        if m_ucb:
            bc = m_ucb.groups()[0]

            if dorcbc:
                bc = str(Bio.Seq.Seq(bc).reverse_complement())

            fout_bc.write(bcl1.strip()+' BC='+bc)
            if doumi:
                fout_bc.write('_UMI='+umi)
            fout_bc.write('\n'+bcl2)
            fout_bc.write('+\n')
            fout_bc.write(bcl4)

            fout_splice.write(insl1.strip()+' BC='+bc)
            if doumi:
                fout_splice.write('_UMI='+umi)
            fout_splice.write('\n'+insl2)
            fout_splice.write('+\n')
            fout_splice.write(insl4)




if __name__=='__main__':
    main()
