import sys

import os
import os.path
import argparse
from collections import defaultdict, OrderedDict
import re
import pysam

def main():
    opts = argparse.ArgumentParser()

    opts.add_argument('--in_bam', dest='in_bam' )
    opts.add_argument('--clone_id_tag', default='BC', dest='clone_id_tag' )
    opts.add_argument('--umi_tag', default='BX', dest='umi_tag' )
    opts.add_argument('--out_bam', dest='out_bam' )
    
    o = opts.parse_args()
    
    bam = pysam.AlignmentFile( o.in_bam, 'rb' )
    bamout = pysam.AlignmentFile( o.out_bam, 'wb', template=bam )

    rebc=re.compile(r'_BC=([ACTGN]+)')
    reumi=re.compile(r'_UMI=([ACTGN]+)')

    for l in bam:

        rn = l.query_name

        mbc = rebc.search(rn)
        mumi = reumi.search(rn)

        if mbc:
            l.tags = l.tags + [(o.clone_id_tag, mbc.groups(1)[0])]

        if mumi:
            l.tags = l.tags + [(o.umi_tag, mumi.groups(1)[0])]

        bamout.write(l)


if __name__ == '__main__':                
    main()

