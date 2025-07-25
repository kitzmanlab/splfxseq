#!/usr/bin/env python3

import pandas as pd
import argparse
from splshared.ssshared import zopen

import re

bc_re = re.compile(r'BC=([ATCGN]+)')

def filter_fq_write_unk(
        bc_tbl: pd.DataFrame,
        fn_fq_in_fwd: str,  
        fn_fq_in_rev: str,  
        fn_fq_out_known_fwd: str, 
        fn_fq_out_known_rev: str, 
        fn_fq_out_unk_fwd: str = None, 
        fn_fq_out_unk_rev: str = None, 
):
    fq_in_fwd = zopen(fn_fq_in_fwd, 'rt')
    fq_in_rev = zopen(fn_fq_in_rev, 'rt')

    filout_known_fwd = zopen(fn_fq_out_known_fwd, 'wt')
    filout_known_rev = zopen(fn_fq_out_known_rev, 'wt')
    filout_unk_fwd = zopen(fn_fq_out_unk_fwd, 'wt') if fn_fq_out_unk_fwd else None
    filout_unk_rev = zopen(fn_fq_out_unk_rev, 'wt') if fn_fq_out_unk_rev else None

    with zopen(fn_fq_in_fwd, 'rt') as fq_in_fwd:

        lfwd = fq_in_fwd.readline()
        lrev = fq_in_rev.readline()

        while lfwd: 
            m = bc_re.search(lfwd)

            bc = m.group(1)
            if bc in bc_tbl.index:
                filout_known_fwd.write(lfwd)
                filout_known_fwd.write(fq_in_fwd.readline())
                filout_known_fwd.write(fq_in_fwd.readline())
                filout_known_fwd.write(fq_in_fwd.readline())

                filout_known_rev.write(lrev)
                filout_known_rev.write(fq_in_rev.readline())
                filout_known_rev.write(fq_in_rev.readline())
                filout_known_rev.write(fq_in_rev.readline())
            elif filout_unk_fwd:
                filout_unk_fwd.write(lfwd)
                filout_unk_fwd.write(fq_in_fwd.readline())
                filout_unk_fwd.write(fq_in_fwd.readline())
                filout_unk_fwd.write(fq_in_fwd.readline())

                filout_unk_rev.write(lrev)
                filout_unk_rev.write(fq_in_rev.readline())
                filout_unk_rev.write(fq_in_rev.readline())
                filout_unk_rev.write(fq_in_rev.readline())
            
            lrev = fq_in_rev.readline()
            lfwd = fq_in_fwd.readline()
    
    filout_known_fwd.close()
    filout_known_rev.close()
    if filout_unk_fwd:
        filout_unk_fwd.close()
    if filout_unk_rev:
        filout_unk_rev.close()


def main():
    parser = argparse.ArgumentParser(description='Compute overlap of barcodes between a barcode pairing table and a set of MPSA RNA-seq barcode histograms')

    parser.add_argument('--pairing_tbl', type=str, required=True,
                      help='path to pairing table')
    
    parser.add_argument('--fq_fwd_in', type=str, required=True,
                      help='rnaseq forward trimmed fq path, read names should have BC=... filled out')
    
    parser.add_argument('--fq_rev_in', type=str, required=True,
                      help='rnaseq reverse trimmed fq path, read names should have BC=... filled out')

    parser.add_argument('--known_fq_fwd_out', type=str, required=True,
                      help='path to output fastq, forward reads from known barcodes')
    
    parser.add_argument('--known_fq_rev_out', type=str, required=True,
                      help='path to output fastq, forward reads from known barcodes')

    parser.add_argument('--unk_fq_fwd_out', type=str, default=None, 
                      help='path to output fastq, forward reads from unknown barcodes')
    
    parser.add_argument('--unk_fq_rev_out', type=str, default=None, 
                      help='path to output fastq, reverse reads from unknown barcodes')

    args = parser.parse_args()
    
    # Read input data
    pairtbl = pd.read_table(args.pairing_tbl).set_index( 'readgroupid' )
    print(f'{args.pairing_tbl} has {pairtbl.shape[0]} rows')
    samp_tbl = pd.read_table(args.rnaseq_samp_tbl)
    
    filter_fq_write_unk( 
        pairtbl,
        args.fq_fwd_in,
        args.fq_rev_in,
        args.known_fq_fwd_out,
        args.known_fq_rev_out,
        args.unk_fq_fwd_out,
        args.unk_fq_rev_out
    )


if __name__ == '__main__':
    main() 