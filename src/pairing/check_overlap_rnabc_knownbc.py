#!/usr/bin/env python3

import pandas as pd
import altair as alt
from typing import List, Tuple, Dict, Set, Callable
from itertools import product
import argparse
from pathlib import Path
import pysam
from splshared.var_mapping import GeneModelMapper
from splshared.ssshared import revComp
import numpy as np
import matplotlib.pyplot as plt

def read_starcode_histo(fpath: str) -> pd.DataFrame:
    bchisto = pd.read_table( fpath, header=None, usecols=[0,1] )
    bchisto.columns=['barcode','count']
    bchisto = bchisto.set_index('barcode')['count']
    return bchisto

def read_mutli_rnaseq_histos(
    samp_tbl: pd.DataFrame,
    col_histo = 'bc_histo',
) -> pd.DataFrame:
    
    jointhisto = []
    for _, row in samp_tbl.iterrows():
        fn=row[col_histo]
        bchisto = read_starcode_histo(fn)
        print(f'{fn} has {bchisto.shape[0]} barcodes')
        bchisto.name = row['libname']
        jointhisto.append( bchisto )

    jointhisto = pd.concat( jointhisto, axis=1 )
    
    return jointhisto

def indivsamp_overlap_knownbc(
    jointhisto: pd.DataFrame,
    pairtbl: pd.DataFrame,
    lmin_read: List[int],
    # revcomp: bool = False,
) -> pd.DataFrame:
    
    pairbc = set(pairtbl.index)
    # if revcomp:
    #     pairbc = set( [ revComp(bc) for bc in pairbc ] )
    
    tbl_out = { cn:[] for cn in 'libname,min_read,n_bc_sample,n_bc_pairing,n_bc_overlap,n_read_sample,n_read_overlap'.split(',') }

    for samp_name in jointhisto.columns:
        for min_read in lmin_read:
            li_sampbc = jointhisto[samp_name]>=min_read 
            li_bcknown = jointhisto.index.isin( pairbc )
            
            tbl_out['libname'].append( samp_name )
            tbl_out['min_read'].append( min_read )
            tbl_out['n_bc_sample'].append( (li_sampbc).sum() )
            tbl_out['n_bc_pairing'].append( len(pairbc) )
            tbl_out['n_bc_overlap'].append( (li_sampbc & li_bcknown).sum() )

            tbl_out['n_read_sample'].append( jointhisto.loc[ li_sampbc, samp_name ].sum() )
            tbl_out['n_read_overlap'].append( jointhisto.loc[ (li_sampbc & li_bcknown), samp_name ].sum() )

    tbl_out = pd.DataFrame(tbl_out)
    tbl_out['n_read_overlap']=tbl_out['n_read_overlap'].astype(int)
    tbl_out['n_read_sample']=tbl_out['n_read_sample'].astype(int)
    return tbl_out

def plot_upset_overlaps(
    jointhisto: pd.DataFrame,
    pairtbl: pd.DataFrame,
    min_read: int,
    pairing_name: str = 'pairing',
) -> pd.DataFrame:
    
    pairbc = set(pairtbl.index)
    msamp_bc = { sn:set(jointhisto.index[ jointhisto[sn]>=min_read ]) for sn in jointhisto.columns }
    msamp_bc[pairing_name] = pairbc

    import upsetplot

    f, ax = plt.subplots(1,1)
    uptbl = upsetplot.from_contents( msamp_bc )
    up = upsetplot.UpSet(uptbl, sort_by='cardinality', show_counts='{:0.3e}', orientation='vertical')
    up.plot(fig=f)

    return f

def plot_indivsamp_waterfall_overlaps(
    jointhisto: pd.DataFrame,
    pairtbl: pd.DataFrame,
    sampname: str,
    **kwargs,
) -> pd.DataFrame:
   
    f, ax = plt.subplots(1,1,**kwargs)

    cursamp = jointhisto[sampname]
    cursamp = cursamp[ cursamp>0 ]

    cts_cursamp_known = np.array(cursamp[ cursamp.index.isin( pairtbl.index ) ])
    cts_cursamp_unk = np.array(cursamp[ ~cursamp.index.isin( pairtbl.index ) ])

    plt.scatter( y=cts_cursamp_known[np.argsort(-cts_cursamp_known)],
                 x=np.arange( cts_cursamp_known.shape[0] ),
                 label='known barcodes',
                 s=4,
    )

    plt.scatter( y=cts_cursamp_unk[np.argsort(-cts_cursamp_unk)],
                 x=np.arange( cts_cursamp_unk.shape[0] ),
                 label='unknown barcodes',
                 s=4,
    )

    plt.yscale('log')
    plt.xscale('log')

    plt.suptitle(sampname)
    plt.ylabel('read count')
    plt.xlabel('barcode rank')
    plt.legend()

    return plt.gcf()

def main():
    parser = argparse.ArgumentParser(description='Compute overlap of barcodes between a barcode pairing table and a set of MPSA RNA-seq barcode histograms')

    parser.add_argument('--pairing_tbl', type=str, required=True,
                      help='path to pairing table')
    
    parser.add_argument('--rnaseq_samp_tbl', type=str, required=True,
                      help='path to sample table of RNA-seq tables; requied columns = libname, bc_histo (containing path to barcode histogram in starcode format)')

    parser.add_argument('--out_base', type=str, default=None,
                      help='Output file path base')

    parser.add_argument('--lmin_barcode', type=str, default='1',
                    help='list of minimum barcode count thresholds, e.g., 1,5,10')

    args = parser.parse_args()
    
    # Read input data
    pairtbl = pd.read_table(args.pairing_tbl).set_index( 'readgroupid' )
    print(f'{args.pairing_tbl} has {pairtbl.shape[0]} rows')
    samp_tbl = pd.read_table(args.rnaseq_samp_tbl)
    
    lmin_read = [ int(x) for x in args.lmin_barcode.split(',') ]
    
    jointhisto = read_mutli_rnaseq_histos( samp_tbl, col_histo='bc_histo' )

    print(jointhisto.shape)
    tbl_sampoverlap = indivsamp_overlap_knownbc( jointhisto, pairtbl, lmin_read )
    tbl_sampoverlap.to_csv( f'{args.out_base}_indivsamp_overlap_knownbc.tsv', sep='\t', index=False )

    upfig = plot_upset_overlaps( jointhisto, pairtbl, min_read=1 )
    upfig.savefig( f'{args.out_base}_upset_overlaps.png' )

    for _,rsamp in samp_tbl.iterrows():
        f = plot_indivsamp_waterfall_overlaps( jointhisto, pairtbl, rsamp['libname'], figsize=(10,4) )
        f.savefig( f'{args.out_base}_knownbc_wfall_{rsamp["libname"]}.png' )
        plt.close(f)

if __name__ == '__main__':
    main() 