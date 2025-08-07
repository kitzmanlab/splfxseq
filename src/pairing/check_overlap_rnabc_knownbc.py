#!/usr/bin/env python3

import pandas as pd
import altair as alt
from typing import List, Tuple, Dict, Set, Callable
from itertools import product
import argparse
# from pathlib import Path
# import pysam
# from splshared.var_mapping import GeneModelMapper
from splshared.ssshared import revComp
import numpy as np
import matplotlib.pyplot as plt

def read_starcode_histo(fpath: str) -> pd.DataFrame:
    bchisto = pd.read_table( fpath, header=None, usecols=[0,1] )
    bchisto.columns=['barcode','count']
    bchisto = bchisto.set_index('barcode')['count']
    return bchisto

# def read_mutli_rnaseq_histos(
#     samp_tbl: pd.DataFrame,
#     col_histo = 'bc_histo',
# ) -> pd.DataFrame:
    
#     jointhisto = []
#     for _, row in samp_tbl.iterrows():
#         fn=row[col_histo]
#         bchisto = read_starcode_histo(fn)
#         print(f'{fn} has {bchisto.shape[0]} barcodes')
#         bchisto.name = row['libname']
#         jointhisto.append( bchisto )

#     jointhisto = pd.concat( jointhisto, axis=1 )
    
#     return jointhisto

def indivsamp_overlap_knownbc(
    cts: pd.Series,
    samp_name: str,
    pairtbl: pd.DataFrame,
    lmin_read: List[int],
    lqile: List[float],
    # revcomp: bool = False,
) -> pd.DataFrame:
    
    # 1. how many barcodes are there at a given read cutoff, and how many of those are in the pairing table?
    # 2. either with or without filtering to pairing table, how many barcodes are there at given quantile of cumulative distribution of barcodes?

    pairbc = set(pairtbl.index)
    # if revcomp:
    #     pairbc = set( [ revComp(bc) for bc in pairbc ] )
    
    tbl_out_atmincts = { cn:[] for cn in 'libname,min_read,n_bc_sample,n_bc_pairing,n_bc_overlap,n_read_sample,n_read_overlap'.split(',') }
    tbl_out_atqiles = { cn:[] for cn in 'libname,qile,n_bc_sample,n_bc_pairing,n_bc_overlap'.split(',') }

    li_bcknown = cts.index.isin( pairbc )

    for min_read in lmin_read:
        li_sampbc = cts>=min_read 
        
        tbl_out_atmincts['libname'].append( samp_name )
        tbl_out_atmincts['min_read'].append( min_read )
        tbl_out_atmincts['n_bc_sample'].append( (li_sampbc).sum() )
        tbl_out_atmincts['n_bc_pairing'].append( len(pairbc) )
        tbl_out_atmincts['n_bc_overlap'].append( (li_sampbc & li_bcknown).sum() )

        tbl_out_atmincts['n_read_sample'].append( cts[ li_sampbc ].sum() )
        tbl_out_atmincts['n_read_overlap'].append( cts[ (li_sampbc & li_bcknown) ].sum() )

    cts_all = np.array(cts,dtype=int)

    li_sort = (-cts_all).argsort()
    cts_all = cts_all[li_sort]
    ccts = cts_all.cumsum()
    ccts = ccts/ccts[-1]

    for qile in lqile:
        nbc = np.where(ccts>qile)[0].min()
        nbc_known = li_bcknown[li_sort[:nbc]].sum()

        tbl_out_atqiles['libname'].append(samp_name)
        tbl_out_atqiles['qile'].append(qile)
        tbl_out_atqiles['n_bc_sample'].append( nbc )
        tbl_out_atqiles['n_bc_pairing'].append( len(pairbc) )
        tbl_out_atqiles['n_bc_overlap'].append( nbc_known )


    tbl_out_atmincts = pd.DataFrame(tbl_out_atmincts)
    tbl_out_atqiles = pd.DataFrame(tbl_out_atqiles)
    
    return tbl_out_atmincts, tbl_out_atqiles


# def plot_upset_overlaps(
#     jointhisto: pd.DataFrame,
#     pairtbl: pd.DataFrame,
#     min_read: int,
#     pairing_name: str = 'pairing',
# ) -> pd.DataFrame:
    
#     pairbc = set(pairtbl.index)
#     msamp_bc = { sn:set(jointhisto.index[ jointhisto[sn]>=min_read ]) for sn in jointhisto.columns }
#     msamp_bc[pairing_name] = pairbc

#     import upsetplot

#     f, ax = plt.subplots(1,1)
#     uptbl = upsetplot.from_contents( msamp_bc )
#     up = upsetplot.UpSet(uptbl, sort_by='cardinality', show_counts='{:0.3e}', orientation='vertical')
#     up.plot(fig=f)

#     return f

def plot_indivsamp_waterfall_overlaps(
    cts: pd.Series,
    samp_name: str,
    pairtbl: pd.DataFrame,
    **kwargs,
) -> pd.DataFrame:
   
    f, ax = plt.subplots(1,1,**kwargs)

    cursamp = cts
    cursamp = cursamp[ cursamp>0 ]

    cts_cursamp_known = np.array(cursamp[ cursamp.index.isin( pairtbl.index ) ])
    cts_cursamp_unk = np.array(cursamp[ ~cursamp.index.isin( pairtbl.index ) ])

    plt.scatter( y=cursamp[np.argsort(-cursamp)],
                 x=np.arange( cursamp.shape[0] ),
                 label='all barcodes',
                 s=4,
    )
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

    plt.suptitle(samp_name)
    plt.ylabel('read count')
    plt.xlabel('barcode rank')
    plt.legend()

    return plt.gcf()

def plot_overlap_summary_charts(df: pd.DataFrame) -> list:
    """Create bar charts showing barcode overlap statistics across samples
    
    Args:
        df: DataFrame with columns:
            - libname: Sample name
            - n_bc_overlap: Number of overlapping barcodes
            - n_bc_sample: Total barcodes in sample
            - min_read: Minimum read count threshold
            - pairing_tbl_name: Name of pairing table used
            
    Returns:
        List of Altair charts showing:
        1. Raw overlap counts
        2. Overlap fractions
    """
    
    # Calculate fraction overlap
    df = df.copy()
    df['overlap_fraction'] = df['n_bc_overlap'] / df['n_bc_sample']
    
    # Base chart settings
    base = alt.Chart(df).encode(
        alt.Y('libname:N', title='Sample'),
        alt.Fill('min_read:N', title='Minimum read count threshold'),
        alt.Tooltip(['libname', 'min_read', 'n_bc_overlap', 'n_bc_sample', 'overlap_fraction'])
    ).properties(
        width=300,
        height=alt.Step(20)  # Controls bar height
    )
    
    loc=list(df.columns)

    # Raw counts chart
    count_chart = base.mark_bar().encode(
        x=alt.X('n_bc_overlap:Q', title='Number of overlapping barcodes'),
        tooltip=alt.Tooltip(loc)
    ).facet(row='pairing_tbl_name', column='min_read', title='Pairing table'
    ).resolve_scale(
        y='independent'
    )
    
    # Fraction chart  
    frac_chart = base.mark_bar().encode(
        x=alt.X('overlap_fraction:Q', title='Fraction of overlapping barcodes'),
        tooltip=alt.Tooltip(loc)
    ).facet(row='pairing_tbl_name', column='min_read', title='Pairing table'
    ).resolve_scale(
        y='independent'
    )
    
    return [count_chart, frac_chart]


def main():
    parser = argparse.ArgumentParser(description='Compute overlap of barcodes between a barcode pairing table and a set of MPSA RNA-seq barcode histograms')

    opts_sps = parser.add_subparsers( dest='cmd' )

    opts_1samp = opts_sps.add_parser('singlesamp')

    opts_1samp.add_argument('--pairing_tbl', type=str, required=True,
                      help='path to pairing table')

    opts_1samp.add_argument('--bc_histo', type=str, required=True,
                      help='path to barcode histogram in starcode format')

    opts_1samp.add_argument('--samp_name', type=str, required=True,
                            help='name of sample')

    opts_1samp.add_argument('--out_base', type=str, default=None,
                      help='Output file path base')

    opts_1samp.add_argument('--lmin_barcode', type=str, default='1',
                    help='list of minimum barcode count thresholds, e.g., "1,5,10"')

    opts_1samp.add_argument('--lqile', type=str, default='0.75,0.90',
                    help='barcode cumulative quantile thresholds, e.g., "0.75,0.90"')
    
    opts_mutlisampsumy = opts_sps.add_parser('summarize_samps')

    opts_mutlisampsumy.add_argument('--overlap_table', type=str, default='overlap_table')
    opts_mutlisampsumy.add_argument('--overlap_qile_table', type=str, default='overlap_qile_table')                      
   
    opts_mutlisampsumy.add_argument('--out_base', type=str, default=None,
                      help='Output file path base')

    args = parser.parse_args()


    if args.cmd == 'singlesamp':
        cts = read_starcode_histo(args.bc_histo)
        pairtbl = pd.read_table(args.pairing_tbl).set_index( 'readgroupid' )
        lmin_read = [ int(x) for x in args.lmin_barcode.split(',') ]
        lqile = [ float(x) for x in args.lqile.split(',') ]
        tbl_out_atmincts, tbl_out_atqiles = indivsamp_overlap_knownbc( cts, args.samp_name, pairtbl, lmin_read, lqile )
        tbl_out_atmincts.to_csv( f'{args.out_base}_indivsamp_overlap_knownbc.tsv', sep='\t', index=False )
        tbl_out_atqiles.to_csv( f'{args.out_base}_indivsamp_overlap_knownbc_qiles.tsv', sep='\t', index=False )
        plot_indivsamp_waterfall_overlaps( cts, args.samp_name, pairtbl, figsize=(10,4) )
        plt.savefig( f'{args.out_base}_knownbc_wfall_{args.samp_name}.png' )

    elif args.cmd == 'summarize_samps':
        tbl_overlap = pd.read_table(args.overlap_table)
        tbl_overlap_qile = pd.read_table(args.overlap_qile_table)
        lcharts = plot_overlap_summary_charts(tbl_overlap)
        lcharts[0].save(f'{args.out_base}_overlap_counts.html')
        lcharts[1].save(f'{args.out_base}_overlap_fractions.html')
        

    # Read input data
    # pairtbl = pd.read_table(args.pairing_tbl).set_index( 'readgroupid' )
    # print(f'{args.pairing_tbl} has {pairtbl.shape[0]} rows')
    # samp_tbl = pd.read_table(args.rnaseq_samp_tbl)
    
    # lmin_read = [ int(x) for x in args.lmin_barcode.split(',') ]
    
    # jointhisto = read_mutli_rnaseq_histos( samp_tbl, col_histo='bc_histo' )

    # print(jointhisto.shape)
    # tbl_sampoverlap = indivsamp_overlap_knownbc( jointhisto, pairtbl, lmin_read )
    # tbl_sampoverlap.to_csv( f'{args.out_base}_indivsamp_overlap_knownbc.tsv', sep='\t', index=False )

    # upfig = plot_upset_overlaps( jointhisto, pairtbl, min_read=1 )
    # upfig.savefig( f'{args.out_base}_upset_overlaps.png' )

    # for _,rsamp in samp_tbl.iterrows():
    #     f = plot_indivsamp_waterfall_overlaps( jointhisto, pairtbl, rsamp['libname'], figsize=(10,4) )
    #     f.savefig( f'{args.out_base}_knownbc_wfall_{rsamp["libname"]}.png' )
    #     plt.close(f)

if __name__ == '__main__':
    main() 