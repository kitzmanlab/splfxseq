from collections import defaultdict
from re import L
import sys
from typing import Dict,List,Set,Callable
import pandas as pd
from splshared.var_mapping import VectorExonTable
from splshared.ssshared import *
import itertools

def var_effects_single_sample(
    libname: str,
    ref_seq_name: str,
    isogrptbl: pd.DataFrame,
    pairtbl: pd.DataFrame,
    fn_bcstatus: str,
    fn_out_perbc: str = None,
    col_var_annot: str = 'variant_list_genome',
    max_var_per_bc: int = 3,
    var_sort_fn: Callable = lambda v:(v.split(':')[0],int(v.split(':')[1])),
):

    # libname - library name
    # seq_name - target reference sequence name
    # isogrptbl - isoform group table
    # pairtbl - pairing table
    # fn_bcstatus - per-barcode isoform counts table

    isogrptbl = isogrptbl.query( f'seq_name == "{ref_seq_name}"' )
    lisogrp = list(isogrptbl['isogrp_name'].unique())

    out_rpt_isogrp = {k:[] for k in 'libname,var,pairing_nbc_varany,pairing_nbc_varsingleton,rna_nbc_varany,rna_nrd_varany,rna_nbc_varsingleton,rna_nrd_varsingleton'.split(',')+
        ['rna_nrd_ok','rna_nrd_bad_ends','rna_nrd_bad_starts','rna_nrd_secondary','rna_nrd_unpaired','rna_nrd_unmapped','rna_nrd_soft_clipped']+
        [f'psi_{isogrp}_singleton_wmean' for isogrp in lisogrp]+
        [f'psi_{isogrp}_singleton_mean' for isogrp in lisogrp]+
        [f'psi_{isogrp}_varany_wmean' for isogrp in lisogrp]+
        [f'psi_{isogrp}_varany_mean' for isogrp in lisogrp] }


    # pairing_nbc_varany - how many barcodes in the pairing table carrying that variant, at all
    # pairing_nbc_varsingleton - how many barcodes in the pairing table carrying that variant, as a singleton
    
    # rna_nbc_varany - how many barcodes in the rna table carrying that variant, at all
    # rna_nrd_varany - how many accepted* reads in the rna table from those barcodes  (*=from named isoforms or OTHER)

    # rna_nbc_varsingleton - how many barcodes in the rna table carrying that variant, as a singleton
    # rna_nrd_varsingleton - how many accepted* reads in the rna table from those barcodes  (*=from named isoforms or OTHER)

    # psi_{isogrp}_varany_wmean - mean psi across barcodes with this variant (allowing other vars), weighted by #reads per barcode
    # psi_{isogrp}_varany_mean - mean psi across barcodes with this variant (allowing other vars), NOT weighted by #reads per barcode
    # psi_{isogrp}_singleton_wmean - mean psi across barcodes with this variant (only), weighted by #reads per barcode
    # psi_{isogrp}_singleton_mean - mean psi across barcodes with this variant (only), NOT weighted by #reads per barcode
    

    lcol_outperbc = ('var,is_singleton,bc,nrd_total_counted'.split(',')+
        [f'nrd_{isogrp}' for isogrp in lisogrp]+
        [f'psi_{isogrp}' for isogrp in lisogrp]+
        [cn for cn in pairtbl.columns])

    if fn_out_perbc:
        out_perbc = {k:[] for k in lcol_outperbc}
        
    pairtbl_nnvar = pairtbl.loc[ (pairtbl[col_var_annot]!="") & ~pd.isnull(pairtbl[col_var_annot]) ]
    pairtbl_nnvar = pairtbl_nnvar.set_index('readgroupid', drop=False)

    # map variants to barcodes for faster lookup
    li_singlevar = ~ pairtbl_nnvar[col_var_annot].str.contains(',')
    m_varsingle_sbc = {}
    for vid, pairbcs_vid in pairtbl_nnvar.groupby(col_var_annot):
        m_varsingle_sbc[vid] = set(pairbcs_vid['readgroupid'])
    
    li_withinmax_var = ( pairtbl_nnvar[col_var_annot].isnull() ) | (pairtbl_nnvar[col_var_annot].str.count(',') < max_var_per_bc)

    m_varmulti_sbc = defaultdict(list)
    for _,r in pairtbl_nnvar[li_withinmax_var &~li_singlevar].iterrows():
        for vid in r[col_var_annot].split(','):
            m_varmulti_sbc[vid].append( r['readgroupid'] )
    for vid in m_varmulti_sbc:
        m_varmulti_sbc[vid] = set(m_varmulti_sbc[vid])

    # find all vars in pairing table    
    allvars = [ lv.split(',') for lv in pairtbl_nnvar[col_var_annot] ]
    allvars = list(itertools.chain(*allvars))
    allvars = list(set(allvars))
    if var_sort_fn is not None:
        allvars = sorted( allvars, key=var_sort_fn )
    else:
        allvars = sorted( allvars )

    rna_bc_iso = pd.read_table( fn_bcstatus ).set_index('bc')

    sys.stderr.write('total %d vars\n'%(len(allvars)))
    sys.stderr.flush()

    first_out = True

    vidctr = 0
    # loop over all vars in the pairing table
    for vid in allvars:

        vidctr+=1
        if vidctr%50==0: 
            sys.stderr.write('%d..'%(vidctr))
            sys.stderr.flush()

        # find barcodes in the pairing table carrying that variant, at all, or as a singleton
        sbc_pair_var_sing = m_varsingle_sbc.get(vid,[])
        sbc_pair_var_multi = m_varmulti_sbc.get(vid,[])

        # find those bcs in the rna table
        rna_curvar_multi = rna_bc_iso.loc[ rna_bc_iso.index.isin(sbc_pair_var_multi) ].reset_index(drop=False)
        rna_curvar_sing = rna_bc_iso.loc[ rna_bc_iso.index.isin(sbc_pair_var_sing) ].reset_index(drop=False)

        lbc_pair_var_multi = list(rna_curvar_multi['bc'])

        cts_bc_x_iso_sing = rna_curvar_sing.groupby( ['bc','isogrp_name'] ).agg(  {'ok_readcount':'sum'} )
        cts_bc_x_iso_multi = rna_curvar_multi.groupby( ['bc','isogrp_name'] ).agg(  {'ok_readcount':'sum'} )
        cts_bc_x_iso = pd.concat( [cts_bc_x_iso_sing, cts_bc_x_iso_multi] )
        cts_bc_x_iso = cts_bc_x_iso.reset_index().pivot( index='bc', columns='isogrp_name', values='ok_readcount' )
        cts_bc_x_iso = cts_bc_x_iso.fillna(0).astype(int)

        for ig in lisogrp:
            if ig not in cts_bc_x_iso.columns:
                cts_bc_x_iso[ig] = 0

        sum_ctd_by_bc = cts_bc_x_iso.sum(axis=1)

        cts_bc_x_iso.columns = ['cts_'+cn for cn in cts_bc_x_iso.columns]
        cts_bc_x_iso['is_singleton'] = True
        cts_bc_x_iso.loc[ lbc_pair_var_multi, 'is_singleton' ] = False
        cts_bc_x_iso['varid'] = vid        
        
        bybc_tbl = cts_bc_x_iso.copy()
        bybc_tbl['rds_counted'] = sum_ctd_by_bc

        for ig in lisogrp:
            bybc_tbl[f'psi_{ig}'] = bybc_tbl[f'cts_{ig}'] / bybc_tbl['rds_counted']

        lcol_overall_stats = 'totalrd_ok,totalrd_bad_ends,totalrd_bad_starts,totalrd_secondary,totalrd_unpaired,totalrd_unmapped,totalrd_soft_clipped'.split(',')
        bc_overall_stats = pd.concat( [ 
            rna_curvar_sing.groupby( 'bc' ).agg( 'first' )[lcol_overall_stats],
            rna_curvar_multi.groupby( 'bc' ).agg( 'first' )[lcol_overall_stats]
        ])

        bybc_tbl = pd.concat( [ bybc_tbl, bc_overall_stats ], axis=1 )

        if fn_out_perbc:
            bybc_tbl_withpairinfo = pd.merge( bybc_tbl, pairtbl_nnvar, left_index=True, right_index=True, how='left', suffixes=('','_pair'), indicator='vennpart' )

            if first_out:
                bybc_tbl_withpairinfo.to_csv( fn_out_perbc, index=True, sep='\t', compression='gzip' )
                first_out = False
            else:
                bybc_tbl_withpairinfo.to_csv( fn_out_perbc, index=True, sep='\t', mode='a', compression='gzip', header=False )
            
        bybc_sing = bybc_tbl.query('is_singleton')

        out_rpt_isogrp['libname'] +=  [libname]
        out_rpt_isogrp['var'] +=  [vid]
        out_rpt_isogrp['pairing_nbc_varany'] +=  [len(sbc_pair_var_multi)+len(sbc_pair_var_sing)]
        out_rpt_isogrp['pairing_nbc_varsingleton'] +=  [len(sbc_pair_var_sing)]
        out_rpt_isogrp['rna_nbc_varany'] +=  [bybc_tbl['rds_counted'].shape[0]]
        out_rpt_isogrp['rna_nrd_varany'] +=  [bybc_tbl['rds_counted'].sum()]
        out_rpt_isogrp['rna_nbc_varsingleton'] +=  [bybc_sing.shape[0]]
        out_rpt_isogrp['rna_nrd_varsingleton'] +=  [bybc_sing['rds_counted'].sum()]

        for cn in lcol_overall_stats:
            out_rpt_isogrp[cn.replace('totalrd_','rna_nrd_')] += [ bc_overall_stats[cn].sum() ]

        for isogrp in lisogrp:
            if bybc_tbl[f'psi_{isogrp}'].shape[0]>0:
                rdsc = bybc_tbl['rds_counted'].sum()
                if rdsc > 0:
                    out_rpt_isogrp[f'psi_{isogrp}_varany_wmean'].append( ( bybc_tbl['rds_counted'] * bybc_tbl[f'psi_{isogrp}']).sum() / rdsc )
                else:
                    out_rpt_isogrp[f'psi_{isogrp}_varany_wmean'].append( None )

                out_rpt_isogrp[f'psi_{isogrp}_varany_mean'].append( bybc_tbl[f'psi_{isogrp}'].mean() )
                
                rdsc = bybc_sing['rds_counted'].sum()
                if rdsc > 0:
                    out_rpt_isogrp[f'psi_{isogrp}_singleton_wmean'].append( ( bybc_sing['rds_counted'] * bybc_sing[f'psi_{isogrp}']).sum() / rdsc )
                else:
                    out_rpt_isogrp[f'psi_{isogrp}_singleton_wmean'].append( None )

                out_rpt_isogrp[f'psi_{isogrp}_singleton_mean'].append( bybc_sing[f'psi_{isogrp}'].mean() )
            else:
                for cn in f'psi_{isogrp}_varany_wmean,psi_{isogrp}_varany_mean,psi_{isogrp}_singleton_wmean,psi_{isogrp}_singleton_mean'.split(','):
                    out_rpt_isogrp[cn].append( None )

    out_rpt_isogrp = pd.DataFrame( out_rpt_isogrp )
    return out_rpt_isogrp


def main():
    import argparse

    parser = argparse.ArgumentParser(description='combine processed rna count tables with vairant-barcode table')
    
    parser.add_argument('--isogrp_tbl', help='input table with per-barcode isoform group counts', dest='isogrp_tbl' )
    parser.add_argument('--bc_x_iso_status_tbl', help='input table with per-barcode isoform group counts', dest='bc_x_iso_status_tbl' )
    parser.add_argument('--seq_name', help='name of target reference sequence', dest='seq_name' )
    parser.add_argument('--bc_pairing_tbl', help='input pairing table', dest='bc_pairing_tbl' )
    parser.add_argument('--libname', help='library name', dest='libname' )
    parser.add_argument('--col_var_annot', help='column in pairing table with variant name (genomic or cDNA)', default='variant_list_genome', dest='col_var_annot' )

    parser.add_argument('--out_var_rpt', help='output variant report', dest='out_var_rpt' )
    parser.add_argument('--out_perbc', default=None, help='output per-barcode report', dest='out_perbc' )

    args = parser.parse_args()

    isogrptbl = pd.read_table( args.isogrp_tbl )
    pairtbl = pd.read_table( args.bc_pairing_tbl )

    out_var_rpt =var_effects_single_sample(
        libname=args.libname,
        ref_seq_name=args.seq_name,
        isogrptbl=isogrptbl,
        pairtbl=pairtbl,
        fn_out_perbc=args.out_perbc,
        fn_bcstatus=args.bc_x_iso_status_tbl,
    )

    out_var_rpt.to_csv( args.out_var_rpt, index=False, sep='\t' )
    
if __name__ == '__main__':
    main()