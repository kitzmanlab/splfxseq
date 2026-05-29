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
    var_sort_fn = lambda v:(v.split(':')[0], int(v.split(':')[1])),
):
    # iso groups for this reference
    isogrp_sub = isogrptbl.loc[isogrptbl['seq_name'].eq(ref_seq_name)]
    lisogrp = isogrp_sub['isogrp_name'].unique().tolist()

    # ---------- Pairing: explode variants to (bc,var,is_singleton) ----------
    pair = pairtbl.loc[pairtbl[col_var_annot].notna() & pairtbl[col_var_annot].ne("")].copy()

    # count variants per bc
    pair['nvar'] = pair[col_var_annot].str.count(',') + 1
    pair['is_singleton'] = pair['nvar'].eq(1)

    single = pair.loc[pair['is_singleton'], ['readgroupid', 'is_singleton', col_var_annot]].rename(
        columns={'readgroupid': 'bc', col_var_annot: 'var'}
    )

    multi = pair.loc[~pair['is_singleton'] & pair['nvar'].le(max_var_per_bc), ['readgroupid', 'is_singleton', col_var_annot]].copy()
    multi['var'] = multi[col_var_annot].str.split(',')
    multi = multi.explode('var').rename(columns={'readgroupid': 'bc'})[['bc', 'is_singleton', 'var']]
    multi = multi.loc[multi['var'].ne("")]  # drop blanks

    pair_long = pd.concat([single[['bc','is_singleton','var']], multi], ignore_index=True)
    pair_long = pair_long.drop_duplicates(['bc', 'var'])

    # pairing-table barcode counts (independent of RNA presence)
    pairing_any = pair_long.groupby('var', sort=False)['bc'].nunique()
    pairing_sing = pair_long.loc[pair_long['is_singleton']].groupby('var', sort=False)['bc'].nunique()

    # ---------- RNA: aggregate once per bc ----------
    lcol_overall_stats = [
        'totalrd_ok','totalrd_bad_ends','totalrd_bad_starts','totalrd_secondary',
        'totalrd_unpaired','totalrd_unmapped','totalrd_soft_clipped'
    ]
    usecols = ['bc','isogrp_name','ok_readcount', *lcol_overall_stats]
    rna = pd.read_table(fn_bcstatus, usecols=usecols)

    # make grouping faster if isogrp_name repeats a lot
    rna['isogrp_name'] = pd.Categorical(rna['isogrp_name'], categories=lisogrp)

    bc_overall = rna.groupby('bc', sort=False)[lcol_overall_stats].first()

    bc_iso = (
        rna.groupby(['bc','isogrp_name'], sort=False, observed=True)['ok_readcount']
           .sum()
           .unstack(fill_value=0)
           .reindex(columns=lisogrp, fill_value=0)
    )
    bc_total = bc_iso.sum(axis=1).rename('rds_counted')

    bc_psi = bc_iso.div(bc_total.replace(0, np.nan), axis=0)  # per-bc psi (NaN when no reads)

    bc_summary = pd.concat(
        [
            bc_total,
            bc_iso.add_prefix('cts_'),
            bc_psi.add_prefix('psi_'),
            bc_overall
        ],
        axis=1
    )

    # ---------- Merge pairing edges with per-bc RNA summary ----------
    pair_rna = pair_long.merge(bc_summary, left_on='bc', right_index=True, how='inner')

    # group “any” and “singleton”
    grp_any = pair_rna.groupby('var', sort=False)
    pair_rna_sing = pair_rna.loc[pair_rna['is_singleton']]
    grp_sing = pair_rna_sing.groupby('var', sort=False)

    # rna_nbc / rna_nrd
    rna_nbc_any = grp_any['bc'].nunique()
    rna_nrd_any = grp_any['rds_counted'].sum()

    rna_nbc_sing = grp_sing['bc'].nunique()
    rna_nrd_sing = grp_sing['rds_counted'].sum()

    # overall read status sums (varany only, matching your current output)
    overall_any = grp_any[lcol_overall_stats].sum()
    overall_any.columns = [c.replace('totalrd_', 'rna_nrd_') for c in overall_any.columns]

    # PSI means
    psi_cols = [f'psi_{ig}' for ig in lisogrp]
    cts_cols = [f'cts_{ig}' for ig in lisogrp]

    psi_any_mean = grp_any[psi_cols].mean()
    psi_sing_mean = grp_sing[psi_cols].mean()

    # PSI weighted means: sum(cts_ig)/sum(total)
    sum_cts_any = grp_any[cts_cols].sum()
    sum_cts_sing = grp_sing[cts_cols].sum()

    denom_any = rna_nrd_any.replace(0, np.nan)
    denom_sing = rna_nrd_sing.replace(0, np.nan)

    psi_any_wmean = sum_cts_any.div(denom_any, axis=0)
    psi_sing_wmean = sum_cts_sing.div(denom_sing, axis=0)

    # ---------- Assemble output ----------
    allvars = pairing_any.index.union(rna_nbc_any.index)
    if var_sort_fn is not None:
        allvars = sorted(allvars, key=var_sort_fn)
    else:
        allvars = sorted(allvars)

    out = pd.DataFrame(index=pd.Index(allvars, name='var'))
    out['libname'] = libname
    out['var'] = out.index

    out['pairing_nbc_varany'] = pairing_any.reindex(out.index).fillna(0).astype(int)
    out['pairing_nbc_varsingleton'] = pairing_sing.reindex(out.index).fillna(0).astype(int)

    out['rna_nbc_varany'] = rna_nbc_any.reindex(out.index).fillna(0).astype(int)
    out['rna_nrd_varany'] = rna_nrd_any.reindex(out.index).fillna(0).astype(int)
    out['rna_nbc_varsingleton'] = rna_nbc_sing.reindex(out.index).fillna(0).astype(int)
    out['rna_nrd_varsingleton'] = rna_nrd_sing.reindex(out.index).fillna(0).astype(int)

    out = out.join(overall_any.reindex(out.index, fill_value=0))

    for ig in lisogrp:
        out[f'psi_{ig}_varany_mean'] = psi_any_mean[f'psi_{ig}'].reindex(out.index)
        out[f'psi_{ig}_singleton_mean'] = psi_sing_mean.get(f'psi_{ig}', pd.Series(dtype=float)).reindex(out.index)

        out[f'psi_{ig}_varany_wmean'] = psi_any_wmean[f'cts_{ig}'].reindex(out.index)
        out[f'psi_{ig}_singleton_wmean'] = psi_sing_wmean.get(f'cts_{ig}', pd.Series(dtype=float)).reindex(out.index)

    # ---------- Optional per-bc output (write once) ----------
    if fn_out_perbc:
        # include original pairing columns by joining back to pairtbl on bc/readgroupid
        perbc = pair_rna.merge(pairtbl, left_on='bc', right_on='readgroupid', how='left', suffixes=('','_pair'))
        perbc.to_csv(fn_out_perbc, sep='\t', index=False, compression='gzip')

    return out.reset_index(drop=True)

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