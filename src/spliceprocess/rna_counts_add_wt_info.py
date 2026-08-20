from collections import defaultdict
from re import L
import sys
from typing import Dict,List,Set,Callable
import pandas as pd
from splshared.var_mapping import VectorExonTable
from splshared.ssshared import *
import itertools

def wt_effects_single_sample(
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

    # ---------- Pairing: explode variants to (bc,var,is_wt) ----------

    # count variants per bc
    pairtbl['is_wt'] = pairtbl['n_variants_passing'].eq(0)

    wt_only = pairtbl.loc[pairtbl['is_wt'], ['readgroupid', 'is_wt']].rename(
        columns={'readgroupid': 'bc'}
    )

    # pairing-table barcode counts (independent of RNA presence)
    pairing_wt = len(wt_only)

    # ---------- RNA: aggregate once per bc ----------
    lcol_overall_stats = [
        'totalrd_ok','totalrd_bad_ends','totalrd_bad_starts','totalrd_secondary',
        'totalrd_unpaired','totalrd_unmapped','totalrd_soft_clipped'
    ]
    usecols = ['bc','isogrp_name','ok_readcount', *lcol_overall_stats]
    try:
        rna = pd.read_table(fn_bcstatus, usecols=usecols)
    except ValueError as exc:
        if pd.read_table(fn_bcstatus).empty:
            print(f'RNA status file {fn_bcstatus!r} is missing required columns: {exc}')
            rna = pd.DataFrame(columns=usecols)
        else:
            raise ValueError(f"Nonempty RNA status file {fn_bcstatus!r} is missing required columns: {exc}")


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
    pair_rna = wt_only.merge(bc_summary, left_on='bc', right_index=True, how='inner')

    # rna_nbc / rna_nrd
    rna_nbc_wt = pair_rna['bc'].nunique()
    rna_nrd_wt = pair_rna['rds_counted'].sum()

    # overall read status sums (varany only, matching your current output)
    overall_wt = pair_rna[lcol_overall_stats].sum()
    overall_wt.index = [c.replace('totalrd_', 'rna_nrd_') for c in overall_wt.index]

    # PSI means
    psi_cols = [f'psi_{ig}' for ig in lisogrp]
    cts_cols = [f'cts_{ig}' for ig in lisogrp]

    psi_wt_mean = pair_rna[psi_cols].mean()

    # PSI weighted means: sum(cts_ig)/sum(total)
    sum_cts_wt = pair_rna[cts_cols].sum()

    denom_wt = rna_nrd_wt if rna_nrd_wt != 0 else np.nan
    
    psi_wt_wmean = sum_cts_wt.div(denom_wt, axis=0)
    
    # ---------- Assemble output ----------

    out = pd.DataFrame([{
        'libname': libname,
        'pairing_nbc_wt': pairing_wt,
        'rna_nbc_wt': rna_nbc_wt,
        'rna_nrd_varany': rna_nrd_wt,
        **overall_wt.to_dict(),
    }])

    for ig in lisogrp:
        out[f'psi_{ig}_mean'] = psi_wt_mean[f'psi_{ig}']
        out[f'psi_{ig}_wmean'] = psi_wt_wmean[f'cts_{ig}']

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

    parser.add_argument('--out_wt_rpt', help='output variant report', dest='out_wt_rpt' )
    parser.add_argument('--out_perbc', default=None, help='output per-barcode report', dest='out_perbc' )

    args = parser.parse_args()

    isogrptbl = pd.read_table( args.isogrp_tbl )
    pairtbl = pd.read_table( args.bc_pairing_tbl )

    out_wt_rpt =wt_effects_single_sample(
        libname=args.libname,
        ref_seq_name=args.seq_name,
        isogrptbl=isogrptbl,
        pairtbl=pairtbl,
        fn_out_perbc=args.out_perbc,
        fn_bcstatus=args.bc_x_iso_status_tbl,
    )

    out_wt_rpt.to_csv( args.out_wt_rpt, index=False, sep='\t' )

if __name__ == '__main__':
    main()