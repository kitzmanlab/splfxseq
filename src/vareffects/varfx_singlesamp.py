import pandas as pd
from splshared.ssshared import *
import altair as alt

def is_snv(
    var_rpt: pd.DataFrame,
    col_genomic_var = 'var'
):
    li_snv = var_rpt[col_genomic_var].str.split(':').apply(lambda x: len(x[2]) == 1 and len(x[3]) == 1) 
    return li_snv

def filter_to_snv(
    var_rpt: pd.DataFrame,
    col_genomic_var = 'var'
):
    li_snv = is_snv( var_rpt, col_genomic_var )
    return var_rpt.loc[ li_snv ]

def plot_varfx(
    var_rpt: pd.DataFrame,
    isogrptbl: pd.DataFrame,
    ref_seq_name: str,
    col_genomic_var = 'var',
    min_rna_bc = 25,
):
    if not is_snv( var_rpt, col_genomic_var ).all():
        raise ValueError('Input table must have only SNVs')

    var_rpt['genomic_pos'] = [ int(x.split(':')[1]) for x in var_rpt[col_genomic_var] ]
    var_rpt['ref'] = [ x.split(':')[2] for x in var_rpt[col_genomic_var] ]
    var_rpt['alt'] = [ x.split(':')[3] for x in var_rpt[col_genomic_var] ]

    _var_rpt = var_rpt

    var_rpt = var_rpt.loc[ var_rpt['rna_nbc_varany'] >= min_rna_bc ].copy()

    if var_rpt.shape[0]==0:
        print(f'WARNING - {_var_rpt.shape[0]} variants before filtering for >= {min_rna_bc} rna bcs, but none pass filter')
        ch = alt.Chart(var_rpt,title=f'NO VARIANTS PASS FILTER').mark_circle()
        return ch

    var_rpt['rna_nrd_bad'] = var_rpt[ 'rna_nrd_bad_ends,rna_nrd_bad_starts,rna_nrd_secondary,rna_nrd_unpaired,rna_nrd_unmapped,rna_nrd_soft_clipped'.split(',')].sum(1)


    libname = var_rpt['libname'].unique()[0]
    if len(set(var_rpt['libname'])) > 1:
        raise ValueError('Input table must have only one libname')
    
    lpl= []

    base = alt.Chart(
        var_rpt
    ).encode(
        x=alt.X('genomic_pos:Q', title=None).scale(zero=False),
    ).properties(width=1000,height=100)

    loc = 'var,pairing_nbc_varsingleton,rna_nbc_varsingleton,rna_nrd_ok,rna_nrd_bad'.split(',') + \
        [ f'psi_{isogrp}_singleton_wmean' for isogrp in isogrptbl['isogrp_name'].unique() ]

    for isogrp in isogrptbl['isogrp_name'].unique():
        p4 = base.mark_circle().encode(
            alt.Y(f'psi_{isogrp}_singleton_wmean:Q', title=[f'PSI {isogrp}','(singleton bc, Wmean']).scale(domain=[0,1]),
            alt.Color('alt:N', title='Alternate Base'),
            alt.Tooltip(loc)
        )
        lpl.append(p4)
    

    p1 = base.mark_circle(
        opacity=0.5,
    ).encode(   
        alt.Y('pairing_nbc_varany:Q', title='#pairing bcs' ),
        alt.Color('alt:N', title='Alternate Base'),
        alt.Tooltip(['genomic_pos','ref','alt','pairing_nbc_varany']),
    )

    p2 = base.mark_circle(
    ).encode(   
        alt.Y('rna_nbc_varany:Q', title='#rna bcs' ),
        alt.Color('alt:N', title='Alternate Base'),
        alt.Tooltip(['genomic_pos','ref','alt','pairing_nbc_varany'])
    )
    
    p3 = base.mark_circle(
    ).encode(   
        alt.Y('rna_nrd_varany:Q', title='#rna reads' ),
        alt.Color('alt:N', title='Alternate Base'),
        alt.Tooltip(['genomic_pos','ref','alt','pairing_nbc_varany'])
    )
    
    lpl += [p1,p2,p3]


    ch = alt.vconcat(*lpl).properties(title=f'{libname} / {ref_seq_name}')

    ch = ch.interactive(bind_y=False)

    return ch


def main():
    import argparse

    parser = argparse.ArgumentParser(description='combine processed rna count tables with vairant-barcode table')

    opts_sps = parser.add_subparsers( dest='cmd' )

    opts_filtsnv = opts_sps.add_parser('filt_to_snv')
    opts_filtsnv.add_argument('--var_tbl_in', help='input table with per-barcode isoform group counts', dest='var_tbl_in' )
    opts_filtsnv.add_argument('--var_tbl_out', help='output table with per-barcode isoform group counts', dest='var_tbl_out' )

    opts_mkplot = opts_sps.add_parser('singlerep_plot')
    opts_mkplot.add_argument('--var_tbl',  dest='var_tbl' )
    opts_mkplot.add_argument('--isogrp_tbl',  dest='isogrp_tbl' )
    opts_mkplot.add_argument('--ref_seq_name',  dest='ref_seq_name' )
    opts_mkplot.add_argument('--out_plot',  dest='out_plot' )

    opts_mkplot.add_argument('--min_rna_bc',  dest='min_rna_bc', type=int, default=20 )

    args = parser.parse_args()

    if args.cmd == 'filt_to_snv':
        var_rpt = pd.read_table( args.var_tbl_in )
        var_rpt_snv = filter_to_snv( var_rpt, col_genomic_var='var' )
        var_rpt_snv.to_csv( args.var_tbl_out, index=False, sep='\t' )

    elif args.cmd == 'singlerep_plot':
        var_rpt = pd.read_table( args.var_tbl )
        isogrp_tbl = pd.read_table( args.isogrp_tbl )
        isogrp_tbl = isogrp_tbl.query( f'seq_name == "{args.ref_seq_name}"' )
        ch = plot_varfx( var_rpt, isogrp_tbl, args.ref_seq_name )
        ch.save( args.out_plot )
    
if __name__ == '__main__':
    main()