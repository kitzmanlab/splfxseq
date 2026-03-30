import pandas as pd
import numpy as np
import scipy.stats as ss
from maxentpy import maxent
from collections import Counter
import postprocess.postprocess_utils as pp
import splshared.var_mapping as mapper



def add_cols(
    tbl_: pd.DataFrame,
    maptbl: pd.DataFrame,
    exon: str
    ): 
    tbl = tbl_.copy()
    
    # make chr, pos, ref, alt columns based on var and add as first columns
    if 'var' not in tbl.columns:
        tbl = tbl.rename(columns={'varid':'var'})
    new_cols = tbl['var'].str.split(':', expand=True)
    new_cols.columns = ['chr', 'pos', 'ref', 'alt']
    new_cols['pos'] = new_cols['pos'].astype(int)
    tbl = pd.concat([new_cols, tbl], axis=1)
    
    # add in exon, and hgvs/vector pos columns
    tbl['exon'] = exon
    tbl['hgvs_pos'] = mapper.genome_to_hgvs(maptbl, tbl['pos'], exon )
    tbl['vector_pos'] = mapper.genome_to_vector(maptbl, tbl['pos'], exon)
   
    return tbl

def prep_tbls(
    var_tbl: pd.DataFrame,
    map_tbl: pd.DataFrame,
    exon: str,
    nbc_per_var_min: int = 5
    ):

    var_dfs, bc_dfs = [], []
    for _, r in var_tbl.iterrows():
        var_df = pd.read_table(r['var_rpt_snvonly'])
        # filter to only include bc above provided min
        var_df = var_df[var_df['rna_nbc_varsingleton'] > nbc_per_var_min].copy()

        bc_df = pd.read_table(r['var_bc_rpt'])
        bc_df = bc_df[bc_df['is_singleton'] ==True].copy()
        bc_df['libname'] = r['libname']
        
        
        bc_df = bc_df[bc_df['varid'].isin(var_df['var'])].copy()
        var_df = add_cols(var_df, map_tbl, exon)
        bc_df = add_cols(bc_df, map_tbl, exon)
        var_df['rna_nrd_bad'] = var_df['rna_nrd_bad_ends,rna_nrd_bad_starts,rna_nrd_secondary,rna_nrd_unpaired,rna_nrd_unmapped,rna_nrd_soft_clipped'.split(',')].sum(1)
    
        var_dfs.append(var_df)
        bc_dfs.append(bc_df)


    byvartbl = pd.concat(var_dfs, axis=0, ignore_index=True)
    bybctbl = pd.concat(bc_dfs, axis=0, ignore_index=True)
    
    return byvartbl, bybctbl

# get the first 90% of reads that account for the total sample
def top_cumulative_rows(df, value_col='totalrd_ok', group_col='libname', threshold=0.90):
    dfs = []
    # Process each sample group independently
    for sample, group in df.groupby(group_col):
        # Sort by value descending
        sorted_group = group.sort_values(value_col, ascending=False).copy()
        tot = sorted_group[value_col].sum()
        # Compute running cumulative fraction
        sorted_group['cumulative_frac'] = sorted_group[value_col].cumsum() / tot
        # Select rows needed to reach at least the threshold
        # The mask is True for all rows up to the first one where cumulative_frac >= threshold
        reach_mask = sorted_group['cumulative_frac'] <= threshold
        # We want at least one row covering the threshold, so include also the first row above it
        if not reach_mask.all():
            first_above = reach_mask.idxmin()
            reach_mask.loc[first_above] = True
        # Keep only these rows
        #print(f'{sample}: {sorted_group.loc[first_above].rds_counted}')
        dfs.append(sorted_group[reach_mask])
    # Concatenate result for all samples
    filtered_df = pd.concat(dfs).drop(columns=['cumulative_frac'])
    return filtered_df

# 
def get_intronic_bcs(exon, maptbl_, bctbl_, bp_from_exon = 10):
    map_tbl = maptbl_.copy()
    bc_tbl = bctbl_.copy()
    assert exon in map_tbl['exon'].unique(), f"{exon} is not present in the mapping table"
    start, end = map_tbl.loc[map_tbl['exon'] == exon, ['ex_genome_start', 'ex_genome_end']].values[0]
    return bc_tbl.loc[(bc_tbl.pos < start - bp_from_exon) | (bc_tbl.pos > end + bp_from_exon)].copy()


def get_bs_stats(bctbl_, vartbl_):
    vartbl = vartbl_.copy()
    bctbl = bctbl_.copy()
    isonames = [ col[4:] for col in bctbl if col.startswith( 'psi_' ) ]
    vartbl = pd.concat( [ pp.bootstrap_varsp_null_distribution( bctbl.loc[bctbl[ 'libname' ] == samp ].set_index( 'bc' ),
                                                                      vartbl.loc[ (vartbl[ 'libname' ] == samp ) ],
                                                                      iso_names = isonames)
                              for samp in vartbl[ 'libname' ].unique().tolist() ],
                              ignore_index = True )
    vartbl = pp.compute_null_zscores( vartbl, 'bs_null', isonames ) 
    vartbl = pd.concat( [ pp.compute_fold_change( vartbl.loc[ vartbl[ 'libname' ] == samp ],
                                                        'wmean_bs_null_',
                                                         'psi_') 
                                                        for samp in vartbl[ 'libname' ].unique().tolist() ],
                                                        ignore_index = True ).sort_values( by = 'pos' ) 
    for iso in isonames:
        
        vartbl[ 'wmean_diff_' + iso ] = vartbl[ 'psi_' + iso + '_singleton_wmean' ] - vartbl[ 'wmean_bs_null_' + iso ]                                                         

    return vartbl

def annotate_stat_sig(bsvartbl_, isonames=None, sample_size=None):
    bsvartbl = bsvartbl_.copy()
    if not isonames:
        isonames = [ col[11:] for col in bsvartbl if col.startswith( 'wmean_diff_' ) ]
    if not sample_size:
        sample_size = bsvartbl.groupby( 'libname' ).size().max()
    bonfer = .05 / sample_size
    for iso in isonames:
        bsvartbl[ f'stat_sig_{iso}' ] = ( bsvartbl[ f'zwmean_bs_null_{iso}' ] >= ss.norm.ppf( 1 - bonfer ) )
    return bsvartbl

## make a wide tbl (i.e. replicate data as columns) from a long table, compute bc_weighted psi, get a stouffers z, annotate stat sig
def make_wide_tbl(bsvartbl_, maptbl_, exon_):
    bsvartbl = bsvartbl_.copy()
    libs = list(bsvartbl['libname'].unique())
    isonames = [ col[11:] for col in bsvartbl if col.startswith( 'wmean_diff_' ) ]
    bs_wide = pp.combine_rep_pervartbls_wide( [ bsvartbl.loc[ bsvartbl[ 'libname' ] == samp ].drop(columns=['libname'])
                                                    for samp in libs ],
                                                  libs,
                                                  indexcols=[ 'chr','pos','ref','alt','var','hgvs_pos','vector_pos','exon' ], 
                                                  group_cols_by_samp = True )
    bs_wide = pp.create_variables_across_samples( bs_wide,
                                                     libs,
                                                     mean_cols = [ 'rna_nbc_varsingleton', 'rna_nrd_varsingleton', 'rna_nrd_ok', 'rna_nrd_bad' ],
                                                     median_cols = [  'rna_nbc_varsingleton', 'rna_nrd_varsingleton', 'rna_nrd_ok', 'rna_nrd_bad' ],
                                                     sum_cols = [ 'rna_nrd_varsingleton', 'rna_nrd_ok', 'rna_nrd_bad' ],
                                                     max_cols = [ 'rna_nbc_varsingleton', 'rna_nrd_varsingleton', 'rna_nrd_ok', 'rna_nrd_bad' ] )
    bs_wide = pp.compute_bc_weighted_psi( bs_wide,
                                             libs,
                                             isonames,
                                             'rna_nbc_varsingleton',
                                              )
    assert exon_ in maptbl_['exon'].unique(), f"{exon_} is not present in the mapping table"
    start, end = maptbl_.loc[maptbl_['exon'] == exon_, ['ex_genome_start', 'ex_genome_end']].iloc[0].tolist()
    bs_wide[ 'is_exon' ] = ( bs_wide.pos >= start ) & ( bs_wide.pos <= end ) 
    bs_wide[ 'is_intron' ] = ( ( bs_wide.pos < start ) | ( bs_wide.pos > end ) )
    
    agg_cols = []
    for iso in isonames:
        agg_cols += [f"fc_{iso}", f"wmean_bs_null_{iso}", f"wmean_diff_{iso}", f"psi_{iso}_singleton_wmean" ]
    
    bs_wide = pp.create_variables_across_samples( bs_wide,
                                                     libs,
                                                     mean_cols = agg_cols,
                                                     median_cols = agg_cols)
    
    bs_wide = pp.stouffers_z(bs_wide,
                                  isonames,
                                  zcol = 'zwmean_bs_null_')
    
    bs_wide = annotate_stat_sig(bs_wide, isonames, len(bs_wide))

    bs_wide['perc_samp_with_var'] = bs_wide[[s + '_rna_nbc_varsingleton' for s in libs]].notna().mean(axis=1)*100

    return bs_wide


def add_splai(widetbl_, splai_path):
    widetbl = widetbl_.copy()
    splai = pd.read_table(splai_path,
                       dtype = { 'chrom': object } )
    splai.rename(columns={'chrom':'chr'}, inplace=True)
    
    index_cols = [ 'chr', 'pos', 'ref', 'alt' ]
    widetbl = widetbl.merge(
        splai,
        how='left',
        on=[ 'chr', 'pos', 'ref', 'alt' ]
    )   

    widetbl['DS_maxm_100'] = widetbl['DS_maxm'] * 100 
    return widetbl


def add_dms(widetbl_, dms_path):
    widetbl = widetbl_.copy()
    dms = pd.read_table(dms_path).rename( columns = { 'Variant': 'protein_var',
                             'LOF score': 'dms_lof' } )
    
    widetbl[ 'hgvs_var' ] = [ f"{p}:{ref}>{alt}" for p, ref, alt in zip(widetbl.hgvs_pos, widetbl.ref, widetbl.alt) ]
    
    widetbl = pp.annotate_protein_hgvs(widetbl, 'NM_000251.2')
    widetbl[ 'protein_var' ] = widetbl.protein_var.apply( lambda x: x.replace( 'p.?', '' ).replace( 'p.(', '' ).replace( ')', '' ) )
    widetbl = widetbl.set_index( 'protein_var' ).merge( dms.set_index( 'protein_var' )[ [ 'dms_lof' ] ],  
                                                                how = 'left',
                                                                left_index = True,
                                                                right_index = True ).reset_index().sort_values( by = [ 'pos', 'alt' ] )
    
    widetbl.protein_var = widetbl.protein_var.fillna( '' )
    widetbl[ 'synon' ] = widetbl.protein_var.str.endswith( '=' )
    widetbl[ 'stop_gain' ] = widetbl.protein_var.str.endswith( '*' )
    widetbl[ 'var_consequence' ] = [ 'synonymous' if sy else 'stop_gain' if sg else 'intronic' if i else 'missense'
                                       for sy, sg, i in zip( widetbl.synon, widetbl.stop_gain, widetbl.is_intron)]
    return widetbl


def add_clinvar(widetbl_, clinvar_path):
    widetbl = widetbl_.copy()
    clinvar = pd.read_table(clinvar_path)
    clinvar['chr'] = clinvar['chr'].apply(lambda x: pp.normalize_chr(x, widetbl['chr'].iloc[0]))
    
    widetbl = widetbl.merge(clinvar,
                     how = 'left',
                     on = [ 'chr', 'pos', 'ref', 'alt' ])
    
    widetbl[ 'clinvar' ] = widetbl[ 'clinvar_interp' ].notnull()
    
    widetbl[ 'clinvar_interp' ] = widetbl.clinvar_interp.fillna( '' )
    widetbl[ 'lit_interp' ] = [ 'Conflicting' if clin == 'Conflicting_classifications_of_pathogenicity' \
                                    else 'VUS' if clin == 'Uncertain_significance' \
                                    else 'B/LB' if 'nign' in clin \
                                    else 'P/LP' if 'ogenic' in clin \
                                    else ''
                                    for clin in widetbl.clinvar_interp ]
    widetbl[ 'clinvar_rev_status' ] = widetbl.clinvar_rev_status.fillna( '' )
    widetbl[ 'clinvar_rev_status' ] = [ rev.split( ',' )[ 1 ].replace( '\'', '' )[ 2: -1 ] if '(' in rev \
                                            else rev for rev in widetbl.clinvar_rev_status ]    
    
    return widetbl


def add_gnomad(widetbl_, gnomad_path):
    widetbl = widetbl_.copy()
    gnomad = pd.read_table(gnomad_path)
    gnomad['CHROM'] = gnomad['CHROM'].apply(lambda x: pp.normalize_chr(x, widetbl['chr'].iloc[0]))
    widetbl = widetbl.merge(gnomad,
                     how = 'left',
                     left_on = [ 'chr', 'pos', 'ref', 'alt' ],
                     right_on = ['CHROM', 'POS', 'REF','ALT']
    )
    return widetbl



def main():
    import argparse

    parser = argparse.ArgumentParser(description='load in varfx tbl and compute statistical signifance for each exon by sample')

    parser.add_argument('--varfx_tbl',  dest='var_tbl' )
    parser.add_argument('--map_tbl',  dest='map_tbl' )
    parser.add_argument('--splai_tbl',  dest='splai_tbl' )
    parser.add_argument('--clinvar_tbl',  dest='clinvar_tbl' )
    parser.add_argument('--dms_tbl',  dest='dms_tbl' )
    parser.add_argument('--gnomad_tbl',  dest='gnomad_tbl' )
    parser.add_argument('--exon', required=True)
    parser.add_argument('--num_bc_per_var',  dest='nbc_per_var', type=int, default=5)
    parser.add_argument('--cumd_cutoff',  dest='cumd_cut' , type=float, default=0.90)
    parser.add_argument('--intronic_bp', help='number of basepair from exon boundary to consider "null"', dest='int_bp', type=int, default=10 )
    parser.add_argument('--out_long', dest='out_long', required=True)
    parser.add_argument('--out_wide', dest='out_wide', required=True)
    parser.add_argument('--out_summary', dest='out_summary', required=True)

    args = parser.parse_args()

    varfx = pd.read_table(args.var_tbl)
    varfx = varfx[varfx['exon'] == args.exon].copy()
    maptbl = pd.read_table(args.map_tbl)
    print('starting')
    vartbl, bctbl = prep_tbls(varfx, maptbl, args.exon, nbc_per_var_min=args.nbc_per_var)
    print('tables prepped')
    bctbl_90 = top_cumulative_rows(bctbl, threshold=args.cumd_cut)
    print('cum_rows done')
    bctbl_intronic = get_intronic_bcs(args.exon, maptbl, bctbl_90, args.int_bp)
    print('got intronic bcs')
    vartbl_bs = get_bs_stats(bctbl_intronic, vartbl)
    print('bs_stats done')
    vartbl_bs_sig = annotate_stat_sig(vartbl_bs)
    print('annotated_stats')
    wide_vartbl = make_wide_tbl(vartbl_bs_sig, maptbl, args.exon)
    
    libs = list(varfx['libname'].unique())
    summ_tbl = wide_vartbl.drop(columns=wide_vartbl.columns[wide_vartbl.columns.str.startswith(tuple(f"{s}_" for s in libs))])
    
    
    if args.splai_tbl:
        print('merging splai')
        summ_tbl = add_splai(summ_tbl, args.splai_tbl)
    if args.clinvar_tbl:
        print('merging clinvar')
        summ_tbl = add_clinvar(summ_tbl, args.clinvar_tbl)
    if args.dms_tbl:
        print('merging dms')
        summ_tbl = add_dms(summ_tbl, args.dms_tbl)
    if args.gnomad_tbl:
        print('merging gnomad')
        summ_tbl = add_gnomad(summ_tbl, args.gnomad_tbl)

    vartbl_bs_sig.to_csv(args.out_long, index=False, sep='\t')
    wide_vartbl.to_csv(args.out_wide, index=False, sep='\t')
    summ_tbl.to_csv(args.out_summary, index=False, sep='\t')


if __name__ == '__main__':
    main()
