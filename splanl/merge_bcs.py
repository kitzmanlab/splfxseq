import os
import pandas as pd
import numpy as np

from collections import OrderedDict as odict  # default is for keys to come back out in order after I think python 3.7

from .coords import pos_to_hgvspos
from .coords import pos_to_gDNA


# make (ordered) dict of lists
def blanktbl(colnames):
    """Make a blank dictionary of column names --> lists
    Can append entries to these lists, then convert it to pandas DataFrame

    Args:
        colnames (list of str): column names, will become keys of resulting dict

    Returns:
        dict of column names, each associated w/ a blank list

    """
    return odict( [ (k,[]) for k in colnames] )


def merge_subasm_and_rna_tbls(
    subasm_tbl,
    rna_tbl,
  ):

    """
    Merge a subassembly and an RNA table; each should be indexed by barcode

    Args:
        subasm_tbl (pd.DataFrame): subassembly results, indexed by barcode seq
        rna_tbl (pd.DataFrame):  RNAseq results (e.g., psi values), indexed by barcode seq

    Returns:
        merged table, in DataFrame, containing the barcodes in both

    """

    # get rid of extra columns
    loc = subasm_tbl
    loc_remove = [ 'ref_target_length','minbq' ] + [c for c in loc if c.startswith('nbp_ge_') ]
    loc = [ c for c in loc if c not in loc_remove ]
    subasm_tbl2 = subasm_tbl[ loc ]

    join_tbl = pd.merge( subasm_tbl2,
                         rna_tbl,
                         how='inner',
                         left_index=True,
                         right_index=True,
                         suffixes=('_subasm','_rna') )

    return join_tbl

def summarize_byvar_singlevaronly(
    subasm_tbl,
    rna_tbl,
    exon_coords,
    min_usable_reads_per_bc,
    isonames=None ):
    """
    Summarize per-variant effects across associated barcodes.
    Considers only single-variant clones; barcodes w/ ≥1 variants are ignored.

    Args:
        subasm_tbl (pd.DataFrame): subassembly results, indexed by barcode seq
        rna_tbl (pd.DataFrame):  RNAseq results (e.g., psi values), indexed by barcode seq
        exon_coords (list of int tuples): coordinates of cloned exons
        min_usable_reads_per_bc (int): min # reads associated with barcode to be considered
        isonames (list of str): names of isoforms; for each entry 'x', a column 'x_psi' should exist in rna_tbl

    Returns:
        pd.DataFrame with per-variant summary values;  mean_x, wmean_x, and median_x are the
        across barcodes mean, read-count-weighted mean, and median psi values for each
        isoform x
    """


    sa_filt = subasm_tbl.query( 'n_variants_passing==1' ).copy()

    li_rna = rna_tbl.index.intersection( sa_filt.index )

    rna_isect = rna_tbl.loc[ li_rna ].copy()

    if isonames is None:
        isonames = [ cn[ :cn.rindex('_') ] for cn in rna_isect.columns if cn.endswith('psi') ]
        assert len(isonames)>0, 'cant infer the isoform name columns; please specify them in parameter isonames'

    rna_isect['varlist'] = sa_filt.loc[ li_rna, 'variant_list' ]

    out_tbl = blanktbl(
        ['chrom','pos','ref','alt','varlist','var_type','n_bc','n_bc_passfilt',
         'sum_reads',
         'sum_usable_reads',
         'sum_unmapped_reads',
         'sum_badstart_reads',
         'sum_otheriso'] +
         [ 'mean_{}'.format(cn) for cn in isonames ] +
         [ 'wmean_{}'.format(cn) for cn in isonames ] +
         [ 'median_{}'.format(cn) for cn in isonames ]
     )

    for singlevar, subtbl in rna_isect.groupby( 'varlist' ):

        subtbl_filt = subtbl.loc[ subtbl.usable_reads >= min_usable_reads_per_bc ].copy()

        out_tbl['varlist'].append(singlevar)
        out_tbl['chrom'].append(singlevar.split(':')[0])
        out_tbl['pos'].append(int(singlevar.split(':')[1]))
        out_tbl['ref'].append(singlevar.split(':')[2])
        out_tbl['alt'].append(singlevar.split(':')[3])

        out_tbl['var_type'].append(None)

        out_tbl['n_bc'].append( subtbl.shape[0] )
        out_tbl['n_bc_passfilt'].append( subtbl_filt.shape[0] )

        if subtbl_filt.shape[0]==0:
            out_tbl['sum_reads'].append( 0 )
            out_tbl['sum_usable_reads'].append( 0 )
            out_tbl['sum_unmapped_reads'].append( 0 )
            out_tbl['sum_badstart_reads'].append( 0 )
            out_tbl['sum_otheriso'].append( 0 )

            for iso in isonames:
                out_tbl[ f'mean_{iso}' ].append( None )
                out_tbl[ f'wmean_{iso}' ].append( None )
                out_tbl[ f'median_{iso}' ].append( None )

            continue
        else:
            out_tbl['sum_reads'].append( subtbl_filt['num_reads'].sum() )
            out_tbl['sum_usable_reads'].append( subtbl_filt['usable_reads'].sum()  )
            out_tbl['sum_unmapped_reads'].append( subtbl_filt['unmapped_reads'].sum()  )
            out_tbl['sum_badstart_reads'].append( subtbl_filt['bad_starts'].sum()  )
            out_tbl['sum_otheriso'].append( subtbl_filt['other_isoform'].sum()  )

            for iso in isonames:
                # mean psi
                out_tbl[ f'mean_{iso}' ].append( subtbl_filt[ f'{iso}_psi' ].mean() )
                # mean psi, weighted by #usable reads
                if subtbl_filt['usable_reads'].sum() != 0:
                    out_tbl[ f'wmean_{iso}' ].append( ( subtbl_filt[ f'{iso}_psi' ] * subtbl_filt['usable_reads'] ).sum() / subtbl_filt['usable_reads'].sum() )
                else:
                    out_tbl[ f'wmean_{iso}' ].append( np.nan )
                # median psi
                out_tbl[ f'median_{iso}' ].append( subtbl_filt[ f'{iso}_psi' ].median() )

    out_tbl = pd.DataFrame( out_tbl )

    return out_tbl

def summarize_byvar_WT(
    subasm_tbl,
    rna_tbl,
    min_usable_reads_per_bc,
    isonames=None ):

    """
    Summarize per-variant effects across associated barcodes.
    Considers only wild type clones; other barcodes are ignored.

    Args:
        subasm_tbl (pd.DataFrame): subassembly results, indexed by barcode seq
        rna_tbl (pd.DataFrame):  RNAseq results (e.g., psi values), indexed by barcode seq
        exon_coords (list of int tuples): coordinates of cloned exons
        min_usable_reads_per_bc (int): min # reads associated with barcode to be considered
        isonames (list of str): names of isoforms; for each entry 'x', a column 'x_psi' should exist in rna_tbl

    Returns:
        pd.DataFrame with per-variant summary values;  mean_x, wmean_x, and median_x are the
        across barcodes mean, read-count-weighted mean, and median psi values for each
        isoform x
    """


    sa_filt = subasm_tbl.query( 'status=="no_variants_input" & n_variants_passing==0' ).copy()

    li_rna = rna_tbl.index.intersection( sa_filt.index )

    rna_isect = rna_tbl.loc[ li_rna ].copy()

    if isonames is None:
        isonames = [ cn for cn in rna_isect.columns if cn.endswith('psi') ]
        assert len(isonames)>0, 'cant infer the isoform name columns; please specify them in parameter isonames'
    out_tbl = blanktbl(
        ['n_bc','n_bc_passfilt',
         'sum_reads',
         'sum_usable_reads',
         'sum_unmapped_reads',
         'sum_badstart_reads',
         'sum_otheriso'] +
         [ 'mean_{}'.format(cn) for cn in isonames ] +
         [ 'wmean_{}'.format(cn) for cn in isonames ] +
         [ 'median_{}'.format(cn) for cn in isonames ] +
         [ 'std_{}'.format(cn) for cn in isonames ]
    )

    rna_isect_filt = rna_isect.loc[ rna_isect.usable_reads >= min_usable_reads_per_bc ].copy()

    out_tbl['n_bc'].append( rna_isect.shape[0] )
    out_tbl['n_bc_passfilt'].append( rna_isect_filt.shape[0] )

    out_tbl['sum_reads'].append( rna_isect_filt['num_reads'].sum() )
    out_tbl['sum_usable_reads'].append( rna_isect_filt['usable_reads'].sum()  )
    out_tbl['sum_unmapped_reads'].append( rna_isect_filt['unmapped_reads'].sum()  )
    out_tbl['sum_badstart_reads'].append( rna_isect_filt['bad_starts'].sum()  )
    out_tbl['sum_otheriso'].append( rna_isect_filt['other_isoform'].sum()  )

    for iso in isonames:
        # mean psi
        out_tbl[ f'mean_{iso}' ].append( rna_isect_filt[ f'{iso}_psi' ].mean() )
        # mean psi, weighted by #usable reads
        out_tbl[ f'wmean_{iso}' ].append( ( rna_isect_filt[ f'{iso}_psi' ] * rna_isect_filt['usable_reads'] ).sum() / rna_isect_filt['usable_reads'].sum() )
        # median psi
        out_tbl[ f'median_{iso}' ].append( rna_isect_filt[ f'{iso}_psi' ].median() )
        # standard deviation psi
        out_tbl[ f'std_{iso}' ].append( rna_isect_filt[ f'{iso}_psi' ].std() )

    out_tbl = pd.DataFrame( out_tbl )

    return out_tbl


####################################
#
# routines to combine replicate per-variant tables
#

def combine_rep_pervartbls_wide(
    ltbls,
    lsampnames,
    indexcols=['chrom','pos','ref','alt','varlist'],
    group_cols_by_samp=False ):

    """
    Combine replicate variant effect tables in wide format

    Args:
        ltbls (list of pd.DataFrame): list of per-variant effect tables, one per replicate or condition
        lsampnames (list of str): list of respective names for those replciates or conditions
        indexcols (list of str): what columns to use to index each variant table
        group_cols_by_samp (bool): should columns from each sample by grouped together

    Returns:
        New pd.DataFrame with by variant effect tables merged together. There may be NAs for variants that are absent from some of the reps/conditions.
    """

    ltbls_ixd = [ tbl.set_index(indexcols).copy() for tbl in ltbls ]

    lcolnames = list(ltbls_ixd[0].columns)

    # all tables must have the same set of columns
    for t in ltbls_ixd:
        assert list(t.columns)==lcolnames

    for (tbl,sampname) in zip(ltbls_ixd,lsampnames):
        tbl.columns = [ '{}_{}'.format( sampname,col ) for col in tbl.columns ]

    tblout = pd.concat( ltbls_ixd, axis=1 )

    if group_cols_by_samp:
        loc = []
        for col in lcolnames:
            loc += [ '{}_{}'.format(sampname,col) for sampname in lsampnames  ]
        tblout = tblout[loc]


    tblout = tblout.reset_index()

    return tblout


def combine_rep_pervartbls_long(
    ltbls,
    lsampnames,
    indexcols=['chrom','pos','ref','alt','varlist'],
):

    """
    Combine replicate variant effect tables in long format

    Args:
        ltbls (list of pd.DataFrame): list of per-variant effect tables, one per replicate or condition
        lsampnames (list of str): list of respective names for those replciates or conditions
        indexcols (list of str): what columns to use to index each variant table

    Returns:
        New pd.DataFrame with by variant effect tables merged together, with each replicate appearing as a separate row
    """

    ltbls_ixd = [ tbl.set_index(indexcols).copy() for tbl in ltbls ]

    lcolnames = list(ltbls_ixd[0].columns)

    # all tables must have the same set of columns
    for t in ltbls_ixd:
        assert list(t.columns)==lcolnames

    for (tbl,sampname) in zip(ltbls_ixd,lsampnames):
        tbl['sample']=sampname

    tblout = pd.concat( ltbls_ixd, axis=0 )

    tblout = tblout[ ['sample']+[cn for cn in lcolnames] ]

    tblout = tblout.reset_index()

    return tblout

def filter_byvartbl_snvonly(
    byvar_tbl
):
    """
    Filter by variant table to only SNVs

    Args:
        byvar_tbl (pd.DataFrame): per-variant effect table
    Returns:
        Copy of per-variant effect table with only SNV lines included
    """

    byvar_snvonly = byvar_tbl.loc[ (byvar_tbl.ref.str.len() == 1) & (byvar_tbl.alt.str.len() == 1) ].copy()
    return byvar_snvonly

def across_sample_stats(ltbls,
                        lsampnames,
                        med_col_names):

    out_tbl = { 'sample_group':[],
                'sample':[],
                'n_reads':[],
                'n_bcs':[],
                'n_bcs_passfilt':[]}

    for col in med_col_names:
        out_tbl['med_'+col] = []

    i=0
    for grp, _lsamp in lsampnames.items():
        for lsamp in _lsamp:
            out_tbl['sample_group'].append(grp)
            out_tbl['sample'].append(grp+'_'+lsamp)
            out_tbl['n_reads'].append( int( ltbls[ i ].query( 'sample=="%s"' % lsamp ).sum_reads.sum() ) )
            out_tbl['n_bcs'].append( int( ltbls[ i ].query( 'sample=="%s"' % lsamp ).n_bc.sum() ) )
            out_tbl['n_bcs_passfilt'].append( int( ltbls[ i ].query( 'sample=="%s"' % lsamp ).n_bc_passfilt.sum() ) )
            for col in med_col_names:
                out_tbl['med_'+col].append( float( ltbls[ i ].query( 'sample=="%s"' % lsamp )[ col ].median() ) )
        i+=1

    out_tbl = pd.DataFrame( out_tbl )

    return out_tbl
