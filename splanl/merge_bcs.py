import pandas as pd
import numpy as np
import pysam
import time
import scipy.stats as ss
from collections import OrderedDict as odict  # default is for keys to come back out in order after I think python 3.7
from collections import Counter
import splanl.plots as sp
import splanl.junction_scorer as jn

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

def compute_psi_values(
    in_df,
    iso_col = None,
    read_count_col = 'usable_reads'
):

    out_df = in_df.copy()

    if not iso_col:
        iso_col = [ col for col in out_df.columns if col.startswith( 'iso' ) ]

    assert iso_col, 'Please specify columns for PSI'

    assert read_count_col in out_df.columns, '%s is not within %s' % ( read_count_col, str( in_df ) )

    for col in iso_col:
        out_df[ col + '_psi' ] = out_df[ col ] / out_df[ read_count_col ]

    return out_df

def summarize_byvar_singlevaronly(
    subasm_tbl,
    rna_tbl,
    min_usable_reads_per_bc,
    isonames=None ):
    """
    Summarize per-variant effects across associated barcodes.
    Considers only single-variant clones; barcodes w/ ≥1 variants are ignored.

    Args:
        subasm_tbl (pd.DataFrame): subassembly results, indexed by barcode seq
        rna_tbl (pd.DataFrame):  RNAseq results (e.g., psi values), indexed by barcode seq
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

    rna_isect_psi = compute_psi_values( rna_isect, iso_col = isonames )

    if isonames is None:
        isonames = [ cn[ :cn.rindex('_') ] for cn in rna_isect_psi.columns if cn.endswith('psi') ]
        assert len(isonames)>0, 'cant infer the isoform name columns; please specify them in parameter isonames'

    rna_isect_psi['varlist'] = sa_filt.loc[ li_rna, 'variant_list' ]

    out_tbl = blanktbl(
        ['chrom','pos','ref','alt','varlist','n_bc','n_bc_passfilt',
         'sum_reads',
         'sum_reads_passfilt',
         'sum_usable_reads',
         'sum_unmapped_reads',
         'sum_badstart_reads',
         'sum_badend_reads',
         'sum_softclipped_reads',
         'sum_otheriso'] +
         [ 'mean_{}'.format(cn) for cn in isonames ] +
         [ 'wmean_{}'.format(cn) for cn in isonames ] +
         [ 'median_{}'.format(cn) for cn in isonames ]
     )

    for singlevar, subtbl in rna_isect_psi.groupby( 'varlist' ):

        subtbl_filt = subtbl.loc[ subtbl.usable_reads > min_usable_reads_per_bc ].copy()

        out_tbl['varlist'].append(singlevar)
        out_tbl['chrom'].append(singlevar.split(':')[0])
        out_tbl['pos'].append(int(singlevar.split(':')[1]))
        out_tbl['ref'].append(singlevar.split(':')[2])
        out_tbl['alt'].append(singlevar.split(':')[3])

        out_tbl['n_bc'].append( subtbl.shape[0] )
        out_tbl['n_bc_passfilt'].append( subtbl_filt.shape[0] )

        out_tbl['sum_reads'].append( subtbl['num_reads'].sum() )

        if subtbl_filt.shape[0]==0:
            out_tbl['sum_reads_passfilt'].append( 0 )
            out_tbl['sum_usable_reads'].append( 0 )
            out_tbl['sum_unmapped_reads'].append( 0 )
            out_tbl['sum_badstart_reads'].append( 0 )
            out_tbl['sum_badend_reads'].append( 0 )
            out_tbl['sum_softclipped_reads'].append( 0 )
            out_tbl['sum_otheriso'].append( 0 )

            for iso in isonames:
                out_tbl[ f'mean_{iso}' ].append( None )
                out_tbl[ f'wmean_{iso}' ].append( None )
                out_tbl[ f'median_{iso}' ].append( None )

            continue
        #this is tricky to think about
        #currently set so that its counting reads after removing the barcodes not passing the filter
        else:
            out_tbl['sum_reads_passfilt'].append( subtbl_filt['num_reads'].sum() )
            out_tbl['sum_usable_reads'].append( subtbl_filt['usable_reads'].sum()  )
            out_tbl['sum_unmapped_reads'].append( subtbl_filt['unmapped_reads'].sum()  )
            out_tbl['sum_badstart_reads'].append( subtbl_filt['bad_starts'].sum()  )
            out_tbl['sum_badend_reads'].append( subtbl_filt['bad_ends'].sum()  )
            out_tbl['sum_softclipped_reads'].append( subtbl_filt['soft_clipped'].sum()  )
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

    out_tbl = count_bcs_per_var_sa( out_tbl,
                                    sa_filt )

    #these two are have the total barcode/read count in the denominator
    out_tbl['per_bc_passfilt'] = 100*( out_tbl.n_bc_passfilt / out_tbl.n_bc )
    out_tbl['per_reads_passfilt'] = 100*( out_tbl.sum_reads_passfilt / out_tbl.sum_reads )

    #these columns are based of barcodes which are passing the filter
    #so only reads from barcodes passing the filter are used in the denominator
    out_tbl['per_reads_usable'] = 100*( out_tbl.sum_usable_reads / out_tbl.sum_reads_passfilt )
    out_tbl['per_unmapped'] = 100*( out_tbl.sum_unmapped_reads / out_tbl.sum_reads_passfilt )
    out_tbl['per_badend'] = 100*( out_tbl.sum_badend_reads / out_tbl.sum_reads_passfilt )
    out_tbl['per_badstart'] = 100*( out_tbl.sum_badstart_reads / out_tbl.sum_reads_passfilt )
    out_tbl['per_softclipped'] = 100*( out_tbl.sum_softclipped_reads / out_tbl.sum_reads_passfilt )
    out_tbl['per_otheriso'] = 100*( out_tbl.sum_otheriso / out_tbl.sum_reads_passfilt )


    return out_tbl

def summarize_byvar_singlevaronly_pe( subasm_tbl,
                                      rna_tbl,
                                      min_usable_reads_per_bc,
                                      summary_cols,
                                      isonames = None, ):
    """
    Summarize per-variant effects across associated barcodes.
    Considers only single-variant clones; barcodes w/ ≥1 variants are ignored.

    Args:
        subasm_tbl (pd.DataFrame): subassembly results, indexed by barcode seq
        rna_tbl (pd.DataFrame):  RNAseq results (e.g., psi values), indexed by barcode seq
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

    rna_isect_psi = compute_psi_values( rna_isect, iso_col = isonames )

    if isonames is None:
        isonames = [ cn[ :cn.rindex( '_' ) ] for cn in rna_isect_psi.columns if cn.endswith( 'psi' ) ]
        assert len( isonames ) > 0, 'cant infer the isoform name columns; please specify them in parameter isonames'

    rna_isect_psi[ 'varlist' ] = sa_filt.loc[ li_rna, 'variant_list' ]

    out_tbl = blanktbl( ['chrom','pos','ref','alt','varlist','n_bc','n_bc_passfilt','sum_reads','sum_reads_passfilt', ] +
                        [ 'sum_{}'.format( cn ) for cn in summary_cols ] +
                        [ 'mean_{}'.format(cn) for cn in isonames ] +
                        [ 'wmean_{}'.format(cn) for cn in isonames ] +
                        [ 'median_{}'.format(cn) for cn in isonames ] )

    for singlevar, subtbl in rna_isect_psi.groupby( 'varlist' ):

        subtbl_filt = subtbl.loc[ subtbl.usable_reads > min_usable_reads_per_bc ].copy()

        out_tbl['varlist'].append(singlevar)
        out_tbl['chrom'].append(singlevar.split(':')[0])
        out_tbl['pos'].append(int(singlevar.split(':')[1]))
        out_tbl['ref'].append(singlevar.split(':')[2])
        out_tbl['alt'].append(singlevar.split(':')[3])

        out_tbl['n_bc'].append( subtbl.shape[0] )
        out_tbl['n_bc_passfilt'].append( subtbl_filt.shape[0] )

        out_tbl['sum_reads'].append( subtbl['num_reads'].sum() )

        if subtbl_filt.shape[0]==0:
            out_tbl['sum_reads_passfilt'].append( 0 )

            for col in summary_cols:
                out_tbl[ f'sum_{col}' ].append( 0 )

            for iso in isonames:
                out_tbl[ f'mean_{iso}' ].append( None )
                out_tbl[ f'wmean_{iso}' ].append( None )
                out_tbl[ f'median_{iso}' ].append( None )

            continue
        #this is tricky to think about
        #currently set so that its counting reads after removing the barcodes not passing the filter
        else:
            out_tbl['sum_reads_passfilt'].append( subtbl_filt['num_reads'].sum() )

            for col in summary_cols:
                out_tbl[ f'sum_{col}' ].append( subtbl_filt[ col ].sum()  )

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

    out_tbl = count_bcs_per_var_sa( out_tbl,
                                    sa_filt )

    #these two are have the total barcode/read count in the denominator
    out_tbl['per_bc_passfilt'] = 100*( out_tbl.n_bc_passfilt / out_tbl.n_bc )
    out_tbl['per_reads_passfilt'] = 100*( out_tbl.sum_reads_passfilt / out_tbl.sum_reads )

    #these columns are based of barcodes which are passing the filter
    #so only reads from barcodes passing the filter are used in the denominator
    for col in summary_cols:
        out_tbl[ f'per_{col}' ] = 100*( out_tbl[ f'sum_{col}' ] / out_tbl.sum_reads_passfilt )

    return out_tbl

def count_bcs_per_var_sa( rna_tbl,
                          satbl ):

    #count the number of rows for each variant
    #this is the number of barcodes for that variant
    count_tbl = satbl.copy().groupby( [ 'variant_list' ] )[ 'variant_list' ].count().rename('n_bc_sa')

    rna_tbl = rna_tbl.copy().set_index( 'varlist' )

    out_tbl = pd.merge( rna_tbl, count_tbl, left_index = True, right_index = True )

    out_tbl.index.name = 'varlist'

    out_tbl = out_tbl.reset_index()

    assert rna_tbl.shape[0] == out_tbl.shape[0], 'RNA table rows were lost in the merge'

    return out_tbl

def summarize_byvar_WT( rna_tbl,
                        exon_coords,
                        min_usable_reads_per_bc,
                        chrom,
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

    rna_psi = compute_psi_values( rna_tbl, iso_col = isonames )

    if isonames is None:
        isonames = [ cn[ :cn.rindex('_') ] for cn in rna_psi.columns if cn.endswith('psi') ]
        assert len(isonames)>0, 'cant infer the isoform name columns; please specify them in parameter isonames'

    rna_psi_filt = rna_psi.loc[ rna_psi.usable_reads > min_usable_reads_per_bc ].copy()

    out_tbl = {}

    out_tbl[ 'chrom' ] = [ chrom ]
    out_tbl[ 'varlist' ] = [ 'WT' ]

    out_tbl[ 'n_bc' ] = [ rna_psi.shape[ 0 ] ]
    out_tbl[ 'n_bc_passfilt' ] = [ rna_psi_filt.shape[ 0 ] ]

    out_tbl[ 'sum_reads' ] = [ rna_psi.num_reads.sum() ]

    if rna_psi_filt.shape[ 0 ] == 0:

        out_tbl[ 'sum_reads_passfilt' ] = [ 0 ]
        out_tbl[ 'sum_usable_reads' ] = [ 0 ]
        out_tbl[ 'sum_unmapped_reads' ] = [ 0 ]
        out_tbl[ 'sum_badstart_reads' ] = [ 0 ]
        out_tbl[ 'sum_badend_reads' ] = [ 0 ]
        out_tbl[ 'sum_softclipped_reads' ] = [ 0 ]
        out_tbl[ 'sum_otheriso' ] = [ 0 ]

        for iso in isonames:
            out_tbl[ f'mean_{iso}' ] = [ None ]
            out_tbl[ f'wmean_{iso}' ] = [ None ]
            out_tbl[ f'median_{iso}' ] = [ None ]
            out_tbl[ f'stdev_{iso}' ] = [ None ]

    else:

        out_tbl[ 'sum_reads_passfilt' ] = [ rna_psi_filt.num_reads.sum() ]
        out_tbl[ 'sum_usable_reads' ] = [ rna_psi_filt.usable_reads.sum() ]
        out_tbl[ 'sum_unmapped_reads' ] = [ rna_psi_filt.unmapped_reads.sum() ]
        out_tbl[ 'sum_badstart_reads' ] = [ rna_psi_filt.bad_starts.sum() ]
        out_tbl[ 'sum_badend_reads' ] = [ rna_psi_filt.bad_ends.sum() ]
        out_tbl[ 'sum_softclipped_reads' ] = [ rna_psi_filt.soft_clipped.sum() ]
        out_tbl[ 'sum_otheriso' ] = [ rna_psi_filt.other_isoform.sum() ]

        for iso in isonames:
            # mean psi
            out_tbl[ f'mean_{iso}' ] = [ rna_psi_filt[ f'{iso}_psi' ].mean() ]
            # mean psi, weighted by #usable reads
            if rna_psi_filt.usable_reads.sum() != 0:
                out_tbl[ f'wmean_{iso}' ] = [ ( rna_psi_filt[ f'{iso}_psi' ] * rna_psi_filt.usable_reads ).sum() / rna_psi_filt.usable_reads.sum() ]
            else:
                out_tbl[ f'wmean_{iso}' ] = [ np.nan ]
                # median psi
            out_tbl[ f'median_{iso}' ] = [ rna_psi_filt[ f'{iso}_psi' ].median() ]
            out_tbl[ f'stdev_{iso}' ] = [ rna_psi_filt[ f'{iso}_psi' ].std() ]
            out_tbl[ f'wstdev_{iso}' ] = [ np.sqrt( ( rna_psi_filt.usable_reads  * ( rna_psi_filt[ f'{iso}_psi' ] - out_tbl[ f'wmean_{iso}' ][ 0 ] )**2 ).sum() \
                                                    / ( rna_psi_filt.usable_reads.sum() - 1 ) ) ]

    out_tbl = pd.DataFrame( out_tbl )

    #these two are have the total barcode/read count in the denominator
    out_tbl['per_bc_passfilt'] = 100*( out_tbl.n_bc_passfilt / out_tbl.n_bc )
    out_tbl['per_reads_passfilt'] = 100*( out_tbl.sum_reads_passfilt / out_tbl.sum_reads )

    #these columns are based of barcodes which are passing the filter
    #so only reads from barcodes passing the filter are used in the denominator
    out_tbl['per_reads_usable'] = 100*( out_tbl.sum_usable_reads / out_tbl.sum_reads_passfilt )
    out_tbl['per_unmapped'] = 100*( out_tbl.sum_unmapped_reads / out_tbl.sum_reads_passfilt )
    out_tbl['per_badend'] = 100*( out_tbl.sum_badend_reads / out_tbl.sum_reads_passfilt )
    out_tbl['per_badstart'] = 100*( out_tbl.sum_badstart_reads / out_tbl.sum_reads_passfilt )
    out_tbl['per_softclipped'] = 100*( out_tbl.sum_softclipped_reads / out_tbl.sum_reads_passfilt )
    out_tbl['per_otheriso'] = 100*( out_tbl.sum_otheriso / out_tbl.sum_reads_passfilt )

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

def create_variables_across_samples( wide_tbl,
                                     lsampnames,
                                     median_cols = [],
                                     mean_cols = [],
                                     sum_cols = [],
                                     max_cols = [] ):

    wide = wide_tbl.copy()

    if median_cols:

        for col in median_cols:

            samp_cols = [ '_'.join( [ samp, col ] ) for samp in lsampnames ]

            wide[ col + '_med' ] = wide[ samp_cols ].median( axis = 1 )

    if mean_cols:

        for col in mean_cols:

            samp_cols = [ '_'.join( [ samp, col ] ) for samp in lsampnames ]

            wide[ col + '_mean' ] = wide[ samp_cols ].mean( axis = 1 )

    if sum_cols:

        for col in sum_cols:

            samp_cols = [ '_'.join( [ samp, col ] ) for samp in lsampnames ]

            wide[ col + '_sum' ] = wide[ samp_cols ].sum( axis = 1 )

    if max_cols:

        for col in max_cols:

            samp_cols = [ '_'.join( [ samp, col ] ) for samp in lsampnames ]

            wide[ col + '_max' ] = wide[ samp_cols ].max( axis = 1 )

    return wide

def compute_bc_weighted_psi( wide_tbl,
                             lsampnames,
                             isonames,
                             bccount,
                              ):

    wide = wide_tbl.copy()

    samp_bc_cols = [ '_'.join( [ samp, bccount ] ) for samp in lsampnames ]

    for col in isonames:

        #we tend to use the wmean more often so this is intentionally a mean of the wmeans
        wide[ 'mean_' + col ] = wide[ [ '_'.join( [ samp, 'wmean', col ] ) for samp in lsampnames ] ].mean( axis = 1 )

        #this probably would look better as a numpy array dot product but we survived
        wide[ 'wmean_' + col ] = pd.DataFrame( ( wide[ '_'.join( [ samp, 'wmean', col ] ) ] * wide[ '_'.join( [ samp, bccount ] ) ]
                                                   for samp in lsampnames ) ).T.sum( axis = 1) \
                                             / wide[ samp_bc_cols ].sum( axis = 1 )

    return wide

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

        #allows us to create long tables from long tables that already include a sample column
        if 'sample' in lcolnames:
            tbl['sample_grp']=sampname
        else:
            tbl['sample']=sampname

    tblout = pd.concat( ltbls_ixd, axis=0 )

    tblout = tblout[ ['sample_grp' if 'sample' in lcolnames else 'sample']+[cn for cn in lcolnames] ]

    tblout = tblout.reset_index()

    return tblout

def combine_rep_perbctbls_long(
    ltbls,
    lsampnames
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

    ltbls_ixd = [ tbl.copy() for tbl in ltbls ]

    lcolnames = list(ltbls_ixd[0].columns)

    # all tables must have the same set of columns
    for t in ltbls_ixd:
        assert list(t.columns)==lcolnames

    for (tbl,sampname) in zip(ltbls_ixd,lsampnames):

        #allows us to create long tables from long tables that already include a sample column
        if 'sample' in lcolnames:
            tbl['sample_grp']=sampname
        else:
            tbl['sample']=sampname

    tblout = pd.concat( ltbls_ixd, axis=0 )

    tblout = tblout[ ['sample_grp' if 'sample' in lcolnames else 'sample']+[cn for cn in lcolnames] ]

    tblout.index.name = 'barcode'

    return tblout

def combine_allisos_pervartbls_long(
                                    ltbls,
                                    lsampnames,
                                    indexcols=['chrom','pos','ref','alt','varlist'],
                                    ):

    """
    Combine replicate variant effect tables with all isoforms (not necessarily matching column names) in long format

    Args:
        ltbls (list of pd.DataFrame): list of per-variant effect tables with all isoforms, one per replicate or condition
        lsampnames (list of str): list of respective names for those replciates or conditions
        indexcols (list of str): what columns to use to index each variant table

    Returns:
        New pd.DataFrame with by variant effect tables merged together, with each replicate appearing as a separate row
        Columns not represented in one input dataframe compared to other input dataframes will contain nan values
    """

    lcolnames = list( set( [ col for tbl in ltbls for col in tbl ] ) )

    ltbls_ixd = [ tbl.set_index( indexcols ).copy() for tbl in ltbls ]

    for (tbl,sampname) in zip(ltbls_ixd,lsampnames):

        #allows us to create long tables from long tables that already include a sample column
        if 'sample' in lcolnames:
            tbl['sample_grp']=sampname
        else:
            tbl['sample']=sampname

    tblout = ltbls_ixd[0].append( ltbls_ixd[1:], sort = True)

    tblout = tblout[ ['sample_grp' if 'sample' in lcolnames else 'sample']+[cn for cn in lcolnames if cn not in indexcols ] ].reset_index()

    return tblout

def combine_allisos_perbctbls_long(
    ltbls,
    lsampnames
):

    """
    Combine replicate barcode effect tables with all isoforms (not necessarily matching column names) in long format

    Args:
        ltbls (list of pd.DataFrame): list of per-variant effect tables with all isoforms, one per replicate or condition
        lsampnames (list of str): list of respective names for those replciates or conditions
        indexcols (list of str): what columns to use to index each variant table

    Returns:
        New pd.DataFrame with by variant effect tables merged together, with each replicate appearing as a separate row
        Columns not represented in one input dataframe compared to other input dataframes will contain nan values
    """

    lcolnames = list( set( [ col for tbl in ltbls for col in tbl ] ) )

    ltbls_ixd = [ tbl.copy() for tbl in ltbls ]

    for (tbl,sampname) in zip(ltbls_ixd,lsampnames):

        #allows us to create long tables from long tables that already include a sample column
        if 'sample' in lcolnames:
            tbl['sample_grp']=sampname
        else:
            tbl['sample']=sampname

    tblout = ltbls_ixd[0].append( ltbls_ixd[1:], sort = True)

    tblout = tblout[ ['sample_grp' if 'sample' in lcolnames else 'sample']+[cn for cn in lcolnames] ]

    #makes all missing columns 0's
    tblout = tblout.fillna(0)

    tblout.index.name = 'barcode'

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

def create_sparse_df( large_df,
                        sparse_val = 'nan' ):

    """Make sure the large_df has no meaningful index columns - ie enter as df.reset_index()"""

    in_df = large_df.copy()

    non_num_col = list( in_df.select_dtypes(exclude=np.number).columns )

    #the sparse dataframe can't have any non-numerical columns so lets set them all as indices
    out_df = in_df.set_index( non_num_col ).select_dtypes(include=np.number)

    if sparse_val == 'nan':
        out_df = out_df.astype( pd.SparseDtype( "float", np.nan ) )
    elif isinstance(sparse_val, int):
        out_df = out_df.astype( pd.SparseDtype( "float", sparse_val ) )

    return out_df.reset_index()

def process_bcs_wrapper( sample,
                         pysam_align,
                         all_isogrpdict,
                         named_isogrpdict,
                         cnst_exons,
                         satbl,
                         spl_tol = 3,
                         indel_tol = 20,
                         min_matches_for = 70,
                         min_matches_rev = 50,
                         tol_first_last = 0,
                         min_reads = 1,
                         other_isos_usable = False,
                         bc_tag = 'RX',
                         max_bc_len = 30,
                         unmapped_pysam = None,
                         unmap_bc_split_char = '_BC=',
                         unmap_col = 'unmapped_reads',
                         read_tot_col = 'num_reads',
                         usable_read_col = 'usable_reads',
                         other_cols = [ 'secondary_reads', 'unpaired_reads', 'bad_starts', 'bad_ends', 'soft_clipped', 'other_isoform' ],
                         waterfall_thresh = [ 75, 95 ],
                         bc_read_cut_thresh = 95,
                         verbose = False,
                         cml = False,
                         plot_out = None, ):

    assert bc_read_cut_thresh in waterfall_thresh, 'BC read cut off thresholds must be a subset of waterfall thresholds!'

    t0 = time.time()

    print( 'Processing sample %s' % sample )

    msamp_bcrnatbl = jn.compute_isoform_counts_pe( pysam_align,
                                                all_isogrpdict,
                                                cnst_exons,
                                                spl_tol = spl_tol,
                                                indel_tol = indel_tol,
                                                min_matches_for = min_matches_for,
                                                min_matches_rev = min_matches_rev,
                                                bc_tag = bc_tag,
                                                verbose = verbose )

    pysam_align.close()

    msamp_bcrnatbl_flen = filter_on_barc_len( msamp_bcrnatbl,
                                              max_len = max_bc_len )

    if unmapped_pysam:

        msamp_bcrnatbl_flen = merge_unmapped_bcs( msamp_bcrnatbl_flen,
                                                  unmapped_pysam,
                                                  bc_split_char = unmap_bc_split_char,
                                                  unmap_col = unmap_col,
                                                  read_tot_col = read_tot_col )

    msamp_bcrnatbl_rename = combine_isogrps_pe( named_isogrpdict,
                                                msamp_bcrnatbl_flen,
                                                keep_cols = [ read_tot_col, usable_read_col, unmap_col ] + other_cols  )

    msamp_varbcrnatbl_flen_allisos = merge_subasm_and_rna_tbls( satbl,
                                                                msamp_bcrnatbl_flen )

    x_cuts,y_cuts = sp.waterfall_plot( msamp_bcrnatbl_flen,
                                    usable_read_col,
                                    waterfall_thresh,
                                    verbose = verbose,
                                    cml = cml,
                                    savefig = plot_out )

    cut_d = { 'x': x_cuts, 'y': y_cuts }

    msamp_byvartbl_allisos = summarize_byvar_singlevaronly_pe( satbl,
                                                               msamp_bcrnatbl_flen,
                                                               cut_d[ 'y' ][ bc_read_cut_thresh ],
                                                               [ usable_read_col, unmap_col ] + other_cols )

    t1 = time.time()

    print( 'Finished processing sample %s in %.2f minutes' % ( sample, ( t1 - t0 ) / 60 ) )

    return { 'msamp_bcrnatbl_rename': msamp_bcrnatbl_rename,
             'msamp_varbcrnatbl_flen_allisos': msamp_varbcrnatbl_flen_allisos,
             'msamp_byvartbl_allisos': msamp_byvartbl_allisos,
             'bc_read_cutoffs': cut_d }

def filter_on_barc_len( jxnbybctbl, max_len=35 ):
    li = jxnbybctbl.index.str.len() <= max_len
    subtbl = jxnbybctbl.loc[li].copy()
    return subtbl

def merge_unmapped_bcs( tbl_by_bc,
                        unmapped_m1,
                        bc_split_char = '_BC=',
                        unmap_col = 'unmapped_reads',
                        read_tot_col = 'num_reads'
                        ):

    tbb = tbl_by_bc.drop( columns = unmap_col ).copy()

    bc_cnt_m1 = pd.Series( Counter( [ entry.name.split( bc_split_char )[ -1 ] for entry in unmapped_m1 ] ),
                           name = unmap_col )

    unmapped_m1.close()

    bc_intsect = bc_cnt_m1.index.intersection( tbb.index )

    bc_cnt_m1 = bc_cnt_m1.loc[ bc_intsect ].copy()

    tbb[ unmap_col ] = bc_cnt_m1

    tbb[ unmap_col ] = tbb[ unmap_col ].fillna( value = 0 ).astype( int )

    tbb[ read_tot_col ] += tbb[ unmap_col ]

    return tbb

def combine_isogrps_pe( new_grpnames_to_old_grpnames,
                        jxnbybctbl,
                        keep_cols = ['num_reads','secondary_reads', 'unpaired_reads', 'unmapped_reads','bad_starts','bad_ends','soft_clipped','other_isoform','usable_reads'] ):
    """ combine barcode groups
    Arguments:
        new_grpnames_to_old_grpnames {[type]} -- [description]
    """

    chg_cols = [ isocol for isocol_l in new_grpnames_to_old_grpnames.values()
                        for isocol in isocol_l ]

    other_cols = [ col for col in jxnbybctbl.columns
                   if col not in chg_cols ]

    if keep_cols != other_cols:

        missing_cols = list( set( other_cols ).difference( set( keep_cols ) ) )

        print( '%i columns will be removed during this process.\nColumns: ' % len( missing_cols ), missing_cols )

    newtbl = pd.DataFrame()
    for c in keep_cols:
        newtbl[c] = jxnbybctbl[c]

    for newgrp in new_grpnames_to_old_grpnames:
        oldgrp = new_grpnames_to_old_grpnames[ newgrp ]

        if type(oldgrp)==str:
            oldgrp=[oldgrp]
        else:
            assert type(oldgrp)==list

        if len(oldgrp)==1:
            newtbl[newgrp] = jxnbybctbl[oldgrp]
        else:
            newtbl[newgrp] = jxnbybctbl[oldgrp].sum(axis=1)

    for newgrp in new_grpnames_to_old_grpnames:
        newtbl[newgrp+'_psi'] = newtbl[newgrp] / newtbl['usable_reads']

    #fills in missing PSI values with 0 so we can create sparse datasets
    newtbl = newtbl.fillna(0)

    return newtbl

def combine_isogrps( new_grpnames_to_old_grpnames,
                     jxnbybctbl,
                     keep_cols = ['num_reads','unmapped_reads','bad_starts','bad_ends','soft_clipped','other_isoform','usable_reads'],
                     ):
    """ combine barcode groups
    Arguments:
        new_grpnames_to_old_grpnames {[type]} -- [description]
    """

    chg_cols = [ isocol for isocol_l in new_grpnames_to_old_grpnames.values()
                        for isocol in isocol_l ]

    other_cols = [ col for col in jxnbybctbl.columns
                   if col not in chg_cols ]

    if keep_cols != other_cols:

        missing_cols = list( set( other_cols ).difference( set( keep_cols ) ) )

        print( '%i columns will be removed during this process.\nColumns: ' % len( missing_cols ), missing_cols )

    newtbl = pd.DataFrame()
    for c in keep_cols:
        newtbl[c] = jxnbybctbl[c]

    for newgrp in new_grpnames_to_old_grpnames:
        oldgrp = new_grpnames_to_old_grpnames[ newgrp ]

        if type(oldgrp)==str:
            oldgrp=[oldgrp]
        else:
            assert type(oldgrp)==list

        if len(oldgrp)==1:
            newtbl[newgrp] = jxnbybctbl[oldgrp]
        else:
            newtbl[newgrp] = jxnbybctbl[oldgrp].sum(axis=1)

    for newgrp in new_grpnames_to_old_grpnames:
        newtbl[newgrp+'_psi'] = newtbl[newgrp] / newtbl['usable_reads']

    #fills in missing PSI values with 0 so we can create sparse datasets
    newtbl = newtbl.fillna(0)

    return newtbl

def create_read_count_df( bysamp_bc_dict,
                          thresholds,
                          read_cutoff_key = 'bc_read_cutoffs' ):

    out_dict = { 'sample': [] }

    for thresh in thresholds:

        out_dict[ str( thresh ) + '_x' ] = []
        out_dict[ str( thresh ) + '_y' ] = []

    for samp in bysamp_bc_dict:

        out_dict[ 'sample' ].append( samp )

        for thresh in thresholds:

            out_dict[ str( thresh ) + '_x' ].append( 10**( bysamp_bc_dict[ samp ][ read_cutoff_key ][ 'x' ][ thresh ] ) )
            out_dict[ str( thresh ) + '_y' ].append( bysamp_bc_dict[ samp ][ read_cutoff_key ][ 'y' ][ thresh ] )

    outdf = pd.DataFrame( out_dict )

    for col in outdf.columns:

        if col.endswith( '_x' ):

            outdf[ col + '_log10' ] = np.log10( outdf[ col ].tolist() )

    return outdf

def across_sample_stats(ltbls,
                        lsampnames,
                        med_col_names):

    out_tbl = { 'sample_group':[],
                'sample':[],
                'psbl_var':[],
                'n_var':[],
                'n_var_ex':[],
                'n_var_in':[],
                'n_reads':[],
                'n_reads_passfilt':[],
                'n_usable_reads':[],
                'n_bc':[],
                'n_bc_passfilt':[],
                'n_unmapped':[],
                'n_badstart':[],
                'n_badend':[],
                'n_softclip':[],
                'n_otheriso':[],
                }

    for col in med_col_names:
        out_tbl['med_'+col] = []

    i=0
    for grp, _lsamp in lsampnames.items():
        for lsamp in _lsamp:

            lsamp_df = ltbls[ i ].query( 'sample=="%s"' % lsamp ).copy()
            lsamp_filt_df = lsamp_df.query( 'n_bc_passfilt > 0' )

            out_tbl['sample_group'].append( grp )
            out_tbl['sample'].append( grp+'_'+lsamp )
            out_tbl['psbl_var'].append( 3*( lsamp_df.pos.max() - lsamp_df.pos.min() ) )
            out_tbl['n_var'].append( int( lsamp_filt_df.shape[0] ) )
            out_tbl['n_var_ex'].append( int( lsamp_filt_df.loc[ lsamp_filt_df.exon ].shape[0] ) )
            out_tbl['n_var_ex'].append( int( lsamp_filt_df.loc[ ~lsamp_filt_df.exon ].shape[0] ) )
            out_tbl['n_reads'].append( int( lsamp_df.sum_reads.sum() ) )
            out_tbl['n_reads_passfilt'].append( int( lsamp_df.sum_reads_passfilt.sum() ) )
            out_tbl['n_usable_reads'].append( int( lsamp_df.sum_usable_reads.sum() ) )
            out_tbl['n_bc'].append( int( lsamp_df.n_bc.sum() ) )
            out_tbl['n_bc_passfilt'].append( int( lsamp_df.n_bc_passfilt.sum() ) )
            out_tbl['n_unmapped'].append( int( lsamp_df.sum_unmapped_reads.sum() ) )
            out_tbl['n_badstart'].append( int( lsamp_df.sum_bad_starts.sum() ) )
            out_tbl['n_badend'].append( int( lsamp_df.sum_bad_ends.sum() ) )
            out_tbl['n_softclip'].append( int( lsamp_df.sum_soft_clipped.sum() ) )
            out_tbl['n_otheriso'].append( int( lsamp_df.sum_other_isoform.sum() ) )

            for col in med_col_names:
                out_tbl['med_'+col].append( float( lsamp_df[ col ].median() ) )

        i+=1

    out_tbl = pd.DataFrame( out_tbl )

    out_tbl['per_var_seen'] = 100*( out_tbl.n_var / out_tbl.psbl_var )
    out_tbl['per_var_ex'] = 100*( out_tbl.n_var_ex / out_tbl.n_var )
    out_tbl['per_var_int'] = 100*( out_tbl.n_var_int / out_tbl.n_var )
    out_tbl['per_reads_passfilt'] = 100*( out_tbl.n_reads_passfilt / out_tbl.n_reads )
    out_tbl['per_bc_passfilt'] = 100*( out_tbl.n_bc_passfilt / out_tbl.n_bc )
    out_tbl['per_usable'] = 100*( out_tbl.n_usable_reads / out_tbl.n_reads_passfilt )
    out_tbl['per_unmapped'] = 100*( out_tbl.n_unmapped / out_tbl.n_reads_passfilt )
    out_tbl['per_badstart'] = 100*( out_tbl.n_badstart / out_tbl.n_reads_passfilt )
    out_tbl['per_badend'] = 100*( out_tbl.n_badend / out_tbl.n_reads_passfilt )
    out_tbl['per_softclip'] = 100*( out_tbl.n_softclip / out_tbl.n_reads_passfilt )
    out_tbl['per_otheriso'] = 100*( out_tbl.n_otheriso / out_tbl.n_reads_passfilt )

    return out_tbl
