import numpy as np
import pandas as pd
import hgvs.parser
import hgvs.dataproviders.uta
import hgvs.assemblymapper
import os

def bootstrap_varsp_null_distribution( null_bc_df,
                                       byvar_df,
                                       seed = 1687,
                                       iso_names = None,
                                       bootstraps = 1000, ):

    null_psi = null_bc_df.copy()
    tbv = byvar_df.copy()

    n_bcs = tbv.rna_nbc_varsingleton.tolist()

    assert all( n > 0 for n in n_bcs ), 'Inputting a sample with 0 barcodes will lead to division by 0 errors'

    null_bcs = null_psi.index.tolist()

    if not iso_names:

        iso_names = [ col[4:] for col in null_psi if col.startswith( 'psi_' ) ]

        assert len( iso_names ) > 0, 'Could not infer isoform names - please input names to use'

    usable_reads = null_psi.totalrd_ok.to_numpy()

    sample_tbl = {}
    null_iso = {}

    for iso in iso_names:

        sample_tbl[ 'wmean_bs_null_' + iso ] = []
        sample_tbl[ 'wstdev_bs_null_' + iso ] = []

        null_iso[ iso ] = null_psi[ 'psi_' + iso ].to_numpy()

    bcs_sampled = {}
    
    np.random.seed( seed )
    
    for i,n_bc in enumerate( n_bcs ):

        if n_bc in bcs_sampled:

            idx = bcs_sampled[ n_bc ]

            for iso in iso_names:

                sample_tbl[ 'wmean_bs_null_' + iso ].append( sample_tbl[ 'wmean_bs_null_' + iso ][ idx ] )
                sample_tbl[ 'wstdev_bs_null_' + iso ].append( sample_tbl[ 'wstdev_bs_null_' + iso ][ idx ] )

            continue

        bcs_sampled[ n_bc ] = i

        null_idx = np.random.randint( len( null_psi ), size = ( bootstraps, int( n_bc ) ) )

        for iso in iso_names:

            mus = ( usable_reads[ null_idx ] * null_iso[ iso ][ null_idx ] ).sum( axis = 1 ) / usable_reads[ null_idx ].sum( axis = 1 )
            
            sample_tbl[ 'wmean_bs_null_' + iso ].append( np.mean( mus ) )
            sample_tbl[ 'wstdev_bs_null_' + iso ].append( np.std( mus ) )

    #samp_df = pd.DataFrame( sample_tbl )

    for iso in iso_names:

        tbv[ 'wmean_bs_null_' + iso ] = sample_tbl[ 'wmean_bs_null_' + iso ]
        tbv[ 'wstdev_bs_null_' + iso ] = sample_tbl[ 'wstdev_bs_null_' + iso ]

    #tbv = pd.concat( [ tbv.reset_index(), samp_df.reset_index() ],
                      #axis = 1, )

    print( 'done' )

    return tbv

def compute_null_zscores( tbl_byvar,
                        null_stem,
                        iso_names ):

    tbv = tbl_byvar.copy()

    for iso in iso_names:
    # edited name here in first line
        tbv[ '_'.join( [ 'zwmean', null_stem, iso ] ) ] = ( tbv[ 'psi_' + iso + '_singleton_wmean'] \
                                                            - tbv[ '_'.join( [ 'wmean', null_stem, iso ] ) ] ) \
                                                            / tbv[ '_'.join( [ 'wstdev', null_stem, iso ] ) ]

    return tbv

def compute_fold_change( tbl_byvar,
                         null_col_stem,
                         test_col_stem,
                         iso_names = None,
                         out_col = 'fc_' ):

    tbv = tbl_byvar.copy()

    if not iso_names:

        iso_names = [ col[ len( null_col_stem ): ] for col in tbv if col.startswith( null_col_stem ) ]

        assert len( iso_names ) > 0, 'Cannot infer isoform names - please provide them directly'

    for iso in iso_names:

        tbv[ out_col + iso ] = tbv[ test_col_stem + iso + '_singleton_wmean' ] / tbv[ null_col_stem + iso ]

    return tbv

def stouffers_z( tbl_byvar_wide,
                 iso_names,
                 zcol = 'zmean_',
                 weight = False ):

    tbv = tbl_byvar_wide.copy()

    for iso in iso_names:

        if not weight:

            tbv[ zcol + iso ] = tbv[ [ col for col in tbv if col.endswith( zcol + iso ) ] ].sum( axis = 1 ) \
                                / np.sqrt( tbv[ [ col for col in tbv if col.endswith( zcol + iso ) ] ].notnull().sum( axis = 1 ) )

        else:

            tbv[ zcol[ 1: ] + iso ] = ( tbv[ [ col for col in tbv if col.endswith( zcol[ :-2 ] + '_' + iso ) ] ] * weight[ iso ] ).sum( axis = 1 ) \
                                / np.sqrt( ( weight[ iso ]**2 ).sum() )

    return tbv

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
        wide[ 'mean_' + col ] = wide[ [ '_'.join( [ samp, 'psi', col, 'singleton', 'wmean' ] ) for samp in lsampnames ] ].mean( axis = 1 )

        #this probably would look better as a numpy array dot product but we survived
        wide[ 'bcw_wmean_' + col ] = pd.DataFrame( ( wide[ '_'.join( [ samp, 'psi', col, 'singleton', 'wmean' ] ) ] * wide[ '_'.join( [ samp, bccount ] ) ]
                                                   for samp in lsampnames ) ).T.sum( axis = 1) \
                                             / wide[ samp_bc_cols ].sum( axis = 1 )

    return wide

# t is the transcript accession, nc is chrom number
def annotate_protein_hgvs(widetbl_, t):
    widetbl = widetbl_.copy()

    # this means hgvs will use a locally installed seqrepo, much faster and do not run into the connection issues you get with using the web server
    # also allows for jobs to be run into parallel without the server thinking we are a bot
    os.environ["HGVS_SEQREPO_DIR"] = "/nfs/turbo/umms-kitzmanj/oldvol2/conbward/msh2_splicing/refs/seqrepo/2024-12-20"
    hp = hgvs.parser.Parser()
    hdp = hgvs.dataproviders.uta.connect()
    am = hgvs.assemblymapper.AssemblyMapper( hdp,
                                            assembly_name = 'GRCh38', 
                                            alt_aln_method = 'splign',
                                            replace_reference = True )

    widetbl[ 'protein_var' ] = [ str( am.c_to_p( hp.parse_hgvs_variant( t + ':' + str( p ) + ra[ 0 ] + '>' + ra[ 1 ] ) ).format( conf = { 'p_3_letter': False } ) ).split( ':' )[ 1 ]
                                    for p, ra in zip ( widetbl.hgvs_pos, zip( widetbl.ref, widetbl.alt ) ) ]

    return widetbl



def normalize_chr(source, target):
    """
    Normalize chromosome identifiers to ensure consistency in format (with or without 'chr' prefix).
    
    Args:
        source (str): Source chromosome identifier.
        target (str): Target chromosome identifier for comparison.
    
    Returns:
        str: Normalized chromosome identifier.
    """

    def has_prefix(x):
        return str(x).startswith('chr')  # Check if a chromosome name has 'chr' prefix

    if has_prefix(source) and not has_prefix(target):
        return source.replace('chr', '', 1)  # Remove 'chr' prefix if target doesn't have it
    elif not has_prefix(source) and has_prefix(target):
        return 'chr' + source  # Add 'chr' prefix if target has it

    return source  # Return source as is if both or neither have 'chr' prefix