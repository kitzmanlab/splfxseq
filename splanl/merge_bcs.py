import os
import pandas as pd
import numpy as np

from collections import OrderedDict as odict  # default is for keys to come back out in order after I think python 3.7

from .coords import pos_to_hgvspos


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
        isonames = [ cn for cn in rna_isect.columns if cn.endswith('psi') ]
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

    aminoCode = {
     'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
     'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
     'TAT': 'Y', 'TAC': 'Y', 'TAA': 'X', 'TAG': 'X',
     'TGT': 'C', 'TGC': 'C',  'TGA': 'X', 'TGG': 'W',
     'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
     'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
     'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
     'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
     'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
     'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
     'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
     'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
     'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
     'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
     'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
     'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
     }

    #dictionary of position: reference to classify variants
    seqDict = {int(singlevar.split(':')[1]) : singlevar.split(':')[2]
                for singlevar, subtbl in rna_isect.groupby( 'varlist' )
                if len( singlevar.split(':')[2] ) == 1 }

    #dictionary of first position of codon: amino acid to classify variants
    aminoDict={}
    shift={}
    for c in exon_coords:
        for j in range(c[0],c[1]+1,3):
            aminoDict[j]=aminoCode[''.join(seqDict[k] for k in range(j,j+3))]
            shift.update({j:0,j+1:1,j+2:2})

    for singlevar, subtbl in rna_isect.groupby( 'varlist' ):

        subtbl_filt = subtbl.loc[ subtbl.usable_reads >= min_usable_reads_per_bc ].copy()

        out_tbl['varlist'].append(singlevar)
        out_tbl['chrom'].append(singlevar.split(':')[0])
        out_tbl['pos'].append(int(singlevar.split(':')[1]))
        out_tbl['ref'].append(singlevar.split(':')[2])
        out_tbl['alt'].append(singlevar.split(':')[3])

        #classifies variants as synonymous, nonsense, misssense, or intronic
        pos=int(singlevar.split(':')[1])
        #lets only deal with SNVs right now - later can be more precise
        if ( len(singlevar.split(':')[2]) > 1 or len(singlevar.split(':')[3]) > 1 ):
            out_tbl['var_type'].append('not SNV')
        elif any([pos in range(coord[0],coord[1]+1) for coord in exon_coords]):
            newAmino=aminoCode[''.join(seqDict[k] if k!= pos else singlevar.split(':')[3]
                                for k in range(pos-shift[pos],pos+3-shift[pos]))]
            if newAmino == aminoDict[pos-shift[pos]]:
                out_tbl['var_type'].append('synonymous')
            elif newAmino == 'X':
                out_tbl['var_type'].append('nonsense')
            else:
                out_tbl['var_type'].append('missense')
        else:
            out_tbl['var_type'].append('intronic SNV')

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
                out_tbl[ f'wmean_{iso}' ].append( ( subtbl_filt[ f'{iso}_psi' ] * subtbl_filt['usable_reads'] ).sum() / subtbl_filt['usable_reads'].sum() )
                # median psi
                out_tbl[ f'median_{iso}' ].append( subtbl_filt[ f'{iso}_psi' ].median() )

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
