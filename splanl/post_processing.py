import pandas as pd
import numpy as np
import scipy.stats as ss
import pysam

def get_refseq( fa_file ):

    refseq = []

    with pysam.FastxFile(fa_file) as fa:
        for entry in fa:
            refseq.append( entry.sequence )

    return refseq

def get_kmers(seq,
               k):

    return( [ seq[ i : i+k ] for i in range( len( seq ) - k + 1 ) ] )

def acceptors_donors(refseq,
                    byvartbl,
                    use_col,
                    use_thresh
                             ):

    tbv = byvartbl.copy()

    assert use_col in tbv, '%s is not in the dataframe columns'

    min_pos = tbv.pos.min()
    max_pos = tbv.pos.max()

    #accounts for 1 based numbering
    vec_seq = refseq.upper()[ min_pos-1: max_pos ]

    di_nts = get_kmers(vec_seq, 2 )

    for col in ['wt_acc','wt_don','psbl_snv_acc','psbl_snv_don']:
        tbv[ col ] = False

    for p, dnt in enumerate( di_nts ):

        #first check for cryptic exon locations
        if dnt == 'AG':
            tbv.loc[ ( tbv.pos == ( p + min_pos ) ) | ( tbv.pos == ( p + min_pos + 1 ) ), 'wt_acc' ] = True
        elif dnt == 'GT':
            tbv.loc[ ( tbv.pos == ( p + min_pos ) ) | ( tbv.pos == ( p + min_pos + 1 ) ), 'wt_don' ] = True
        #ok now look for possible acceptors and donors
        elif dnt.startswith( 'A' ):
            tbv.loc[ ( tbv.pos == ( p + min_pos + 1 ) ) & ( tbv.alt == 'G' ), 'psbl_snv_acc' ] = True
        elif dnt.endswith( 'G' ):
            tbv.loc[ ( tbv.pos == ( p + min_pos ) ) & ( tbv.alt == 'A' ), 'psbl_snv_acc' ] = True
        elif dnt.startswith( 'G' ):
            tbv.loc[ ( tbv.pos == ( p + min_pos + 1 ) ) & ( tbv.alt == 'T' ), 'psbl_snv_don' ] = True
        elif dnt.endswith( 'T' ):
            tbv.loc[ ( tbv.pos == ( p + min_pos ) ) & ( tbv.alt == 'G' ), 'psbl_snv_don' ] = True

    tbv[ 'snv_acc' ] = ( tbv[ use_col ] >= use_thresh ) & ( tbv.psbl_snv_acc == True )
    tbv[ 'snv_don' ] = ( tbv[ use_col ] >= use_thresh ) & ( tbv.psbl_snv_don == True )

    return tbv

def sdvs(byvartbl,
        sdv_col,
        sdv_thresh,
        abs_vals = True):

    tbv = byvartbl.copy()

    if abs_vals:
        tbv['sdv'] = np.abs( tbv[ sdv_col ] ) >= sdv_thresh
    else:
        tbv['sdv'] = tbv[ sdv_col ] >= sdv_thresh

    return tbv

def print_summary_info( byvartbl ):

    tbv = byvartbl.loc[ byvartbl.n_bc_passfilt > 0  ].copy()

    pos_var = 3*( tbv.pos.max() - tbv.pos.min() )
    seen_var = tbv.shape[0]
    seen_per = 100*( seen_var / pos_var )

    print( 'Out of %i possible variants, we see %i (%.2f%%).' % ( pos_var, seen_var, seen_per ) )

    sdvs = tbv.sdv.sum()
    sdv_per = 100*( sdvs / seen_var )

    print( 'Out of %i variants, %i (%.2f%%) are splice disrupting.' % ( seen_var, sdvs, sdv_per ) )

    if 'var_type' in tbv.columns:

        syn_tbv = tbv.query( 'var_type == "Synonymous"' ).copy()
        syn = syn_tbv.shape[0]
        syn_sdvs = syn_tbv.sdv.sum()
        syn_sdv_per = 100*( syn_sdvs / syn )

        print( 'Out of %i synononymous variants, %i (%.2f%%) are splice disrupting.' % ( syn, syn_sdvs, syn_sdv_per ) )

    if 'sdv_exon' in tbv.columns:

        ex_tbv = tbv.query( 'var_type != "Intronic"' ).copy()
        ex = ex_tbv.shape[0]
        ex_sdvs = ex_tbv.sdv_exon.sum()
        ex_sdv_per = 100*( ex_sdvs / ex )

        print( 'Out of %i exonic variants, %i (%.2f%%) are splice disrupting.' % ( ex, ex_sdvs, ex_sdv_per ) )

    if 'sdv_intron' in tbv.columns:

        intron_tbv = tbv.query( 'var_type == "Intronic"' ).copy()
        intron = intron_tbv.shape[0]
        intron_sdvs = intron_tbv.sdv_intron.sum()
        intron_sdv_per = 100*( intron_sdvs / intron )

        print( 'Out of %i intronic variants, %i (%.2f%%) are splice disrupting.' % ( intron, intron_sdvs, intron_sdv_per ) )

    pos_acc = tbv.psbl_snv_acc.sum()
    used_acc = tbv.snv_acc.sum()
    acc_per = 100*( used_acc / pos_acc )

    print( 'Out of %i possible alternate acceptors, %i (%.2f%%) have a high OTHER value.' % ( pos_acc, used_acc, acc_per ) )

    pos_don = tbv.psbl_snv_don.sum()
    used_don = tbv.snv_don.sum()
    don_per = 100*( used_don / pos_don )

    print( 'Out of %i possible alternate donors, %i (%.2f%%) have a high OTHER value.\n' % ( pos_don, used_don, don_per ) )

def stdize_cols_by_sample(tbv,
                            std_cols):

    out_tbl = tbv.copy()

    #if this is a long dataset
    if 'sample' in out_tbl.columns:

        for col in std_cols:
            #creates z score by sample while ignoring any missing values
            out_tbl[ 'z' + col ] = out_tbl.groupby( [ 'sample' ] )[ col ].transform( lambda x : ss.zscore( x, nan_policy='omit' ) )

    #if the data contains only one sample
    else:

        for col in std_cols:
            print(col)
            out_tbl[ 'z' + col ] = ss.zscore( out_tbl[ col ], nan_policy='omit' )

    return( out_tbl )

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
                'n_otheriso':[],
                'n_sdv':[],
                'n_sdv_ex':[],
                'n_sdv_in':[],
                'psbl_alt_acc':[],
                'psbl_alt_don':[],
                'n_alt_acc':[],
                'n_alt_don':[]
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
            out_tbl['n_reads'].append( int( lsamp_df.sum_reads.sum() ) )
            out_tbl['n_reads_passfilt'].append( int( lsamp_df.sum_reads_passfilt.sum() ) )
            out_tbl['n_usable_reads'].append( int( lsamp_df.sum_usable_reads.sum() ) )
            out_tbl['n_bc'].append( int( lsamp_df.n_bc.sum() ) )
            out_tbl['n_bc_passfilt'].append( int( lsamp_df.n_bc_passfilt.sum() ) )
            out_tbl['n_unmapped'].append( int( lsamp_df.sum_unmapped_reads.sum() ) )
            out_tbl['n_badstart'].append( int( lsamp_df.sum_badstart_reads.sum() ) )
            out_tbl['n_otheriso'].append( int( lsamp_df.sum_otheriso.sum() ) )
            out_tbl['n_sdv'].append( int( lsamp_filt_df.sdv.sum() ) )
            out_tbl['psbl_alt_acc'].append( int( lsamp_filt_df.psbl_snv_acc.sum() ) )
            out_tbl['psbl_alt_don'].append( int( lsamp_filt_df.psbl_snv_don.sum() ) )
            out_tbl['n_alt_acc'].append( int( lsamp_filt_df.snv_acc.sum() ) )
            out_tbl['n_alt_don'].append( int( lsamp_filt_df.snv_don.sum() ) )
            for col in med_col_names:
                out_tbl['med_'+col].append( float( lsamp_df[ col ].median() ) )

            if 'var_type' in lsamp_filt_df.columns:
                out_tbl['n_var_ex'].append( int( lsamp_filt_df.query( 'var_type != "Intronic"' ).shape[0] ) )
                out_tbl['n_var_in'].append( int( lsamp_filt_df.query( 'var_type == "Intronic"' ).shape[0] ) )
                out_tbl['n_sdv_ex'].append( int( lsamp_filt_df.sdv_exon.sum() ) )
                out_tbl['n_sdv_in'].append( int( lsamp_filt_df.sdv_intron.sum() ) )
            else:
                out_tbl['n_var_ex'].append( None )
                out_tbl['n_var_in'].append( None )
                out_tbl['n_sdv_ex'].append( None )
                out_tbl['n_sdv_in'].append( None )
        i+=1

    out_tbl = pd.DataFrame( out_tbl )

    out_tbl['per_var_seen'] = 100*( out_tbl.n_var / out_tbl.psbl_var )
    out_tbl['per_reads_passfilt'] = 100*( out_tbl.n_reads_passfilt / out_tbl.n_reads )
    out_tbl['per_bc_passfilt'] = 100*( out_tbl.n_bc_passfilt / out_tbl.n_bc )
    out_tbl['per_usable'] = 100*( out_tbl.n_usable_reads / out_tbl.n_reads_passfilt )
    out_tbl['per_unmapped'] = 100*( out_tbl.n_unmapped / out_tbl.n_reads_passfilt )
    out_tbl['per_badstart'] = 100*( out_tbl.n_badstart / out_tbl.n_reads_passfilt )
    out_tbl['per_otheriso'] = 100*( out_tbl.n_otheriso / out_tbl.n_reads_passfilt )
    out_tbl['per_sdv'] = 100*( out_tbl.n_sdv / out_tbl.n_var )
    out_tbl['per_sdv_ex'] = 100*( out_tbl.n_sdv_ex / out_tbl.n_var_ex )
    out_tbl['per_sdv_in'] = 100*( out_tbl.n_sdv_in / out_tbl.n_var_in )
    out_tbl['per_acc_used'] = 100*( out_tbl.n_alt_acc / out_tbl.psbl_alt_acc )
    out_tbl['per_don_used'] = 100*( out_tbl.n_alt_don / out_tbl.psbl_alt_don )

    return out_tbl

#This is the one letter universal translation table
#It handles cases of DNA ambiguity where the encoded amino acid is unambiguous.
#You need to deal with the missing cases where ambiguity codes would result in
#an ambiguous amino acid assignment. It is suggested that you use 'X' in these
#cases as this is the standard character for an unknown amino acid.
#Only Y (pyrimidine), R (purine) and N (any) degeneracy symbols are handled at
#this time. (need to add M,K,W,S,B,D,H,V where appropriate)
#Stop codons are symbolized as X
#Reassign TAA, TAG, TAR and TGA to change the stop codon symbol if desired.

transTab1L = {
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

def get_ref_amino( vartbl,
                    refseq,
                    exon_coords,
                    frame_shift = 0 ):

    assert isinstance( frame_shift, int) and frame_shift < 3, 'Frameshift must be a non-negative integer less than 3'

    tbv = vartbl.set_index( 'pos' ).copy()

    exon_seq = refseq[ exon_coords[0] - 1 + frame_shift : exon_coords[1] ].upper()

    daminos = { ( i + exon_coords[0] + frame_shift ) : transTab1L[ exon_seq[ i: i+3 ] ]
                for i in range( 0, len( exon_seq ), 3 )
                if len( exon_seq[ i: i+3 ] ) == 3  }

    #fills in the other positions
    for pos in list( daminos.keys() ):
        amino = daminos[pos]
        daminos[ pos + 1 ] = amino
        daminos[ pos + 2 ] = amino

    aa_df = pd.DataFrame.from_dict( daminos, orient='index' ).reset_index()
    aa_df = aa_df.rename( columns={'index':'pos', 0:'ref_aa'}).sort_values(by='pos').set_index('pos')

    out_tbl = pd.merge( tbv, aa_df, left_index=True, right_index=True, how = 'outer' ).reset_index()

    #drop rows that aren't associated with a variant, can't do an inner merge bc we'll lose any intronic variants
    out_tbl = out_tbl.dropna( subset=['varlist'] )

    return out_tbl

def get_snv_alt_amino( vartbl,
                        refseq,
                        exon_coords,
                        frame_shift = 0 ):

    assert isinstance( frame_shift, int) and frame_shift < 3, 'Frameshift must be a non-negative integer less than 3'

    tbv = vartbl.set_index( [ 'pos', 'alt' ] ).copy()

    exon_seq = refseq[ exon_coords[0] - 1 + frame_shift : exon_coords[1] ].upper()

    lcodons = [ exon_seq[ i: i+3 ] for i in range( 0, len( exon_seq ), 3 ) ]

    nt_sub = ['A', 'C', 'G', 'T']

    daminos = {}

    #handles the case in which the initial bases go to a codon in the upstream exon
    for i in range( frame_shift ):

        #adjust for 1 and 0 based coordinates
        ref = refseq[ exon_coords[0] - 1 + i ].upper()

        if ( i + exon_coords[0] ) not in daminos:
            daminos[ i + exon_coords[0] ] = {}

        for snv in nt_sub:

            if snv == ref:
                continue

            daminos[ i + exon_coords[0] ][ snv ] = 'Exonic - out of frame'

    for i in range( exon_coords[1] - exon_coords[0] - frame_shift ):

        #adjust for 1 and 0 based coordinates
        ref = refseq[ exon_coords[0] - 1 + frame_shift + i ].upper()
        codon = lcodons[ i // 3 ]
        snv_pos = i % 3

        if ( i + exon_coords[0] + frame_shift ) not in daminos:
            daminos[ i + exon_coords[0] + frame_shift ] = {}

        for snv in nt_sub:

            if snv == ref:
                continue

            if len( codon ) == 3:
                daminos[ i + exon_coords[0] + frame_shift ][ snv ] = transTab1L[ codon[ :snv_pos ] + snv + codon[ snv_pos + 1: ] ]
            #handles the case in which the final codon goes into the next downstream exon
            else:
                daminos[ i + exon_coords[0] + frame_shift ][ snv ] = 'Exonic - out of frame'

    #creates dataframe indexed on position and alt with alternate aminos as the only column
    aa_df = pd.concat( { k: pd.DataFrame.from_dict( v, 'index', columns = [ 'alt_aa' ] ) for k, v in daminos.items() }, axis=0 ).reset_index()
    aa_df = aa_df.rename( columns={ 'level_0':'pos', 'level_1':'alt' } ).set_index( [ 'pos', 'alt' ] )

    out_tbl = pd.merge( tbv, aa_df, left_index=True, right_index=True, how = 'outer' ).reset_index()

    #drop rows that aren't associated with a variant, can't do an inner merge bc we'll lose any intronic variants
    out_tbl = out_tbl.dropna( subset=['varlist'] )

    return out_tbl

def extract_var_type( vartbl,
                      refseq,
                        exon_coords,
                        frame_shift = 0 ):

    assert isinstance( frame_shift, int) and frame_shift < 3, 'Frameshift must be a non-negative integer less than 3'

    tbv = vartbl.sort_values( by = 'pos' ).copy()

    if 'ref_aa' not in tbv.columns:
        tbv = get_ref_amino( tbv, refseq, exon_coords, frame_shift )

    if 'alt_aa' not in tbv.columns:
        tbv = get_snv_alt_amino( tbv, refseq, exon_coords, frame_shift )

    var_type = []

    for ref, alt in zip( tbv.ref_aa.values, tbv.alt_aa.values ):

        #checking if alt is missing but with a work around since np.isnan fails with strings
        if not isinstance( alt, str ):
            var_type.append( 'Intronic' )
        elif alt == 'Exonic - out of frame':
            var_type.append('Exonic - out of frame')
        elif ref == alt:
            var_type.append( 'Synonymous' )
        elif alt == 'X':
            var_type.append( 'Nonsense' )
        else:
            var_type.append( 'Missense' )

    tbv['var_type'] = var_type

    return tbv

def sdv_by_var_type( vartbl ):

    tbv = vartbl.copy()

    #check that variant type is in the data and not all missing
    assert ( 'var_type' in tbv.columns ) and ( tbv.shape[0] != tbv.var_type.isnull().sum() ), \
    'var_type must be a non-empty column in the dataframe'

    assert ( 'sdv' in tbv.columns ), 'sdv must be a column in the dataframe'

    tbv['sdv_exon'] = ( ( tbv.sdv ) & ( tbv.var_type != 'Intronic' ) )

    tbv['sdv_intron'] = ( ( tbv.sdv ) & ( tbv.var_type == 'Intronic' ) )

    return tbv
