import pandas as pd
import numpy as np
from keras.models import load_model
from pkg_resources import resource_filename
from spliceai.utils import one_hot_encode

def rev_complement( seq ):
    """
    Creates reverse complement of DNA string (not case sensitive)

    Args: seq (str)

    Returns: comp (str) reverse complement of input string
    """

    trans_tbl = str.maketrans( 'ACGTNacgtn', 'TGCANtgcan' )

    rev = seq[::-1]

    comp = rev.translate( trans_tbl )

    return comp

def get_gene_bds( annots_df,
                  chrom,
                  position,
                  scored_context,
                  unscored_context = 5000,
                ):
    """
    Gets the number of bases which are within the sequence context but outside the gene boundaries on each side.

    Args: annots_df (pandas df) - columns: #NAME, CHROM, STRAND, TX_START, TX_END, EXON_START, EXON_END
                                    EXON_START, EXON_END are both comma separated strings of all exon bds
          chrom - (str) chromosome of center variant ( format should match your annots file (ie chr3 v 3) )
          position - (int) position (genomic coords) of center variant
          scored_context - (int) number of bases to score on each side of the variant
          unscored_context - (int) number of flanking unscored bases on each side of the variant

    Returns: gene_bds - (tuple of ints) ( upstream bases outside of gene bds, downstream bases outside of gene bds )
                        if sequence is contained entirely within the gene, returns ( 0, 0 )
                        These are the number of bases that will be replaced with N's on either side of the sequence
    """

    ann_df = annots_df.copy()

    idx = ann_df.index[ ( ann_df.CHROM == chrom )
                      & ( ann_df.TX_START <= position )
                      & ( ann_df.TX_END >= position ) ]

    assert len( idx ) == 1, \
    'The chromosome and position is not matching exactly one gene!'

    tx_startsd = position - ( ann_df.at[ idx[ 0 ], 'TX_START' ] + 1 )

    tx_endsd = ann_df.at[ idx[ 0 ], 'TX_END' ] - position

    flanks = unscored_context + scored_context

    gene_bds = ( max( flanks - tx_startsd, 0 ), max( flanks - tx_endsd, 0 ) )

    return gene_bds

def get_2exon_bds( annots_df,
                  chrom,
                  position,
                  rev_strand = False,
                ):
    """
    Gets the distance from the center variant to the nearest acceptor and donor.

    Args: annots_df (pandas df) - columns: #NAME, CHROM, STRAND, TX_START, TX_END, EXON_START, EXON_END
                                    EXON_START, EXON_END are both comma separated strings of all exon bds
          chrom - (str) chromosome of center variant ( format should match your annots file (ie chr3 v 3) )
          position - (int) position (genomic coords) of center variant
          rev_strand - (bool) is the variant on the reverse strand?

    Returns: 2exon_bds - (tuple of ints) ( distance to nearest acceptor, distance to nearest donor )
                         SpliceAI default scoring would only return distance to nearest donor OR acceptor
                         These values are used for masking..
    """

    ann_df = annots_df.copy()

    idx = ann_df.index[ ( ann_df.CHROM == chrom )
                      & ( ann_df.TX_START <= position )
                      & ( ann_df.TX_END >= position ) ]

    assert len( idx ) == 1, \
    'The chromosome and position is not matching exactly one gene!'

    #add 1 to adjust to 0-based coords
    #compute distance to center variant
    exon_startd = [ ( int( start ) + 1 ) - position
                        for start in ann_df.at[ idx[ 0 ], 'EXON_START' ].split( ',' )
                        if start != '' ]

    start = exon_startd[ np.argmin( np.abs( exon_startd ) ) ]

    #compute distance to center variant
    exon_endd = [ int( end ) - position
                      for end in ann_df.at[ idx[ 0 ], 'EXON_END' ].split( ',' )
                      if end != '' ]

    end = exon_endd[ np.argmin( np.abs( exon_endd ) ) ]

    #if reverse strand, flip acceptors/donors
    exon_bds = ( start, end ) if not rev_strand else ( end, start )

    return exon_bds

def create_input_seq( refseq,
                      center_var,
                      haplotype,
                      ref_var,
                      gene_bds,
                      scored_context,
                      rev_strand = False,
                      unscored_context = 5000,
                    ):
    """
    Creates the reference and variant sequences to input into SpliceAI models

    Args: refseq (str) - fasta file for an entire chromosome
          center_var (tuple) - center_variant to be scored:
                                ( position - (int) position (genomic coords) of center variant,
                                 reference base(s) - (str) reference base(s) relative to forward strand,
                                 alternate base(s) - (str) alternate base(s) relative to forward strand, )
          haplotype ( list of tuples ) - other variants to be added to the variant sequences
                                        ( position - (int) position (genomic coords) of variant,
                                        reference base(s) - (str) reference base(s) relative to forward strand,
                                        alternate base(s) - (str) alternate base(s) relative to forward strand, )
                                        Empty list adds no additional variants
          ref_var ( list of tuples ) - other variants to be added to the reference AND variant sequences
                                       ( position - (int) position (genomic coords) of variant,
                                       reference base(s) - (str) reference base(s) relative to forward strand,
                                       alternate base(s) - (str) alternate base(s) relative to forward strand, )
                                       Empty list adds no additional variants
          gene_bds ( tuple of ints ) - ( number of upstream bases outside of gene, number of downstream bases outside of gene)
                                       Adds N's to sequence outside of gene
                                       ( 0, 0 ) adds no N bases
          rev_strand - (bool) is the center variant on the reverse strand?
          scored_context - (int) number of bases to score on each side of the variant
          unscored_context - (int) number of flanking unscored bases on each side of the variant

    Returns: refvar_seq - (tuple of str) ( reference sequence, variant sequence )
                          Both sequences will be 2*( scored_context + unscored_context ) + 1 long if none of the variants are indels
    """

    flanks = unscored_context + scored_context

    #subtracting one to adjust from 1 based to 0 based coords
    #gene bds is adding N's for locations outside the gene
    refseq = 'N'*gene_bds[ 0 ] \
             + refseq[ center_var[ 0 ] - flanks - 1 + gene_bds[ 0 ] : center_var[ 0 ] + flanks - gene_bds[ 1 ] ] \
             + 'N'*gene_bds[ 1 ]

    #print( len( refseq ) )
    assert len( refseq ) == flanks*2 + 1, 'Debug your math'

    refvar_seqpos = [ ( flanks + ( p - center_var[ 0 ] ), r, a ) for p,r,a in ref_var ]

    #start from the back so the coords stay the same despite indels
    #refvar_seqpos.sort( reverse = True )

    for pos,ref,alt in refvar_seqpos:

        assert len( ref ) == len( alt ), 'Code is not ready for indels within the reference variants yet'

        assert refseq[ pos: pos + len( ref ) ].upper() == ref, 'Reference base within reference variants does not match given position'

        refseq = refseq[ : pos ] + alt + varseq[ pos + len( ref ): ]

    hap_seqpos = [ ( flanks + ( p - center_var[ 0 ] ), r, a ) for p,r,a in haplotype + [ center_var ] ]

    #start from the back so the coords stay the same despite indels
    hap_seqpos.sort( reverse = True )

    varseq = refseq

    for pos,ref,alt in hap_seqpos:

        assert refseq[ pos: pos + len( ref ) ].upper() == ref, 'Reference base within haplotype does not match given position'

        varseq = varseq[ : pos ] + alt + varseq[ pos + len( ref ): ]

    if rev_strand:

        refseq = rev_complement( refseq )
        varseq = rev_complement( varseq )

    return ( refseq, varseq )

def splai_score_variants( annots_df,
                           models,
                           refvarseqs,
                           ref_name,
                           chrom,
                           center_var,
                           haplotypes,
                           scored_context = 50,
                           rev_strand = False,
                           mask_value = 0,
                         ):
    """
    Uses SpliceAI default scoring to compute SDV probabilities.

    Args: annots_df (pandas df) - columns: #NAME, CHROM, STRAND, TX_START, TX_END, EXON_START, EXON_END
                                    EXON_START, EXON_END are both comma separated strings of all exon bds
          models (list) - SpliceAI models
          refvarseqs ( list of tuples of strings ) - list of all of the reference and variant sequence pairs to score
          refname (str) - reference and/or annotation name to label the entry
          center_var (list of tuples) - center_variants to be scored:
                                ( position - (int) position (genomic coords) of center variant,
                                 reference base(s) - (str) reference base(s) relative to forward strand,
                                 alternate base(s) - (str) alternate base(s) relative to forward strand, )
          haplotypes ( list of list of tuples ) - other variants added to the variant sequences
                                        ( position - (int) position (genomic coords) of variant,
                                        reference base(s) - (str) reference base(s) relative to forward strand,
                                        alternate base(s) - (str) alternate base(s) relative to forward strand, )
                                        Empty list indicates no additional variants
          scored_context - (int) number of bases to score on each side of the variant
          rev_strand - (bool) is the center variant on the reverse strand?
          mask_value - (int) value to used to mask scores

    Returns: outdf - (pandas df) Pandas dataframe with probabilities for each type of splice site event
                                 and separate probabilities after masking
    """

    outtbl = { 'ref_name': [ ref_name ]*len( refvarseqs ),
               'chrom': [ chrom ]*len( refvarseqs ),
               'pos': [],
               'ref': [],
               'alt': [],
               'other_var': [],
               'acc_abs_chg': [],
               'don_abs_chg': [],
               'DS_AG': [],
               'DS_AL': [],
               'DS_DG': [],
               'DS_DL': [],
               'DP_AG': [],
               'DP_AL': [],
               'DP_DG': [],
               'DP_DL': [],
               'DS_max': [],
               'DS_max_type': [],
               'POS_max': [],
               'DS_AGm': [],
               'DS_ALm': [],
               'DS_DGm': [],
               'DS_DLm': [],
               'DS_maxm': [],
               'DS_maxm_type': [],
             }

    for idx, refvarseq in enumerate( refvarseqs ):

        refseq, varseq = refvarseq

        x_ref = one_hot_encode( refseq )[ None, : ]
        y_ref = np.mean( [ models[ m ].predict( x_ref ) for m in range( 5 ) ], axis=0 )

        x_var = one_hot_encode( varseq )[ None, : ]
        y_var = np.mean( [ models[ m ].predict( x_var ) for m in range( 5 ) ], axis=0 )

        #flips the results so the positions are relative to the forward strand
        if rev_strand:
                y_ref = y_ref[:, ::-1]
                y_var = y_var[:, ::-1]

        ref_acc = y_ref[0, :, 1]
        ref_don = y_ref[0, :, 2]
        var_acc = y_var[0, :, 1]
        var_don = y_var[0, :, 2]

        #transforms variants into sequence position coords
        allvars_seqpos = [ ( scored_context + ( p - center_var[ idx ][ 0 ] ), r, a )
                           for p,r,a in [ center_var[ idx ] ] + haplotypes[ idx ] ]

        var_acc, var_don = adjust_for_indels( var_acc,
                                              var_don,
                                              allvars_seqpos,
                                             )

        diff_acc = ref_acc - var_acc
        diff_don = ref_don - var_don

        outtbl[ 'pos' ].append( center_var[ idx ][ 0 ] )
        outtbl[ 'ref' ].append( center_var[ idx ][ 1 ] )
        outtbl[ 'alt' ].append( center_var[ idx ][ 2 ] )
        outtbl[ 'other_var' ].append( ';'.join( [ ':'.join( [ str( p ), '>'.join( [ r, a ] ) ] )
                                      for p,r,a in haplotypes[ idx ] ] ) )

        outtbl[ 'acc_abs_chg' ].append( sum( np.abs( diff_acc ) ) )
        outtbl[ 'don_abs_chg' ].append( sum( np.abs( diff_don ) ) )

        outtbl[ 'DS_AG' ].append( np.abs( np.min( [ 0, np.min( diff_acc ) ] ) ) )
        outtbl[ 'DS_AL' ].append( np.max( [ 0, np.max( diff_acc ) ] ) )
        outtbl[ 'DS_DG' ].append( np.abs( np.min( [ 0, np.min( diff_don ) ] ) ) )
        outtbl[ 'DS_DL' ].append( np.max( [ 0, np.max( diff_don ) ] ) )

        outtbl[ 'DP_AG' ].append( np.argmin( diff_acc ) - scored_context )
        outtbl[ 'DP_AL' ].append( np.argmax( diff_acc ) - scored_context )
        outtbl[ 'DP_DG' ].append( np.argmin( diff_don ) - scored_context )
        outtbl[ 'DP_DL' ].append( np.argmax( diff_don ) - scored_context )

        score_keys = [ 'DS_AG', 'DS_AL', 'DS_DG', 'DS_DL' ]

        #first get the maximum probability across the difference scores
        outtbl[ 'DS_max' ].append( max( outtbl[ key ][ -1 ] for key in score_keys ) )
        #then get the type of event that represents the maximum probability
        outtbl[ 'DS_max_type' ].append( [ key for key in score_keys
                                          if outtbl[ key ][ -1 ] == outtbl[ 'DS_max' ][ -1 ] ][ 0 ] )
        #finally, get the location of the event associated with the highest difference score
        outtbl[ 'POS_max' ].append( outtbl[ outtbl[ 'DS_max_type' ][ -1 ].replace( 'DS', 'DP' ) ][ -1 ] )

        outtbl[ 'DS_AGm' ].append( outtbl[ 'DS_AG' ][ -1 ] )
        outtbl[ 'DS_ALm' ].append( outtbl[ 'DS_AL' ][ -1 ] )
        outtbl[ 'DS_DGm' ].append( outtbl[ 'DS_DG' ][ -1 ] )
        outtbl[ 'DS_DLm' ].append( outtbl[ 'DS_DL' ][ -1 ] )

        exon_bds = get_2exon_bds( annots_df,
                                  chrom,
                                  center_var[ idx ][ 0 ],
                                  rev_strand = rev_strand,
                                  )

        if outtbl[ 'DP_AL' ][ -1 ] != exon_bds[ 0 ]:
            outtbl[ 'DS_ALm' ][ -1 ] = mask_value

        if outtbl[ 'DP_AG' ][ -1 ] == exon_bds[ 0 ]:
            outtbl[ 'DS_AGm' ][ -1 ] = mask_value

        if outtbl[ 'DP_DL' ][ -1 ] != exon_bds[ 1 ]:
            outtbl[ 'DS_DLm' ][ -1 ] = mask_value

        if outtbl[ 'DP_DG' ][ -1 ] == exon_bds[ 1 ]:
            outtbl[ 'DS_DGm' ][ -1 ] = mask_value

        score_keys = [ key + 'm' for key in score_keys ]

        #first get the maximum probability across the difference scores
        outtbl[ 'DS_maxm' ].append( max( outtbl[ key ][ -1 ] for key in score_keys ) )
        #then get the type of event that represents the maximum probability
        outtbl[ 'DS_maxm_type' ].append( [ key for key in score_keys
                                              if outtbl[ key ][ -1 ] == outtbl[ 'DS_maxm' ][ -1 ] ][ 0 ] )

    outdf = pd.DataFrame( outtbl )

    return outdf

def adjust_for_indels( var_accp,
                       var_donp,
                       variants,
                     ):
    """
    Adjusts the length of the variant sequence probabilities when there are indels.
    Specifically, for deletions, fills deleted bases with a probability of 0 and for insertions,
    take the maximum probability across the inserted bases.

    Args:
          var_accp - (np array) acceptor probabilities for variant sequence
          var_donp - (np array) donor probabilities for variant sequence
          variants - ( list of tuples ) all center and haplotype variants on variant sequence
                     ( position - (int) distance from center variant,
                      reference base(s) - (str) reference base(s) relative to forward strand,
                      alternate base(s) - (str) alternate base(s) relative to forward strand, )

    Returns: variantp - (tuple of np arrays) ( acceptor probabilities, donor probabilities )
    """

    #make sure the variants are sorted
    variants.sort()

    for pos,ref,alt in variants:

        #if there's a deletion, fill the missing locations with 0's
        if len( ref ) > len( alt ):

            var_accp = np.concatenate( [ var_accp[ : pos + len( alt ) ],
                                       np.zeros( len( ref ) - len( alt ) ),
                                       var_accp[ pos + len( alt ): ] ] )

            var_donp = np.concatenate( [ var_donp[ : pos + len( alt ) ],
                                       np.zeros( len( ref ) - len( alt ) ),
                                       var_donp[ pos + len( alt ): ] ] )

        #if there's an insertion, fill use the maximum across the insertion as the variant prob
        #need to add 1 here since the final bd in python is exclusive..
        elif len( alt ) > len( ref ):

            var_accp = np.concatenate( [ var_accp[ : pos ],
                                       [ np.max( var_accp[ pos: pos + ( len( alt ) - len( ref ) ) + 1 ] ) ],
                                       var_accp[ pos + ( len( alt ) - len( ref ) ) + 1: ] ] )

            var_donp = np.concatenate( [ var_donp[ : pos ],
                                       [ np.max( var_donp[ pos: pos + ( len( alt ) - len( ref ) ) + 1 ] ) ],
                                       var_donp[ pos + ( len( alt ) - len( ref ) ) + 1: ] ] )

    return ( var_accp, var_donp )

def splai_score_mult_variants_onegene( annots_df,
                                      models,
                                      refseq,
                                      ref_name,
                                      chrom,
                                      center_var,
                                      haplotypes = None,
                                      ref_vars = None,
                                      mask_value = 0,
                                      scored_context = 50,
                                      unscored_context = 5000,
                                      rev_strand = False ):
    """
    Wrapper function to compute SpliceAI default probabilities across one gene.

    Args: annots_df (pandas df) - columns: #NAME, CHROM, STRAND, TX_START, TX_END, EXON_START, EXON_END
                                    EXON_START, EXON_END are both comma separated strings of all exon bds
          models (list) - SpliceAI models
          refseq (str) - fasta file for an entire chromosome
          ref_name (str) - reference and/or annotation name to label the entries
          chrom - (str) chromosome of center variant ( format should match your annots file (ie chr3 v 3) )
          center_var (list of tuples) - center_variants to be scored:
                                ( position - (int) position (genomic coords) of center variant,
                                 reference base(s) - (str) reference base(s) relative to forward strand,
                                 alternate base(s) - (str) alternate base(s) relative to forward strand, )
          haplotypes ( list of list of tuples ) - other variants added to the variant sequences
                                                  each ith list corresponds to the ith center_var
                                        ( position - (int) position (genomic coords) of variant,
                                        reference base(s) - (str) reference base(s) relative to forward strand,
                                        alternate base(s) - (str) alternate base(s) relative to forward strand, )
                                        None indicates no additional variants
          ref_vars ( list of list of tuples ) - other variants to be added to the reference AND variant sequences
                                                each ith list corresponds to the ith center_var
                                             ( position - (int) position (genomic coords) of variant,
                                             reference base(s) - (str) reference base(s) relative to forward strand,
                                             alternate base(s) - (str) alternate base(s) relative to forward strand, )
                                             None adds no additional variants
                                             NOTE: currently only substitutions can be handled - indels will raise error
          mask_value - (int) value to used to mask scores
          scored_context - (int) number of bases to score on each side of the center variant
          unscored_context - (int) number of flanking unscored bases on each side of the center variant
          rev_strand - (bool) is the gene on the reverse strand?

    Returns: outdf - (pandas df) Pandas dataframe with probabilities for each type of splice site event for each center_var
                                 and separate probabilities after masking
    """

    flanks = scored_context + unscored_context

    if not haplotypes:

        haplotypes = [ [] for var in center_var ]

    if not ref_vars:

        ref_vars = [ [] for var in center_var ]

    else:

        for p,r,a in ref_vars:

            ref_name += '+' + str( p ) + ':' + r + '>' + a

    assert all( len(i) == len( center_var ) for i in [ haplotypes, ref_vars ] ), \
    'Haplotypes and ref_vars input must be either missing or the same length as the center_var list'

    #creates a giant list of ( reference seq, variant seq ) tuples
    refvarseqs = [ create_input_seq( refseq,
                                     center,
                                     hapref[ 0 ],
                                     hapref[ 1 ],
                                     get_gene_bds( annots_df,
                                                   chrom,
                                                   center[ 0 ],
                                                   scored_context = scored_context,
                                                   unscored_context = unscored_context ),
                                     scored_context,
                                     rev_strand = rev_strand,
                                     unscored_context = unscored_context )
                  for center,hapref in zip( center_var, zip( haplotypes, ref_vars ) ) ]

    #this will fail if any of the other variants are indels...
    #maybe I should add some functionality to the splai_score_variant fn to handle this...
    outdf = splai_score_variants( annots_df,
                                  models,
                                  refvarseqs,
                                  ref_name,
                                  chrom,
                                  center_var,
                                  haplotypes,
                                  scored_context = scored_context,
                                  rev_strand = rev_strand,
                                  mask_value = mask_value
                          )

    return outdf

def score_mult_variants_multgene( annots_df,
                                  models,
                                  refseqs,
                                  ref_name,
                                  center_var,
                                  haplotypes = {},
                                  ref_vars = {},
                                  mask_value = 0,
                                  scored_context = 50,
                                  unscored_context = 5000,
                                ):
    """
    Wrapper function to compute SpliceAI default probabilities across multiple genes.

    Args: annots_df (pandas df) - columns: #NAME, CHROM, STRAND, TX_START, TX_END, EXON_START, EXON_END
                                    EXON_START, EXON_END are both comma separated strings of all exon bds
          models (list) - SpliceAI models
          refseqs (dict of str) - { chrom: fasta file for chromosome } ( format should match your annots file (ie chr3 v 3) )
          ref_name (str) - reference and/or annotation name to label the entries
          chrom - (str) chromosome of center variant ( format should match your annots file (ie chr3 v 3) )
          center_var (dict of list of tuples) - center_variants to be scored:
                            { chrom: [ ( position - (int) position (genomic coords) of center variant,
                                         reference base(s) - (str) reference base(s) relative to forward strand,
                                         alternate base(s) - (str) alternate base(s) relative to forward strand,
                                         rev_strand - (bool) is the variant on the reverse strand? ) ] }
          haplotypes ( dict of list of list of tuples ) - other variants added to the variant sequences
                                                          each ith list corresponds to the ith center_var
                                        { chrom: [ [ ( position - (int) position (genomic coords) of center variant,
                                                     reference base(s) - (str) reference base(s) relative to forward strand,
                                                     alternate base(s) - (str) alternate base(s) relative to forward strand,
                                                     rev_strand - (bool) is the variant on the reverse strand? ) ] ] }
                                        Empty dictionary indicates no additional variants
          ref_vars ( list of list of tuples ) - other variants to be added to the reference AND variant sequences
                                                each ith list corresponds to the ith center_var
                                             { chrom: [ [ ( position - (int) position (genomic coords) of center variant,
                                                          reference base(s) - (str) reference base(s) relative to forward strand,
                                                          alternate base(s) - (str) alternate base(s) relative to forward strand,
                                                          rev_strand - (bool) is the variant on the reverse strand? ) ] ] }
                                             Empty dictionary indicates no additional variants
                                             NOTE: currently only substitutions can be handled - indels will raise error
          mask_value - (int) value to used to mask scores
          scored_context - (int) number of bases to score on each side of the center variant
          unscored_context - (int) number of flanking unscored bases on each side of the center variant


    Returns: outdf - (pandas df) Pandas dataframe with probabilities for each type of splice site event for each center_var
                                 and separate probabilities after masking
    """

    outdfs = []

    for chrom in center_var.keys():

        chrseq = refseqs[ chrom ]

        assert len( chrseq ) > 0

        for_center = []
        for_indices = []
        rev_center = []
        rev_indices = []

        for i,var in enumerate( center_var[ chrom ] ):

            #have to remove the last entry in the tuple for it to fit in the other fn
            if var[ -1 ] == True:
                rev_center.append( ( var[ 0 ], var[ 1 ], var[ 2 ] ) )
                rev_indices.append( i )

            else:
                for_center.append( ( var[ 0 ], var[ 1 ], var[ 2 ] ) )
                for_indices.append( i )

        if chrom not in haplotypes:
            for_hap = None
            rev_hap = None
        else:
            for_hap = [ haplotypes[ chrom ][ i ] for i in for_indices ]
            rev_hap = [ haplotypes[ chrom ][ i ] for i in rev_indices ]

        if chrom not in ref_vars:
            for_rv = None
            rev_rv = None
        else:
            for_rv = [ ref_vars[ chrom ][ i ] for i in for_indices ]
            rev_rv = [ ref_vars[ chrom ][ i ] for i in rev_indices ]

        if len( for_center ) > 0:

            outdfs.append( score_mult_variants_onegene( annot,
                                                        models,
                                                        chrseq,
                                                        ref_name,
                                                        chrom,
                                                        for_center,
                                                        haplotypes = for_hap,
                                                        ref_vars = for_rv,
                                                        mask_value = mask_value,
                                                        scored_context = scored_context,
                                                        unscored_context = unscored_context,
                                                     )
                         )

        if len( rev_center ) > 0:

            outdfs.append( score_mult_variants_onegene( annot,
                                                        models,
                                                        chrseq,
                                                        refname,
                                                        chrom,
                                                        for_center,
                                                        haplotypes = for_hap,
                                                        ref_vars = for_rv,
                                                        mask_value = mask_value,
                                                        scored_context = scored_context,
                                                        unscored_context = unscored_context,
                                                        rev_strand = True
                                                     )
                         )

    outdf = pd.concat( outdfs, ignore_index = True )

    return outdf

def get_allexon_bds( annots_df,
                     chrom,
                     position,
                     scored_context,
                     rev_strand = False,
                    ):
    """
    Gets the distance from the center variant to all acceptors and donors within range

    Args: annots_df (pandas df) - columns: #NAME, CHROM, STRAND, TX_START, TX_END, EXON_START, EXON_END
                                    EXON_START, EXON_END are both comma separated strings of all exon bds
          chrom - (str) chromosome of center variant ( format should match your annots file (ie chr3 v 3) )
          position - (int) position (genomic coords) of center variant
          scored_context - (int) number of bases on each side to look for donors and acceptors
          rev_strand - (bool) is the variant on the reverse strand?

    Returns: allexon_bds - (list of tuples of ints)
                         [ ( distance to acceptor1,...,distance to acceptorn ),
                           ( distance to donor1,...,distance to donorn ) ]
                         SpliceAI default scoring would only return distance to nearest donor OR acceptor
                         get_2exon_bds returns only two values
                         This fn returns all values within the scored context range
                         These values are used for masking and getting relative jn use
    """

    ann_df = annots_df.copy()

    idx = ann_df.index[ ( ann_df.CHROM == chrom )
                      & ( ann_df.TX_START <= position )
                      & ( ann_df.TX_END >= position ) ]

    assert len( idx ) == 1, \
    'The chromosome and position is not matching exactly one gene - the code is not ready for that!'

    exon_startd = [ ( int( start ) + 1 ) - position
                            for start in ann_df.at[ idx[ 0 ], 'EXON_START' ].split( ',' )
                            if start != '' ]

    starts = tuple( p for p in exon_startd if np.abs( p ) <= scored_context )

    exon_endd = [ int( end ) - position
                          for end in ann_df.at[ idx[ 0 ], 'EXON_END' ].split( ',' )
                          if end != '' ]

    ends = tuple( p for p in exon_endd if np.abs( p ) <= scored_context )

    exon_bds = [ starts, ends ] if not rev_strand else [ ends, starts ]

    return exon_bds

def get_relative_jn_use( gtex_df,
                         chrom,
                         position,
                         exon_bds,
                         scored_context,
                         rev_strand = False, ):

    """
        Gets the distance from the center variant to all acceptors and donors within range

        Args: gtex_df (pandas df) - columns: chrom, jn (1-based), reads
              chrom - (str) chromosome of center variant ( format should match your gtex file (ie chr3 v 3) )
              position - (int) position (genomic coords) of center variant
              scored_context - (int) number of bases on each side to look for donors and acceptors
              rev_strand - (bool) is the variant on the reverse strand?

        Returns: acceptor and donor relative use for each acceptor and donor (list of dictionaries)
                 [ { acceptor1_pos: acceptor1_rel_use, ... acceptorn_pos: acceptorn_rel_use },
                   { donor1_pos: donor1_rel_use, ... donorn_pos: donorn_rel_use } ]
                 The values in both dictionaries should sum to 1
                 The positions are the relative distances to the center variant
    """

    gtex = gtex_df.copy()

    acc, don = exon_bds.copy()

    #gets jns back into gdna coords
    #adjusts back to 0 based coords
    if rev_strand:
        acc_jn = [ position + acc_pos - scored_context for acc_pos in acc ]
        don_jn = [ position + don_pos - scored_context - 1 for don_pos in don ]

    else:
        acc_jn = [ position + acc_pos - scored_context - 1 for acc_pos in acc ]
        don_jn = [ position + don_pos - scored_context for don_pos in don ]

    if len( acc_jn ) > 1:
        #create dictionaries to hold sequence position jn: relative expression (proportion of total reads) pairs
        acc_exp = { acc_spos: gtex.loc[ ( gtex.chrom == chrom ) & ( gtex.jn == acc_gpos ) ].n_rds.values[ 0 ]
                    for acc_spos, acc_gpos in zip( acc, acc_jn ) }
        tot_reads = sum( acc_exp.values() )
        acc_exp = { jn: rds / tot_reads for jn, rds in acc_exp.items() }

    else:
        acc_exp = { acc[ 0 ]: 1 }

    if len( don_jn ) > 1:
        don_exp = { don_spos: gtex.loc[ ( gtex.chrom == chrom ) & ( gtex.jn == don_gpos ) ].n_rds.values[ 0 ]
                    for don_spos, don_gpos in zip( don, don_jn ) }
        tot_reads = sum( don_exp.values() )
        don_exp = { jn: rds / tot_reads for jn, rds in don_exp.items() }

    else:
        don_exp = { don[ 0 ]: 1 }

    return [ acc_exp, don_exp ]

def jnuse_score_variants(  models,
                           refvarseqs,
                           ref_name,
                           ss_jn_use,
                           chrom,
                           center_var,
                           haplotypes,
                           scored_context,
                           rev_strand = False,
                 ):
    """
    Uses relative junction use to compute SDV probabilities.

    Args:
          models (list) - SpliceAI models
          refvarseqs ( list of tuples of strings ) - list of all of the reference and variant sequence pairs to score
          refname (str) - reference and/or annotation name to label the entry
          ss_jn_use (list of dicts) - acceptor and donor relative use for each acceptor and donor (list of dictionaries)
                   [ { acceptor1_pos: acceptor1_rel_use, ... acceptorn_pos: acceptorn_rel_use },
                     { donor1_pos: donor1_rel_use, ... donorn_pos: donorn_rel_use } ]
                   The values in both dictionaries should sum to 1
                   The positions are the relative distances to the center variant
          center_var (list of tuples) - center_variants to be scored:
                                ( position - (int) position (genomic coords) of center variant,
                                 reference base(s) - (str) reference base(s) relative to forward strand,
                                 alternate base(s) - (str) alternate base(s) relative to forward strand, )
          haplotypes ( list of list of tuples ) - other variants added to the variant sequences
                                        ( position - (int) position (genomic coords) of variant,
                                        reference base(s) - (str) reference base(s) relative to forward strand,
                                        alternate base(s) - (str) alternate base(s) relative to forward strand, )
                                        Empty list indicates no additional variants
          scored_context - (int) number of bases to score on each side of the variant
          rev_strand - (bool) is the center variant on the reverse strand?

    Returns: outdf - (pandas df) Pandas dataframe with probabilities for each type of splice site event
                                 and separate probabilities after masking
    """

    outtbl = { 'ref_name': [ ref_name ]*len( refvarseqs ),
               'chrom': [ chrom ]*len( refvarseqs ),
               'pos': [],
               'ref': [],
               'alt': [],
               'other_var': [],
               'AL_chg': [],
               'AG_chg': [],
               'DL_chg': [],
               'DG_chg': [],
               'DS_AGrw': [],
               'DS_ALrw': [],
               'DS_DGrw': [],
               'DS_DLrw': [],
               'DP_AGrw': [],
               'DP_ALrw': [],
               'DP_DGrw': [],
               'DP_DLrw': [],
               'DS_maxrw': [],
               'DS_maxrw_type': [],
               'POS_maxrw': [],
             }

    for idx, refvarseq in enumerate( refvarseqs ):

        refseq, varseq = refvarseq

        x_ref = one_hot_encode( refseq )[ None, : ]
        y_ref = np.mean( [ models[ m ].predict( x_ref ) for m in range( 5 ) ], axis=0 )

        x_var = one_hot_encode( varseq )[ None, : ]
        y_var = np.mean( [ models[ m ].predict( x_var ) for m in range( 5 ) ], axis=0 )

        #flips the results so the positions are relative to the forward strand
        if rev_strand:
                y_ref = y_ref[:, ::-1]
                y_var = y_var[:, ::-1]

        ref_acc = y_ref[0, :, 1]
        ref_don = y_ref[0, :, 2]
        var_acc = y_var[0, :, 1]
        var_don = y_var[0, :, 2]

        #transforms variants into sequence position coords
        allvars_seqpos = [ ( scored_context + ( p - center_var[ idx ][ 0 ] ), r, a )
                           for p,r,a in [ center_var[ idx ] ] + haplotypes[ idx ] ]

        var_acc, var_don = adjust_for_indels( var_acc,
                                              var_don,
                                              allvars_seqpos,
                                             )

        diff_acc = ref_acc - var_acc
        diff_don = ref_don - var_don

        outtbl[ 'AL_chg' ].append( sum( ( diff_acc > 0 ) * diff_acc ) )
        outtbl[ 'AG_chg' ].append( sum( ( diff_acc < 0 ) * diff_acc ) )
        outtbl[ 'DL_chg' ].append( sum( ( diff_don > 0 ) * diff_don ) )
        outtbl[ 'DG_chg' ].append( sum( ( diff_don < 0 ) * diff_don ) )

        acc_jn_use, don_jn_use = ss_jn_use[ idx ]

        acc_wt = np.zeros( len( diff_acc ) )

        for acc_jn in acc_jn_use.keys():

            acc_wt[ acc_jn ] = acc_jn_use[ acc_jn ]

        diff_acc = ( diff_acc > 0 ) * ( acc_wt ) * ( diff_acc ) \
                     + ( diff_acc < 0 ) * ( 1 - acc_wt ) * ( diff_acc )

        don_wt = np.zeros( len( diff_acc ) )

        for don_jn in don_jn_use.keys():

            don_wt[ don_jn ] = don_jn_use[ don_jn ]

        diff_don = ( diff_don > 0 ) * ( don_wt ) * ( diff_don ) \
                     + ( diff_don < 0 ) * ( 1 - don_wt ) * ( diff_don )

        outtbl[ 'pos' ].append( center_var[ idx ][ 0 ] )
        outtbl[ 'ref' ].append( center_var[ idx ][ 1 ] )
        outtbl[ 'alt' ].append( center_var[ idx ][ 2 ] )
        outtbl[ 'other_var' ].append( ';'.join( [ ':'.join( [ str( p ), '>'.join( [ r, a ] ) ] )
                                      for p,r,a in haplotypes[ idx ] ] ) )

        outtbl[ 'DS_AGrw' ].append( np.abs( np.min( [ 0, np.min( diff_acc ) ] ) ) )
        outtbl[ 'DS_ALrw' ].append( np.max( [ 0, np.max( diff_acc ) ] ) )
        outtbl[ 'DS_DGrw' ].append( np.abs( np.min( [ 0, np.min( diff_don ) ] ) ) )
        outtbl[ 'DS_DLrw' ].append( np.max( [ 0, np.max( diff_don ) ] ) )

        outtbl[ 'DP_AGrw' ].append( np.argmin( diff_acc ) - scored_context )
        outtbl[ 'DP_ALrw' ].append( np.argmax( diff_acc ) - scored_context )
        outtbl[ 'DP_DGrw' ].append( np.argmin( diff_don ) - scored_context )
        outtbl[ 'DP_DLrw' ].append( np.argmax( diff_don ) - scored_context )

        score_keys = [ 'DS_AGrw', 'DS_ALrw', 'DS_DGrw', 'DS_DLrw' ]

        #first get the maximum probability across the difference scores
        outtbl[ 'DS_maxrw' ].append( max( outtbl[ key ][ -1 ] for key in score_keys ) )
        #then get the type of event that represents the maximum probability
        outtbl[ 'DS_maxrw_type' ].append( [ key for key in score_keys
                                          if outtbl[ key ][ -1 ] == outtbl[ 'DS_maxrw' ][ -1 ] ][ 0 ] )
        #finally, get the location of the event associated with the highest difference score
        outtbl[ 'POS_maxrw' ].append( outtbl[ outtbl[ 'DS_maxrw_type' ][ -1 ].replace( 'DS', 'DP' ) ][ -1 ] )

    outdf = pd.DataFrame( outtbl )

    return outdf

def jnuse_score_mult_variants_oneexon( annots_df,
                                        gtex_df,
                                        models,
                                        refseq,
                                        ref_name,
                                        exon_cds,
                                        chrom,
                                        center_var,
                                        haplotypes = None,
                                        ref_vars = None,
                                        scored_context_pad = 10,
                                        unscored_context = 5000,
                                        rev_strand = False ):
    """
    Wrapper function to compute junction use weighted probabilities across one exon.

    Args: annots_df (pandas df) - columns: #NAME, CHROM, STRAND, TX_START, TX_END, EXON_START, EXON_END
                                    EXON_START, EXON_END are both comma separated strings of all exon bds
          gtex_df (pandas df) - columns: chrom, jn (1-based), reads
          models (list) - SpliceAI models
          refseq (str) - fasta file for an entire chromosome
          exon_cds (tuple of ints) - exon bds in genomic coords to determined scored context length
          ref_name (str) - reference and/or annotation name to label the entries
          chrom - (str) chromosome of center variant ( format should match your annots file (ie chr3 v 3) )
          center_var (list of tuples) - center_variants to be scored:
                                ( position - (int) position (genomic coords) of center variant,
                                 reference base(s) - (str) reference base(s) relative to forward strand,
                                 alternate base(s) - (str) alternate base(s) relative to forward strand, )
          haplotypes ( list of list of tuples ) - other variants added to the variant sequences
                                                  each ith list corresponds to the ith center_var
                                        ( position - (int) position (genomic coords) of variant,
                                        reference base(s) - (str) reference base(s) relative to forward strand,
                                        alternate base(s) - (str) alternate base(s) relative to forward strand, )
                                        None indicates no additional variants
          ref_vars ( list of list of tuples ) - other variants to be added to the reference AND variant sequences
                                                each ith list corresponds to the ith center_var
                                             ( position - (int) position (genomic coords) of variant,
                                             reference base(s) - (str) reference base(s) relative to forward strand,
                                             alternate base(s) - (str) alternate base(s) relative to forward strand, )
                                             None adds no additional variants
                                             NOTE: currently only substitutions can be handled - indels will raise error
          scored_context_pad - (int) number of additional bases above exon length to score on each side of the center variant
          unscored_context - (int) number of flanking unscored bases on each side of the center variant
          rev_strand - (bool) is the exon on the reverse strand?

    Returns: outdf - (pandas df) Pandas dataframe with probabilities for each type of splice site event for each center_var
                                 and separate probabilities after masking
    """

    exon_len = exon_cds[ 1 ] - exon_cds[ 0 ]

    scored_context = exon_len + scored_context_pad

    flanks = scored_context + unscored_context

    if not haplotypes:

        haplotypes = [ [] for var in center_var ]

    if not ref_vars:

        ref_vars = [ [] for var in center_var ]

    else:

        for p,r,a in ref_vars:

            ref_name += '+' + str( p ) + ':' + r + '>' + a

    assert all( len(i) == len( center_var ) for i in [ haplotypes, ref_vars ] ), \
    'Haplotypes and ref_vars input must be either missing or the same length as the center_var list'

    #creates a giant list of ( reference seq, variant seq ) tuples
    refvarseqs = [ create_input_seq( refseq,
                                     center,
                                     hapref[ 0 ],
                                     hapref[ 1 ],
                                     get_gene_bds( annots_df,
                                                   chrom,
                                                   center[ 0 ],
                                                   scored_context,
                                                   unscored_context = unscored_context ),
                                     scored_context,
                                     rev_strand = rev_strand,
                                     unscored_context = unscored_context )
                  for center,hapref in zip( center_var, zip( haplotypes, ref_vars ) ) ]

    #creates a giant list of lists of relative acceptor/donor use
    rel_jn_use = [ get_relative_jn_use( gtex_df,
                                        chrom,
                                        center[ 0 ],
                                        get_allexon_bds( annots_df,
                                                         chrom,
                                                         center[ 0 ],
                                                         scored_context,
                                                         rev_strand = rev_strand
                                                        ),
                                        scored_context,
                                        rev_strand = rev_strand )
                   for center in center_var ]

    #this will fail if any of the other variants are indels...
    #maybe I should add some functionality to the score_variant fn to handle this...
    outdf = custom_score_variants( annots_df,
                                    models,
                                    refvarseqs,
                                    ref_name,
                                    rel_jn_use,
                                    chrom,
                                    center_var,
                                    haplotypes,
                                    scored_context,
                                    rev_strand = rev_strand,
                                  )
    return outdf
