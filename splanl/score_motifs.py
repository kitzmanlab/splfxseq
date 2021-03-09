import pandas as pd
import numpy as np
import zipfile as zp
from maxentpy import maxent

def ImportFastA(fastafile):
    """Opens a fasta file, discards the header, and returns the sequence.

    Args:
        fasta file (str): /path/and/name/to.fa
            (file must be unzipped)

    Returns:
        sequence (string): sequence of nucleotides with the case of the base retained.

    Usage:
        seq = ImportFastA('/path/and/name/to.fa')
    """
    seq=''
    with open(fastafile,'rt') as f:
        for line in f:
            if line.startswith('>'):
                next
            else:
                seq+=line.strip()
    return(seq)

def RNA_to_DNA( RNA_motif_tbl ):

    DNA_mtf_tbl = RNA_motif_tbl.copy()

    DNA_mtf_tbl.index = DNA_mtf_tbl.index.str.replace( 'U', 'T' )
    DNA_mtf_tbl.index = DNA_mtf_tbl.index.str.replace( 'u', 't' )

    return( DNA_mtf_tbl )

def ExtractKmers( seq,
                  k ):
    """Takes a sequence of length 2*k-1 and extracts all kmers of length k from that sequence.

    Args:
        seq (str): sequence to extract kmers from
        k (int): length of kmers

    Returns:
        kmers (list of str): list of k kmers from the sequence

    Usage:
        kmers=ExtractKmers('GGCATGTAACT',6)
        print(kmers)
        ['GGCATG', 'GCATGT', 'CATGTA', 'ATGTAA', 'TGTAAC', 'GTAACT']
    """
    assert len( seq ) == 2*k - 1
    return( [ seq[ i: i + k ] for i in range( k ) ] )

def score_motifs_max( refseq,
                motif_df,
               score_col,
                vardf,
                k,
                col_stem):

    tbv = vardf.copy()
    mdf = motif_df.copy()

    pos_min = tbv.pos.min()

    cloneseq = refseq[ pos_min - k: tbv.pos.max() + k ]

    wt = { i + pos_min: np.max( [ float( mdf.loc[ kseq, score_col ] ) if kseq in mdf.index else 0
          for kseq in ExtractKmers( cloneseq[ i: i + ( 2*k - 1) ], k ) ] )
          for i in range( len( cloneseq ) - 2*( k - 1 ) ) }

    wtdf = pd.DataFrame( list( wt.values() ), index = list( wt.keys() ), columns = [ col_stem + '_wtMAX' ] )

    mut = { ( i + pos_min, alt ) : np.max( [ float( mdf.loc[ kseq, score_col ] ) if kseq in mdf.index else 0
             for kseq in ExtractKmers( cloneseq[ i: i + ( k - 1) ] + alt + cloneseq[ i + k : i + ( 2*k - 1) ], k ) ] )
             for i in range( len( cloneseq ) - 2*( k - 1 ) )
             for alt in [ 'A', 'C', 'G', 'T' ] if cloneseq[ i + k - 1 ] != alt
           }

    mutdf = pd.DataFrame( list( mut.values() ), index=pd.MultiIndex.from_tuples( mut ), columns=[ col_stem+'_snvMAX' ] )
    mutdf.index = mutdf.index.set_names ( [ 'pos', 'alt' ] )

    change = { ( p, a ): mut[ ( p, a ) ] - wt[ p ] for p,a in mutdf.index }

    changedf = pd.DataFrame( list( change.values() ), index=pd.MultiIndex.from_tuples( change ), columns=[ col_stem+'_chgMAX' ] )
    changedf.index = changedf.index.set_names ( [ 'pos', 'alt' ] )

    outdf = pd.merge( tbv, wtdf, how = 'left', left_on = 'pos', right_on = wtdf.index ).set_index( [ 'pos', 'alt' ] )
    outdf = pd.merge( outdf, mutdf, how = 'left', left_index = True, right_index = True )
    outdf = pd.merge( outdf, changedf, how = 'left', left_index = True, right_index = True ).reset_index()
    firstCol = [ 'chrom', 'pos', 'ref', 'alt' ]
    outdf = outdf[ firstCol + [ c for c in outdf.columns if c not in firstCol ] ]

    return outdf

def score_motifs_mean( refseq,
                motif_df,
               score_col,
                vardf,
                k,
                col_stem):

    tbv = vardf.copy()
    mdf = motif_df.copy()

    pos_min = tbv.pos.min()

    cloneseq = refseq[ pos_min - k: tbv.pos.max() + k ]

    wt = { i + pos_min: np.mean( [ float( mdf.loc[ kseq, score_col ] ) if kseq in mdf.index else 0
          for kseq in ExtractKmers( cloneseq[ i: i + ( 2*k - 1) ], k ) ] )
          for i in range( len( cloneseq ) - 2*( k - 1 ) ) }

    wtdf = pd.DataFrame( list( wt.values() ), index = list( wt.keys() ), columns = [ col_stem + '_wtMEAN' ] )

    mut = { ( i + pos_min, alt ) : np.mean( [ float( mdf.loc[ kseq, score_col ] ) if kseq in mdf.index else 0
             for kseq in ExtractKmers( cloneseq[ i: i + ( k - 1) ] + alt + cloneseq[ i + k : i + ( 2*k - 1) ], k ) ] )
             for i in range( len( cloneseq ) - 2*( k - 1 ) )
             for alt in [ 'A', 'C', 'G', 'T' ] if cloneseq[ i + k - 1 ] != alt
           }

    mutdf = pd.DataFrame( list( mut.values() ), index=pd.MultiIndex.from_tuples( mut ), columns=[ col_stem+'_snvMEAN' ] )
    mutdf.index = mutdf.index.set_names ( [ 'pos', 'alt' ] )

    change = { ( p, a ): mut[ ( p, a ) ] - wt[ p ] for p,a in mutdf.index }

    changedf = pd.DataFrame( list( change.values() ), index=pd.MultiIndex.from_tuples( change ), columns=[ col_stem+'_chgMEAN' ] )
    changedf.index = changedf.index.set_names ( [ 'pos', 'alt' ] )

    outdf = pd.merge( tbv, wtdf, how = 'left', left_on = 'pos', right_on = wtdf.index ).set_index( [ 'pos', 'alt' ] )
    outdf = pd.merge( outdf, mutdf, how = 'left', left_index = True, right_index = True )
    outdf = pd.merge( outdf, changedf, how = 'left', left_index = True, right_index = True ).reset_index()
    firstCol = [ 'chrom', 'pos', 'ref', 'alt' ]
    outdf = outdf[ firstCol + [ c for c in outdf.columns if c not in firstCol ] ]

    return outdf

def WT_binding_df( refseq,
                    motifdf,
                    score_col,
                    vardf,
                    k,
                    col_stem):

    mdf = motifdf.copy()
    vdf = vardf.copy()

    wt_dict = { pos + 1 : [ np.mean( [ mdf.loc[ s ][ score_col ]
                for s in [ refseq[ i: i+k ]
                          for i in range( pos-k, pos ) ] ] ),
                       vdf.loc[ vdf.pos == pos + 1 ][ 'ref' ].values[0],
                       vdf.loc[ vdf.pos == pos + 1 ][ 'hgvs_pos' ].values[0]
                      ]
             for pos in range( vdf.pos.min() + k, vdf.pos.max() - k ) }

    wt_df = pd.DataFrame.from_dict( wt_dict, orient='index', columns=[col_stem+'_WT','ref','hgvs_pos' ] )
    wt_df.index.name = 'pos'
    wt_df.reset_index(inplace=True)

    return wt_df

def create_hal_input_df( varlist,
                         refseq,
                         exon_coords,
                         wt_psi ):

    out_dict = { 'var_name': [],
                 'refseq': [],
                 'altseq': [],
                 'psi': [] }

    for var in varlist:

        pos = int( var.split( ':' )[ 1 ] )

        if exon_coords[ 1 ]>= pos >= exon_coords[ 0 ]:

            ref = var.split( ':' )[ 2 ]
            alt = var.split( ':' )[ 3 ]

            #adjust for 0 and 1 based numbering
            assert refseq[ pos - 1 ] == ref, \
            'Expected reference allele does not match sequence - check numbering'

            out_dict[ 'var_name' ].append( var )
            #hal only works on exonic sequences and wants 6 downstream intronic bases in lowercase
            out_dict[ 'refseq' ].append( refseq[ exon_coords[ 0 ] - 1: exon_coords[ 1 ] ] \
                                         + refseq[ exon_coords[ 1 ]: exon_coords[ 1 ] + 6 ].lower() )
            #adds in the alternate allele at the reference spot - only tested for SNVS
            out_dict[ 'altseq' ].append( refseq[ exon_coords[ 0 ] - 1: pos - 1 ] \
                                        + alt
                                        + refseq[ pos: exon_coords[ 1 ] ]
                                         + refseq[ exon_coords[ 1 ]: exon_coords[ 1 ] + 6 ].lower() )
            out_dict[ 'psi' ] = wt_psi

    out_tbl = pd.DataFrame( out_dict )

    return out_tbl

def merge_hal( var_df,
               hal_df,
               index = [ 'varlist' ]):

    tbv = var_df.copy().set_index( index )

    thal = hal_df.copy()
    thal = thal.rename( columns = { 'VARIANT_NAME': 'varlist', 'DELTA_PSI': 'hal_chgPER' } )
    thal = thal[ index + [ 'hal_chgPER' ] ].set_index( index )

    out_tbl = pd.merge( tbv, thal, how = 'left', on = index ).reset_index()

    return out_tbl

def get_spidex_scores( chrom,
                       coords,
                       spidex_dir,
                       rev_strand = False ):

    if len( str( chrom ) ) > 2:
        chrom = int( chrom[3:] )

    if rev_strand:
        trans_tbl = str.maketrans( 'ACGT', 'TGCA' )

    out_dict = { 'chrom': [],
                 'gdna_pos_hg19': [],
                 'ref': [],
                 'alt': [],
                 'spanr_chgPER_tissue': [],
                 'spanr_chgZPER': []
               }

    zf = zp.ZipFile( spidex_dir )
    with zf.open( 'hg19_spidex.txt' ) as f:

        first_line = True

        for line in f:

            if first_line:

                first_line = False
                continue

            else:

                row = line.decode('utf-8').strip().split( '\t' )

                if int( row[ 0 ] ) != chrom:
                    continue

                pos = int( row[1] )

                if coords[ 0 ] <= pos <= coords[ 1 ]:

                    out_dict[ 'chrom' ].append( int( row[ 0 ] ) )
                    out_dict[ 'gdna_pos_hg19' ].append( pos )

                    if rev_strand:
                        out_dict[ 'ref' ].append( row[ 3 ].translate( trans_tbl ) )
                        out_dict[ 'alt' ].append( row[ 4 ].translate( trans_tbl ) )
                    else:
                        out_dict[ 'ref' ].append( row[ 3 ] )
                        out_dict[ 'alt' ].append( row[ 4 ] )
                    out_dict[ 'spanr_chgPER_tissue' ].append( float( row[ 5 ] ) )
                    out_dict[ 'spanr_chgZPER' ].append( float( row[ 6 ] ) )

                #trying to save some time since I'm not iterating efficiently here
                elif pos > coords[ 1 ]:

                    break

    out_tbl = pd.DataFrame( out_dict )

    return out_tbl

def merge_spidex( var_df,
                  spidex_df,
                  index = [ 'gdna_pos_hg19', 'ref', 'alt' ] ):

    tbv = var_df.copy().set_index( index )

    tsp = spidex_df.copy().set_index( index )
    tsp = tsp.drop( columns = 'chrom' )

    out_tbl = pd.merge( tbv, tsp, how = 'left', on = index ).reset_index()

    return out_tbl

def compute_maxent_scores( byvartbl,
                           refseq,
                           wt_accept_col,
                           wt_donor_col,
                           var_accept_col,
                           var_donor_col ):

    tbv = byvartbl.copy()

    wt_a = list( tbv.loc[ tbv[ wt_accept_col ] ].pos )
    wt_a_scores = { pos: score_acceptor( pos, refseq ) for pos in wt_a }
    wt_a_df = pd.DataFrame( wt_a_scores.items(), columns = [ 'pos', 'wt_acc_maxent'] ).set_index( 'pos' )

    wt_d = list( tbv.loc[ tbv[ wt_donor_col ] ].pos )
    wt_d_scores = { pos: score_donor( pos, refseq ) for pos in wt_d }
    wt_d_df = pd.DataFrame( wt_d_scores.items(), columns = [ 'pos', 'wt_don_maxent'] ).set_index( 'pos' )

    var_a = list( zip( tbv.loc[ tbv[ var_accept_col ] ].pos,
                       tbv.loc[ tbv[ var_accept_col ] ].alt ) )
    var_a_scores = { ( pos, alt ): score_acceptor( pos, refseq, alt ) for pos, alt in var_a }
    var_a_df = pd.Series( var_a_scores ).reset_index()
    var_a_df.columns = [ 'pos', 'alt', 'snv_acc_maxent' ]
    var_a_df = var_a_df.set_index( [ 'pos', 'alt' ] )

    var_d = list( zip( tbv.loc[ tbv[ var_donor_col ] ].pos,
                       tbv.loc[ tbv[ var_donor_col ] ].alt ) )
    var_d_scores = { ( pos, alt ): score_donor( pos, refseq, alt ) for pos, alt in var_d }
    var_d_df = pd.Series( var_d_scores ).reset_index()
    var_d_df.columns = [ 'pos', 'alt', 'snv_don_maxent' ]
    var_d_df = var_d_df.set_index( [ 'pos', 'alt' ] )

    out_tbl = tbv.set_index( [ 'pos' ] )

    out_tbl = pd.merge( out_tbl, wt_a_df, left_index = True, right_index = True, how = 'outer')
    out_tbl = pd.merge( out_tbl, wt_d_df, left_index = True, right_index = True, how = 'outer')

    out_tbl = out_tbl.reset_index().set_index( [ 'pos', 'alt' ] )

    out_tbl = pd.merge( out_tbl, var_a_df, left_index = True, right_index = True, how = 'outer')
    out_tbl = pd.merge( out_tbl, var_d_df, left_index = True, right_index = True, how = 'outer')

    out_tbl = out_tbl.reset_index()

    return out_tbl

def score_acceptor( pos,
                    refseq,
                    alt_allele = False ):

    if not alt_allele:

        assert refseq[ pos - 1 ] == 'A' or refseq[ pos - 1 ] == 'G', \
        'Reference does not contain A or G at position %i' % pos

        if refseq[ pos - 1 ] == 'A':
            score = maxent.score3( refseq[ pos - 19: pos + 4 ] )
        else:
            score = maxent.score3( refseq[ pos - 20: pos + 3 ] )

    else:

        assert alt_allele == 'A' or alt_allele == 'G', \
        'Alternate allele is not A or G'

        if alt_allele == 'A':
            score = maxent.score3( refseq[ pos - 19: pos - 1 ] + alt_allele + refseq[ pos: pos + 4 ] )
        else:
            score = maxent.score3( refseq[ pos - 20: pos - 1 ] + alt_allele + refseq[ pos: pos + 3 ] )

    return score

def score_donor( pos,
                 refseq,
                 alt_allele = False ):

    if not alt_allele:

        assert refseq[ pos - 1 ] == 'G' or refseq[ pos - 1 ] == 'T', \
        'Reference does not contain G or T at position %i' % pos

        if refseq[ pos - 1 ] == 'G':
            score = maxent.score5( refseq[ pos - 4: pos + 5 ] )
        else:
            score = maxent.score5( refseq[ pos - 5: pos + 4 ] )

    else:

        assert alt_allele == 'G' or alt_allele == 'T', \
        'Alternate allele is not G or T'

        if alt_allele == 'G':
            score = maxent.score5( refseq[ pos - 4: pos - 1 ] + alt_allele + refseq[ pos: pos + 5 ] )
        else:
            score = maxent.score5( refseq[ pos - 5: pos - 1 ] + alt_allele + refseq[ pos: pos + 4 ] )

    return score

def score_RBP_motifs(refseq,
                     motifdf,
                     vardf,
                     k,
                     col_stem ):
    """Computes change scores from WT for motifs using an existing database

    Args:
        fasta file (str): /path/and/name/to.fa (file must be unzipped)
        motifdf (pandas df) - pandas df of motif scores with the motif as the index
        vardf (pandas df) - pandas df of splicing scores by variants
        k (int) - size of the kmer
        col_stem (str) - desired stem for the columns containing the new scores

    Returns:
        vardf (pandas df) - same dataframe with the change in motif score (mean) appended
    """
    mdf = motifdf.copy()
    vdf = vardf.copy().set_index( [ 'pos', 'alt' ] )

    #score kmers for wt and mutant - assigns 0 if not in motifDict
    #numbering adjusts for 0-based and 1-based coordinates
    wt = { p: np.mean( [ mdf.loc[ kmer ] if kmer in mdf.index else 0
          for kmer in sm.ExtractKmers( refseq[ p - k: p + ( k - 1 ) ], k ) ] )
          for p in vdf.index.unique( level = 'pos' ) }

    wtdf = pd.DataFrame( list( wt.values() ), index = list( wt.keys() ), columns = [ col_stem+'_wtMEAN' ] )

    #same as above but substitute the mutated base
    mut = { ( p, a ): np.mean( [ mdf.loc[ kmer ] if kmer in mdf.index else 0
            for kmer in sm.ExtractKmers( refseq[ p - k : p - 1 ] + a + refseq[ p : p + ( k - 1 ) ], k ) ] )
            for p,a in vdf.index }

    mutdf = pd.DataFrame( list( mut.values() ), index=pd.MultiIndex.from_tuples( mut ), columns=[ col_stem+'_snvMEAN' ] )

    #outdf = pd.concat( [ vdf, changedf ],axis=1 )
    #outdf = pd.merge( vdf, wtdf, how = 'left', left_on = 'pos', right_on = wtdf.index ).set_index( [ 'pos', 'alt' ] )
    outdf = pd.merge( vdf.reset_index().set_index( 'pos' ),
                      wtdf,
                      how = 'left',
                      left_index = True,
                      right_index = True )
    outdf.index.name = 'pos'

    outdf = pd.merge( outdf.reset_index().set_index( [ 'pos', 'alt' ] ),
                      ( mutdf.reset_index()
                        .rename( columns = { 'level_0': 'pos', 'level_1': 'alt' } )
                        .set_index( [ 'pos', 'alt' ] ) ),
                      how = 'left',
                      left_index = True,
                      right_index = True ).reset_index()

    outdf[ col_stem + '_chgMEAN' ] = outdf[ col_stem + '_snvMEAN' ] - outdf[ col_stem + '_wtMEAN' ]


    firstCol = [ 'chrom', 'pos', 'ref', 'alt' ]
    outdf = outdf[ firstCol + [ c for c in outdf.columns if c not in firstCol ] ]

    return outdf

def merge_mmsplice( var_df,
                    mmsplice_df,
                    index = [ 'gdna_pos_hg19', 'ref', 'alt' ] ):

    tbv = var_df.set_index( index ).copy()

    tmm = mmsplice_df.copy()
    tmm = tmm.rename( columns = { 'delta_logit_psi': 'mmsplice_chgPERlogit' } )
    tmm = tmm[ index + [ 'mmsplice_chgPERlogit' ] ].set_index( index )

    out_tbl = pd.merge( tbv, tmm, how = 'left', on = index ).reset_index()

    return out_tbl
