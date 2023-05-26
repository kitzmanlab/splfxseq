import pandas as pd
import numpy as np
import zipfile as zp
from maxentpy import maxent
import scipy.stats as ss
import splanl.coords as cd

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
                       exon_coords,
                       motif_df,
                       vardf,
                       k,
                      col_stem = 'esrseq',
                      score_col = 'score',
                      rev_strand = False,
                      var_merge_cols = [ 'pos', 'ref', 'alt' ] ):

    tbv = vardf.copy()
    mdf = motif_df.copy()

    if not rev_strand:

        score_seq = refseq[ exon_coords[ 0 ] : exon_coords[ 1 ] ]

    else:

        score_seq = cd.rev_complement( refseq[ exon_coords[ 0 ] - 1 : exon_coords[ 1 ] ] )

    wt = { ( i, score_seq[ i + k - 1 ] ) :
            np.mean( [ float( mdf.loc[ kseq, score_col ] ) if kseq in mdf.index else 0
          for kseq in ExtractKmers( score_seq[ i: i + ( 2*k - 1) ], k ) ] )
          for i in range( len( score_seq ) - 2*( k - 1 ) ) }

    wtdf = pd.DataFrame( list( wt.values() ),
                         index = pd.MultiIndex.from_tuples( wt ),
                         columns = [ col_stem + '_wtMEAN' ] )
    wtdf.index = wtdf.index.set_names( [ var_merge_cols[ 0 ], var_merge_cols[ 1 ] ] )
    wtdf = wtdf.reset_index()

    mut = { ( i, score_seq[ i + k - 1 ], alt ) : np.mean( [ float( mdf.loc[ kseq, score_col ] ) if kseq in mdf.index else 0
             for kseq in ExtractKmers( score_seq[ i: i + ( k - 1) ] + alt + score_seq[ i + k : i + ( 2*k - 1) ], k ) ] )
             for i in range( len( score_seq ) - 2*( k - 1 ) )
             for alt in [ 'A', 'C', 'G', 'T' ] if score_seq[ i + k - 1 ] != alt
           }

    mutdf = pd.DataFrame( list( mut.values() ),
                          index=pd.MultiIndex.from_tuples( mut ),
                          columns=[ col_stem+'_snvMEAN' ] )
    mutdf.index = mutdf.index.set_names ( var_merge_cols )
    mutdf = mutdf.reset_index()

    change = { ( p, score_seq[ p + k - 1 ], a ): mut[ ( p, score_seq[ p + k - 1 ],  a ) ] - wt[ ( p, score_seq[ p + k - 1 ] ) ]
                for p, a in zip( mutdf[ var_merge_cols [ 0 ] ], mutdf[ var_merge_cols [ 2 ] ] ) }

    changedf = pd.DataFrame( list( change.values() ), index=pd.MultiIndex.from_tuples( change ), columns=[ col_stem+'_chgMEAN' ] )
    changedf.index = changedf.index.set_names ( var_merge_cols )
    changedf = changedf.reset_index()

    outdf = wtdf.set_index( [ var_merge_cols[ 0 ], var_merge_cols[ 1 ] ] ).merge( mutdf.set_index( [ var_merge_cols[ 0 ], var_merge_cols[ 1 ] ] ),
                        how = 'outer',
                        left_index = True,
                        right_index = True ).reset_index()

    outdf = outdf.set_index( var_merge_cols ).merge( changedf.set_index( var_merge_cols ),
                                                     how = 'outer',
                                                     left_index = True,
                                                     right_index = True ).reset_index()

    if not rev_strand:

        outdf[ var_merge_cols[ 0 ] ] += exon_coords[ 0 ] + ( k - 1 )

    else:

        outdf[ var_merge_cols[ 0 ] ] = [ exon_coords[ 1 ] - ( k - 1 ) - p for p in outdf[ var_merge_cols[ 0 ] ] ]

    return outdf

def merge_ke( tbl_by_var,
              ke_scores,
              index_cols = [ 'pos', 'ref', 'alt' ] ):

    tbv = tbl_by_var.set_index( index_cols ).copy()

    ke = ke_scores.set_index( index_cols ).copy()

    outdf = tbv.merge( ke,
                       how = 'left',
                       left_index = True,
                       right_index = True ).reset_index()

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

def create_hal_skip_input_df( varlist,
                              refseq,
                              exon_coords,
                              wt_psi,
                              rev_strand = False ):
#this was formerly known as create_hal_input_df - previous notebooks will use that function

    out_dict = { 'var_name': [],
                 'altseq': [], }

    ex_start, ex_end = exon_coords

    if not rev_strand:
        #hal only works on exonic sequences and wants 6 downstream intronic bases in lowercase
        score_ref = refseq[ ex_start: ex_end ].upper() + refseq[ ex_end: ex_end + 6 ].lower()
    else:
        score_ref = cd.rev_complement( refseq[ ex_start - 6: ex_start ].lower() + refseq[ ex_start : ex_end ].upper() )

    for var in varlist:

        var_l = var.split( ':' )

        pos = int( var_l[ 1 ] )

        if ex_end >= pos > ex_start:

            ref = var_l[ 2 ].upper()
            alt = var_l[ 3 ].upper()

            #adjust for 0 and 1 based numbering
            assert refseq[ pos - 1: pos - 1 + len( ref ) ].upper() == ref, \
            'Expected reference allele does not match sequence - check numbering'

            out_dict[ 'var_name' ].append( var )
            #out_dict[ 'refseq' ].append( score_ref )

            if not rev_strand:
                #adds in the alternate allele at the reference spot - only tested for SNVS
                out_dict[ 'altseq' ].append( refseq[ ex_start : pos - 1 ].upper() \
                                            + alt
                                            + refseq[ pos - 1 + len( ref ): exon_coords[ 1 ] ].upper()
                                            + refseq[ ex_end: ex_end + 6 ].lower() )

            else:
                out_dict[ 'altseq' ].append( cd.rev_complement( refseq[ ex_start - 6: ex_start ].lower() \
                                                                 + refseq[ ex_start : pos - 1 ].upper() \
                                                                 + alt \
                                                                 + refseq[ pos - 1 + len( ref ): ex_end ].upper() )  )

    out_tbl = pd.DataFrame( out_dict )

    out_tbl[ 'refseq' ] = score_ref
    out_tbl[ 'psi' ] = wt_psi

    return out_tbl

def create_hal_5ss_input_df( varlist,
                             refseq,
                             acc_coords,
                             don_coords,
                             wt_psi,
                             rev_strand = False ):

    out_dict = { 'var_name': [],
                 'refseq': [],
                 'altseq': [],
                 'psi': [] }

    assert len( don_coords ) == 2, 'HAL only scores one alt donor at a time!'

    assert don_coords[ 0 ] < don_coords[ 1 ], 'Please enter your donor coordinates from least to greatest value'

    if not rev_strand:

        for var in varlist:

            pos = int( var.split( ':' )[ 1 ] )

            if acc_coords >= pos > don_coords[ 1 ]:

                ref = var.split( ':' )[ 2 ].upper()
                alt = var.split( ':' )[ 3 ].upper()

                #adjust for 0 and 1 based numbering
                assert refseq[ pos - 1: pos - 1 + len( ref ) ].upper() == ref, \
                'Expected reference allele does not match sequence - check numbering'

                out_dict[ 'var_name' ].append( var )

                out_dict[ 'refseq' ].append( refseq[ acc_coords : don_coords[ 0 ] ].upper() \
                                            + refseq[ don_coords[ 0 ]: don_coords[ 1 ] ].lower() )

                if pos <= don_coords[ 0 ]:
                    out_dict[ 'altseq' ].append( refseq[ acc_coords[ 0 ] : pos - 1 ].upper()
                                                 + alt
                                                 + refseq[ pos - 1 + len( ref ): don_coords[ 0 ] ].upper()
                                                 + refseq[ don_coords[ 0 ]: don_coords[ 1 ] ].lower() )
                else:
                    out_dict[ 'altseq' ].append( refseq[ acc_coords[ 0 ] : don_coords[ 0 ] ].upper()
                                                 + refseq[ don_coords[ 0 ]: pos - 1 ].lower()
                                                 + alt.lower()
                                                 + refseq[ pos - 1 + len( ref ): don_coords[ 1 ] ].lower() )

                out_dict[ 'psi' ].append( wt_psi )

    else:

        for var in varlist:

            pos = int( var.split( ':' )[ 1 ] )

            if don_coords[ 1 ] >= pos > acc_coords:

                ref = var.split( ':' )[ 2 ].upper()
                alt = var.split( ':' )[ 3 ].upper()

                #adjust for 0 and 1 based numbering
                assert refseq[ pos - 1: pos - 1 + len( ref ) ].upper() == ref, \
                'Expected reference allele does not match sequence - check numbering'

                out_dict[ 'var_name' ].append( var )

                out_dict[ 'refseq' ].append( cd.rev_complement( refseq[ don_coords[ 0 ]: don_coords[ 1 ] ].lower()
                                                                 + refseq[ don_coords[ 1 ] : acc_coords ].upper() ) )

                if pos <= don_coords[ 1 ]:
                    out_dict[ 'altseq' ].append( cd.rev_complement( refseq[ don_coords[ 0 ]: pos - 1 ].lower()
                                                                     + alt.lower()
                                                                     + refseq[ pos - 1 + len( ref ): don_coords[ 1 ] ].lower()
                                                                     + refseq[ don_coords[ 1 ] : acc_coords ].upper() ) )
                else:
                    out_dict[ 'altseq' ].append( cd.rev_complement( refseq[ don_coords[ 0 ]: don_coords[ 1 ] ].lower()
                                                                     + refseq[ don_coords[ 1 ] : pos - 1 ].upper()
                                                                     + alt
                                                                     + refseq[ pos - 1 + len( ref ): acc_coords ].upper() ) )

                out_dict[ 'psi' ].append( wt_psi )

    out_tbl = pd.DataFrame( out_dict )

    return out_tbl

def merge_hal( var_df,
               hal_df,
               index = [ 'varlist' ],
               out_raw_col = 'hal_PSI',
               out_chg_col = 'hal_chgPER',
             ):

    tbv = var_df.copy().set_index( index )

    thal = hal_df.copy()
    thal = thal.rename( columns = { 'VARIANT_NAME': 'varlist', 'MUT_PSI': out_raw_col, 'DELTA_PSI': out_chg_col } )
    thal = thal[ index + [ out_raw_col, out_chg_col ] ].set_index( index )

    out_tbl = pd.merge( tbv, thal, how = 'left', on = index ).reset_index()

    return out_tbl

def merge_squirls( tbl_by_var,
                   squirls_df,
                   out_score_col = 'squirls_score',
                   out_sdv_col = 'squirls_sdv',
                   index_cols = [ 'transcript_id', 'chrom', 'hg19_pos', 'ref', 'alt' ] ):

    tbv = tbl_by_var.set_index( index_cols ).copy()

    #if transcript:

        #squirls_df = squirls_df.loc[ squirls_df.tx_accession.str.startswith( transcript ) ].copy()

    squirls = squirls_df.rename( columns = { 'squirls_score': out_score_col,
                                              'pos': [ col for col in index_cols if 'pos' in col ][ 0 ],
                                              'tx_accession': 'squirls_transcript_id' } ).set_index( index_cols ).copy()

    assert len( squirls ) > 0, 'No data selected to merge - is the transcript in the table?'

    squirls[ out_sdv_col ] = squirls.interpretation == 'pathogenic'

    tbv = tbv.merge( squirls[ [ out_score_col, out_sdv_col, 'squirls_transcript_id' ] ],
                     how = 'left',
                     left_index = True,
                     right_index = True ).reset_index()

    return tbv

def get_spidex_scores_by_region( chrom,
                                 coords,
                                 spidex_dir,
                                 rev_strand = False ):

    print( 'Sit back and relax - it takes ~4 minutes to unzip the spidex file...' )

    if len( str( chrom ) ) > 2:
        chrom = chrom[ 3: ]

    zf = zp.ZipFile( spidex_dir )

    spidex = pd.read_csv( zf.open( 'hg19_spidex.txt' ),
                          sep = '\t',
                          dtype = { '#Chr': str } )

    zf.close()

    out_df = spidex.loc[ ( spidex[ '#Chr' ] == chrom ) &
                         ( spidex.Start > coords[ 0 ] ) &
                         ( coords[ 1 ] >= spidex.End ) ].copy()

    out_df = out_df.drop( columns = 'End' )

    out_df = out_df.rename( columns = { '#Chr': 'chrom',
                                        'Start': 'hg19_pos',
                                        'Ref': 'ref',
                                        'Alt': 'alt',
                                        'dpsi_max_tissue': 'spanr_chgPER_tissue',
                                        'dpsi_zscore': 'spanr_chgZPER' } )

    if rev_strand:

        trans_tbl = str.maketrans( 'ACGTNacgtn', 'TGCANtgcan' )

        out_df[ 'ref' ] = [ r.translate( trans_tbl ) for r in out_df.ref ]
        out_df[ 'alt' ] = [ a.translate( trans_tbl ) for a in out_df.alt ]

    return out_df

def get_spidex_scores_by_variant( varlist,
                                  spidex_dir,
                                  rev_strand = False ):

    print( 'Sit back and relax - it takes ~4 minutes to unzip the spidex file...' )

    zf = zp.ZipFile( spidex_dir )

    spidex = pd.read_csv( zf.open( 'hg19_spidex.txt' ),
                          sep = '\t',
                          dtype = { '#Chr': str } )

    zf.close()

    if rev_strand:

        trans_tbl = str.maketrans( 'ACGTNacgtn', 'TGCANtgcan' )

        spidex[ 'Ref' ] = [ r.translate( trans_tbl ) for r in spidex.Ref ]
        spidex[ 'Alt' ] = [ a.translate( trans_tbl ) for a in spidex.Alt ]

    spidex = spidex.set_index( [ '#Chr', 'Start', 'Ref', 'Alt' ] )

    out_df = spidex.loc[ varlist ].reset_index().copy()

    out_df = out_df.drop( columns = 'End' )

    out_df = out_df.rename( columns = { '#Chr': 'chrom',
                                        'Start': 'hg19_pos',
                                        'Ref': 'ref',
                                        'Alt': 'alt',
                                        'dpsi_max_tissue': 'spanr_chgPER_tissue',
                                        'dpsi_zscore': 'spanr_chgZPER' } )

    return out_df

def clean_tbx_spidex( spidex,
                      rev_strand = False ):

    spid = spidex.copy()

    if rev_strand:

        trans_tbl = str.maketrans( 'ACGTNacgtn', 'TGCANtgcan' )

        spid[ 'Ref' ] = [ r.translate( trans_tbl ) for r in spid.Ref ]
        spid[ 'Alt' ] = [ a.translate( trans_tbl ) for a in spid.Alt ]

    spid = spid.rename( columns = { '#Chr': 'chrom',
                                    'Start': 'hg19_pos',
                                    'Ref': 'ref',
                                    'Alt': 'alt',
                                    'dpsi_max_tissue': 'spanr_chgPER_tissue',
                                    'dpsi_zscore': 'spanr_chgZPER' } )

    return spid

def merge_spidex( var_df,
                  spidex_df,
                  index = [ 'chrom', 'hg19_pos', 'ref', 'alt' ] ):

    tbv = var_df.copy().set_index( index )

    tsp = spidex_df.copy().set_index( index )
    tsp = tsp.drop( columns = [ 'End' ] )
    #tsp = tsp.drop( columns = [ 'chrom', 'End' ] )

    out_tbl = pd.merge( tbv, tsp, how = 'left', on = index ).reset_index()

    return out_tbl

def get_scap_scores( scap_vcf,
                     replace_common = np.nan ):

    out_d = { 'chrom': [],
              'hg19_pos': [],
              'ref': [],
              'alt': [],
              'scap_type': [],
              'scap_raw': [],
              'scap_sens': [],
              'scap_raw_dom': [],
              'scap_sens_dom': [],
              'scap_raw_rec': [],
              'scap_sens_rec': [],}

    for var in scap_vcf:

        info_tup = var.info[ 'SCAP' ]

        for idx, alt in enumerate( var.alts ):

            out_d[ 'chrom' ].append( str( var.chrom ) )
            out_d[ 'hg19_pos' ].append( var.pos )
            out_d[ 'ref' ].append( var.ref )
            out_d[ 'alt' ].append( alt )

            info_cols = info_tup[ idx ].split( ':' )

            out_d[ 'scap_type' ].append( info_cols[ 1 ] )
            out_d[ 'scap_raw' ].append( info_cols[ 2 ] )
            out_d[ 'scap_sens' ].append( info_cols[ 3 ] )
            out_d[ 'scap_raw_dom' ].append( info_cols[ 4 ] )
            out_d[ 'scap_sens_dom' ].append( info_cols[ 5 ] )
            out_d[ 'scap_raw_rec' ].append( info_cols[ 6 ] )
            out_d[ 'scap_sens_rec' ].append( info_cols[ 7 ] )

    outdf = pd.DataFrame( out_d )

    outdf = outdf.replace( '.', np.nan )
    outdf = outdf.replace( 'COMMON', replace_common )

    try:

        score_cols = [ 'scap_raw', 'scap_sens', 'scap_raw_dom', 'scap_sens_dom', 'scap_raw_rec', 'scap_sens_rec', ]
        outdf[ score_cols ] = outdf[ score_cols ].astype( float )

    except:

        print( 'Cannot coerce scap scores to float' )

    outdf[ 'scap_raw_max' ] = outdf[ [ col for col in outdf if 'scap_raw' in col ] ].max( axis = 1 )
    outdf[ 'scap_sens_min' ] = outdf[ [ col for col in outdf if 'scap_sens' in col ] ].min( axis = 1 )

    return outdf

def merge_scap( tbl_by_var,
                scap_scores,
              index_cols = [ 'chrom', 'hg19_pos', 'ref', 'alt' ] ):

    tbv = tbl_by_var.set_index( index_cols ).copy()

    scap = scap_scores.set_index( index_cols ).copy()

    outdf = tbv.merge( scap,
                       how = 'left',
                       left_index = True,
                       right_index = True ).reset_index()

    return outdf

def merge_splai( tbl_by_var,
                 splai_scores,
                 index_cols = [ 'chrom', 'hg19_pos', 'ref', 'alt' ],
                 pos_col_idx = 1, ):

    tbv = tbl_by_var.set_index( index_cols ).copy()

    if index_cols[ pos_col_idx ] not in splai_scores:
        splai_scores = splai_scores.rename( columns = { 'pos': index_cols[ pos_col_idx ] } )

    splai = splai_scores.set_index( index_cols ).copy()

    outdf = tbv.merge( splai,
                       how = 'left',
                       left_index = True,
                       right_index = True ).reset_index()

    return outdf

def remove_training( tbl_by_var,
                     train_df,
                     score_cols,
                     train_col,
                     out_col_suff = '_not',
                     index = [ 'chrom', 'hg19_pos', 'ref', 'alt' ] ):

    tbv = tbl_by_var.set_index( index ).copy()

    train = train_df.set_index( index ).copy()

    tbv[ train_col ] = False

    in_train = train.index.intersection( tbv.index )

    print( '%i of training data variants intersect with measured data' % len( in_train ) )

    if len( in_train ) > 0:

        tbv.loc[ in_train, train_col ] = True

        for col in score_cols:

            tbv[ col + out_col_suff ] = [ score if not in_t else np.nan
                                                for score,in_t in zip( tbv[ col ], tbv[ train_col ] ) ]

    tbv = tbv.reset_index()

    return tbv

def maxent_score_wt( refseq,
                     pos_l,
                     rev_strand = False,
                     pos_out_col = 'pos' ):

    outtbl = {}

    outtbl[ pos_out_col  ] = pos_l

    outtbl[ 'ref' ] = [ refseq[ p - 1 ] for p in pos_l ]

    if not rev_strand:

        outtbl[ 'maxent_wt_acc' ] = [ maxent.score3( refseq[ pos - 1 - 20: pos - 1 + 3 ] ) for pos in pos_l ]
        outtbl[ 'maxent_wt_don' ] = [ maxent.score5( refseq[ pos - 1 - 2: pos - 1 + 7 ] ) for pos in pos_l ]

    else:

        outtbl[ 'maxent_wt_acc' ] = [ maxent.score3( cd.rev_complement( refseq[ pos - 1 - 2: pos - 1 + 21 ] ) )
                                  for pos in pos_l ]
        outtbl[ 'maxent_wt_don' ] = [ maxent.score5( cd.rev_complement( refseq[ pos - 1 - 6: pos - 1 + 3 ] ) )
                                  for pos in pos_l ]

    outdf = pd.DataFrame( outtbl )

    return outdf

def maxent_score_donors( refseq,
                         tbl_by_var,
                         pos_col,
                         ref_col,
                         alt_col,
                         donor_pos_l,
                         rev_strand = False ):

    #NOW in score_motifs!

    tbv = tbl_by_var.copy()

    if not rev_strand:

        m_pos = { pos: ( pos - 1 - 2, pos - 1 + 7 ) for pos in donor_pos_l }

    else:

        m_pos = { pos: ( pos - 1 - 6, pos - 1 + 3 ) for pos in donor_pos_l }

    maxent_scores = { pos: [] for pos in donor_pos_l }

    for dpos in donor_pos_l:

        for p,ra in zip( tbv[ pos_col ], zip( tbv[ ref_col ], tbv[ alt_col ] ) ):

            p = p - 1
            ref,alt = ra

            if p < m_pos[ dpos ][ 0 ] or p >= m_pos[ dpos ][ 1 ]:

                maxent_scores[ dpos ].append( np.nan )

            else:

                if not rev_strand:

                    assert ref.upper() == refseq[ p ].upper(), 'Refseq does not match for %i:%s>%s: refseq has %s as reference' % ( p + 1, ref, alt, refseq[ p ].upper() )
                    maxent_scores[ dpos ].append( maxent.score5( refseq[ m_pos[ dpos ][ 0 ]: p ] + alt + refseq[ p + 1: m_pos[ dpos ][ 1 ] ] ) )

                else:

                    assert cd.rev_complement( ref.upper() ) == cd.rev_complement( refseq[ p ] ).upper(), 'Refseq does not match for %i:%s>%s: refseq has %s as reference' % ( p + 1, ref, alt, refseq[ p ].upper() )
                    maxent_scores[ dpos ].append( maxent.score5( cd.rev_complement( refseq[ m_pos[ dpos ][ 0 ]: p ] + alt + refseq[ p + 1: m_pos[ dpos ][ 1 ] ] ) ) )

    return maxent_scores

def maxent_score_acceptors( refseq,
                            tbl_by_var,
                            pos_col,
                            ref_col,
                            alt_col,
                            acceptor_pos_l,
                            rev_strand = False ):

    tbv = tbl_by_var.copy()

    if not rev_strand:

        m_pos = { pos: ( pos - 1 - 20, pos - 1 + 3 ) for pos in acceptor_pos_l }

    else:

        m_pos = { pos: ( pos - 1 - 2, pos - 1 + 21 ) for pos in acceptor_pos_l }

    maxent_scores = { pos: [] for pos in acceptor_pos_l }

    for apos in acceptor_pos_l:

        for p,ra in zip( tbv[ pos_col ], zip( tbv[ ref_col ], tbv[ alt_col ] ) ):

            p = p - 1
            ref,alt = ra

            if p < m_pos[ apos ][ 0 ] or p >= m_pos[ apos ][ 1 ]:

                maxent_scores[ apos ].append( np.nan )

            else:

                if not rev_strand:

                    assert ref.upper() == refseq[ p ].upper(), 'Refseq does not match for %i:%s>%s: refseq has %s as reference' % ( p + 1, ref, alt, refseq[ p ].upper() )
                    maxent_scores[ apos ].append( maxent.score3( refseq[ m_pos[ apos ][ 0 ]: p ] + alt + refseq[ p + 1: m_pos[ apos ][ 1 ] ] ) )

                else:

                    assert cd.rev_complement( ref ).upper() == cd.rev_complement( refseq[ p ] ).upper(), 'Refseq does not match for %i:%s>%s: refseq has %s as reference' % ( p + 1, ref, alt, refseq[ p ].upper() )
                    maxent_scores[ apos ].append( maxent.score3( cd.rev_complement( refseq[ m_pos[ apos ][ 0 ]: p ] + alt + refseq[ p + 1: m_pos[ apos ][ 1 ] ] ) ) )

    return maxent_scores

def get_pangolin_scores( pangolin_vcf ):

    out_d = { 'chrom': [],
              'hg19_pos': [],
              'ref': [],
              'alt': [],
              'gene_id': [],
              'pang_incr': [],
              'pang_incr_pos': [],
              'pang_decr': [],
              'pang_decr_pos': [],
              'pang_max': [],
              'pang_max_type': [],
              'pang_max_pos': [] }

    for var in pangolin_vcf:

        if 'Pangolin' not in var.info:

            print( 'Variant %s:%s>%s was not scored' % ( var.pos, var.ref, var.alts[ 0 ] ) )
            continue

        info_str = var.info[ 'Pangolin' ][ 0 ]

        rows_l = info_str.split( 'Warnings:' )

        assert rows_l[ -1 ] == '', 'Careful - you might be cutting out a transcript score!'

        rows_l = rows_l[ : -1 ]

        for row in rows_l:

            out_d[ 'chrom' ].append( str( var.chrom ) )
            out_d[ 'hg19_pos' ].append( var.pos )
            out_d[ 'ref' ].append( var.ref )
            out_d[ 'alt' ].append( var.alts[ 0 ] )
            out_d[ 'gene_id' ].append( row.split( '|' )[ 0 ] )
            out_d[ 'pang_incr' ].append( float( row.split( '|' )[ 1 ].split( ':' )[ 1 ] ) )
            out_d[ 'pang_incr_pos' ].append( int( row.split( '|' )[ 1 ].split( ':' )[ 0 ] ) )
            out_d[ 'pang_decr' ].append( float( row.split( '|' )[ 2 ].split( ':' )[ 1 ] ) )
            out_d[ 'pang_decr_pos' ].append( int( row.split( '|' )[ 2 ].split( ':' )[ 0 ] ) )

            score_keys = [ 'pang_incr', 'pang_decr' ]

            #first get the maximum probability across the difference scores
            out_d[ 'pang_max' ].append( max( [ out_d[ key ][ -1 ] for key in score_keys ], key = abs ) )

            #then get the type of event that represents the maximum probability
            out_d[ 'pang_max_type' ].append( [ key for key in score_keys
                                          if out_d[ key ][ -1 ] == out_d[ 'pang_max' ][ -1 ] ][ 0 ] )

            #finally, get the location of the event associated with the highest difference score
            out_d[ 'pang_max_pos' ].append( out_d[ out_d[ 'pang_max_type' ][ -1 ] + '_pos' ][ -1 ] )

    outdf = pd.DataFrame( out_d )

    return outdf

def merge_pangolin( tbl_by_var,
                    pangolin_scores,
                    index_cols = [ 'gene_id', 'chrom', 'hg19_pos', 'ref', 'alt' ] ):

    tbv = tbl_by_var.set_index( index_cols ).copy()

    pang = pangolin_scores.set_index( index_cols ).copy()

    outdf = tbv.merge( pang,
                       how = 'left',
                       left_index = True,
                       right_index = True ).reset_index()

    return outdf

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
                    index = [ 'transcript_id', 'chrom', 'hg19_pos', 'ref', 'alt' ] ):

    tbv = var_df.set_index( index ).copy()

    tmm = mmsplice_df.copy()

    #if transcript:

        #tmm = tmm.loc[ tmm.transcript_id.str.startswith( transcript ) ].copy()

    tmm = tmm.rename( columns = { 'delta_logit_psi': 'mmsplice_chgPERlogit' } )
    tmm = tmm.rename( columns = { 'pathogenicity': 'mmsplice_pathogen' } )
    tmm = tmm.rename( columns = { 'efficiency': 'mmsplice_eff' } )
    if 'transcript_id' in index:
        tmm[ 'mmsplice_trans' ] = tmm.transcript_id
    else:
        tmm = tmm.rename( columns = { 'transcript_id': 'mmsplice_trans' } )
    tmm = tmm[ index + [ 'mmsplice_chgPERlogit', 'mmsplice_pathogen', 'mmsplice_eff', 'mmsplice_trans' ] ].set_index( index )

    if not tmm.index.is_unique:

        tmm = tmm.reset_index().groupby( index + [ 'mmsplice_trans' ] )[ [ 'mmsplice_chgPERlogit', 'mmsplice_pathogen', 'mmsplice_eff' ] ].max()

    out_tbl = pd.merge( tbv, tmm, how = 'left', on = index ).reset_index()

    return out_tbl

def nom_RBPs( tbl_byvar,
              rbp_info,
              sdv_col,
              region_cds,
              wt_bind_thresh,
              chg_thresh,
              p_thresh,
              rbp_suff = '_wtMAX', ):

    tbv = tbl_byvar.loc[ ( tbl_byvar.pos >= region_cds[ 0 ] ) & ( tbl_byvar.pos <= region_cds[ 1 ] ) ].copy()
    info = rbp_info.copy()

    rbp_cols = [ col for col in tbv if col.endswith( rbp_suff ) ]

    #counts number of sdvs
    sdv = tbv[ sdv_col ].sum()

    #counts total number of variants in region
    total = len( tbv )

    out_dict = {}

    for motif_id in rbp_cols:

        rbp_name = info[ ( info.Motif_ID == motif_id[ : -len( rbp_suff ) ] ) \
                         & ( info.RBP_Species == 'Homo_sapiens') ].RBP_Name.tolist()

        #some of the RBP aren't from humans - skip these
        if len( rbp_name ) == 0:
            continue

        #gets the maximum WT binding 'affinity' to SDVs in the defined region
        wtmax = ( tbv.loc[ ( tbv[ sdv_col ] ) ][ [ 'pos',  motif_id ] ]
                .groupby( [ 'pos' ] )
                .agg( np.nanmean ) ).max()[ 0 ]

        #if it is below the set threshold, the protein doesn't bind there - move on
        if wtmax < wt_bind_thresh:
            continue

        #counts the number of variants (not necessarily SDV) which break the motif
        #counting the number of variants in the defined region with a change score lower than the threshold
        snv_below =  ( tbv[ motif_id.replace( 'wt', 'chg' ) ] <= chg_thresh ).sum()

        #if no variants break the motif, move on
        if snv_below == 0:
            continue

        #counts number of sdv positions with WT binding above threshold
        wt_above = ( ( tbv.loc[ ( tbv[ sdv_col ] ) ][ [ 'pos', motif_id ] ]
                    .groupby( [ 'pos', ] )
                    .agg( np.nanmean ) ) >= wt_bind_thresh ).sum()

        #counts number of SDVs with a change score below threshold
        sdv_snv_below = ( tbv.loc[ tbv[ sdv_col ] ]
                         [ motif_id.replace( 'wt', 'chg' ) ] <= chg_thresh ).sum()

        assert sdv_snv_below + ( sdv - sdv_snv_below ) + ( snv_below - sdv_snv_below ) + ( total - sdv ) - ( snv_below - sdv_snv_below ) == len( tbv )
        #compares the number of variants breaking the motif between SDVs and neutral variants
        odds, p = ss.fisher_exact( [ [ sdv_snv_below, sdv - sdv_snv_below ],
                                       [ snv_below - sdv_snv_below, ( total - sdv ) - ( snv_below - sdv_snv_below ) ] ]
                                   )

        #checks that the relationship is significant, and in the right direction (odds ratio > 1)
        if p < p_thresh and odds > 1:
            out_dict[ '_'.join( rbp_name + [ motif_id[ : -len( rbp_suff ) ] ] ) ] = \
                           ( wtmax,
                             ss.fisher_exact( [ [ sdv_snv_below,
                                                  sdv - sdv_snv_below ],
                                                [ snv_below - sdv_snv_below,
                                                ( total - sdv ) - ( snv_below - sdv_snv_below ) ] ] )
                            )

    return out_dict

def merge_cadd( var_df,
                cadd_df,
                index = [ 'pos', 'ref', 'alt' ] ):

    tbv = var_df.copy().set_index( index )

    tcadd = cadd_df.copy()
    tcadd = tcadd.rename( columns = { 'Pos': 'pos',
                                      'Ref': 'ref',
                                      'Alt': 'alt',
                                      'RawScore': 'cadd_raw',
                                      'PHRED': 'cadd_scaled' } )

    #the annotations are creating duplicate rows - scores appear to be equal
    tcadd = tcadd[ index + [ 'cadd_raw', 'cadd_scaled' ] ].groupby( index ).mean()

    out_tbl = pd.merge( tbv, tcadd, how = 'left', on = index ).reset_index()

    return out_tbl

def branchpointer_ds_scores( branchpointer_df ):

    bp = branchpointer_df.copy()

    bp_ref = bp.loc[ bp.status == 'REF' ].copy()

    bp_alt = bp.loc[ bp.status == 'ALT' ].copy()

    #to3primepoint is the key column to compare the probabilities
    idx_cols = [ 'seqnames', 'start', 'end', 'width', 'strand', 'id', 'ref_allele', 'alt_allele', 'to_3prime',
                 'to_5prime', 'same_gene', 'exon_3prime', 'exon_5prime', 'to_3prime_point', 'to_5prime_point',
                  'test_site', ]
    bp_m = bp_ref.set_index( idx_cols ).merge( bp_alt.set_index( idx_cols ),
                                               how = 'outer',
                                               left_index = True,
                                               right_index = True,
                                               suffixes = ( '_ref', '_alt' ) ).reset_index()


    bp_m[ 'DBP' ] = bp_m.branchpoint_prob_alt - bp_m.branchpoint_prob_ref

    return bp_m

def branchpointer_ds_max( branchpointer_ds_df ):

    bp = branchpointer_ds_df.copy()

    bp[ 'DBP_abs' ] = bp.DBP.abs()

    ds_max = bp.groupby( [ 'id', 'start', 'ref_allele', 'alt_allele' ] ).DBP_abs.max().reset_index().copy()

    ds_sign = bp.loc[ bp.groupby( [ 'id', 'start', 'ref_allele', 'alt_allele' ] ).DBP.apply( lambda x: x.abs().idxmax() ) ][ [ 'id', 'start', 'ref_allele', 'alt_allele', 'DBP' ] ].copy()

    ds_sign = ds_sign.rename( columns = { 'DBP': 'DBP_max_signed' } )
    ds_sign[ 'DBP_event' ] = [ 'GAIN' if dbp > 0 else 'LOSS' for dbp in ds_sign.DBP_max_signed ]

    pos_max = bp.loc[ bp.groupby( [ 'id', 'start', 'ref_allele', 'alt_allele' ] ).DBP.apply( lambda x: x.abs().idxmax() ) ][ [ 'id', 'start', 'ref_allele', 'alt_allele', 'test_site' ] ].copy()

    pos_max[ 'DBP_POS_max' ] = pos_max.test_site - pos_max.start

    idx_cols = [ 'id', 'start', 'ref_allele', 'alt_allele' ]
    outdf = ds_max.set_index( idx_cols )[ [ 'DBP_abs' ] ].merge( ds_sign.set_index( idx_cols )[ [ 'DBP_event' ] ],
                                                                 how = 'outer',
                                                                 left_index = True,
                                                                 right_index = True )
    outdf = outdf.merge( pos_max.set_index( idx_cols )[ [ 'DBP_POS_max' ] ],
                         how = 'outer',
                         left_index = True,
                         right_index = True ).reset_index()

    outdf = outdf.rename( columns = { 'start': 'hg19_pos',
                                      'ref_allele': 'ref',
                                      'alt_allele': 'alt',
                                      'DBP_abs': 'DBP_max' } )

    return outdf

def merge_branchpointer( tbl_by_var,
                         branchpointer_max_df,
                         idx_cols = [ 'hg19_pos', 'ref', 'alt' ] ):

    tbv = tbl_by_var.copy()
    bp = branchpointer_max_df.copy()

    outdf = tbv.set_index( idx_cols ).merge( bp.set_index( idx_cols )[ [ 'DBP_max', 'DBP_event', 'DBP_POS_max' ] ],
                                             how = 'outer',
                                             left_index = True,
                                             right_index = True ).reset_index()

    return outdf

def comb_bpp_scores( bpp_fn,
                     bpp_wt,
                     rev_strand = False ):

    merged_files = []

    idx_cols = [ 'bp_pos' ]

    for bpp_df in bpp_fn.values():

        bpp_df[ 'hg19_pos' ] = bpp_df[ '#id' ].apply( lambda x: int( x[ 1: ].split( ':' )[ 0 ] ) )
        bpp_df[ 'ref' + rev_strand*'_c' ] = bpp_df[ '#id' ].apply( lambda x: x.split( ':' )[ 1 ].split( '>' )[ 0 ] )
        bpp_df[ 'alt' + rev_strand*'_c' ] = bpp_df[ '#id' ].apply( lambda x: x.split( ':' )[ 1 ].split( '>' )[ 1 ] )

        merged_files.append( bpp_df.set_index( idx_cols ).merge( bpp_wt.set_index( idx_cols )[ [ 'sc_bps', 'sc_ppt', 'sc', 'zsc_bps', 'zsc_ppt', 'zsc' ] ],
                                                                 left_index = True,
                                                                 right_index = True,
                                                                 how = 'outer',
                                                                 suffixes = ( '', '_wt' ) ).reset_index() )

    allvar_bpp = pd.concat( merged_files,
                             ignore_index = True )

    for col in [ 'zsc_bps', 'zsc_ppt', 'zsc' ]:

        allvar_bpp[ 'DS_' + col ] = allvar_bpp[ col ] - allvar_bpp[ col + '_wt' ]

    return allvar_bpp

def bpp_ds_max( bpp_ds_df,
                first_ex_bp,
                zcols = [ 'DS_zsc_bps', 'DS_zsc_ppt', 'DS_zsc' ],
                rev_strand = False ):

    bpp = bpp_ds_df.copy()

    idx_cols = [ '#id', 'ref' + rev_strand*'_c', 'alt' + rev_strand*'_c' ]

    for col in zcols:

        bpp[ col + '_abs' ] = bpp[ col ].abs()

    ds_max = bpp.groupby( idx_cols + [ 'hg19_pos' ] )[ [ col + '_abs' for col in zcols ] ].max().reset_index().copy()

    ds_sign = bpp.groupby( idx_cols )[ zcols ].apply( lambda x: x.abs().idxmax() ).reset_index().copy()

    for col in zcols:

        ds_sign[ col + '_max' ] = bpp.loc[ ds_sign[ col ], col ].tolist()

    pos_max = bpp.groupby( idx_cols )[ zcols ].apply( lambda x: x.abs().idxmax() ).reset_index().copy()

    for col in zcols:

        pos_max[ col + 'bp_pos' ] = bpp.loc[ pos_max[ col ] ].bp_pos.tolist()
        pos_max[ col + 'hg19_pos' ] = bpp.loc[ pos_max[ col ] ].hg19_pos.tolist()

        pos_max[ col + '_POS_max' ] = ( first_ex_bp - ( 1 + -2*rev_strand )*pos_max[ col + 'bp_pos' ] ) - pos_max[ col + 'hg19_pos' ]


    outdf = ds_max.set_index( idx_cols )[ [ 'hg19_pos' ] ].merge( ds_sign.set_index( idx_cols )[ [ col + '_max' for col in zcols ] ],
                                                                  how = 'outer',
                                                                  left_index = True,
                                                                  right_index = True )

    outdf = outdf.merge( pos_max.set_index( idx_cols )[ [ col + '_POS_max' for col in zcols ] ],
                         how = 'outer',
                         left_index = True,
                         right_index = True ).reset_index()

    return outdf

def merge_bpp( tbl_by_var,
               bpp_max_df,
               rev_strand = False,
               idx_cols = [ 'hg19_pos' ] ):

    tbv = tbl_by_var.copy()
    bpp = bpp_max_df.copy()

    idx_cols = [ col for col in idx_cols ] + [ 'ref' + rev_strand*'_c', 'alt' + rev_strand*'_c' ]

    outdf = tbv.set_index( idx_cols ).merge( bpp.set_index( idx_cols )[ [ col for col in bpp if col.startswith( 'DS_' ) ] ],
                                             how = 'outer',
                                             left_index = True,
                                             right_index = True ).reset_index()

    return outdf

def merge_wt_branchpointer( tbl_by_var,
                            bp_wt,
                            tbl_by_var_merge_cols = [ 'hg19_pos' ] ):

    tbv = tbl_by_var.copy()

    bp = bp_wt.copy()

    bp = bp.rename( columns = { 'test_site': 'hg19_pos',
                                'branchpoint_prob': 'bp_wt_prob' } )

    tbv = tbv.set_index( tbl_by_var_merge_cols ).merge( bp.set_index( tbl_by_var_merge_cols )[ [ 'bp_wt_prob' ] ],
                                                        left_index = True,
                                                        right_index = True,
                                                        how = 'left' ).reset_index()

    return tbv

def merge_bpp_wt( tbl_by_var,
                  bpp_wt,
                  first_ex_bp,
                  merge_pos_col = 'hg19_pos',
                  rev_strand = False ):

    tbv = tbl_by_var.copy()

    bpp = bpp_wt.copy()

    bpp[ merge_pos_col ] = first_ex_bp - ( 1 + -2*rev_strand )*bpp.bp_pos

    bpp = bpp.rename( columns = { 'zsc_bps': 'zbpp_wt_bps',
                                  'zsc_ppt': 'zbpp_wt_ppt',
                                  'zsc': 'zbpp_wt' } )

    tbv = tbv.set_index( merge_pos_col ).merge( bpp.set_index( merge_pos_col )[ [ 'zbpp_wt_bps','zbpp_wt_ppt','zbpp_wt' ] ],
                                                left_index = True,
                                                right_index = True,
                                                how = 'outer' ).reset_index()

    return tbv

def get_consplice_scores( pysam_vcf,
                          liftover_bed, ):

    lift = liftover_bed.copy()

    hg38tohg19_lift = { chrom: { hg38:hg19
                                 for hg38, hg19 in zip( chrom_df.hg38_pos, chrom_df.hg19_pos ) }
                         for chrom, chrom_df in lift.groupby( 'chrom' ) }

    outd = { 'chrom': [],
             'hg38_pos': [],
             'hg19_pos': [],
             'ref': [],
             'alt': [],
             'gene_name': [],
             'consplice': [],
             'conspliceml': [] }

    for rec in pysam_vcf:

        if rec.chrom not in hg38tohg19_lift:
            if rec.chrom.startswith( 'chr' ):
                chrom = rec.chrom[ 3: ]
            else:
                chrom = 'chr' + str( rec.chrom )
        else:
            chrom = rec.chrom

        if rec.pos not in hg38tohg19_lift[ chrom ]:
            continue

        outd[ 'chrom' ].append( chrom )
        outd[ 'hg38_pos' ].append( rec.pos )
        outd[ 'hg19_pos' ].append( hg38tohg19_lift[ outd[ 'chrom' ][ -1 ] ][ outd[ 'hg38_pos' ][ -1 ] ] )
        outd[ 'ref' ].append( rec.ref )
        outd[ 'alt' ].append( rec.alts[ 0 ] )
        assert len( rec.info[ 'ConSpliceML' ][ 0 ].split( '|' ) ) == 2, 'You might be missing info - check your code!'
        outd[ 'gene_name' ].append( rec.info[ 'ConSpliceML' ][ 0 ].split( '|' )[ 0 ] )
        outd[ 'conspliceml' ].append( float( rec.info[ 'ConSpliceML' ][ 0 ].split( '|' )[ -1 ] ) )
        assert any( outd[ 'gene_name' ][ -1 ] in tup for tup in rec.info[ 'ConSplice' ] ), \
        'Your gene name is not matching for %s:%i%s>%s- code will fail!' % ( outd[ 'chrom' ][ -1 ], outd[ 'hg38_pos' ][ -1 ], outd[ 'ref' ][ -1 ], outd[ 'alt' ][ -1 ] )
        consplice_info = rec.info[ 'ConSplice' ][ 0 ].split( ',' )
        for info in rec.info[ 'ConSplice' ]:
            if info.split( '|' )[ 0 ] == outd[ 'gene_name' ][ -1 ]:
                outd[ 'consplice' ].append( float( info.split( '|' )[ -1 ] ) )

    pysam_vcf.close()

    outdf = pd.DataFrame( outd )

    return outdf

def merge_consplice( tbl_by_var,
                     consplice_df,
                     idx_cols = [ 'gene_name', 'chrom', 'hg19_pos', 'ref', 'alt' ] ):

    tbv = tbl_by_var.copy()
    conspl = consplice_df.copy()

    assert 'hg38_pos' not in tbv.columns or 'hg38_pos' in idx_cols, 'Add hg38_pos to idx_cols if its already in tbl_by_var!'

    outdf = tbv.set_index( idx_cols ).merge( conspl.set_index( idx_cols ),
                                             how = 'left',
                                             left_index = True,
                                             right_index = True ).reset_index()

    return outdf

def get_splai_precanned_scores( splai_vcf,
                                masked = True ):

    outtbl = { 'chrom': [],
               'pos': [],
               'ref': [],
               'alt': [],
               'gene_name': [],
               'DS_AG' + masked*'m': [],
               'DS_AL' + masked*'m': [],
               'DS_DG' + masked*'m': [],
               'DS_DL' + masked*'m': [],
               'DP_AG' + masked*'m': [],
               'DP_AL' + masked*'m': [],
               'DP_DG' + masked*'m': [],
               'DP_DL' + masked*'m': [],
               'DS_max' + masked*'m': [],
               'DS_max' + masked*'m' + '_type': [],
               'POS_max' + masked*'m': [],
             }

    DS_scores = [ 'DS_AG' + masked*'m', 'DS_AL' + masked*'m', 'DS_DG' + masked*'m', 'DS_DL' + masked*'m' ]
    DP_scores = [ 'DP_AG' + masked*'m', 'DP_AL' + masked*'m', 'DP_DG' + masked*'m', 'DP_DL' + masked*'m' ]

    for var in splai_vcf:

        outtbl[ 'chrom' ].append( var.chrom )
        outtbl[ 'pos' ].append( int( var.pos ) )
        outtbl[ 'ref' ].append( var.ref )
        assert len( var.alts ) == 1, 'Alt longer than one character!'
        outtbl[ 'alt' ].append( var.alts[ 0 ] )

        info_l = var.info[ 'SpliceAI'][ 0 ].split( '|' )

        outtbl[ 'gene_name' ].append( info_l[ 1 ] )

        for i,col in enumerate( DS_scores ):
            outtbl[ col ].append( info_l[ i + 2 ] )

        for i,col in enumerate( DP_scores ):
            outtbl[ col ].append( info_l[ i + 6 ] )

        #first get the maximum probability across the difference scores
        outtbl[ 'DS_max' + masked*'m' ].append( max( outtbl[ col ][ -1 ] for col in DS_scores ) )
        #then get the type of event that represents the maximum probability
        outtbl[ 'DS_max' + masked*'m' + '_type' ].append( [ col for col in DS_scores
                                                          if outtbl[ col ][ -1 ] == outtbl[ 'DS_max' + masked*'m' ][ -1 ] ][ 0 ] )
        #finally, get the location of the event associated with the highest difference score
        outtbl[ 'POS_max' + masked*'m' ].append( outtbl[ outtbl[ 'DS_max' + masked*'m' + '_type' ][ -1 ].replace( 'DS', 'DP' ) ][ -1 ] )

    outdf = pd.DataFrame( outtbl )

    return outdf
