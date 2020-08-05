import pandas as pd
import numpy as np

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

def score_motifs(refseq,
                motifdf,
                vardf,
                k,
                col_stem):
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
    wt = { p: np.mean( [ mdf.loc[kmer].values[ 0 ] if kmer in mdf.index else 0
          for kmer in ExtractKmers( refseq[ p - k: p + ( k - 1 ) ], k ) ] )
          for p in vdf.index.unique( level = 'pos' ) }

    #same as above but substitute the mutated base
    mut = { ( p, a ): np.mean( [ mdf.loc[kmer].values[ 0 ] if kmer in mdf.index else 0
            for kmer in ExtractKmers( refseq[ p - k : p - 1 ] + a + refseq[ p : p + ( k - 1 ) ], k ) ] )
            for p,a in vdf.index }

    #compute difference: mut-wt and take the mean
    change = { ( p, a ): mut[ ( p, a ) ] - wt[ p ] for p,a in vdf.index }

    changedf = pd.DataFrame( list( change.values() ), index=pd.MultiIndex.from_tuples( change ), columns=[ col_stem+'_chgMEAN' ] )

    outdf = pd.concat( [ vdf, changedf ],axis=1 )
    outdf = outdf.reset_index( level = [ 'pos', 'alt' ] )
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
