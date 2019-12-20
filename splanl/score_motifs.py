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

def ExtractKmers(seq,k):
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
    assert len(seq)==2*k-1
    return([seq[i:i+k] for i in range(k)])

def ScoreMotifs(fastafile,
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
    seq=ImportFastA(fastafile)
    motifDict=motifdf.T.to_dict('records')[0]
    vardf=vardf.set_index(['pos','alt'])

    #score kmers for wt and mutant - assigns 0 if not in motifDict
    #numbering adjusts for 0-based and 1-based coordinates
    wt={ p:[motifDict[kmer] if kmer in motifDict else 0
        for kmer in ExtractKmers(seq[p-6:p+5],k)]
        for p in vardf.index.unique(level='pos') }

    #same as above but substitute the mutated base
    mut={ (p,a):[ motifDict[kmer] if kmer in motifDict else 0
        for kmer in ExtractKmers( seq[ p-6:p-1 ]+a+seq[ p:p+5 ],k ) ]
         for p,a in vardf.index }

    #compute difference: mut-wt and take the mean
    change={(p,a):np.mean([i-j for i,j in zip(mut[(p,a)],wt[p])])
            for p,a in vardf.index}

    changedf = pd.DataFrame(list(change.values()), index=pd.MultiIndex.from_tuples(change), columns=[col_stem+'_chgMEAN'])

    vardf = pd.concat([vardf,changedf],axis=1)
    vardf = vardf.reset_index(level=['pos','alt'])
    firstCol=['chrom','pos','ref','alt']
    vardf=vardf[ firstCol + [c for c in vardf if c not in firstCol ] ]

    return (vardf)
