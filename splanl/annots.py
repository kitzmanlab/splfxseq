import pandas as pd

def gene_specific_gtf( gtf_file,
                       gene_name,
                       transcript_name ):

    gtf = gtf_file.copy()

    gtf[ 'gene_id' ] = [ a.split( ' ' )[ 1 ].replace( '"', '' )
                         for att in gtf.attribute
                         for a in att.split( '; ' )
                         if a.split( ' ' )[ 0 ] == 'gene_id' ]

    gene = gtf.loc[ gtf.gene_id == gene_name ].copy()

    gene.loc[ gtf.feature != 'gene', 'transcript_id' ] = [ a.split( ' ' )[ 1 ].replace( '"', '' )
                                                              for att in gene.attribute
                                                              for a in att.split( '; ' )
                                                              if a.split( ' ' )[ 0 ] == 'transcript_id' ]

    gene.transcript_id = gene.transcript_id.fillna( '' )

    assert transcript_name in gene.transcript_id.unique().tolist(), \
    'Transcript %s not in table - only transcript %s' % ( transcript_name, gene.loc[ gene.transcript_id != '' ].transcript_id.unique().tolist() )

    #if there's more than one transcript for the gene just grab those rows
    if len( gene.loc[ gene.transcript_id != '' ].transcript_id.unique().tolist() ) > 1:
        trans = gene.loc[ ( gene.transcript_id == '' ) | ( gene.transcript_id == transcript_name ) ].copy()

    else:
        trans = gene.copy()

    trans = trans.drop( columns = [ 'gene_id', 'transcript_id' ] )

    return trans

def gtf_to_splai_annot( gtf_file ):

    gtf = gtf_file.copy()

    outtbl = { '#NAME': [],
               'CHROM': [],
               'STRAND': [],
               'TX_START': [],
               'TX_END': [],
               'EXON_START': [],
               'EXON_END': [],
               'TRANSCRIPT_ID': [] }

    gtf[ 'gene_id' ] = [ a.split( ' ' )[ 1 ].replace( '"', '' )
                         for att in gtf.attribute
                         for a in att.split( '; ' )
                         if a.split( ' ' )[ 0 ] == 'gene_id' ]

    gtf.loc[ gtf.feature != 'gene', 'transcript_id' ] = [ a.split( ' ' )[ 1 ].replace( '"', '' )
                                                              for att in gtf.attribute
                                                              for a in att.split( '; ' )
                                                              if a.split( ' ' )[ 0 ] == 'transcript_id' ]

    gtf[ 'transcript_id' ] = gtf[ 'transcript_id' ].fillna( '' )

    gtf[ 'start' ] -= 1

    for tran in list( set( gtf.transcript_id ) ):

        #gene rows are missing the transcript id - we want to skip these obviously
        if not tran:
            continue

        trans = gtf.loc[ gtf.transcript_id == tran ].copy()

        assert len( trans.gene_id.unique() ) == 1, 'Transcript %s is mapping to more than one gene' % tran

        outtbl[ '#NAME' ].append( trans.gene_id.unique()[ 0 ] )
        #this is also removing the 'chr' at the begginning of each chromosome
        outtbl[ 'CHROM' ].append( trans.chrom.unique()[ 0 ][ 3: ] )
        outtbl[ 'STRAND' ].append( trans.strand.unique()[ 0 ] )

        tx_tbl = trans.loc[ ( trans.feature == 'transcript' ) ].copy()
        outtbl[ 'TX_START' ].append( int( tx_tbl.start ) )
        outtbl[ 'TX_END' ].append( int( tx_tbl.end ) )

        ex_tbl = trans.loc[ trans.feature == 'exon' ].copy()
        outtbl[ 'EXON_START' ].append( ','.join( ex_tbl.start.apply( str ).tolist() ) )
        outtbl[ 'EXON_END' ].append( ','.join( ex_tbl.end.apply( str ).tolist() ) )

        outtbl[ 'TRANSCRIPT_ID' ].append( tran )

    outdf = pd.DataFrame( outtbl )

    return outdf
