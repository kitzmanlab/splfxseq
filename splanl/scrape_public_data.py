import pysam
import pandas as pd
import numpy as np
import subprocess as subp
import splanl.coords as cd

def create_gnomad_df( gnomad_tbx,
                      chrom,
                      coords,
                      genome, ):

    out_tbl = { 'chrom': [],
                genome + '_pos': [],
                'ref': [],
                'alt': [],
                'n_alt': [],
                'n_allele': [],
                'n_homo': [], }

    for row in gnomad_tbx.fetch( chrom ,
                                 coords[0],
                                 coords[1], ):

        if 'PASS' not in row.filter.keys() and '.' not in row.filter.keys():
            continue

        out_tbl[ 'chrom' ].append( row.chrom )
        out_tbl[ genome + '_pos' ].append( row.pos )

        out_tbl[ 'ref' ].append( row.ref )
        assert len( row.alts ) == 1, 'Unexpectedly, more than one alt - check data carefully and alter code if needed.'
        out_tbl[ 'alt' ].append( row.alts[ 0 ] )

        #keep int counts to avoid any loss of info for small allele freq
        assert len( row.info['AC'] ) == 1, 'Unexpectedly, more than one alt allele count - check data carefully and alter code if needed.'
        out_tbl[ 'n_alt' ].append( int( row.info['AC'][ 0 ] ) )
        out_tbl[ 'n_allele' ].append( int( row.info['AN'] ) )
        assert len( row.info['nhomalt'] ) == 1, 'Unexpectedly, more than one homozygotes count - check data carefully and alter code if needed.'
        out_tbl[ 'n_homo' ].append( int( row.info['nhomalt'][ 0 ] ) )

    out_df = pd.DataFrame( out_tbl )

    out_df[ 'af' ] = out_df.n_alt / out_df.n_allele

    return out_df

def merge_gnomad2_data( gnomad2_genomes,
                        gnomad2_exomes,
                        merge_cols = [ 'chrom', 'hg19_pos', 'ref', 'alt' ] ):

    gen = gnomad2_genomes.set_index( merge_cols ).copy()
    ex = gnomad2_exomes.set_index( merge_cols ).copy()

    #so no overlapping entries - just merge normally!
    if len( gen.index.intersection( ex.index ) ) == 0:
        outdf = gen.merge( ex,
                           how = 'outer',
                           left_index = True,
                           right_index = True ).reset_index()
    else:
        out_d = { col: [] for col in merge_cols + gen.columns.tolist() if col != 'af' }

        idx_all = gen.index.union( ex.index )

        for idx in idx_all:

            out_d[ 'chrom' ].append( idx[ 0 ] )
            out_d[ 'hg19_pos' ].append( idx[ 1 ] )
            out_d[ 'ref' ].append( idx[ 2 ] )
            out_d[ 'alt' ].append( idx[ 3 ] )

            if idx in gen.index and idx in ex.index:
                out_d[ 'n_alt' ].append( int( gen.loc[ idx ].n_alt ) + int( ex.loc[ idx ].n_alt ) )
                out_d[ 'n_allele' ].append( int( gen.loc[ idx ].n_allele ) + int( ex.loc[ idx ].n_allele ) )
                out_d[ 'n_homo' ].append( int( gen.loc[ idx ].n_homo ) + int( ex.loc[ idx ].n_homo ) )
            elif idx in gen.index:
                out_d[ 'n_alt' ].append( int( gen.loc[ idx ].n_alt ) )
                out_d[ 'n_allele' ].append( int( gen.loc[ idx ].n_allele ) )
                out_d[ 'n_homo' ].append( int( gen.loc[ idx ].n_homo ) )
            elif idx in ex.index:
                out_d[ 'n_alt' ].append( int( ex.loc[ idx ].n_alt ) )
                out_d[ 'n_allele' ].append( int( ex.loc[ idx ].n_allele ) )
                out_d[ 'n_homo' ].append( int( ex.loc[ idx ].n_homo ) )

        outdf = pd.DataFrame( out_d )

        outdf[ 'af' ] = outdf.n_alt / outdf.n_allele

    return outdf

def is_int( str ):
  try:
    int( str )
    return True
  except ValueError:
    return False

def merge_data_gnomad( psi_df,
                       gnomad_df,
                       indexcols = [ 'hg19_pos', 'ref', 'alt' ],
                       suffix = None ):

    psi = psi_df.set_index( indexcols ).copy()
    gnomad = gnomad_df.set_index( indexcols ).copy()

    if 'chrom' not in indexcols:
        gnomad = gnomad.drop( columns = [ 'chrom' ] )

    if suffix:
        renamed = { col: col + '_' + suffix for col in gnomad.columns if 'pos' not in col }
        gnomad = gnomad.rename( columns = renamed )

    out = psi.join( gnomad, how = 'left' ).reset_index()

    return out

def gnomad_var( vartbl ):

    tbv = vartbl.copy()

    tbv[ 'gnomad_var' ] = False

    tbv.loc[ ~( tbv.n_allele_v2.isnull() ) | ~( tbv.n_allele_v3.isnull() ), 'gnomad_var' ] = True

    return( tbv )

def download_phylop( chrom,
                     start,
                     end,
                     outdir ):

    if is_int( chrom ):
        chrom = 'chr' + str( chrom )

    assert end > start, 'End coordinates must be larger than then start coordinates'

    out = subp.run( '/home/smithcat/bigWigToBedGraph  \
                                   -chrom=%s \
                                   -start=%i \
                                   -end=%i \
                                   http://hgdownload.soe.ucsc.edu/goldenPath/hg19/phyloP100way/hg19.100way.phyloP100way.bw \
                                  %s%s_%s_%s.phyloP100way.bed' % ( chrom, start, end, outdir, chrom, str( start ), str( end ) ),
                    stdout = subp.PIPE,
                    stderr = subp.PIPE,
                    shell = True,
                )

def is_int( str ):
  try:
    int( str )
    return True
  except ValueError:
    return False

def read_phylop( phylop_file ):

    outtbl = pd.read_table( phylop_file,
                        names = [ 'chrom', 'start', 'end', 'phylop' ] )

    outtbl_nomiss = enforce_adj_starts( outtbl )

    return( outtbl_nomiss )


def enforce_adj_starts( bedgraph_tbl ):

    #if you rewrite this function to take in a 'value' column name instead of phylop
    #it would convert any bedgraph into a bed

    bg = bedgraph_tbl.sort_values( by = 'start' ).copy()

    begin = bg.iloc[ 0 ].start
    fin = bg.iloc[ bg.shape[0] - 1 ].end

    outtbl = {
               'chrom': [],
               'gdna_pos': [],
               'phylop': []
             }


    for pos in range( begin, fin ):

        #add one to make 1 based
        outtbl[ 'gdna_pos' ].append( pos + 1 )

        if pos in bg.start.tolist():

            outtbl[ 'chrom' ].append( bg.loc[ bg.start == pos ].chrom.values[0] )
            outtbl[ 'phylop' ].append( bg.loc[ bg.start == pos ].phylop.values[0] )

        #position isn't in dataframe - lets use previous positions values instead
        else:

            #first make sure the position is still within the previous interval
            assert bg.loc[ bg.start < pos ].iloc[ -1 ].end > pos, \
            ' Start coordinate %i is not contained in any interval' % pos

            #since everything is position sorted, we're grabbing the lowest row of all rows that start before
            #this given position
            outtbl[ 'chrom' ].append( bg.loc[ bg.start < pos ].iloc[ -1 ].chrom )
            outtbl[ 'phylop' ].append( bg.loc[ bg.start < pos ].iloc[ -1 ].phylop )

    outtbl = pd.DataFrame( outtbl )

    return outtbl

def merge_data_phylop( byvartbl,
                       phylop_tbl,
                       index_cols = [ 'gdna_pos' ] ):

    tbv = byvartbl.copy()
    p = phylop_tbl.copy()

    merge_tbl = merge_data_gnomad( tbv,
                                      p,
                                      indexcols = index_cols )

    return merge_tbl

def create_clinvar_df( clinvar_tbx,
                       chrom,
                       coords,
                       genome ):

    if genome != 'hg19':
        print( 'Processing genome other than hg19 - MPSA driver script assumes merge on hg19_pos column!' )
        print( 'You will need to do a liftover to merge with MPSA data - could download hg19 VCF instead!' )

    all_info_keys = list( set( [ key for rec in clinvar_tbx.fetch( chrom ,
                                                               coords[0],
                                                               coords[1], )
                                 for key in rec.info.keys() ] ) )

    outd = { col: [] for col in [ 'chrom', genome + '_pos', 'ref', 'alt' ] + all_info_keys }

    for rec in clinvar_tbx.fetch( chrom ,
                              coords[0],
                              coords[1], ):

        ref = rec.ref

        if rec.alts is None or len( rec.alts ) > 1:
            continue

        alt = rec.alts[ 0 ]

        #only get SNVs
        if len( ref ) != 1 or len( alt ) != 1:
            continue

        info_keys = list( rec.info.keys() )

        if 'CLNSIG' not in info_keys or rec.info[ 'CLNSIG' ] == '':
            continue

        outd[ 'chrom' ].append( int( rec.chrom ) )
        outd[ 'hg19_pos' ].append( int( rec.pos ) )
        outd[ 'ref' ].append( ref )
        outd[ 'alt' ].append( alt )

        for key in all_info_keys:

            if key in info_keys:

                if isinstance( rec.info[ key ], tuple ) and len( rec.info[ key ] ) == 1:
                    outd[ key ].append( rec.info[ key ][ 0 ] )
                else:
                    outd[ key ].append( rec.info[ key ] )

            else:
                outd[ key ].append( '' )

    outdf = pd.DataFrame( outd ).rename( columns = { 'CLNSIG': 'clinvar_interp' } )

    if 'GENEINFO' in outdf:
        outdf[ 'clinvar_gene' ] = outdf.GENEINFO.apply( lambda x: x.split( ':' )[ 0 ] )

    clinvar_possible = [ 'Uncertain_significance', 'Likely_benign',
                         'Conflicting_interpretations_of_pathogenicity',
                         'Benign/Likely_benign', 'Benign', 'Pathogenic',
                         'Pathogenic/Likely_pathogenic', 'Likely_pathogenic',
                         'not_provided' ]

    if len( set( outdf.clinvar_interp.unique() ).difference( set( clinvar_possible ) ) ) > 0:
        print( 'Some ClinVar interpretations not getting processed into the clinvar_abbrev column properly!', set( outdf.clinvar_interp.unique() ).difference( set( clinvar_possible ) ) )

    outdf[ 'clinvar_abbrev' ] = [ 'LBB' if interp == 'Likely_benign' or interp == 'Benign/Likely_benign' or interp == 'Benign'
                                  else 'LPP' if interp == 'Pathogenic' or interp == 'Pathogenic/Likely_pathogenic' or interp == 'Likely_pathogenic'
                                  else 'VUS' if interp == 'Uncertain_significance' or interp == 'Conflicting_interpretations_of_pathogenicity'
                                  else ''
                                  for interp in outdf.clinvar_interp ]

    return outdf

def merge_data_clinvar( psi_df,
                       clinvar_df,
                       indexcols = [ 'hg19_pos', 'ref', 'alt' ],
                       keep_cols = [ 'clinvar_interp', 'clinvar_abbrev', 'clinvar_gene' ] ):

    psi = psi_df.set_index( indexcols ).copy()
    clinvar = clinvar_df.set_index( indexcols ).copy()

    if 'chrom' not in indexcols and 'chrom' in clinvar:
        clinvar = clinvar.drop( columns = [ 'chrom' ] )

    out = psi.join( clinvar[ keep_cols ], how = 'left' ).reset_index()

    out[ 'clinvar' ] = out.clinvar_interp.notnull()

    return out
