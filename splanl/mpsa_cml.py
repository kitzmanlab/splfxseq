import pybedtools as pbt
import pysam
import pandas as pd
import numpy as np
import argparse
import os
import sys
import time
import splanl.plots as sp
import splanl.coords as cd
import splanl.post_processing as pp
import splanl.junction_scorer as jn
import splanl.merge_bcs as mbcs
import splanl.score_motifs as sm
import splanl.scrape_public_data as spd

def main():

    parser = argparse.ArgumentParser( description = 'Collect isoforms, merge BCs, and create variant isoform use percentages!',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument( '-rn', '--refseqname',
                         help = 'Name of sequence in vector fasta file (str) ie jkp618' )
    parser.add_argument( '-crn', '--chrom_refseqname',
                         help = 'Name of sequence in genomic fasta file (str) ie chr11' )
    parser.add_argument( '-sn', '--samp_name_col',
                         default = 'libname',
                         help = 'Sample name column in reference table (str)' )
    parser.add_argument( '-sbn', '--samp_name_bam_col',
                         default = 'libname',
                         help = 'Sample name column in sample to bam table (str)' )
    parser.add_argument( '-s', '--splicing_err_tol',
                         default = 3,
                         type = int,
                         help = 'Buffer of error allowed at splice junctions (int)' )
    parser.add_argument( '-i', '--indel_tol',
                         default = 20,
                         type = int,
                         help = 'Collapses indels within reads (int)' )
    parser.add_argument( '-rl', '--read_length',
                         default = 150,
                         type = int,
                         help = 'Length of the reads (int)' )
    parser.add_argument( '-bcl', '--max_bc_length',
                         default = 30,
                         type = int,
                         help = 'Maximum barcode length - longer barcodes are filtered out (int)' )
    parser.add_argument( '-bct', '--bc_read_cut_thresh',
                         default = 95,
                         type = int,
                         help = 'Cumulative read count percentile from waterfall plot used to discard barcodes with few reads (int)' )
    parser.add_argument( '-wft', '--waterfall_thresh',
                         default = 90,
                         type = int,
                         help = 'Cumulative read count percentile to explore on waterfall plot - default settings with show cutoffs at 75th, 90th, and 95th percentiles (int)' )
    parser.add_argument( '-b', '--barcode_tag',
                         default = 'BC',
                         help = 'Barcode tag in bamfiles (str)' )
    parser.add_argument( '-ubc', '--unmapped_bc_split',
                         default = '_BC=',
                         help = 'Split character to extract barcode from unmapped reads in STAR  (str)' )
    parser.add_argument( '-u', '--unmapped_bam_file_dir',
                         help = 'Directory of unmapped bams - Required to count unmapped reads with STAR' )
    parser.add_argument( '-g2', '--gnomad2_file',
                         help = 'Directory + file of gnomAD2 processed by scrape_public_data module' )
    parser.add_argument( '-g3', '--gnomad3_file',
                         help = 'Directory + file of gnomAD3 processed by scrape_public_data module - MUST HAVE hg19 coords as hg19_pos column!' )
    parser.add_argument( '-cl', '--clinvar_file',
                         help = 'Directory + file of ClinVar processed by scrape_public_data module - MUST HAVE hg19 coords as hg19_pos column!' )
    parser.add_argument( '-spl', '--spliceai_file',
                         help = 'Directory + file of SpliceAI file computed by custom_spliceai_scores module' )
    parser.add_argument( '-splcol', '--merge_splai_basic_cols_only',
                         action = 'store_false',
                         help = 'Merge the SpliceAI DS columns only - turn on if the SpliceAI file has additional columns to drop')
    parser.add_argument( '-me', '--maxentscan',
                         action = 'store_true',
                         help = 'Include MaxEntScan scores?' )
    parser.add_argument( '-v', '--verbose',
                         action = 'store_true',
                         help = 'Increase number of printouts!' )
    parser.add_argument( 'vector_fafile',
                         help = 'Fasta file of plasmid sequence: dir + filename (str)' )
    parser.add_argument( 'genomic_fafile',
                         help = 'Fasta file of genomic sequence: dir + filename (str)' )
    parser.add_argument( 'subasm',
                         help = 'Subassembly file: dir + filename (str)' )
    parser.add_argument( 'exonbed',
                         help = 'Exon bedfile: dir + filename (str) - used to create list of known isoforms' )
    parser.add_argument( 'reference_table',
                         help = 'File containing sample names by exon (needs column named exon!!!): dir + filename (str)' )
    parser.add_argument( 'exon_name',
                         help = 'Exon name used in exon columns of the reference table (str)' )
    parser.add_argument( 'samp_to_bam_table',
                         help = 'File containing bam files by sample names: dir + filename (str)' )
    parser.add_argument( 'bam_column',
                         help = 'Column to collect bam file locations in sample to bam table (str)' )
    parser.add_argument( 'bam_dir',
                         help = 'Directory of bam files (str)' )
    parser.add_argument( 'exon_number',
                         type = int,
                         help = 'Exon number (int)' )
    parser.add_argument( 'ds_const_start',
                         type = int,
                         help = 'Vector start position of downstream constant exon - 0 based (likely 649) (int)' )
    parser.add_argument( 'ds_const_end',
                         type = int,
                         help = 'Vector end position of downstream constant exon (int)' )
    parser.add_argument( 'us_const_start',
                         type = int,
                         help = 'Vector start position of upstream constant exon - 0 based (int)' )
    parser.add_argument( 'us_const_end',
                         type = int,
                         help = 'Vector end position of upstream constant exon (int)' )
    parser.add_argument( 'cloned_vstart',
                         type = int,
                         help = 'Vector start position of cloned region - 1 based (int)' )
    parser.add_argument( 'cloned_vend',
                         type = int,
                         help = 'Vector end position of cloned region (int)' )
    parser.add_argument( 'exon_vstart',
                         type = int,
                         help = 'Vector start position of exon - 1 based (int)' )
    parser.add_argument( 'exon_vend',
                         type = int,
                         help = 'Vector end position of exon (int)' )
    parser.add_argument( 'exon_hgvs_start',
                         type = int,
                         help = 'cDNA start position of exon (int)' )
    parser.add_argument( 'exon_hgvs_end',
                         type = int,
                         help = 'cDNA end position of exon (int)' )
    parser.add_argument( 'exon_hg19_start',
                         type = int,
                         help = 'hg19 start position of exon (int) - 1 based' )
    parser.add_argument( 'exon_hg19_end',
                         type = int,
                         help = 'hg19 end position of exon (int)' )
    parser.add_argument( 'chrom',
                         help = 'Chromosome - can be with or without chr (str)' )
    parser.add_argument( 'strand',
                         choices = [ '+', '-' ],
                         help = 'Strand - either + or - (str)' )
    parser.add_argument( 'min_for_matches',
                         type = int,
                         help = 'Minimum matches required for forward reads to not be "softclipped" (int)' )
    parser.add_argument( 'min_rev_matches',
                         type = int,
                         help = 'Minimum matches required for reverse reads to not be "softclipped" (int)' )
    parser.add_argument( 'softclip_out_dir',
                         help = 'Directory to store softclipped reads bams (str)' )
    parser.add_argument( 'plots_out_dir',
                         help = 'Directory to store plots (str)' )
    parser.add_argument( 'data_out_dir',
                         help = 'Directory to store datasets (str)' )
    args = parser.parse_args()
    config = vars(args)

    T0 = time.time()

    if not config[ 'softclip_out_dir' ].endswith( '/' ):
        config[ 'softclip_out_dir' ] = config[ 'softclip_out_dir' ] + '/'
    if not config[ 'plots_out_dir' ].endswith( '/' ):
        config[ 'plots_out_dir' ] = config[ 'plots_out_dir' ] + '/'
    if not config[ 'data_out_dir' ].endswith( '/' ):
        config[ 'data_out_dir' ] = config[ 'data_out_dir' ] + '/'
    if not config[ 'bam_dir' ].endswith( '/' ):
        config[ 'bam_dir' ] = config[ 'bam_dir' ] + '/'
    if config[ 'unmapped_bam_file_dir' ] and not config[ 'unmapped_bam_file_dir' ].endswith( '/' ):
        config[ 'unmapped_bam_file_dir' ] = config[ 'unmapped_bam_file_dir' ] + '/'

    if config[ 'gnomad2_file' ]:
        assert config[ 'gnomad3_file' ], 'This script assumes you entered both versions of gnomAD or neither! Add gnomAD3!'
        gnomad2 = pd.read_table( config[ 'gnomad2_file' ] )
    if config[ 'gnomad3_file' ]:
        assert config[ 'gnomad2_file' ], 'This script assumes you entered both versions of gnomAD or neither! Add gnomAD2!'
        gnomad3 = pd.read_table( config[ 'gnomad3_file' ] )
        assert 'hg19_pos' in gnomad3, 'Your gnomAD3 dataframe must have hg19_pos as one of the columns!'

    if config[ 'clinvar_file' ]:
        clinvar = pd.read_table( config[ 'clinvar_file' ] )
        assert 'hg19_pos' in clinvar, 'Your ClinVar dataframe must have hg19_pos as one of the columns!'

    if config[ 'spliceai_file' ]:
        splai = pd.read_table( config[ 'spliceai_file' ],
                               dtype = { 'chrom': object } )
        if config[ 'merge_splai_basic_cols_only' ]:
            splai = splai[ [ 'chrom', 'pos', 'ref', 'alt' ] + [ col for col in splai if col.startswith( 'DS_' ) or col.startswith( 'DP_' ) ] ].copy()

    chrom = config[ 'chrom' ] if not config[ 'chrom' ].startswith( 'chr' ) else config[ 'chrom' ][ 3: ]
    chrom_chr = config[ 'chrom' ] if config[ 'chrom' ].startswith( 'chr' ) else 'chr' + config[ 'chrom' ]

    date_string = time.strftime( '%Y-%m%d', time.localtime() )

    refseq_d = pp.get_refseq( config[ 'vector_fafile' ] )

    if config[ 'refseqname' ]:
        if config[ 'refseqname' ] in refseq_d:
            refseq = refseq_d[ config[ 'refseqname' ] ]
        else:
            if config[ 'verbose' ]:
                print( 'Refseqname (-rn) %s not in fasta - trying %s' % ( config[ 'refseqname' ], chrom_chr  ) )
            with open( config[ 'data_out_dir' ] + 'mpsa_cml_log.' + date_string + '.txt', 'a' ) as f:
                f.write( 'Refseqname (-rn) %s not in fasta - trying %s\n' % ( config[ 'refseqname' ], chrom_chr  ) )
            assert chrom_chr in refseq_d, 'Refseqname (-rn) %s and chromosome %s not in provided fasta file - only %s!' % ( config[ 'refseqname' ], chrom_chr, ', '.join( list( refseq_d.keys() ) ) )
            refseq = refseq_d[ chrom_chr ]
    elif len( refseq_d ) == 1:
        if config[ 'verbose' ]:
            print( 'No refseqname (-rn) provided - using %s' % list( refseq_d.keys() )[ 0 ] )
        with open( config[ 'data_out_dir' ] + 'mpsa_cml_log.' + date_string + '.txt', 'a' ) as f:
            f.write( 'No refseqname (-rn) provided - using %s\n' % list( refseq_d.keys() )[ 0 ] )
        refseq = list( refseq_d.values() )[ 0 ]
    else:
        if config[ 'verbose' ]:
            print( 'No refseqname (-rn) provided - trying %s' % chrom_chr  )
        with open( config[ 'data_out_dir' ] + 'mpsa_cml_log.' + date_string + '.txt', 'a' ) as f:
            f.write( 'No refseqname (-rn) provided - trying %s\n' % chrom_chr )
        assert chrom_chr in refseq_d, 'Chromosome %s not in provided fasta file - only %s! Try specifying refseqname (-rn)...' % ( chrom_chr, ', '.join( list( refseq_d.keys() ) ) )
        refseq = refseq_d[ chrom_chr ]

    refseq_d = pp.get_refseq( config[ 'genomic_fafile' ] )

    if config[ 'chrom_refseqname' ]:
        assert config[ 'chrom_refseqname' ] in refseq_d, 'Chromosome refseqname (-crn) not in provided fasta file!'
        chrom_refseq = refseq_d[ config[ 'chrom_refseqname' ] ]
    elif len( refseq_d ) == 1:
        if config[ 'verbose' ]:
            print( 'No refseqname (-crn) provided - using ' + list( refseq_d.keys() )[ 0 ] )
        with open( config[ 'data_out_dir' ] + 'mpsa_cml_log.' + date_string + '.txt', 'a' ) as f:
            f.write( 'No refseqname (-crn) provided - using ' + list( refseq_d.keys() )[ 0 ] + '\n' )
        chrom_refseq = list( refseq_d.values() )[ 0 ]
    else:
        sys.exit( 'Multiple sequences in chromosome fasta file! Please provide the name of the sequence in the fasta file using the -crn argument' )

    #garbage collect the memory in case a lot of sequences provided
    refseq_d = {}

    ref_tbl = pd.read_table( config[ 'reference_table' ] )
    samp_names = ref_tbl.loc[ ref_tbl.exon == config[ 'exon_name' ] ][ config[ 'samp_name_col' ] ].tolist()

    samp_tbl = pd.read_table( config[ 'samp_to_bam_table' ] ).set_index( config[ 'samp_name_bam_col' ] )

    bams = { samp: samp_tbl.loc[ samp ][ config[ 'bam_column' ] ] for samp in samp_names }

    msamp_fn = { samp: bam if config[ 'bam_dir' ] in bam else config[ 'bam_dir' ] + bam
                 for samp,bam in bams.items() }

    msamp_rnabam = { samp: pysam.AlignmentFile( msamp_fn[ samp ], 'rb' ) for samp in msamp_fn }

    cnst_exons = [ ( config[ 'ds_const_start' ], config[ 'ds_const_end' ] ),
                   ( config[ 'us_const_start' ], config[ 'us_const_end' ] ) ]

    t0 = time.time()

    isos_dfs = { samp: jn.get_all_isoforms_pe( msamp_rnabam[ samp ],
                                            cnst_exons,
                                            spl_tol = config[ 'splicing_err_tol' ],
                                            indel_tol = config[ 'indel_tol' ],
                                            min_matches_for = config[ 'min_for_matches' ],
                                            min_matches_rev = config[ 'min_rev_matches' ],
                                            softclip_outbam = config[ 'softclip_out_dir' ] + '/%s_iso_collect_softclip_%i_%i.bam' % ( samp,
                                                                                                                                     config[ 'min_for_matches' ],
                                                                                                                                     config[ 'min_rev_matches' ] ) )
                for samp in msamp_rnabam }

    t1 = time.time()

    if config[ 'verbose' ]:
        print( 'Collected isoforms across samples in %.2f minutes!' % ( ( t1 - t0 ) / 60 ) )
    with open( config[ 'data_out_dir' ] + 'mpsa_cml_log.' + date_string + '.txt', 'a' ) as f:
        f.write( 'Collected isoforms across samples in %.2f minutes!\n' % ( ( t1 - t0 ) / 60 ) )

    if config[ 'unmapped_bam_file_dir' ]:

        mate1_files = [ f for f in os.listdir( config[ 'unmapped_bam_file_dir' ] ) if f.endswith( '.mate1') ]
        mate2_files = [ f for f in os.listdir( config[ 'unmapped_bam_file_dir' ] ) if f.endswith( '.mate2' ) ]

        unmap_mate1_fn = { samp: pysam.FastxFile( config[ 'unmapped_bam_file_dir' ] + f )
                           for samp in samp_names
                           for f in mate1_files if samp in f }

        unmap_mate2_fn = { samp: pysam.FastxFile( config[ 'unmapped_bam_file_dir' ] + f )
                           for samp in samp_names
                           for f in mate2_files if samp in f }

        isos_dfs = { samp: jn.append_unmapped( isos_dfs[ samp ],
                                            unmap_mate1_fn[ samp ],
                                            unmap_mate2_fn[ samp ] )
                     for samp in isos_dfs }

    isogrp_df = jn.number_and_merge_isoforms( isos_dfs )

    satbl = pd.read_table( config[ 'subasm' ],
                            sep='\t' )
    satbl = satbl.set_index( 'readgroupid' )

    exonbed = pbt.BedTool( config[ 'exonbed' ] )

    isos = jn.make_junction_graph( exonbed )

    if config[ 'verbose' ]:
        print( 'Known isoforms from provided exon bedfile: ')
        print( isos )
    with open( config[ 'data_out_dir' ] + 'mpsa_cml_log.' + date_string + '.txt', 'a' ) as f:
        f.write( 'Known isoforms from provided exon bedfile:\n' )
        print( isos, file = f )
        f.write( '\n' )

    unique_jns = list( { jn for grp,jn_tups in isos.items() for jn_tup in jn_tups for jn in jn_tup
                         if cnst_exons[ 0 ][ 1 ] < jn < cnst_exons[ 1 ][ 0 ] } )

    msamp_rnabam = { samp : pysam.AlignmentFile( msamp_fn[ samp ], 'rb' ) for samp in msamp_fn }

    canonical_isos = list( set( [ iso[ 1:-1 ] for iso in isos.values() if len( iso ) > 2 ] \
                                + [ (), ( ( config[ 'exon_vstart' ], config[ 'exon_vend' ] ), ) ] ) )

    print( 'Getting isoform statistics - this takes some time...' )

    t0 = time.time()

    iso_df_stats = jn.summarize_isos_by_var_bc_pe( msamp_rnabam,
                                                cnst_exons,
                                                satbl,
                                                isogrp_df,
                                                unique_jns,
                                                canonical_isos,
                                                config[ 'data_out_dir' ],
                                                spl_tol = config[ 'splicing_err_tol' ],
                                                indel_tol = config[ 'indel_tol' ],
                                                min_matches_for = config[ 'min_for_matches' ],
                                                min_matches_rev = config[ 'min_rev_matches' ],
                                                bc_tag = config[ 'barcode_tag' ],
                                                verbose = config[ 'verbose' ]
                                          )

    t1 = time.time()

    if config[ 'verbose' ]:
        print( 'Collected isoform statistics in %.2f minutes!' % ( ( t1 - t0 ) / 60 ) )
    with open( config[ 'data_out_dir' ] + 'mpsa_cml_log.' + date_string + '.txt', 'a' ) as f:
        f.write( 'Collected isoform statistics in %.2f minutes!\n' % ( ( t1 - t0 ) / 60 ) )

    iso_df_stats.reset_index().to_csv( config[ 'data_out_dir' ] + config[ 'exon_name' ] + '_isoforms.' + date_string + '.txt',
                                        sep = '\t',
                                        index = False )

    iso_rows = [ 'secondary', 'unmapped', 'unpaired', 'bad_starts', 'bad_ends', 'soft_clipped' ] + canonical_isos

    reads_df = sp.plot_iso_stats( iso_df_stats,
                               iso_rows = iso_rows,
                               cml = True,
                               plot_out = config[ 'plots_out_dir' ], )

    reads_df.to_csv( config[ 'data_out_dir' ] + config[ 'exon_name' ] + '_read_counts.' + date_string + '.txt',
                     sep = '\t',
                     index = False )

    if config[ 'verbose' ]:
        print( 'Total isoforms collected: %i' % iso_df_stats.shape[ 0 ] )
        print( 'Total isoforms passing filters in at least one sample: %i' % iso_df_stats.query( 'total_passfilt > 0' ).shape[ 0 ] )
        print( 'Percentage of reads in isoforms NOT passing filters: %.2f' % ( iso_df_stats.query( 'total_passfilt == 0' ).total_read_count.sum() / iso_df_stats.total_read_count.sum() )*100 )
        print( 'Percentage of barcodes in isoforms NOT passing filters: %.2f' % ( iso_df_stats.query( 'total_passfilt == 0' ).total_num_bcs.sum() / iso_df_stats.total_num_bcs.sum() )*100 )
        print( 'Total isoforms passing filters in ALL samples: %i' % len( iso_df_stats.loc[ iso_df_stats.total_passfilt == len( msamp_fn ) ] ) )
        print( 'Percentage of reads in isoforms passing filters in ALL samples %.2f' % ( iso_df_stats.loc[ iso_df_stats.total_passfilt == len( msamp_fn ) ].total_read_count.sum() / iso_df_stats.total_read_count.sum() )*100 )
        print( 'Percentage of barcodes in isoforms passing filters in ALL samples %.2f' % ( iso_df_stats.loc[ iso_df_stats.total_passfilt == len( msamp_fn ) ].total_num_bcs.sum() / iso_df_stats.total_num_bcs.sum() )*100 )

    with open( config[ 'data_out_dir' ] + 'mpsa_cml_log.' + date_string + '.txt', 'a' ) as f:
        f.write( 'Total isoforms collected: %i\n' % iso_df_stats.shape[ 0 ] )
        f.write( 'Total isoforms passing filters in at least one sample: %i\n' % iso_df_stats.query( 'total_passfilt > 0' ).shape[ 0 ] )
        f.write( 'Percentage of reads in isoforms NOT passing filters: %.2f\n' % ( ( iso_df_stats.query( 'total_passfilt == 0' ).total_read_count.sum() / iso_df_stats.total_read_count.sum() )*100 ) )
        f.write( 'Percentage of barcodes in isoforms NOT passing filters: %.2f\n' % ( ( iso_df_stats.query( 'total_passfilt == 0' ).total_num_bcs.sum() / iso_df_stats.total_num_bcs.sum() )*100 ) )
        f.write( 'Total isoforms passing filters in ALL samples: %i\n' % len( iso_df_stats.loc[ iso_df_stats.total_passfilt == len( msamp_fn ) ] ) )
        f.write( 'Percentage of reads in isoforms passing filters in ALL samples %.2f\n' % ( ( iso_df_stats.loc[ iso_df_stats.total_passfilt == len( msamp_fn ) ].total_read_count.sum() / iso_df_stats.total_read_count.sum() )*100 ) )
        f.write( 'Percentage of barcodes in isoforms passing filters in ALL samples %.2f\n' % ( ( iso_df_stats.loc[ iso_df_stats.total_passfilt == len( msamp_fn ) ].total_num_bcs.sum() / iso_df_stats.total_num_bcs.sum() )*100 ) )

    isogrpdict = { samp: jn.create_iso_dict_no_cnst( iso_df_stats.query( samp + '_filter!=0' ),
                                                     string_isos = False )
                   for samp in msamp_fn }

    incl_iso = ( ( config[ 'exon_vstart' ], config[ 'exon_vend' ] ), )
    isonamedict = { 'CAN_ISO' + str( i ): [ iso[ 0 ] ]
                    for i,iso in enumerate( canonical_isos ) if iso != () and iso != incl_iso }
    isonamedict[ 'SKIP' ] = [ ]
    isonamedict[ 'INCL' ] = [ incl_iso[ 0 ] ]

    if config[ 'verbose' ]:
        print( 'Named isoforms are: ', isonamedict )
    with open( config[ 'data_out_dir' ] + 'mpsa_cml_log.' + date_string + '.txt', 'a' ) as f:
        f.write( 'Named isoforms are:\n' )
        print( isonamedict, file = f )
        f.write( '\n' )

    namedisogrps = { samp: jn.create_named_isogrps( isogrpdict[ samp ],
                                                    isonamedict,
                                                    [],
                                                    cnst_exons[ 1 ][ 1 ] - cnst_exons[ 1 ][ 0 ],
                                                    config[ 'read_length' ],
                                                    config[ 'splicing_err_tol' ] )
                     for samp in msamp_fn }

    msamp_rnabam = { samp : pysam.AlignmentFile( msamp_fn[ samp ], 'rb' ) for samp in msamp_fn }

    if config[ 'unmapped_bam_file_dir' ]:

        unmap_mate1_fn = { samp: pysam.FastxFile( config[ 'unmapped_bam_file_dir' ] + f )
                           for samp in samp_names
                           for f in mate1_files if samp in f }

    else:

        unmap_mate1_fn = { samp: None for samp in samp_names }

    waterfall_thresh = list( set( [ 75, 95, config[ 'bc_read_cut_thresh' ], config[ 'waterfall_thresh' ] ] ) )
    waterfall_thresh.sort()

    print( 'Merging barcodes now...' )

    msamp_bcs_processed = { samp: mbcs.process_bcs_wrapper( samp,
                                                       msamp_rnabam[ samp ],
                                                       isogrpdict[ samp ],
                                                       namedisogrps[ samp ],
                                                       cnst_exons,
                                                       satbl,
                                                       spl_tol = config[ 'splicing_err_tol' ],
                                                       indel_tol = config[ 'indel_tol' ],
                                                       min_matches_for = config[ 'min_for_matches' ],
                                                       min_matches_rev = config[ 'min_rev_matches' ],
                                                       bc_tag = config[ 'barcode_tag' ],
                                                       max_bc_len = config[ 'max_bc_length' ],
                                                       unmap_bc_split_char = config[ 'unmapped_bc_split' ],
                                                       unmapped_pysam = unmap_mate1_fn[ samp ],
                                                       waterfall_thresh = waterfall_thresh,
                                                       bc_read_cut_thresh = config[ 'bc_read_cut_thresh' ],
                                                       verbose = config[ 'verbose' ],
                                                       cml = True,
                                                       plot_out =  config[ 'plots_out_dir' ] + samp + '_waterfall.pdf'  )
                           for samp in msamp_rnabam }

    read_cut_unfilt = mbcs.create_read_count_df( msamp_bcs_processed,
                                            waterfall_thresh )

    read_cut_unfilt.to_csv( config[ 'data_out_dir' ] + config[ 'exon_name' ] + '.n_bcs_bysamp.' + date_string + '.txt',
                            sep = '\t',
                            index = False )

    cut_dict = { samp: cut for samp,cut in zip( read_cut_unfilt[ 'sample' ],
                                                read_cut_unfilt[ str( config[ 'bc_read_cut_thresh' ] ) + '_y' ] ) }

    sp.plot_waterfall_bysamp( read_cut_unfilt,
                           cutoffs = tuple( waterfall_thresh ),
                           cml = True,
                           savefig = config[ 'plots_out_dir' ] + config[ 'exon_name' ] + '_waterfall_cutoffs_bysamp.pdf',
                             )

    for samp in msamp_bcs_processed:

        msamp_bcs_processed[ samp ][ 'msamp_varbcrnatbl_flen_allisos' ].to_csv( config[ 'data_out_dir' ] + samp + '_' + config[ 'exon_name' ] + '_by_bc_effects_allvars_allisos-' + date_string + '.txt',
                                                                                sep = '\t' )

    msamp_varbcrnatbl_flen_rename = { samp: mbcs.merge_subasm_and_rna_tbls( satbl,
                                                                       msamp_bcs_processed[ samp ][ 'msamp_bcrnatbl_rename' ] )
                                      for samp in msamp_bcs_processed }

    msamp_varbcrnatbl_flen_rename_filt = { samp: msamp_varbcrnatbl_flen_rename[ samp ].loc[ msamp_varbcrnatbl_flen_rename[ samp ].usable_reads > cut_dict[ samp ] ].copy()
                                           for samp in msamp_varbcrnatbl_flen_rename }

    bybcvartbl_filt_long = mbcs.combine_rep_perbctbls_long( [ msamp_varbcrnatbl_flen_rename_filt[ samp ] for samp in msamp_varbcrnatbl_flen_rename_filt ],
                                                       [ samp for samp in msamp_varbcrnatbl_flen_rename_filt ] )

    bybcvartbl_long = mbcs.combine_rep_perbctbls_long( [ msamp_varbcrnatbl_flen_rename[ samp ] for samp in msamp_varbcrnatbl_flen_rename ],
                                                  [ samp for samp in msamp_varbcrnatbl_flen_rename ] )

    bybcvartbl_long.to_csv( config[ 'data_out_dir' ] + config[ 'exon_name' ] + '_by_bc_effects_allvars-' + date_string + '.txt',
                                sep='\t'
                                )

    msamp_byvartbl_allisos_snvs = { samp: mbcs.filter_byvartbl_snvonly( msamp_bcs_processed[ samp ][ 'msamp_byvartbl_allisos' ] )
                                    for samp in msamp_bcs_processed }

    if min( msamp_byvartbl_allisos_snvs[ samp ].pos.min() for samp in msamp_byvartbl_allisos_snvs ) < config[ 'cloned_vstart' ] \
       or max( msamp_byvartbl_allisos_snvs[ samp ].pos.max() for samp in msamp_byvartbl_allisos_snvs ) > config[ 'cloned_vend' ]:

       if config[ 'verbose' ]:
           print( 'Your data goes outside the specified cloned vector end or cloned vector start positions - truncating or else the plots will be funny' )

       with open( config[ 'data_out_dir' ] + 'mpsa_cml_log.' + date_string + '.txt', 'a' ) as f:
           f.write( 'Your data goes outside the specified cloned vector end or cloned vector start positions - truncating or else the plots will be funny\n' )
           for samp in msamp_byvartbl_allisos_snvs:
               samp_len = len( msamp_byvartbl_allisos_snvs[ samp ] )
               msamp_byvartbl_allisos_snvs[ samp ] = msamp_byvartbl_allisos_snvs[ samp ].loc[ ( msamp_byvartbl_allisos_snvs[ samp ].pos >= config[ 'cloned_vstart' ] ) \
                                                                                            &  ( msamp_byvartbl_allisos_snvs[ samp ].pos <= config[ 'cloned_vend' ] ) ].copy()
               if len( msamp_byvartbl_allisos_snvs[ samp ] ) != samp_len:
                   f.write( 'Truncated %i variants outside cloned region in sample %s\n' % ( samp_len - len( msamp_byvartbl_allisos_snvs[ samp ] ), samp ) )

    for samp in msamp_byvartbl_allisos_snvs:

        msamp_byvartbl_allisos_snvs[ samp ][ 'hgvs_pos' ] = cd.pos_to_hgvspos( msamp_byvartbl_allisos_snvs[ samp ].pos,
                                                                            ( config[ 'cloned_vstart' ], config[ 'cloned_vend' ] ),
                                                                            [ incl_iso[ 0 ] ],
                                                                            [ ( config[ 'exon_hgvs_start' ], config[ 'exon_hgvs_end' ] ), ]
                                                                              )

    byvartbl_allisos_long = mbcs.combine_allisos_pervartbls_long( [ msamp_byvartbl_allisos_snvs[ samp ] for samp in msamp_byvartbl_allisos_snvs ],
                                                             [ samp for samp in msamp_byvartbl_allisos_snvs ] )

    byvartbl_allisos_long.to_csv( config[ 'data_out_dir' ] + config[ 'exon_name' ] + '_by_var_effects_allisos_snvs-' + date_string + '.txt',
                                  sep = '\t',
                                  index = False
                                )

    msamp_byvartbl = { samp: mbcs.summarize_byvar_singlevaronly_pe( satbl,
                                                                msamp_bcs_processed[ samp ][ 'msamp_bcrnatbl_rename' ],
                                                                cut_dict[ samp ],
                                                                [ 'secondary_reads', 'unpaired_reads', 'unmapped_reads', 'bad_starts', 'bad_ends', 'soft_clipped', 'other_isoform', 'usable_reads', ],
                                                                list( isonamedict.keys() ) + [ 'OTHER' ] )
                      for samp in msamp_bcs_processed }

    msamp_byvartbl_snvs = { samp: mbcs.filter_byvartbl_snvonly( msamp_byvartbl[ samp ] )
                            for samp in msamp_byvartbl }

    for samp in msamp_byvartbl_snvs:

        msamp_byvartbl_snvs[ samp ] = msamp_byvartbl_snvs[ samp ].loc[ ( msamp_byvartbl_snvs[ samp ].pos >= config[ 'cloned_vstart' ] ) \
                                                                       &  ( msamp_byvartbl_snvs[ samp ].pos <= config[ 'cloned_vend' ] ) ].copy()

    for samp in msamp_byvartbl_snvs:

        msamp_byvartbl_snvs[ samp ][ 'hgvs_pos' ] = cd.pos_to_hgvspos( msamp_byvartbl_snvs[ samp ].pos,
                                                                    ( config[ 'cloned_vstart' ], config[ 'cloned_vend' ] ),
                                                                    [ incl_iso[ 0 ] ],
                                                                    [ ( config[ 'exon_hgvs_start' ], config[ 'exon_hgvs_end' ] ), ]
                                                                   )

    for samp in msamp_byvartbl_snvs:

        msamp_byvartbl_snvs[ samp ][ 'hg19_pos' ] = cd.vpos_to_gpos( msamp_byvartbl_snvs[ samp ].pos,
                                                                    incl_iso[ 0 ],
                                                                    [ config[ 'exon_hg19_start' ], config[ 'exon_hg19_end' ] ],
                                                                    rev_strand = config[ 'strand' ] == '-' )

    byvartbl_long = mbcs.combine_rep_pervartbls_long( [ msamp_byvartbl_snvs[ samp ] for samp in msamp_byvartbl_snvs ],
                                                  [ samp for samp in msamp_byvartbl_snvs ] )

    byvartbl_long[ 'exon_num' ] = config[ 'exon_number' ]
    byvartbl_long[ 'exon' ] = ( byvartbl_long.pos >= config[ 'exon_vstart' ] ) & ( byvartbl_long.pos <= config[ 'exon_vend' ] )
    #puts ref/alt relative to forward strand - helps with future merges
    if config[ 'strand' ] == '-':
        byvartbl_long[ 'ref_c' ] = byvartbl_long.ref
        byvartbl_long[ 'alt_c' ] = byvartbl_long.alt

        byvartbl_long[ 'ref' ] = [ cd.rev_complement( r ) for r in byvartbl_long.ref_c ]
        byvartbl_long[ 'alt' ] = [ cd.rev_complement( a ) for a in byvartbl_long.alt_c ]

    byvartbl_wide = mbcs.combine_rep_pervartbls_wide( [ byvartbl_long.loc[ byvartbl_long[ 'sample' ] == samp ][ [ col for col in byvartbl_long if col != 'sample' ] ]
                                                    for samp in byvartbl_long[ 'sample' ].unique() ],
                                                    [ samp for samp in byvartbl_long[ 'sample' ].unique() ],
                                                    indexcols=[ 'chrom','hg19_pos','hgvs_pos','pos','ref','alt','varlist','exon_num','exon' ] + ( config[ 'strand' ] == '-')*[ 'ref_c', 'alt_c' ],
                                                    group_cols_by_samp = True )

    byvartbl_wide = mbcs.create_variables_across_samples( byvartbl_wide,
                                                     [ samp for samp in byvartbl_long[ 'sample' ].unique() ],
                                                     mean_cols = [ 'n_bc','n_bc_passfilt', 'sum_reads','sum_reads_passfilt', ],
                                                     median_cols = [ 'n_bc','n_bc_passfilt', 'sum_reads','sum_reads_passfilt', ],
                                                     sum_cols = [ 'sum_reads','sum_reads_passfilt', ],
                                                     max_cols = [ 'n_bc','n_bc_passfilt', 'sum_reads','sum_reads_passfilt', ] )

    byvartbl_wide = mbcs.compute_bc_weighted_psi( byvartbl_wide,
                                             [ samp for samp in byvartbl_long[ 'sample' ].unique() ],
                                             list( isonamedict.keys() ) + [ 'OTHER' ],
                                             'n_bc_passfilt',
                                            )
    if config[ 'strand' ] == '+':
        byvartbl_wide[ 'hgvs_var' ] = [ p + r + '>' + a
                                        for p,( r, a ) in zip( byvartbl_wide.hgvs_pos,
                                                               zip( byvartbl_wide.ref, byvartbl_wide.alt ) ) ]
    elif config[ 'strand' ] == '-':
        byvartbl_wide[ 'hgvs_var' ] = [ p + r + '>' + a
                                        for p,( r, a ) in zip( byvartbl_wide.hgvs_pos,
                                                               zip( byvartbl_wide.ref_c, byvartbl_wide.alt_c ) ) ]
    print( 'Adding in public data and bioinformatics predictions if requested..' )

    if config[ 'gnomad2_file' ]:
        byvartbl_wide = spd.merge_data_gnomad( byvartbl_wide,
                                                gnomad2,
                                                suffix = 'v2' )
        byvartbl_wide = spd.merge_data_gnomad( byvartbl_wide,
                                                gnomad3,
                                                suffix = 'v3' )
        byvartbl_long = spd.merge_data_gnomad( byvartbl_long,
                                                gnomad2,
                                                suffix = 'v2' )
        byvartbl_long = spd.merge_data_gnomad( byvartbl_long,
                                                gnomad3,
                                                suffix = 'v3' )

        byvartbl_wide[ 'gnomad' ] = ( byvartbl_wide.n_alt_v2 > 0 ) | ( byvartbl_wide.n_alt_v3 > 0 )
        byvartbl_long[ 'gnomad' ] = ( byvartbl_long.n_alt_v2 > 0 ) | ( byvartbl_long.n_alt_v3 > 0 )

        byvartbl_long[ 'gnomad_af_max' ] = byvartbl_long[ ['af_v2', 'af_v3'] ].max( axis = 1 )
        byvartbl_wide[ 'gnomad_af_max' ] = byvartbl_wide[ ['af_v2', 'af_v3'] ].max( axis = 1 )

        with open( config[ 'data_out_dir' ] + 'mpsa_cml_log.' + date_string + '.txt', 'a' ) as f:
            f.write( '%i variants intersect with gnomAD entries\n' % byvartbl_wide.gnomad.sum() )
            f.write( '%i variants intersect with common (AF>1%%) gnomAD entries\n' % ( byvartbl_wide.gnomad_af_max >= .01 ).sum() )
            f.write( '%.4f is the maximum AF among gnomAD variants within these data\n' % byvartbl_wide.gnomad_af_max.max() )

    if config[ 'clinvar_file' ]:
        byvartbl_wide = spd.merge_data_clinvar( byvartbl_wide,
                                            clinvar,
                                            keep_cols = [ 'clinvar_interp', 'clinvar_abbrev', 'clinvar_gene' ] )
        byvartbl_long = spd.merge_data_clinvar( byvartbl_long,
                                                clinvar,
                                                keep_cols = [ 'clinvar_interp', 'clinvar_abbrev', 'clinvar_gene' ] )

        byvartbl_wide.loc[ byvartbl_wide.clinvar_interp.isnull(), 'clinvar_abbrev' ] = ''
        byvartbl_long.loc[ byvartbl_long.clinvar_interp.isnull(), 'clinvar_abbrev' ] = ''
        byvartbl_wide.loc[ byvartbl_wide.clinvar_interp.isnull(), 'clinvar_interp' ] = ''
        byvartbl_long.loc[ byvartbl_long.clinvar_interp.isnull(), 'clinvar_interp' ] = ''

        with open( config[ 'data_out_dir' ] + 'mpsa_cml_log.' + date_string + '.txt', 'a' ) as f:
                    f.write( '%i variants intersect with ClinVar entries\n' % ( byvartbl_wide.clinvar_interp != '' ).sum() )
                    f.writelines( [ '%i ClinVar variants classified as %s\n' % ( ( byvartbl_wide.clinvar_interp == interp ).sum(), interp )
                                    for interp in byvartbl_wide.clinvar_interp.unique() if interp != '' ] )

    if config[ 'maxentscan' ]:

        maxent_wt = sm.maxent_score_wt( chrom_refseq,
                                         byvartbl_wide.hg19_pos.unique().tolist(),
                                         rev_strand = config[ 'strand' ] == '-',
                                         pos_out_col = 'hg19_pos' )

        byvartbl_wide = byvartbl_wide.set_index( [ 'hg19_pos', 'ref' ] ).merge( maxent_wt.set_index( [ 'hg19_pos', 'ref' ] ),
                                                                            how = 'left',
                                                                            left_index = True,
                                                                            right_index = True ).reset_index()
        byvartbl_long = byvartbl_long.set_index( [ 'hg19_pos', 'ref' ] ).merge( maxent_wt.set_index( [ 'hg19_pos', 'ref' ] ),
                                                                       how = 'left',
                                                                       left_index = True,
                                                                       right_index = True ).reset_index()

        if config[ 'strand' ] == '+':
            ref_col = 'ref'
            alt_col = 'alt'
        elif config[ 'strand' ] == '-':
            ref_col = 'ref_c'
            alt_col = 'alt_c'

        acceptors = list( set( [ iso[ 0 ][ 0 ] for iso in canonical_isos if iso != () ] ) )
        donors = list( set( [ iso[ 0 ][ 1 ] for iso in canonical_isos if iso != () ] ) )

        acceptor_dw = sm.maxent_score_acceptors( refseq,
                                                byvartbl_wide,
                                                'pos',
                                                ref_col,
                                                alt_col,
                                                acceptors,
                                            )
        acceptor_dl = sm.maxent_score_acceptors( refseq,
                                                byvartbl_long,
                                                'pos',
                                                ref_col,
                                                alt_col,
                                                acceptors,
                                            )
        donor_dw = sm.maxent_score_donors( refseq,
                                          byvartbl_wide,
                                          'pos',
                                          ref_col,
                                          alt_col,
                                          donors,
                                        )
        donor_dl = sm.maxent_score_donors( refseq,
                                          byvartbl_long,
                                          'pos',
                                          ref_col,
                                          alt_col,
                                          donors,
                                        )

        for acc in acceptor_dw:
            byvartbl_wide[ 'maxent_acc_%i' % acc ] = acceptor_dw[ acc ]
            byvartbl_wide[ 'maxent_acc_%i_diff' % acc ] = byvartbl_wide[ 'maxent_acc_%i' % acc ] - byvartbl_wide.loc[ byvartbl_wide.pos == acc ].maxent_wt_acc.mean()
            byvartbl_long[ 'maxent_acc_%i' % acc ] = acceptor_dl[ acc ]
            byvartbl_long[ 'maxent_acc_%i_diff' % acc ] = byvartbl_long[ 'maxent_acc_%i' % acc ] - byvartbl_long.loc[ byvartbl_long.pos == acc ].maxent_wt_acc.mean()
        for don in donor_dw:
            byvartbl_wide[ 'maxent_don_%i' % don ] = donor_dw[ don ]
            byvartbl_wide[ 'maxent_don_%i_diff' % don ] = byvartbl_wide[ 'maxent_don_%i' % don ] - byvartbl_wide.loc[ byvartbl_wide.pos == don ].maxent_wt_don.mean()
            byvartbl_long[ 'maxent_don_%i' % don ] = donor_dl[ don ]
            byvartbl_long[ 'maxent_don_%i_diff' % don ] = byvartbl_long[ 'maxent_don_%i' % don ] - byvartbl_long.loc[ byvartbl_long.pos == don ].maxent_wt_don.mean()

    byvartbl_long.to_csv( config[ 'data_out_dir' ] + config[ 'exon_name' ] + '_by_var_effects_snvs-' + date_string + '.txt',
                          sep = '\t',
                          index = False )

    byvartbl_wide.to_csv( config[ 'data_out_dir' ] + config[ 'exon_name' ] + '_by_var_effects_snvs_wide-' + date_string + '.txt',
                          sep = '\t',
                          index = False )


    if config[ 'spliceai_file' ]:
        if len( splai.chrom.unique() ) == 1 and 'chr' in splai.chrom.unique()[ 0 ]:
            byvartbl_wide[ 'chrom' ] = chrom_chr
            byvartbl_long[ 'chrom' ] = chrom_chr
        elif len( splai.chrom.unique() ) == 1 and 'chr' not in splai.chrom.unique()[ 0 ]:
            byvartbl_wide[ 'chrom' ] = chrom
            byvartbl_long[ 'chrom' ] = chrom

        byvartbl_wide = sm.merge_splai( byvartbl_wide,
                                        splai, )
        byvartbl_long = sm.merge_splai( byvartbl_long,
                                        splai, )

        with open( config[ 'data_out_dir' ] + 'mpsa_cml_log.' + date_string + '.txt', 'a' ) as f:
            f.write( '%i variants classified as SDV at SpliceAI threshold of .2\n' % ( ( byvartbl_wide.DS_maxm >= .2 ).sum() ) )

    byvartbl_long.to_csv( config[ 'data_out_dir' ] + config[ 'exon_name' ] + '_by_var_effects_snvs-' + date_string + '.txt',
                          sep = '\t',
                          index = False )

    byvartbl_wide.to_csv( config[ 'data_out_dir' ] + config[ 'exon_name' ] + '_by_var_effects_snvs_wide-' + date_string + '.txt',
                          sep = '\t',
                          index = False )

    print( 'Data saved - plotting now!' )

    byvartbl_long[ 'n_bc_passfilt_log10' ] = np.log10( byvartbl_long[ 'n_bc_passfilt' ] + .1 )

    #variants that don't match the genomic reference screw up a lot of the plots
    ref_idx = pd.DataFrame( { 'hg19_pos': [ p for p in range( byvartbl_long.hg19_pos.min(), byvartbl_long.hg19_pos.max() + 1 ) ],
                              'ref': [ chrom_refseq[ p - 1 ] for p in range( byvartbl_long.hg19_pos.min(), byvartbl_long.hg19_pos.max() + 1 ) ] } )
    ref_idx = ref_idx.set_index( [ 'hg19_pos', 'ref' ] )
    keep_idx = ref_idx.index.intersection( byvartbl_long.set_index( [ 'hg19_pos', 'ref' ] ).index )
    byvartbl_long_filt = byvartbl_long.set_index( [ 'hg19_pos', 'ref' ] ).loc[ keep_idx ].reset_index().copy()

    sat_by_samp = { samp: sp.saturate_variants( byvartbl_long_filt.loc[ byvartbl_long_filt[ 'sample' ] == samp ],
                                             chrom_refseq,
                                             'hg19_pos',
                                             'exon_num',
                                             rev_strand = config[ 'strand' ] == '-' )[ config[ 'exon_number' ] ]
                    for samp in byvartbl_long_filt[ 'sample' ].unique() }

    #we need the position filled in even if the variant isn't present...
    if config[ 'strand' ] == '+':
        offset = byvartbl_long.hg19_pos - byvartbl_long.pos
    elif config[ 'strand' ] == '-':
        offset = 2*byvartbl_long.hg19_pos.min() - byvartbl_long.hg19_pos - byvartbl_long.pos

    for samp in sat_by_samp:
        if config[ 'strand' ] == '+':
            sat_by_samp[ samp ].pos = sat_by_samp[ samp ].hg19_pos - int( offset.median() )
        elif config[ 'strand' ] == '-':
            sat_by_samp[ samp ].pos = sat_by_samp[ samp ].hg19_pos - int( offset.median() ) - 2*( sat_by_samp[ samp ].hg19_pos - sat_by_samp[ samp ].hg19_pos.min() )

    #check that you didn't alter the positions!
    pos19 = int( byvartbl_long.hg19_pos.sample() )
    posv = int( byvartbl_long.loc[ byvartbl_long.hg19_pos == pos19 ].pos.unique()[ 0 ] )
    for samp in sat_by_samp:
        if pos19 in sat_by_samp[ samp ].hg19_pos.values:
            assert sat_by_samp[ samp ].loc[ sat_by_samp[ samp ].hg19_pos == pos19 ].pos.unique()[ 0 ] == posv, 'Positioning is off in saturated data!'

    #have to put back the hgvs positions for new variants that are missing in the data but included in the saturated data
    for samp in sat_by_samp:
        sat_by_samp[ samp ][ 'hgvs_pos' ] = cd.pos_to_hgvspos( sat_by_samp[ samp ].pos,
                                                               ( config[ 'cloned_vstart' ], config[ 'cloned_vend' ] ),
                                                               [ incl_iso[ 0 ] ],
                                                               [ ( config[ 'exon_hgvs_start' ], config[ 'exon_hgvs_end' ] ), ]
                                                      )

    #want to plot the transcribed base
    if config[ 'strand' ] == '-':
        for samp in sat_by_samp:
            sat_by_samp[ samp ] = sat_by_samp[ samp ].rename( columns = { 'ref': 'r',
                                                                          'alt': 'a',
                                                                          'ref_c': 'ref',
                                                                          'alt_c': 'alt' } )

    github_colors = '3182bd6baed69ecae1c6dbefe6550dfd8d3cfdae6bfdd0a231a35474c476a1d99bc7e9c0756bb19e9ac8bcbddcdadaeb636363969696bdbdbdd9d9d9'
    light_colors = [ '#' + github_colors[i:i+6] for i in range( 0, len( github_colors ), 6 ) ]

    for samp in sat_by_samp:

        sp.sat_subplots_wrapper( sat_by_samp[ samp ],
                                 [ col for col in sat_by_samp[ samp ] if col.startswith( 'wmean_' ) ] + [ 'n_bc_passfilt_log10' ],
                                 'hgvs_pos',
                                 [ l for idx,l in enumerate( light_colors ) if idx%4 == 0 ],
                                 fig_size = ( 40, 3*len( [ col for col in sat_by_samp[ samp ] if col.startswith( 'wmean_' ) ] + [ 'n_bc_passfilt_log10' ] ) ),
                                 share_y = False,
                                 y_ax_lim = [ ( 0, 1 ), ]*len( [ col for col in sat_by_samp[ samp ] if col.startswith( 'wmean_' ) ] ) + [ ( -1, 4 ) ],
                                 y_ax_title = [ col[ 6: ] for col in sat_by_samp[ samp ] if col.startswith( 'wmean_' ) ] + [ 'n_bc_log10' ],
                                 x_ax_title = 'cDNA position',
                                 y_label_rotation = 0,
                                 tick_spacing = 10,
                                 cml = True,
                                 savefile = config[ 'plots_out_dir' ] + config[ 'exon_name' ] + '_' + samp + '_isoform_psis.' + date_string + '.pdf',
                                )

    sample_stats = mbcs.across_sample_stats( [ byvartbl_long ],
                                            { 'all': byvartbl_long[ 'sample' ].unique().tolist() },
                                            [ col for col in byvartbl_long if col.startswith( 'wmean_' ) ]
                                            )

    sample_stats.to_csv( config[ 'data_out_dir' ] + config[ 'exon_name' ] + '_filtered_read_summary_stats.' + date_string + '.txt',
                        sep = '\t',
                        index = False )

    for col in sample_stats:

        if sample_stats[ col ].notnull().sum() == 0:
            continue
        #checks if that the column is numeric - i is int, u is unsigned int, f is float
        elif sample_stats[ col ].dtype.kind not in 'iuf':
            continue
        elif col.startswith( 'per_' ):
            sp.barplot_across_samples( sample_stats,
                                   col,
                                   y_label = col,
                                   y_lim = ( 0, 100 ),
                                   color_col = 'sample_group',
                                   color_dict = { 'all': 'green' },
                                   cml = True,
                                   savefile = config[ 'plots_out_dir' ] + 'filtered_sampsummary_' + col + '.pdf',
                                 )
        elif col.startswith( 'per_' ) and sample_stats[ col ].max() <= 20:
            sp.barplot_across_samples( sample_stats,
                                    col,
                                    y_label = col,
                                    color_col = 'sample_group',
                                    color_dict = { 'all': 'green' },
                                    cml = True,
                                    savefile = config[ 'plots_out_dir' ] + 'filtered_sampsummary_' + col + '_zoomed.pdf',
                                  )
        elif col.startswith( 'med_' ):
            sp.barplot_across_samples( sample_stats,
                                    col,
                                    y_label = col,
                                    y_lim = ( 0, 1 ),
                                    color_col = 'sample_group',
                                    color_dict = { 'all': 'green' },
                                    cml = True,
                                    savefile = config[ 'plots_out_dir' ] + 'filtered_sampsummary_' + col + '.pdf',
                                  )
        else:
            sp.barplot_across_samples( sample_stats,
                                    col,
                                    y_label = col,
                                    color_col = 'sample_group',
                                    color_dict = { 'all': 'green' },
                                    cml = True,
                                    savefile = config[ 'plots_out_dir' ] + 'filtered_sampsummary_' + col + '_counts.pdf',
                                 )

    byvartbl_wide_filt = byvartbl_wide.set_index( [ 'hg19_pos', 'ref' ] ).loc[ keep_idx ].reset_index().copy()

    if len( byvartbl_wide ) != len( byvartbl_wide_filt ):
        if config[ 'verbose' ]:
            print( '%i variants do not match the genomic reference base - removing for plots but they are in the saved files!' % ( len( byvartbl_wide ) - len( byvartbl_wide_filt ) ) )
        with open( config[ 'data_out_dir' ] + 'mpsa_cml_log.' + date_string + '.txt', 'a' ) as f:
            f.write( '%i variants do not match the genomic reference base - removing for plots but they are in the saved files!\n' % ( len( byvartbl_wide ) - len( byvartbl_wide_filt ) ) )

    byvartbl_wide_sat = sp.saturate_variants( byvartbl_wide_filt,
                                              chrom_refseq,
                                              'hg19_pos',
                                              'exon_num',
                                              rev_strand = config[ 'strand' ] == '-' )[ config[ 'exon_number' ] ]

    if config[ 'strand' ] == '+':
        byvartbl_wide_sat.pos = byvartbl_wide_sat.hg19_pos - int( offset.median() )
    elif config[ 'strand' ] == '-':
        byvartbl_wide_sat.pos = byvartbl_wide_sat.hg19_pos - int( offset.median() ) - 2*( byvartbl_wide_sat.hg19_pos - byvartbl_wide_sat.hg19_pos.min() )

    #check that you didn't alter the positions!
    pos19 = int( byvartbl_wide.hg19_pos.sample() )
    posv = int( byvartbl_wide.loc[ byvartbl_wide.hg19_pos == pos19 ].pos.unique()[ 0 ] )
    assert byvartbl_wide_sat.loc[ byvartbl_wide_sat.hg19_pos == pos19 ].pos.unique()[ 0 ] == posv, 'Positioning is off in saturated data!'

    #have to put back the hgvs positions for any variants missing in measured dataset
    byvartbl_wide_sat[ 'hgvs_pos' ] = cd.pos_to_hgvspos( byvartbl_wide_sat.pos,
                                                      ( config[ 'cloned_vstart' ], config[ 'cloned_vend' ] ),
                                                      [ incl_iso[ 0 ] ],
                                                      [ ( config[ 'exon_hgvs_start' ], config[ 'exon_hgvs_end' ] ), ] )

    max_y = byvartbl_wide.n_bc_passfilt_mean.max()
    median_y = byvartbl_wide.n_bc_passfilt_mean.median()

    wmean_cols = [ col for col in byvartbl_wide_sat if col.startswith( 'wmean_' ) ]

    #this makes sure we're plotting base subs on the transcribed strand
    if config[ 'strand' ] == '-':
        byvartbl_wide_sat = byvartbl_wide_sat.rename( columns = { 'ref': 'r',
                                                                  'alt': 'a',
                                                                  'ref_c': 'ref',
                                                                  'alt_c': 'alt' } )

    sp.split_ax_bcs( byvartbl_wide_sat,
                      [ 'n_bc_passfilt_mean' ],
                      'hgvs_pos',
                      [ l for idx,l in enumerate( light_colors ) if idx%4 == 0 ],
                      [ ( median_y + ( ( max_y - median_y ) / 2 ) + .5, max_y ), ( 0, median_y + ( ( max_y - median_y ) / 2 ) ) ],
                      fig_size = ( 40, 7.5 ),
                      legend = False,
                      hratios = [ 1, 4 ],
                      x_ax_title = 'cDNA position',
                      cml = True,
                      savefile = config[ 'plots_out_dir' ] + config[ 'exon_name' ] + '_bcs_filt_by_pos.' + date_string + '.pdf',
                  )

    sp.sat_subplots_wrapper( byvartbl_wide_sat,
                              wmean_cols,
                              'hgvs_pos',
                              [ l for idx,l in enumerate( light_colors ) if idx%4 == 0 ],
                              fig_size = ( 40, 3*len( wmean_cols ) ),
                              share_y = True,
                              legend = False,
                              y_ax_lim = [ ( 0, 1 ), ],
                              y_ax_title = [ col[ 6: ] for col in wmean_cols ],
                              x_ax_title = 'cDNA position',
                              y_label_rotation = 0,
                              tick_spacing = 10,
                              cml = True,
                              savefile = config[ 'plots_out_dir' ] + config[ 'exon_name' ] + '_isoform_psis.' + date_string + '.pdf',
                          )

    if config[ 'gnomad2_file' ]:

        lit_marker_d = { True: ( 'o', 'white', 'black', 3, 200 ), }

        sp.sat_subplots_wrapper( byvartbl_wide_sat,
                                 wmean_cols,
                                 'hgvs_pos',
                                 [ l for idx,l in enumerate( light_colors ) if idx%4 == 0 ],
                                 fig_size = ( 40, 3*len( wmean_cols ) ),
                                 share_y = True,
                                 legend = False,
                                 y_ax_lim = [ ( 0, 1 ), ],
                                 y_ax_title = [ col[ 6: ] for col in wmean_cols ],
                                 x_ax_title = 'cDNA position',
                                 y_label_rotation = 0,
                                 tick_spacing = 10,
                                 bar_labels = [ ( 'gnomad', lit_marker_d, 1.2 ) ],
                                 tight = False,
                                 save_margin = 1,
                                 cml = True,
                                 savefile = config[ 'plots_out_dir' ] + config[ 'exon_name' ] + '_isoform_psis_gnomad.' + date_string + '.pdf',
                                 )

    if config[ 'clinvar_file' ]:

        byvartbl_wide_sat[ 'clinvar_abbrev' ] = byvartbl_wide_sat.clinvar_abbrev.fillna( '' )

        lit_marker_d = { 'LBB': ( 's', 'white', 'black', 3, 100 ),
                             'LPP': ( '^', 'black', 'face', 1.5, 100 ),
                             'VUS': ( 'd', 'black', 'face', 1.5, 100 ) }

        sp.sat_subplots_wrapper( byvartbl_wide_sat,
                                 wmean_cols,
                                 'hgvs_pos',
                                 [ l for idx,l in enumerate( light_colors ) if idx%4 == 0 ],
                                 fig_size = ( 40, 3*len( wmean_cols ) ),
                                 share_y = True,
                                 legend = False,
                                 y_ax_lim = [ ( 0, 1 ), ],
                                 y_ax_title = [ col[ 6: ] for col in wmean_cols ],
                                 x_ax_title = 'cDNA position',
                                 y_label_rotation = 0,
                                 tick_spacing = 10,
                                 bar_labels = [ ( 'clinvar_abbrev', lit_marker_d, 1.2 ) ],
                                 tight = False,
                                 save_margin = 1,
                                 cml = True,
                                 savefile = config[ 'plots_out_dir' ] + config[ 'exon_name' ] + '_isoform_psis_clinvar.' + date_string + '.pdf',
                                 )

    if config[ 'maxentscan' ]:

        sp.sat_lollipop_subplots_wrapper( byvartbl_wide_sat,
                                             wmean_cols + [ col for col in byvartbl_wide_sat if col.endswith( '_diff' ) ],
                                             'hgvs_pos',
                                             [ l for idx,l in enumerate( light_colors ) if idx%4 == 0 ],
                                             fig_size = ( 40, 3*( len( wmean_cols ) + len( [ col for col in byvartbl_wide_sat if col.endswith( '_diff' ) ] ) )),
                                             share_y = False,
                                             legend = False,
                                             lollipop_size = 8,
                                             linewidth = 4,
                                             y_ax_lim = [ ( 0, 1 ) ]*len( wmean_cols ) + [ ( -15, 15 ) ]*len( [ col for col in byvartbl_wide_sat if col.endswith( '_diff' ) ] ),
                                             y_ax_title = [ col[ 6: ] for col in wmean_cols ] + [ col[ 7: ] for col in byvartbl_wide_sat if col.endswith( '_diff' ) ],
                                             x_ax_title = 'hgvs position',
                                             y_label_rotation = 0,
                                             hlines = [ ( 0, 'black', 'solid' ) for i in range( len( wmean_cols ) + len( [ col for col in byvartbl_wide_sat if col.endswith( '_diff' ) ] ) ) ],
                                             tick_spacing = 10,
                                             cml = True,
                                             savefile = config[ 'plots_out_dir' ] + config[ 'exon_name' ] + '_maxent.' + date_string + '.pdf',
                                             )

    if config[ 'spliceai_file' ]:

        sp.sat_subplots_wrapper( byvartbl_wide_sat,
                                  wmean_cols + [ 'DS_maxm' ],
                                  'hgvs_pos',
                                  [ l for idx,l in enumerate( light_colors ) if idx%4 == 0 ],
                                  fig_size = ( 40, 3*( len( wmean_cols ) + 1 ) ),
                                  share_y = True,
                                  legend = False,
                                  y_ax_lim = [ ( 0, 1 ), ],
                                  y_ax_title = [ col[ 6: ] for col in wmean_cols ] + [ 'DS_maxm' ],
                                  x_ax_title = 'cDNA position',
                                  y_label_rotation = 0,
                                  tick_spacing = 10,
                                  cml = True,
                                  savefile = config[ 'plots_out_dir' ] + config[ 'exon_name' ] + '_isoform_psis_spliceai.' + date_string + '.pdf',
                              )

    if config[ 'clinvar_file' ]:

        lit_marker_d = { 'LBB': ( 's', 'white', 'black', 3, 200 ),
                         'LPP': ( '^', 'black', 'face', 1.5, 200 ),
                         'VUS': ( 'd', 'black', 'face', 1.5, 200 ) }

        if config[ 'gnomad2_file' ]:

            sp.plot_clinvar_by_interp( byvartbl_wide.loc[ ( byvartbl_wide.clinvar_abbrev != '' ) ],
                                    wmean_cols,
                                    lit_marker_d,
                                    interp_col = 'clinvar_abbrev',
                                    figsize = ( ( byvartbl_wide.clinvar_abbrev != '' ).sum() / 2, 3*len( wmean_cols ) ),
                                    ylim = ( 0, 1 ),
                                    row_label = [ col[ 6: ] for col in wmean_cols ],
                                    gnomad_labels = [ ( 'gnomad_af_max',
                                                      ( 'o', 'lightgray', 200 ),
                                                      1.3 ) ],
                                    cml = True,
                                    savefig = config[ 'plots_out_dir' ] + config[ 'exon_name' ] + '_clinvar.' + date_string + '.pdf'
                                    )
        else:

            sp.plot_clinvar_by_interp( byvartbl_wide.loc[ ( byvartbl_wide.clinvar_abbrev != '' ) ],
                                    wmean_cols,
                                    lit_marker_d,
                                    interp_col = 'clinvar_abbrev',
                                    figsize = ( ( byvartbl_wide.clinvar_abbrev != '' ).sum() / 2, 3*len( wmean_cols ) ),
                                    ylim = ( 0, 1 ),
                                    row_label = [ col[ 6: ] for col in wmean_cols ],
                                    cml = True,
                                    savefig = config[ 'plots_out_dir' ] + config[ 'exon_name' ] + '_clinvar.' + date_string + '.pdf'
                                    )

    T1 = time.time()

    if config[ 'verbose' ]:
        print( 'Total runtime: %.2f hours!' % ( ( T1 - T0 ) / 60 / 60 ) )

    with open( config[ 'data_out_dir' ] + 'mpsa_cml_log.' + date_string + '.txt', 'a' ) as f:
        f.write( 'Total runtime: %.2f hours\n' % ( ( T1 - T0 ) / 60 / 60 ) )
        f.write( '\nFigure legends:\n' )
        f.write( 'Alt colors: A = Blue, C = Orange, G = Green, T = Purple\n' )
        f.write( 'gnomAD plots: In gnomAD = open circle\n')
        f.write( 'ClinVar plots: LBB = open square, LPP = black triangle, VUS/Conflicting = black diamond\n' )
        f.write( 'Calling variants as SDV takes some time and consideration - do it in a notebook!\n' )
        f.write(  'Check out a post_processing notebook as an example of how to call SDV\n' )

    print( 'Calling variants as SDV takes some time and consideration - do it in a notebook!' )
    print( 'Check out a post_processing notebook as an example of how to call SDV' )
    print( 'See ya next time on The Splice is Right - MPSA edition!' )

    return

if __name__ == "__main__":
    main()
