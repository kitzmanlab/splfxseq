import pandas as pd
import numpy as np
import argparse
import sys
import time
import csv
from keras.models import load_model
from pkg_resources import resource_filename
import splanl.plots as sp
import splanl.post_processing as pp
import splanl.coords as cd
import splanl.annots as ann
import splanl.create_vcf as vcf
import splanl.custom_splai_scores as css
import splanl.scrape_public_data as spd
import splanl.score_motifs as sm

def main():

    T0 = time.time()

    parser = argparse.ArgumentParser( description = 'Score SNVs with SpliceAI!',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument( '-rn', '--refseqname',
                         help = 'Name of sequence in reference fasta file (str) ie chrom11' )
    parser.add_argument( '-ann', '--annotation_file_override',
                         help = 'Annotation file in splai format to use - overrides custom transcript specific file! : dir + filename (str)' )
    parser.add_argument( '-d', '--scoring_distance',
                         type = int,
                         help = 'Scored distance (default is exon length + 10) (int)' )
    parser.add_argument( '-altss', '--score_alt_ss',
                         help = 'Compute WT and SS probabilities for a site other than the exon bds (comma sep ints - ie 1234,)' )
    parser.add_argument( '-vcf', '--vcf_out_dir',
                         help = 'Directory to store vcfs - only enter if you want a vcf of unscored variants (str)' )
    parser.add_argument( '-g2', '--gnomad2_file',
                         help = 'Directory + file of gnomAD2 processed by scrape_public_data module' )
    parser.add_argument( '-g3', '--gnomad3_file',
                         help = 'Directory + file of gnomAD3 processed by scrape_public_data module - MUST HAVE hg19 coords as hg19_pos column!' )
    parser.add_argument( '-cl', '--clinvar_file',
                         help = 'Directory + file of ClinVar processed by scrape_public_data module - MUST HAVE hg19 coords as hg19_pos column!' )
    parser.add_argument( '-mcol', '--merge_pos_col',
                         help = 'Column used to merge gnomAD and/or ClinVar with Splai - use same coord system for all! (str)' )
    parser.add_argument( '-me', '--maxentscan',
                         action = 'store_true',
                         help = 'Include MaxEntScan scores?' )
    parser.add_argument( '-v', '--verbose',
                         action = 'store_true',
                         help = 'Increase number of printouts!' )
    parser.add_argument( 'gencode_basic',
                         help = 'Gencode basic annotation file: dir + filename (str)' )
    parser.add_argument( 'gencode_gene_id',
                         help = 'Gencode gene identifier (starts with ENSG) (str)' )
    parser.add_argument( 'gencode_transcript_id',
                         help = 'Gencode transcript identifier (starts with ENST) (str)' )
    parser.add_argument( 'reference_fafile',
                         help = 'Reference fasta file: dir + filename (str)' )
    parser.add_argument( 'score_start_pos',
                         type = int,
                         help = 'Start position to score - be consistent with reference_fafile coords! (int)' )
    parser.add_argument( 'score_end_pos',
                         type = int,
                         help = 'End position to score - be consistent with reference_fafile coords! (int)' )
    parser.add_argument( 'exon_acceptor_pos',
                         type = int,
                         help = 'First bp of exon (transcribed strand) - used for scoring distance unless overridden by -d! (int)' )
    parser.add_argument( 'exon_donor_pos',
                         type = int,
                         help = 'Last bp of exon (transcribed strand) - used for scoring distance unless overridden by -d! (int)' )
    parser.add_argument( 'chrom',
                         help = 'Chromosome - can be with or without chr (str)' )
    parser.add_argument( 'strand',
                         choices = [ '+', '-' ],
                         help = 'Strand: either + or - (- will score rev comp of your reference_fafile so do not flip for vector coords!) (str)' )
    parser.add_argument( 'annot_out_dir',
                         help = 'Directory to store annotation files (str)' )
    parser.add_argument( 'plots_out_dir',
                         help = 'Directory to store plots (str)' )
    parser.add_argument( 'data_out_dir',
                         help = 'Directory to store datasets (str)' )
    args = parser.parse_args()
    config = vars(args)

    if not config[ 'annot_out_dir' ].endswith( '/' ):
        config[ 'annot_out_dir' ] = config[ 'softclip_out_dir' ] + '/'
    if not config[ 'plots_out_dir' ].endswith( '/' ):
        config[ 'plots_out_dir' ] = config[ 'plots_out_dir' ] + '/'
    if not config[ 'data_out_dir' ].endswith( '/' ):
        config[ 'data_out_dir' ] = config[ 'data_out_dir' ] + '/'
    if config[ 'vcf_out_dir' ] and not config[ 'vcf_out_dir' ].endswith( '/' ):
        config[ 'vcf_out_dir' ] = config[ 'vcf_out_dir' ] + '/'

    if config[ 'strand' ] == '+':
        assert config[ 'exon_acceptor_pos' ] < config[ 'exon_donor_pos' ], \
        'Acceptor should have a LOWER position than the donor for + strand! Check you exon_acceptor and exon_donor positions!'
    elif config[ 'strand' ] == '-':
        assert config[ 'exon_donor_pos' ] < config[ 'exon_acceptor_pos' ], \
        'Acceptor should have a HIGHER position than the donor for - strand! Check you exon_acceptor and exon_donor positions!'

    if config[ 'gnomad2_file' ]:
        assert config[ 'gnomad3_file' ], 'This script assumes you entered both versions of gnomAD or neither! Add gnomAD3!'
        gnomad2 = pd.read_table( config[ 'gnomad2_file' ] )
        assert config[ 'merge_pos_col' ] in gnomad2, 'Your gnomAD2 dataframe must have the merge_pos_col (-mcol) as one of the columns!'
    if config[ 'gnomad3_file' ]:
        assert config[ 'gnomad2_file' ], 'This script assumes you entered both versions of gnomAD or neither! Add gnomAD2!'
        gnomad3 = pd.read_table( config[ 'gnomad3_file' ] )
        assert config[ 'merge_pos_col' ] in gnomad3, 'Your gnomAD3 dataframe must have the merge_pos_col (-mcol) as one of the columns!'

    if config[ 'clinvar_file' ]:
        clinvar = pd.read_table( config[ 'clinvar_file' ] )
        assert config[ 'merge_pos_col' ] in clinvar, 'Your ClinVar dataframe must have the merge_pos_col (-mcol) as one of the columns!'

    chrom = config[ 'chrom' ] if not config[ 'chrom' ].startswith( 'chr' ) else config[ 'chrom' ][ 3: ]
    chrom_chr = config[ 'chrom' ] if config[ 'chrom' ].startswith( 'chr' ) else 'chr' + config[ 'chrom' ]

    date_string = time.strftime( '%Y-%m%d', time.localtime() )

    gencode_annot = pd.read_table(  config[ 'gencode_basic' ],
                                    compression = 'gzip',
                                    sep = '\t',
                                    comment = '#',
                                    names = [ 'chrom', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute' ] )

    trans_gtf = ann.gene_specific_gtf( gencode_annot,
                                       config[ 'gencode_gene_id' ],
                                       config[ 'gencode_transcript_id' ] )

    trans_gtf.to_csv( config[ 'annot_out_dir' ] + 'gencode_' + config[ 'gencode_transcript_id' ] + '.gtf.gz',
                      sep = '\t',
                      index = False,
                      header = False,
                      quoting = csv.QUOTE_NONE,
                      compression = 'gzip' )

    trans_splai = ann.gtf_to_splai_annot( trans_gtf )

    trans_splai.to_csv( config[ 'annot_out_dir' ] + 'splai_' + config[ 'gencode_transcript_id' ] + '.txt',
                        sep = '\t',
                        index = False, )

    if config[ 'verbose' ]:
        print( 'Transcript specific annotation files created and saved at %s' % config[ 'annot_out_dir' ] )
    with open( config[ 'data_out_dir' ] + 'splai_log.' + date_string + '.txt', 'a' ) as f:
        f.write( 'Transcript specific annotation files created and saved at %s\n' % config[ 'annot_out_dir' ] )

    if config[ 'annotation_file_override' ]:
        annot = config[ 'annotation_file_override' ]
        if config[ 'verbose' ]:
            print( 'Scoring against overide annotation file provided: %s' % config[ 'annotation_file_override' ] )
        with open( config[ 'data_out_dir' ] + 'splai_log.' + date_string + '.txt', 'a' ) as f:
            f.write( 'Scoring against overide annotation file provided: %s\n' % config[ 'annotation_file_override' ] )
    else:
        annot = trans_splai

    refseq_d = pp.get_refseq( config[ 'reference_fafile' ] )

    if config[ 'refseqname' ]:
        if config[ 'refseqname' ] in refseq_d:
            refseq = refseq_d[ config[ 'refseqname' ] ]
        else:
            if config[ 'verbose' ]:
                print( 'Refseqname (-rn) %s not in fasta - trying %s' % ( config[ 'refseqname' ], chrom_chr  ) )
            with open( config[ 'data_out_dir' ] + 'splai_log.' + date_string + '.txt', 'a' ) as f:
                f.write( 'Refseqname (-rn) %s not in fasta - trying %s\n' % ( config[ 'refseqname' ], chrom_chr  ) )
            assert chrom_chr in refseq_d, 'Refseqname (-rn) %s and chromosome %s not in provided fasta file - only %s!' % ( config[ 'refseqname' ], chrom_chr, ', '.join( list( refseq_d.keys() ) ) )
            refseq = refseq_d[ chrom_chr ]
    elif len( refseq_d ) == 1:
        if config[ 'verbose' ]:
            print( 'No refseqname (-rn) provided - using %s' % list( refseq_d.keys() )[ 0 ] )
        with open( config[ 'data_out_dir' ] + 'splai_log.' + date_string + '.txt', 'a' ) as f:
            f.write( 'No refseqname (-rn) provided - using %s\n' % list( refseq_d.keys() )[ 0 ] )
        refseq = list( refseq_d.values() )[ 0 ]
    else:
        if config[ 'verbose' ]:
            print( 'No refseqname (-rn) provided - trying %s' % chrom_chr  )
        with open( config[ 'data_out_dir' ] + 'splai_log.' + date_string + '.txt', 'a' ) as f:
            f.write( 'No refseqname (-rn) provided - trying %s\n' % chrom_chr )
        assert chrom_chr in refseq_d, 'Chromosome %s not in provided fasta file - only %s! Try specifying refseqname (-rn)...' % ( chrom_chr, ', '.join( list( refseq_d.keys() ) ) )
        refseq = refseq_d[ chrom_chr ]

    #garbage collect the memory in case a lot of sequences provided
    refseq_d = {}

    assert config[ 'score_start_pos' ] != config[ 'score_end_pos' ], 'Distance between score_start_pos and score_end_pos must be > 0!'

    if config[ 'score_start_pos' ] > config[ 'score_end_pos' ]:
        junk = config[ 'score_start_pos' ]
        config[ 'score_start_pos' ] = config[ 'score_end_pos' ]
        config[ 'score_end_pos' ] = junk

    if config[ 'vcf_out_dir' ]:

        vcf.create_input_vcf( vcf.vcf_header(),
                              vcf.tup_to_vcf( [ ( chrom, pos, refseq[ pos - 1 ].upper(), alt )
                                                for pos in range( config[ 'score_start_pos' ], config[ 'score_end_pos' ] + 1 )
                                                for alt in [ 'A', 'C', 'G', 'T' ]
                                                if alt != refseq[ pos - 1 ].upper() ] ),
                              config[ 'vcf_out_dir' ] + config[ 'gencode_transcript_id' ] + '.vcf' )

        vcf.create_input_vcf( vcf.vcf_header_chr(),
                              vcf.tup_to_vcf( [ ( chrom_chr, pos, refseq[ pos - 1 ].upper(), alt )
                                                for pos in range( config[ 'score_start_pos' ], config[ 'score_end_pos' ] + 1 )
                                                for alt in [ 'A', 'C', 'G', 'T' ]
                                                if alt != refseq[ pos - 1 ].upper() ] ),
                               config[ 'vcf_out_dir' ] + config[ 'gencode_transcript_id' ] + '_chr' + '.vcf' )

        if config[ 'verbose' ]:
            print( 'VCFs created and saved at %s' % config[ 'vcf_out_dir' ] )
        with open( config[ 'data_out_dir' ] + 'splai_log.' + date_string + '.txt', 'a' ) as f:
            f.write( 'VCFs created and saved at %s\n' % config[ 'vcf_out_dir' ] )

    if config[ 'verbose' ]:
        print( 'Whirling up SpliceAI - takes a few minutes..\n' )
        print( 'No need to worry about the training configuration warning!\n' )

    paths = ('models/spliceai{}.h5'.format(x) for x in range(1, 6))
    models = [load_model(resource_filename('spliceai', x)) for x in paths]

    centers = [ ( p, refseq[ p - 1 ].upper(), alt )
                  for p in range( config[ 'score_start_pos' ], config[ 'score_end_pos' ] + 1 )
                  for alt in [ 'A', 'C', 'G', 'T' ]
                  if alt != refseq[ p - 1 ].upper() ]

    if not config[ 'scoring_distance' ]:
        dist = np.abs( config[ 'exon_acceptor_pos' ] - config[ 'exon_donor_pos' ] ) + 10
    else:
        dist = config[ 'scoring_distance' ]

    if config[ 'verbose' ]:
        print( 'Scoring all SNVs!' )

    t0 = time.time()

    #hacky way to suppress the tensorflow messages
    save_stdout = sys.stdout
    sys.stdout = open('trash', 'w')

    splai_scores = css.splai_score_mult_variants_onegene( annot,
                                                            models,
                                                            refseq,
                                                            config[ 'gencode_transcript_id' ],
                                                            str( annot.CHROM.values[ 0 ]  ),
                                                            centers,
                                                            scored_context = dist,
                                                            rev_strand = config[ 'strand' ] == '-' )
    sys.stdout = save_stdout

    if config[ 'strand' ] == '-':
        splai_scores[ 'ref_c' ] = [ cd.rev_complement( r ) for r in splai_scores.ref ]
        splai_scores[ 'alt_c' ] = [ cd.rev_complement( a ) for a in splai_scores.alt ]

    splai_scores.to_csv( config[ 'data_out_dir' ] + config[ 'gencode_transcript_id' ] + '_splai_SNVs_' + date_string + '.txt',
                         sep = '\t',
                         index = False )

    t1 = time.time()

    if config[ 'verbose' ]:
        print( 'Scored %i SNVs in %.2f minutes!' % ( len( centers ), ( t1 - t0 ) / 60 ) )

    with open( config[ 'data_out_dir' ] + 'splai_log.' + date_string + '.txt', 'a' ) as f:
        f.write( 'Scored %i SNVs in %.2f minutes!\n' % ( len( centers ), ( t1 - t0 ) / 60 ) )

    if config[ 'verbose' ]:
        print( 'Scoring WT and SS probabilities now - this is faster..' )

    t0 = time.time()

    ss_list = []
    if config[ 'score_start_pos' ] <= config[ 'exon_acceptor_pos' ] and config[ 'score_end_pos' ] >= config[ 'exon_acceptor_pos' ]:
        ss_list.append( config[ 'exon_acceptor_pos' ] )
    if config[ 'score_start_pos' ] <= config[ 'exon_donor_pos' ] and config[ 'score_end_pos' ] >= config[ 'exon_donor_pos' ]:
        ss_list.append( config[ 'exon_donor_pos' ] )

    if not config[ 'score_alt_ss' ]:
        alt_ss = []
    else:
        alt_ss = [ int( p ) for p in config[ 'score_alt_ss' ].split( ',' )
                    if p != '' and int( p ) >= config[ 'score_start_pos' ] and int( p ) <= config[ 'score_end_pos' ] ]
        ss_list += alt_ss
        if config[ 'verbose' ]:
            print( 'Scoring alternate splice sites at %s' % ( ', '.join( [ str( ss ) for ss in alt_ss ] ) ) )
        with open( config[ 'data_out_dir' ] + 'splai_log.' + date_string + '.txt', 'a' ) as f:
            f.write( 'Scoring alternate splice sites at %s\n' % ( ', '.join( [ str( ss ) for ss in alt_ss ] ) ) )

    save_stdout = sys.stdout
    sys.stdout = open('trash', 'w')

    splai_wt_pr  =  css.splai_score_wt_onegene( annot,
                                                  models,
                                                  refseq,
                                                  config[ 'gencode_transcript_id' ],
                                                  str( annot.CHROM.values[ 0 ]  ),
                                                  ss_list,
                                                  scored_context = dist,
                                                  rev_strand = config[ 'strand' ] == '-' )

    sys.stdout = save_stdout

    splai_wt_pr.to_csv( config[ 'data_out_dir' ] + config[ 'gencode_transcript_id' ] + '_splai_WT_' + date_string + '.txt',
                         sep = '\t',
                         index = False )

    if ss_list:

        save_stdout = sys.stdout
        sys.stdout = open('trash', 'w')

        splai_ss_pr, splai_acc_pr, splai_don_pr = css.splai_ss_prob_mult_variants_onegene( annot,
                                                                                            models,
                                                                                            refseq,
                                                                                            config[ 'gencode_transcript_id' ],
                                                                                            str( annot.CHROM.values[ 0 ] ),
                                                                                            centers,
                                                                                            ss_list,
                                                                                            scored_context = dist,
                                                                                            rev_strand = config[ 'strand' ] == '-' )

        sys.stdout = save_stdout

        if config[ 'strand' ] == '-':
            splai_ss_pr[ 'ref_c' ] = [ cd.rev_complement( r ) for r in splai_ss_pr.ref ]
            splai_ss_pr[ 'alt_c' ] = [ cd.rev_complement( a ) for a in splai_ss_pr.alt ]

        splai_ss_pr.to_csv( config[ 'data_out_dir' ] + config[ 'gencode_transcript_id' ] + '_splai_SS_probs_' + date_string + '.txt',
                            sep = '\t',
                            index = False )

    else:

        if config[ 'verbose' ]:
            print( 'No exon bounds or alternate splice sites within scoring range - splice site probabilities not provided!')
        with open( config[ 'data_out_dir' ] + 'splai_log.' + date_string + '.txt', 'a' ) as f:
            f.write( 'No exon bounds or alternate splice sites within scoring range - splice site probabilities not provided!\n' )

    t1 = time.time()

    if config[ 'verbose' ]:
        print( 'Scored WT and SS probs in %.2f minutes!' % ( ( t1 - t0 ) / 60 ) )

    with open( config[ 'data_out_dir' ] + 'splai_log.' + date_string + '.txt', 'a' ) as f:
        f.write( 'Scored WT and SS probs in %.2f minutes!\n' % ( ( t1 - t0 ) / 60 ) )

    if config[ 'verbose' ]:
        print( 'Merging public databases and MaxEntScan scores, if requested!' )

    if config[ 'gnomad2_file' ]:

        splai_scores[ config[ 'merge_pos_col' ] ] = splai_scores.pos

        splai_scores = spd.merge_data_gnomad( splai_scores,
                                                gnomad2,
                                                suffix = 'v2' )
        splai_scores = spd.merge_data_gnomad( splai_scores,
                                                gnomad3,
                                                suffix = 'v3' )
        splai_scores[ 'gnomad' ] = ( splai_scores.n_alt_v2 > 0 ) | ( splai_scores.n_alt_v3 > 0 )
        splai_scores[ 'gnomad_af_max' ] = splai_scores[ ['af_v2', 'af_v3'] ].max( axis = 1 )

        if ss_list:
            splai_ss_pr[ config[ 'merge_pos_col' ] ] = splai_ss_pr.pos
            splai_ss_pr = spd.merge_data_gnomad( splai_ss_pr,
                                                gnomad2,
                                                suffix = 'v2' )
            splai_ss_pr = spd.merge_data_gnomad( splai_ss_pr,
                                                gnomad3,
                                                suffix = 'v3' )
            splai_ss_pr[ 'gnomad' ] = ( splai_ss_pr.n_alt_v2 > 0 ) | ( splai_ss_pr.n_alt_v3 > 0 )
            splai_ss_pr[ 'gnomad_af_max' ] = splai_ss_pr[ ['af_v2', 'af_v3'] ].max( axis = 1 )

    if config[ 'clinvar_file' ]:

        splai_scores[ config[ 'merge_pos_col' ] ] = splai_scores.pos
        splai_scores = spd.merge_data_clinvar( splai_scores,
                                               clinvar,
                                               keep_cols = [ 'clinvar_interp', 'clinvar_abbrev', 'clinvar_gene' ] )
        splai_scores.loc[ splai_scores.clinvar_interp.isnull(), 'clinvar_abbrev' ] = ''

        if ss_list:
            splai_ss_pr[ config[ 'merge_pos_col' ] ] = splai_ss_pr.pos
            splai_ss_pr = spd.merge_data_clinvar( splai_ss_pr,
                                                  clinvar,
                                                  keep_cols = [ 'clinvar_interp', 'clinvar_abbrev', 'clinvar_gene' ] )
            splai_ss_pr.loc[ splai_ss_pr.clinvar_interp.isnull(), 'clinvar_abbrev' ] = ''

    if config[ 'maxentscan' ]:

        maxent_wt = sm.maxent_score_wt( refseq,
                                         splai_scores.pos.unique().tolist(),
                                         rev_strand = config[ 'strand' ] == '-' )

        acceptors = [ config[ 'exon_acceptor_pos' ] ] + alt_ss
        donors = [ config[ 'exon_donor_pos' ] ] + alt_ss

        acceptor_d = sm.maxent_score_acceptors( refseq,
                                                 splai_scores,
                                                 'pos',
                                                 'ref',
                                                 'alt',
                                                 acceptors,
                                                 rev_strand = config[ 'strand' ] == '-'
                                                 )

        donor_d = sm.maxent_score_donors( refseq,
                                          splai_scores,
                                          'pos',
                                          'ref',
                                          'alt',
                                          donors,
                                          rev_strand = config[ 'strand' ] == '-'
                                            )

        splai_scores = splai_scores.set_index( [ 'pos', 'ref' ] ).merge( maxent_wt.set_index( [ 'pos', 'ref' ] ),
                                                                     how = 'left',
                                                                     left_index = True,
                                                                      right_index = True ).reset_index()
        if ss_list:
            splai_ss_pr = splai_ss_pr.set_index( [ 'pos', 'ref' ] ).merge( maxent_wt.set_index( [ 'pos', 'ref' ] ),
                                                                        how = 'left',
                                                                        left_index = True,
                                                                        right_index = True ).reset_index()

        for acc in acceptor_d:
            splai_scores[ 'maxent_acc_%i' % acc ] = acceptor_d[ acc ]
            splai_scores[ 'maxent_acc_%i_diff' % acc ] = splai_scores[ 'maxent_acc_%i' % acc ] - splai_scores.loc[ splai_scores.pos == acc ].maxent_wt_acc.mean()
            if ss_list:
                splai_ss_pr[ 'maxent_acc_%i' % acc ] = acceptor_d[ acc ]
                splai_ss_pr[ 'maxent_acc_%i_diff' % acc ] = splai_ss_pr[ 'maxent_acc_%i' % acc ] - splai_ss_pr.loc[ splai_ss_pr.pos == acc ].maxent_wt_acc.mean()
        for don in donor_d:
            splai_scores[ 'maxent_don_%i' % don ] = donor_d[ don ]
            splai_scores[ 'maxent_don_%i_diff' % don ] = splai_scores[ 'maxent_don_%i' % don ] - splai_scores.loc[ splai_scores.pos == don ].maxent_wt_don.mean()
            if ss_list:
                splai_ss_pr[ 'maxent_don_%i' % don ] = donor_d[ don ]
                splai_ss_pr[ 'maxent_don_%i_diff' % don ] = splai_ss_pr[ 'maxent_don_%i' % don ] - splai_ss_pr.loc[ splai_ss_pr.pos == don ].maxent_wt_don.mean()

    splai_scores.to_csv( config[ 'data_out_dir' ] + config[ 'gencode_transcript_id' ] + '_splai_SNVs_' + date_string + '.txt',
                         sep = '\t',
                         index = False )
    if ss_list:
        splai_ss_pr.to_csv( config[ 'data_out_dir' ] + config[ 'gencode_transcript_id' ] + '_splai_SS_probs_' + date_string + '.txt',
                            sep = '\t',
                            index = False )

    if config[ 'verbose' ]:
        print( 'Time to plot!' )

    github_colors = '3182bd6baed69ecae1c6dbefe6550dfd8d3cfdae6bfdd0a231a35474c476a1d99bc7e9c0756bb19e9ac8bcbddcdadaeb636363969696bdbdbdd9d9d9'
    light_colors = [ '#' + github_colors[i:i+6] for i in range( 0, len( github_colors ), 6 ) ]

    wt_probs = [ float( splai_wt_pr.loc[ splai_wt_pr.pos == config[ 'exon_acceptor_pos' ] ].wt_acc_pr ),
                 float( splai_wt_pr.loc[ splai_wt_pr.pos == config[ 'exon_donor_pos' ] ].wt_don_pr ), ] + \
               [ float( splai_wt_pr.loc[ splai_wt_pr.pos == ss ].wt_acc_pr ) for ss in alt_ss ] + \
               [ float( splai_wt_pr.loc[ splai_wt_pr.pos == ss ].wt_don_pr ) for ss in alt_ss ]

    #this is just for plotting purposes - won't be saving the data again
    if config[ 'strand' ] == '-':
        splai_scores[ 'pos' ] = -splai_scores.pos
        splai_scores = splai_scores.rename( columns = { 'alt': 'a',
                                                        'ref': 'r',
                                                        'alt_c': 'alt',
                                                        'ref_c': 'ref' } )
        if ss_list:
            splai_ss_pr.pos = -splai_ss_pr.pos
            splai_ss_pr = splai_ss_pr.rename( columns = { 'alt': 'a',
                                                          'ref': 'r',
                                                          'alt_c': 'alt',
                                                          'ref_c': 'ref' } )

    splai_scores[ 'coord_pos' ] = splai_scores.pos
    if ss_list:
        splai_ss_pr[ 'coord_pos' ] = splai_ss_pr.pos

    sp.sat_subplots_wrapper( splai_scores,
                            [ 'DS_maxm', 'DS_ALm', 'DS_AGm', 'DS_DLm', 'DS_DGm', ],
                            'coord_pos',
                            [ l for idx,l in enumerate( light_colors ) if idx%4 == 0 ],
                            fig_size = ( 40, 15 ),
                            y_ax_lim = [ ( 0, 1 ) ],
                            y_ax_title = [ 'DS_maxm', 'DS_ALm', 'DS_AGm', 'DS_DLm', 'DS_DGm', ],
                            x_ax_title = 'Coordinate position',
                            tick_spacing = 10,
                            cml = True,
                            savefile = config[ 'plots_out_dir' ] + config[ 'gencode_transcript_id' ] + '_splai_SNVs_' + date_string + '.pdf',
                            )
    if ss_list:
        sp.sat_subplots_wrapper( splai_ss_pr,
                                [ 'ss_acc_prob_' + str( config[ 'exon_acceptor_pos' ] ), 'ss_don_prob_' + str( config[ 'exon_donor_pos' ] ) ] + \
                                [ 'ss_acc_prob_' + str( ss ) for ss in alt_ss  ] +  [ 'ss_don_prob_' + str( ss ) for ss in alt_ss  ],
                                'coord_pos',
                                [ l for idx,l in enumerate( light_colors ) if idx%4 == 0 ],
                                fig_size = ( 40, 3*( 2 + 2*len( alt_ss ) ) ),
                                y_ax_lim = [ ( 0, 1 ) ],
                                y_ax_title = [ 'ACC', 'DON', ] + [ 'ACC_' + str( ss ) for ss in alt_ss ] + [ 'DON_' + str( ss ) for ss in alt_ss ],
                                x_ax_title = 'Coordinate position',
                                tick_spacing = 10,
                                y_label_rotation = 0,
                                hlines = [ ( wt_pr, 'gray', 'solid' ) for wt_pr in wt_probs ],
                                cml = True,
                                savefile = config[ 'plots_out_dir' ] + config[ 'gencode_transcript_id' ] + '_splai_SS_probs_' + date_string + '.pdf',
                                )

    if config[ 'clinvar_file' ]:

        lit_marker_d = { 'LBB': ( 's', 'white', 'black', 3, 100 ),
                         'LPP': ( '^', 'black', 'face', 1.5, 100 ),
                         'VUS': ( 'd', 'black', 'face', 1.5, 100 ) }

        sp.sat_subplots_wrapper( splai_scores,
                                [ 'DS_maxm', 'DS_ALm', 'DS_AGm', 'DS_DLm', 'DS_DGm', ],
                                config[ 'merge_pos_col' ],
                                [ l for idx,l in enumerate( light_colors ) if idx%4 == 0 ],
                                fig_size = ( 40, 15 ),
                                y_ax_lim = [ ( 0, 1 ) ],
                                y_ax_title = [ 'DS_maxm', 'DS_ALm', 'DS_AGm', 'DS_DLm', 'DS_DGm', ],
                                x_ax_title = config[ 'merge_pos_col' ],
                                tick_spacing = 10,
                                bar_labels = [ ( 'clinvar_abbrev', lit_marker_d, 1.2 ) ],
                                tight = False,
                                save_margin = 1,
                                cml = True,
                                savefile = config[ 'plots_out_dir' ] + config[ 'gencode_transcript_id' ] + '_splai_SNVs_clinvar_' + date_string + '.pdf',
                                 )

        if ss_list:
            sp.sat_subplots_wrapper( splai_ss_pr,
                                    [ 'ss_acc_prob_' + str( config[ 'exon_acceptor_pos' ] ), 'ss_don_prob_' + str( config[ 'exon_donor_pos' ] ) ] + \
                                    [ 'ss_acc_prob_' + str( ss ) for ss in alt_ss  ] +  [ 'ss_don_prob_' + str( ss ) for ss in alt_ss  ],
                                    config[ 'merge_pos_col' ],
                                    [ l for idx,l in enumerate( light_colors ) if idx%4 == 0 ],
                                    fig_size = ( 40, 3*( 2 + 2*len( alt_ss ) ) ),
                                    y_ax_lim = [ ( 0, 1 ) ],
                                    y_ax_title = [ 'ACC', 'DON', ] + [ 'ACC_' + str( ss ) for ss in alt_ss ] + [ 'DON_' + str( ss ) for ss in alt_ss ],
                                    x_ax_title = config[ 'merge_pos_col' ],
                                    tick_spacing = 10,
                                    bar_labels = [ ( 'clinvar_abbrev', lit_marker_d, 1.2 ) ],
                                    tight = False,
                                    save_margin = 1,
                                    y_label_rotation = 0,
                                    hlines = [ ( wt_pr, 'gray', 'solid' ) for wt_pr in wt_probs ],
                                    cml = True,
                                    savefile = config[ 'plots_out_dir' ] + config[ 'gencode_transcript_id' ] + '_splai_SS_probs_clinvar_' + date_string + '.pdf',
                                    )

    if config[ 'gnomad2_file' ]:

        lit_marker_d = { True: ( 'o', 'white', 'black', 3, 200 ), }

        sp.sat_subplots_wrapper( splai_scores,
                                [ 'DS_maxm', 'DS_ALm', 'DS_AGm', 'DS_DLm', 'DS_DGm', ],
                                config[ 'merge_pos_col' ],
                                [ l for idx,l in enumerate( light_colors ) if idx%4 == 0 ],
                                fig_size = ( 40, 15 ),
                                y_ax_lim = [ ( 0, 1 ) ],
                                y_ax_title = [ 'DS_maxm', 'DS_ALm', 'DS_AGm', 'DS_DLm', 'DS_DGm', ],
                                x_ax_title = config[ 'merge_pos_col' ],
                                tick_spacing = 10,
                                bar_labels = [ ( 'gnomad', lit_marker_d, 1.2 ) ],
                                tight = False,
                                save_margin = 1,
                                cml = True,
                                savefile = config[ 'plots_out_dir' ] + config[ 'gencode_transcript_id' ] + '_splai_SNVs_gnomad_' + date_string + '.pdf',
                                 )

        if ss_list:
            sp.sat_subplots_wrapper( splai_ss_pr,
                                    [ 'ss_acc_prob_' + str( config[ 'exon_acceptor_pos' ] ), 'ss_don_prob_' + str( config[ 'exon_donor_pos' ] ) ] + \
                                    [ 'ss_acc_prob_' + str( ss ) for ss in alt_ss  ] +  [ 'ss_don_prob_' + str( ss ) for ss in alt_ss  ],
                                    config[ 'merge_pos_col' ],
                                    [ l for idx,l in enumerate( light_colors ) if idx%4 == 0 ],
                                    fig_size = ( 40, 3*( 2 + 2*len( alt_ss ) ) ),
                                    y_ax_lim = [ ( 0, 1 ) ],
                                    y_ax_title = [ 'ACC', 'DON', ] + [ 'ACC_' + str( ss ) for ss in alt_ss ] + [ 'DON_' + str( ss ) for ss in alt_ss ],
                                    x_ax_title = config[ 'merge_pos_col' ],
                                    tick_spacing = 10,
                                    bar_labels = [ ( 'gnomad', lit_marker_d, 1.2 ) ],
                                    tight = False,
                                    save_margin = 1,
                                    y_label_rotation = 0,
                                    hlines = [ ( wt_pr, 'gray', 'solid' ) for wt_pr in wt_probs ],
                                    cml = True,
                                    savefile = config[ 'plots_out_dir' ] + config[ 'gencode_transcript_id' ] + '_splai_SS_probs_gnomad_' + date_string + '.pdf',
                                    )

    if config[ 'maxentscan' ]:

        sp.sat_lollipop_subplots_wrapper( splai_scores,
                                            [ 'DS_maxm' ] + [ col for col in splai_scores if col.endswith( '_diff' ) ],
                                            'coord_pos',
                                            [ l for idx,l in enumerate( light_colors ) if idx%4 == 0 ],
                                            fig_size = ( 40, 3*( 1 + len( [ col for col in splai_scores if col.endswith( '_diff' ) ] ) ) ),
                                             share_y = False,
                                             lollipop_size = 6,
                                             linewidth = 3,
                                             y_ax_lim = [ ( 0, 1 ) ] + [ ( -15, 15 ) ]*len( [ col for col in splai_scores if col.endswith( '_diff' ) ] ),
                                             y_ax_title = [ 'SplAI' ] + [ col[ 7: ] for col in splai_scores if col.endswith( '_diff' ) ],
                                             x_ax_title = 'Coord pos',
                                             hlines = [ ( 0, 'black', 'solid' ) for i in range( 1 + len( [ col for col in splai_scores if col.endswith( '_diff' ) ] ) ) ],
                                             tick_spacing = 10,
                                             y_label_rotation = 0,
                                             cml = True,
                                             savefile = config[ 'plots_out_dir' ] + config[ 'gencode_transcript_id' ] + '_splai_SNVs_maxent_' + date_string + '.pdf',
                                             )

    with open( config[ 'data_out_dir' ] + 'splai_log.' + date_string + '.txt', 'a' ) as f:
        f.write( '\nFigure legends:\n' )
        f.write( 'Alt colors: A = Blue, C = Orange, G = Green, T = Purple\n' )
        f.write( 'gnomAD plots: In gnomAD = open circle\n')
        f.write( 'ClinVar plots: LBB = open square, LPP = black triangle, VUS = black diamond\n')

    T1 = time.time()

    if config[ 'verbose' ]:
        print( 'Completing scoring and plots in %.2f minutes!' % ( ( T1 - T0 ) / 60 ) )

    with open( config[ 'data_out_dir' ] + 'splai_log.' + date_string + '.txt', 'a' ) as f:
        f.write( '\nCompleting scoring and plots in %.2f minutes!\n\n' % ( ( T1 - T0 ) / 60 ) )

    print( 'See ya next time on The Splice is Right!' )

    return

if __name__ == "__main__":
    main()
