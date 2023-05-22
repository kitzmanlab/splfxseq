import pandas as pd
import numpy as np
import argparse
import sys
import time
import csv
from keras.models import load_model
from pkg_resources import resource_filename
import splanl.post_processing as pp
import splanl.coords as cd
import splanl.annots as ann
import splanl.custom_splai_scores as css

def main():

    T0 = time.time()

    parser = argparse.ArgumentParser( description = 'Score DNVs with SpliceAI!',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument( '-rn', '--refseqname',
                         help = 'Name of sequence in reference fasta file (str) ie chrom11' )
    parser.add_argument( '-ann', '--annotation_file_override',
                         help = 'Annotation file in splai format to use - overrides custom transcript specific file! : dir + filename (str)' )
    parser.add_argument( '-d', '--scoring_distance',
                         type = int,
                         help = 'Scored distance (default is exon length + 10) (int)' )
    parser.add_argument( '-sspr', '--score_ss_probs',
                         action = 'store_false',
                         help = 'Compute SS probabilities?' )
    parser.add_argument( '-altss', '--score_alt_ss',
                         help = 'Compute SS probabilities for a site other than the exon bds (comma sep ints - ie 1234,)' )
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

    if config[ 'strand' ] == '+':
        assert config[ 'exon_acceptor_pos' ] < config[ 'exon_donor_pos' ], \
        'Acceptor should have a LOWER position than the donor for + strand! Check you exon_acceptor and exon_donor positions!'
    elif config[ 'strand' ] == '-':
        assert config[ 'exon_donor_pos' ] < config[ 'exon_acceptor_pos' ], \
        'Acceptor should have a HIGHER position than the donor for - strand! Check you exon_acceptor and exon_donor positions!'

    if config[ 'score_alt_ss' ] is not None:
        with open( config[ 'data_out_dir' ] + 'splai_dnv_log.' + date_string + '.txt', 'a' ) as f:
            f.write( 'Your alternate splice sites (-altss) will not be scored unless you turn on score splice site probs (-sspr)!\n' )

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
    with open( config[ 'data_out_dir' ] + 'splai_dnv_log.' + date_string + '.txt', 'a' ) as f:
        f.write( 'Transcript specific annotation files created and saved at %s\n' % config[ 'annot_out_dir' ] )

    if config[ 'annotation_file_override' ]:
        annot = config[ 'annotation_file_override' ]
        if config[ 'verbose' ]:
            print( 'Scoring against overide annotation file provided: %s' % config[ 'annotation_file_override' ] )
        with open( config[ 'data_out_dir' ] + 'splai_dnv_log.' + date_string + '.txt', 'a' ) as f:
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
            with open( config[ 'data_out_dir' ] + 'splai_dnv_log.' + date_string + '.txt', 'a' ) as f:
                f.write( 'Refseqname (-rn) %s not in fasta - trying %s\n' % ( config[ 'refseqname' ], chrom_chr  ) )
            assert chrom_chr in refseq_d, 'Refseqname (-rn) %s and chromosome %s not in provided fasta file - only %s!' % ( config[ 'refseqname' ], chrom_chr, ', '.join( list( refseq_d.keys() ) ) )
            refseq = refseq_d[ chrom_chr ]
    elif len( refseq_d ) == 1:
        if config[ 'verbose' ]:
            print( 'No refseqname (-rn) provided - using %s' % list( refseq_d.keys() )[ 0 ] )
        with open( config[ 'data_out_dir' ] + 'splai_dnv_log.' + date_string + '.txt', 'a' ) as f:
            f.write( 'No refseqname (-rn) provided - using %s\n' % list( refseq_d.keys() )[ 0 ] )
        refseq = list( refseq_d.values() )[ 0 ]
    else:
        if config[ 'verbose' ]:
            print( 'No refseqname (-rn) provided - trying %s' % chrom_chr  )
        with open( config[ 'data_out_dir' ] + 'splai_dnv_log.' + date_string + '.txt', 'a' ) as f:
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

    if config[ 'verbose' ]:
        print( 'Whirling up SpliceAI - takes a few minutes..' )
        print( '\nNo need to worry about the training configuration warning!\n' )

    paths = ('models/spliceai{}.h5'.format(x) for x in range(1, 6))
    models = [load_model(resource_filename('spliceai', x)) for x in paths]

    snv_centers = [ ( p, refseq[ p - 1 ].upper(), alt )
                  for p in range( config[ 'score_start_pos' ], config[ 'score_end_pos' ] + 1 )
                  for alt in [ 'A', 'C', 'G', 'T' ]
                  if alt != refseq[ p - 1 ].upper() ]

    dnvs_centers = []
    dnvs_haps = []
    for i, snv_i in enumerate( snv_centers ):
        for j, snv_j in enumerate( snv_centers ):
            if j <= i:
                continue

            dnvs_centers.append( snv_i )
            dnvs_haps.append( [ snv_j ] )

    if not config[ 'scoring_distance' ]:
        dist = np.abs( config[ 'exon_acceptor_pos' ] - config[ 'exon_donor_pos' ] ) + 10
    else:
        dist = config[ 'scoring_distance' ]

    print( 'Scoring all DNVs!' )
    print( 'This will use ~13 GB of memory and take ~1.5 seconds per variant' )
    print( 'You are scoring %i DNVs so expect at least %.2f hours runtime! Might want to use tmux!' % ( len( dnvs_centers ), ( len( dnvs_centers )*1.5 ) / 60 / 60  ) )

    t0 = time.time()

    splai_scores = css.splai_score_mult_variants_onegene( annot,
                                                            models,
                                                            refseq,
                                                            config[ 'gencode_transcript_id' ],
                                                            str( annot.CHROM.values[ 0 ]  ),
                                                            dnvs_centers,
                                                            haplotypes = dnvs_haps,
                                                            scored_context = dist,
                                                            rev_strand = config[ 'strand' ] == '-' )

    if config[ 'strand' ] == '-':
        splai_scores[ 'ref_c' ] = [ cd.rev_complement( r ) for r in splai_scores.ref ]
        splai_scores[ 'alt_c' ] = [ cd.rev_complement( a ) for a in splai_scores.alt ]

    splai_scores.to_csv( config[ 'data_out_dir' ] + config[ 'gencode_transcript_id' ] + '_splai_DNVs_' + date_string + '.txt',
                         sep = '\t',
                         index = False )

    t1 = time.time()

    if config[ 'verbose' ]:
        print( 'Scored %i DNVs in %.2f hours!' % ( len( dnvs_centers ), ( t1 - t0 ) / 60 / 60 ) )

    with open( config[ 'data_out_dir' ] + 'splai_dnv_log.' + date_string + '.txt', 'a' ) as f:
        f.write( 'Scored %i DNVs in %.2f hours!\n' % ( len( dnvs_centers ), ( t1 - t0 ) / 60 / 60 ) )

    if config[ 'score_ss_probs' ]:

        if config[ 'verbose' ]:
            print( 'Scoring SS probabilities now - this is faster..' )

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
                        if p != '' and p >= config[ 'score_start_pos' ] and p <= config[ 'score_end_pos' ] ]
            ss_list += alt_ss
            if config[ 'verbose' ]:
                print( 'Scoring alternate splice sites at %s' % ( ', '.join( alt_ss ) ) )
            with open( config[ 'data_out_dir' ] + 'splai_dnv_log.' + date_string + '.txt', 'a' ) as f:
                f.write( 'Scoring alternate splice sites at %s\n' % ( ', '.join( alt_ss ) ) )

        if ss_list:
            splai_ss_pr, splai_acc_pr, splai_don_pr = css.splai_ss_prob_mult_variants_onegene( annot,
                                                                                                models,
                                                                                                refseq,
                                                                                                config[ 'gencode_transcript_id' ],
                                                                                                str( annot.CHROM.values[ 0 ] ),
                                                                                                dnvs_centers,
                                                                                                ss_list,
                                                                                                scored_context = dist,
                                                                                                haplotypes = dnvs_haps,
                                                                                                rev_strand = config[ 'strand' ] == '-' )

            if config[ 'strand' ] == '-':
                splai_ss_pr[ 'ref_c' ] = [ cd.rev_complement( r ) for r in splai_ss_pr.ref ]
                splai_ss_pr[ 'alt_c' ] = [ cd.rev_complement( a ) for a in splai_ss_pr.alt ]

            splai_ss_pr.to_csv( config[ 'data_out_dir' ] + config[ 'gencode_transcript_id' ] + '_splai_DNV_SS_probs_' + date_string + '.txt',
                                sep = '\t',
                                index = False )

        else:
            if config[ 'verbose' ]:
                print( 'No exon bounds or alternate splice sites within scoring range - splice site probabilities not provided!')
            with open( config[ 'data_out_dir' ] + 'splai_dnv_log.' + date_string + '.txt', 'a' ) as f:
                f.write( 'No exon bounds or alternate splice sites within scoring range - splice site probabilities not provided!\n' )

        t1 = time.time()

        if config[ 'verbose' ]:
            print( 'Scored SS probs in %.2f hours!' % ( ( t1 - t0 ) / 60 / 60 ) )

        with open( config[ 'data_out_dir' ] + 'splai_dnv_log.' + date_string + '.txt', 'a' ) as f:
            f.write( 'Scored SS probs in %.2f hours!\n' % ( ( t1 - t0 ) / 60 / 60 ) )

    else:

        if config[ 'verbose' ]:
            print( 'Not scoring splice site probabilities! Turn (-sspr) on if you want them scored!' )
        with open( config[ 'data_out_dir' ] + 'splai_dnv_log.' + date_string + '.txt', 'a' ) as f:
            f.write( 'Not scoring splice site probabilities! Turn (-sspr) on if you want them scored!\n' )

    T1 = time.time()

    if config[ 'verbose' ]:
        print( 'Completing scoring in %.2f hours!' % ( ( T1 - T0 ) / 60 / 60 ) )

    with open( config[ 'data_out_dir' ] + 'splai_dnv_log.' + date_string + '.txt', 'a' ) as f:
        f.write( '\nCompleting scoring in %.2f hours!\n\n' % ( ( T1 - T0 ) / 60 / 60 ) )

    print( 'See ya next time on The Splice is Right - DNV spin off!' )

    return

if __name__ == "__main__":
    main()
