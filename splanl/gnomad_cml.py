import splanl.scrape_public_data as spd
import splanl.coords as cd
import argparse
import pysam
import pandas as pd
import subprocess as subp
import time

def main():

    parser = argparse.ArgumentParser( description = 'Scrape gnomAD to merge with MSPA data',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument( '-l', '--liftover_location',
                         help = 'Location of UCSC liftover software: dir (str)' )
    parser.add_argument( '-lc', '--hg19tohg38_chain_location',
                         help = 'Location of UCSC liftover hg19 to hg38 chain: dir (str)' )
    parser.add_argument( '-g', '--gene_name',
                         default = '',
                         help = 'Gene name to append to outfile names (str)' )
    parser.add_argument( '-ex', '--exon_name',
                         default = '',
                         help = 'Exon name to append to outfile names (str)' )
    parser.add_argument( 'gnomad2_exomes',
                         help = 'Dir + file for tabix indexed gnomAD2 exomes dataset - tabix gene/small region from gnomAD online (str)' )
    parser.add_argument( 'gnomad2_genomes',
                         help = 'Dir + file for tabix indexed gnomAD2 genomes dataset - tabix gene/small region from gnomAD online (str)' )
    parser.add_argument( 'gnomad3',
                         help = 'Dir + file for tabix indexed gnomAD3 dataset - tabix gene/small region from gnomAD online (str)' )
    parser.add_argument( 'chrom',
                         help = 'Chromosome of region to scrape (str)' )
    parser.add_argument( 'hg19_start',
                         type = int,
                         help = 'HG19 start position to scrape (int)' )
    parser.add_argument( 'hg19_end',
                         type = int,
                         help = 'HG19 end position to scrape (int)' )
    parser.add_argument( 'hg38_start',
                         type = int,
                         help = 'HG38 start position to scrape (int)' )
    parser.add_argument( 'hg38_end',
                         type = int,
                         help = 'HG38 end position to scrape (int)' )
    parser.add_argument( 'gnomad_outdir',
                         help = 'Out directory to place scraped gnomAD files (str)' )
    args = parser.parse_args()
    config = vars(args)

    if config[ 'liftover_location' ] or config[ 'hg19tohg38_chain_location' ]:
        assert config[ 'liftover_location' ] and config[ 'hg19tohg38_chain_location' ], \
        'Need both UCSC liftover software location (-l) and UCSC hg19 to hg38 chain location (-lc) to perform liftover!\nCan skip both and liftover yourself outside of this script!'

    if not config[ 'gnomad_outdir' ].endswith( '/' ):
        config[ 'gnomad_outdir' ] = config[ 'gnomad_outdir' ] + '/'

    if config[ 'chrom' ].startswith( 'chr' ):
        chrom_chr = config[ 'chrom' ]
        chrom = config[ 'chrom' ][ 3: ]
    else:
        chrom = config[ 'chrom' ]
        chrom_chr = 'chr' + config[ 'chrom' ]

    gnomad2_exomes_tbx = pysam.VariantFile( config[ 'gnomad2_exomes' ] )
    gnomad2_genomes_tbx = pysam.VariantFile( config[ 'gnomad2_genomes' ] )
    gnomad3_tbx = pysam.VariantFile( config[ 'gnomad3' ] )

    gnomad2_exomes_df = spd.create_gnomad_df( gnomad2_exomes_tbx,
                                              chrom,
                                              ( config[ 'hg19_start' ], config[ 'hg19_end' ] ),
                                              'hg19' )

    gnomad2_genomes_df = spd.create_gnomad_df( gnomad2_genomes_tbx,
                                              chrom,
                                              ( config[ 'hg19_start' ], config[ 'hg19_end' ] ),
                                              'hg19' )

    gnomad3_df = spd.create_gnomad_df( gnomad3_tbx,
                                       chrom_chr,
                                       ( config[ 'hg38_start' ], config[ 'hg38_end' ] ),
                                       'hg38' )

    gnomad2_df = spd.merge_gnomad2_data( gnomad2_genomes_df,
                                         gnomad2_exomes_df )

    if config[ 'gene_name' ] or config[ 'exon_name' ]:
        outstem = '_'.join( [ config[ 'gene_name' ], config[ 'exon_name' ] ] )
    else:
        outstem = ''

    if config[ 'liftover_location' ]:

        date_string = time.strftime( '%Y-%m%d', time.localtime() )

        hg19_lift = cd.create_liftover_bed( chrom_chr,
                                            config[ 'hg19_start' ] - 5,
                                            config[ 'hg19_end' ] + 5 )

        hg19_lift.to_csv( config[ 'gnomad_outdir' ] + 'hg19_lift_' + outstem + '.bed',
                          sep = '\t',
                          index = False,
                          header = False )

        out = subp.run( '%s %s \
                         %s \
                         %s \
                         %s' % ( config[ 'liftover_location' ],
                                 config[ 'gnomad_outdir' ] + 'hg19_lift_' + outstem + '.bed',
                                 config[ 'hg19tohg38_chain_location' ],
                                 config[ 'gnomad_outdir' ] + 'hg38_lift_' + outstem + '.bed',
                                 config[ 'gnomad_outdir' ] + 'unmapped_lift_' + outstem + '.bed' ),
                        stdout = subp.PIPE,
                        stderr = subp.PIPE,
                        shell = True,
                     )

        stdout_nl = out.stdout.decode('utf-8').split( '\n' )
        err_nl = out.stderr.decode('utf-8').split( '\n' )

        with open( config[ 'gnomad_outdir' ] + 'gnomad_cml_log.' + date_string + '.txt', 'a' ) as f:
            f.write( 'Liftover stdout\n' )
            f.writelines( stdout_nl )
            f.write( 'Liftover stderr\n' )
            f.writelines( err_nl )

        hg38_lift = pd.read_table( config[ 'gnomad_outdir' ] + 'hg38_lift_' + outstem + '.bed',
                                    names = [ 'chrom', 'hg38_pos', 'end', 'name' ] )

        lift = hg19_lift.drop( columns = 'end' ).merge( hg38_lift.drop( columns = 'end' ),
                                                        on = [ 'chrom', 'name' ],
                                                        how = 'outer' )

        gnomad3_df = cd.liftover_hg38_to_hg19( gnomad3_df,
                                               lift )

    else:

        print( 'gnomAD3 dataset will not have HG19 coords - please add them before attempting to merge in main MPSA driver script!')

        with open( config[ 'gnomad_outdir' ] + 'gnomad_cml_log.' + date_string + '.txt', 'a' ) as f:
            f.write( 'gnomAD3 dataset will not have HG19 coords - please add them before attempting to merge in main MPSA driver script!\n' )

    gnomad2_df.to_csv( config[ 'gnomad_outdir' ] + 'gnomad_v2_' + outstem + '.txt',
                       sep = '\t',
                       index = False
                       )

    gnomad3_df.to_csv( config[ 'gnomad_outdir' ] + 'gnomad_v3_' + outstem + '.txt',
                       sep = '\t',
                       index = False
                       )

    return

if __name__ == "__main__":
    main()
