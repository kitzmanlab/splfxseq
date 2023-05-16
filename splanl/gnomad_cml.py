import splanl.scrape_public_data as spd
import argparse
import pysam
import pandas as pd

def main():

    parser = argparse.ArgumentParser( description = 'Scrape gnomAD to merge with MSPA data',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument( '-ll', '--lazy_liftover',
                         action = 'store_true',
                         help = 'Perform a lazy liftover (HG38 to HG19) on gnomAD3 based on the difference in HG19 and HG38 starts' )
    parser.add_argument( '-g', '--gene_name',
                         default = '',
                         help = 'Gene name to append to outfile names (str)' )
    parser.add_argument( '-ex', '--exon_name',
                         default = '',
                         help = 'Exon name to append to outfile names (str)' )
    parser.add_argument( 'gnomad2',
                         help = 'Dir + file for tabix indexed gnomAD2 dataset - tabix gene/small region from gnomAD online (str)' )
    parser.add_argument( 'gnomad3',
                         help = 'Dir + file for tabix indexed gnomAD3 dataset - tabix gene/small region from gnomAD online (str)' )
    parser.add_argument( 'chrom',
                         help = 'Chromosome of region to scrape (str)' )
    parser.add_argument( 'hg19_start',
                         type = int,
                         help = 'HG19 start position to scrape - have this be the liftover of hg38 start! (int)' )
    parser.add_argument( 'hg19_end',
                         type = int,
                         help = 'HG19 end position to scrape (int)' )
    parser.add_argument( 'hg38_start',
                         type = int,
                         help = 'HG38 start position to scrape - have this be the liftover of hg19 start! (int)' )
    parser.add_argument( 'hg38_end',
                         type = int,
                         help = 'HG38 end position to scrape (int)' )
    parser.add_argument( 'gnomad_outdir',
                         help = 'Out directory to place scraped gnomAD files (str)' )
    args = parser.parse_args()
    config = vars(args)

    if not config[ 'gnomad_outdir' ].endswith( '/' ):
        config[ 'gnomad_outdir' ] = config[ 'gnomad_outdir' ] + '/'

    if config[ 'chrom' ].startswith( 'chr' ):
        chrom_chr = config[ 'chrom' ]
        chrom = config[ 'chrom' ][ 3: ]
    else:
        chrom = config[ 'chrom' ]
        chrom_chr = 'chr' + config[ 'chrom' ]

    gnomad2_tbx = pysam.TabixFile( config[ 'gnomad2' ] )
    gnomad3_tbx = pysam.TabixFile( config[ 'gnomad3' ] )

    gnomad2_df = spd.create_gnomad_df( gnomad2_tbx,
                                       chrom,
                                       ( config[ 'hg19_start' ], config[ 'hg19_end' ] ),
                                       'hg19' )

    gnomad3_df = spd.create_gnomad_df( gnomad3_tbx,
                                       chrom_chr,
                                       ( config[ 'hg38_start' ], config[ 'hg38_end' ] ),
                                       'hg38' )

    if config[ 'lazy_liftover' ]:

        print( 'CAUTION! This script does a lazy liftover based on the HG19 and HG38 start positions supplied!' )
        print( 'Please be sure the HG19 and HG38 starts are at the same bp!' )

        offset = config[ 'hg19_start' ] - config[ 'hg38_start' ]
        gnomad3_df[ 'hg19_pos' ] = gnomad3_df.hg38_pos + offset

    else:

        print( 'gnomAD3 dataset will not have HG19 coords - please add them before attempting to merge in main MPSA driver script!')

    if config[ 'gene_name' ] or config[ 'exon_name' ]:
        outstem = '_'.join( [ config[ 'gene_name' ], config[ 'exon_name' ] ] )
    else:
        outstem = ''

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
