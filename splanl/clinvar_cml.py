import splanl.scrape_public_data as spd
import argparse
import pysam
import pandas as pd

def main():

    parser = argparse.ArgumentParser( description = 'Scrape ClinVar to merge with MSPA data',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument( '-g', '--gene_name',
                         default = '',
                         help = 'Gene name to append to outfile names (str)' )
    parser.add_argument( '-ex', '--exon_name',
                         default = '',
                         help = 'Exon name to append to outfile names (str)' )
    parser.add_argument( 'clinvar',
                         help = 'Dir + file for tabix indexed ClinVar hg19 dataset (str)' )
    parser.add_argument( 'chrom',
                         help = 'Chromosome of region to scrape (str)' )
    parser.add_argument( 'hg19_start',
                         type = int,
                         help = 'HG19 start position to scrape (int)' )
    parser.add_argument( 'hg19_end',
                         type = int,
                         help = 'HG19 end position to scrape (int)' )
    parser.add_argument( 'clinvar_outdir',
                         help = 'Out directory to place scraped ClinVar file (str)' )
    args = parser.parse_args()
    config = vars(args)

    if not config[ 'clinvar_outdir' ].endswith( '/' ):
        config[ 'clinvar_outdir' ] = config[ 'clinvar_outdir' ] + '/'

    if config[ 'chrom' ].startswith( 'chr' ):
        chrom = config[ 'chrom' ][ 3: ]
    else:
        chrom = config[ 'chrom' ]

    clinvar_tbx = pysam.VariantFile( config[ 'clinvar' ] )

    clinvar_df = spd.create_clinvar_df( clinvar_tbx,
                                        chrom,
                                        ( config[ 'hg19_start' ], config[ 'hg19_end' ] ),
                                        'hg19' )

    if config[ 'gene_name' ] or config[ 'exon_name' ]:
        outstem = '_'.join( [ config[ 'gene_name' ], config[ 'exon_name' ] ] )
    else:
        outstem = ''

    clinvar_df.to_csv( config[ 'clinvar_outdir' ] + 'clinvar_' + outstem + '.txt',
                       sep = '\t',
                       index = False
                       )
    return

if __name__ == "__main__":
    main()
