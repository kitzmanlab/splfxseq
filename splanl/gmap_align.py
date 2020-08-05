import pandas as pd
import numpy as np
import pysam
import subprocess as subp
from os import path

def extract_snv_bcs( satbl ):

    sa = satbl.reset_index().copy()

    sa_singlevar = sa.query( 'n_variants_passing == 1' )[[ 'variant_list', 'readgroupid' ]].copy()

    sa_singlevar[ 'ref' ] = [ str( s ).split( ':' )[ 2 ] for s in sa_singlevar.variant_list ]
    sa_singlevar[ 'alt' ] = [ str( s ).split( ':' )[ 3 ] for s in sa_singlevar.variant_list ]

    sa_snvs = sa_singlevar.loc[ (sa_singlevar.ref.str.len() == 1) & (sa_singlevar.alt.str.len() == 1) ].copy()

    var = sa_snvs.variant_list.tolist()
    bc = sa_snvs.readgroupid.tolist()

    var_to_bc_d = {}
    for v, b in zip(var, bc):
        if v in var_to_bc_d:
            var_to_bc_d[v].append(b)
        else:
            var_to_bc_d[v] = [ b ]

    return var_to_bc_d

def create_varseq( refseq,
                    variant ):

    pos = int( variant.split( ':' )[ 1 ] )
    ref = str( variant.split( ':' )[ 2 ] )
    alt = str( variant.split( ':' )[ 3 ] )

    assert refseq[ pos - 1 ] == ref, 'Reference sequence does not match variant reference allele'

    return refseq[ : ( pos -1 ) ] + alt + refseq[ pos : ]

def write_temp_fa( refseq,
                    variant,
                    tempdir ):

    with open( tempdir + 'temp.fa', 'w' ) as fa:
        fa.write( '>' + variant + '\n' )
        fa.write( create_varseq( refseq,
                                 variant ) )

def extract_read_names( fq_file2,
                        bc_l ):

    with pysam.FastxFile( fq_file2 ) as fq:

        read_names = [ entry.name for entry in fq if any( entry.sequence == bc for bc in bc_l ) ]

    return read_names

def write_temp_fq( fq_file1,
                    read_names_l,
                    tempdir ):

    with pysam.FastxFile( fq_file1 ) as fq_in, open( tempdir + 'temp.fq', 'w' ) as fq_out:

        fq_out.writelines( str( entry ) + '\n' for entry in fq_in if any( entry.name == name for name in read_names_l ) )

def build_index( tempdir,
                reffile = 'temp.fa',
                print_err = False ):

    out = subp.run( '/home/smithcat/bin/gmap_build -d tmp_idx -D %s -k 8 -w 0 %s' % ( tempdir, tempdir + reffile ),
                    stdout = subp.PIPE,
                    stderr = subp.PIPE,
                    shell = True,
                )

    if print_err:
        #the formatting of the output is annoying - trying to make it look nice
        err_nl = out.stderr.decode('utf-8').split( '\n' )
        print(*err_nl, sep='\n')

    #the formatting of the output is annoying - trying to make it look nice
    out_nl = out.stdout.decode('utf-8').split( '\n' )
    print(*out_nl, sep='\n')

def align_reads( tempdir,
                  fqfile = 'temp.fq',
                  print_err = False ):

    out = subp.run( '/home/smithcat/bin/gmap -d tmp_idx -D %s -t 8 -f samse --microexon-spliceprob=1.0 --allow-close-indels=2 %s > %s'
                    % ( tempdir, tempdir + fqfile, tempdir + 'temp.sam' ) ,
                    stderr = subp.PIPE,
                    shell = True,
                    )

    if print_err:
        #the formatting of the output is annoying - trying to make it look nice
        err_nl = out.stderr.decode('utf-8').split( '\n' )
        print(*err_nl, sep='\n')

def coordsort_bam( tempdir,
                    bamfile = 'temp.bam',
                    print_err = False ):

    out = subp.run( 'samtools sort %s > %s ' % ( tempdir + 'temp.sam', tempdir + bamfile ),
                    stderr = subp.PIPE,
                    shell = True,
                    )

    if print_err:
        #the formatting of the output is annoying - trying to make it look nice
        err_nl = out.stderr.decode('utf-8').split( '\n' )
        print(*err_nl, sep='\n')

def merge_bam( tempdir,
               bamfile1,
               bamfile2 = 'temp.bam',
               print_err = False ):

    out = subp.run( 'samtools merge -f --threads 8 %s %s %s' % ( tempdir + bamfile1, tempdir + bamfile1, tempdir + bamfile2 ),
                    stderr = subp.PIPE,
                    shell = True,
                    )

    if print_err:
        #the formatting of the output is annoying - trying to make it look nice
        err_nl = out.stderr.decode('utf-8').split( '\n' )
        print(*err_nl, sep='\n')

def align_sample( satbl,
                  refseq,
                  tempdir,
                  fq_file1,
                  fq_file2,
                  sample_name = 'sample',
                  print_err = False
                   ):

    var_to_bc_d = extract_snv_bcs( satbl )

    for i, var in enumerate( var_to_bc_d.keys() ):

        write_temp_fa( refseq, var, tempdir )

        read_names = extract_read_names( fq_file2, var_to_bc_d[ var ] )

        write_temp_fq( fq_file1, read_names, tempdir )

        build_index( tempdir, print_err = print_err )

        align_reads( tempdir, print_err = print_err )

        #set up the final bam file if its the first variant
        if i == 0:
            coordsort_bam( tempdir, bamfile = sample_name + '.bam', print_err = print_err)
            first_var = False

        #merge latest variant into final bam file after the first variant
        else:
            coordsort_bam( tempdir, print_err = print_err )
            merge_bam( tempdir, sample_name + '.bam', print_err = print_err )

        if i % 100 == 0:
            print( 'Processed %i variants: %.2f%% of total variants in sample' % ( i, 100*( i / len( var_to_bc_d.keys() ) ) ) )
