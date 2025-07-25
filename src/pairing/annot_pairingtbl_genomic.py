#!/usr/bin/env python

import pysam
import os.path as op
import argparse
import pandas as pd
from splshared.var_mapping import GenomePlasmidMapper

def parse_variant(variant_str):
    """Parse a variant string of format vector_chrom:pos:ref:mut"""
    try:
        vector_chrom, pos, ref, mut = variant_str.split(':')
        return vector_chrom, int(pos), ref, mut
    except ValueError:
        raise ValueError(f"Invalid variant format: {variant_str}. Expected format: vector_chrom:pos:ref:mut")

def convert_variant_list(variant_list_str, mapper):
    """Convert a comma-separated list of vector variants to genomic coordinates"""
    if pd.isna(variant_list_str) or not variant_list_str:
        return ""
        
    genomic_variants = []
    for variant in variant_list_str.split(','):
        variant = variant.strip()
        if not variant:
            continue
            
        vector_chrom, vector_pos, ref, mut = parse_variant(variant)
        genomic = mapper.vector_to_genomic(vector_chrom, vector_pos)
        if not genomic:
            genomic_variants.append(f"")
        else:
            chrom, pos, strand = genomic
            genomic_variants.append(f"{chrom}:{pos}:{ref}:{mut}:{strand}")
            
    return ','.join(genomic_variants)

def main():
    parser = argparse.ArgumentParser(description='Convert vector coordinates to genomic coordinates in variant list')
    
    parser.add_argument('--in_table', required=True,
                      help='Input table with variant_list column (tab-delimited)')
    parser.add_argument('--mapping_table', required=True,
                      help='Genome to vector mapping table (tab-delimited)')
    parser.add_argument('--out_table', required=True,
                      help='Output table path')
    parser.add_argument('--chunksize', type=int, default=10000,
                      help='Number of rows to process at a time (default: 10000)')
    
    parser.add_argument('--ref_fasta', required=True,
                      help='Reference fasta file')
    
    parser.add_argument('--plas_fasta', required=True,
                      help='Plasmid fasta file')

    args = parser.parse_args()
    
    # Load mapping table and create mapper
    mapper = GenomePlasmidMapper.from_file(args.mapping_table)

    ref_fasta = pysam.FastaFile(args.ref_fasta)
    plas_fasta = pysam.FastaFile(args.plas_fasta)

    mapper.check_sequences(ref_fasta, plas_fasta)
    
    # Process input table in chunks
    first_chunk = True
    for chunk in pd.read_csv(args.in_table, sep='\t', chunksize=args.chunksize):
        if 'variant_list' not in chunk.columns and first_chunk:
            raise ValueError("Input table must contain 'variant_list' column")
            
        # Convert variants for this chunk
        chunk['variant_list_genome'] = chunk['variant_list'].apply(
            lambda x: convert_variant_list(x, mapper)
        )
        
        # put variant_list_genome column immedaitely after variant_list column, 
        loc = [cn for cn in list(chunk.columns) if cn != 'variant_list_genome']
        loc[ loc.index('variant_list'):loc.index('variant_list')+1 ] = ['variant_list', 'variant_list_genome']
        chunk = chunk[loc]

        # Write chunk to output
        mode = 'w' if first_chunk else 'a'
        chunk.to_csv(args.out_table, sep='\t', index=False, mode=mode, header=first_chunk)
        first_chunk = False

if __name__ == '__main__':
    main() 
