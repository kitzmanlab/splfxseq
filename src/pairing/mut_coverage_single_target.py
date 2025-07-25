#!/usr/bin/env python3

import pandas as pd
import altair as alt
from typing import List, Tuple, Dict, Set, Callable
from itertools import product
import argparse
from pathlib import Path
import pysam
from splshared.var_mapping import GeneModelMapper

import numpy as np
import matplotlib.pyplot as plt

def parse_variant(variant_str: str) -> Tuple[str, int, str, str, str]:
    """Parse a variant string in the format 'chrom:pos:ref:alt:strand'."""
    chrom, pos, ref, alt, strand = variant_str.split(':')
    return chrom, int(pos), ref, alt, strand

def get_variants_from_list(variant_list: str) -> List[Tuple[str, int, str, str, str]]:
    """Parse a comma-separated list of variants."""
    if pd.isna(variant_list) or not variant_list:
        return []
    return [parse_variant(v.strip()) for v in variant_list.split(',') if len(v.strip())>0]

def generate_single_base_variants(fasta: pysam.FastaFile, chrom: str, start_pos: int, end_pos: int) -> List[Tuple[str, int, str, str]]:
    """
    Generate all possible single base variants within a given interval using the reference genome.
    
    Args:
        fasta: pysam.FastaFile object for the reference genome
        chrom: Chromosome to analyze
        start_pos: Start position of interval (1-based)
        end_pos: End position of interval (1-based)
    
    Returns:
        List of tuples (chrom, pos, ref, alt) for all possible SNVs
    """
    bases = set(['A', 'C', 'G', 'T'])
    variants = []
    
    # Get the reference sequence for the interval
    # pysam uses 0-based coordinates, so subtract 1 from start_pos
    ref_seq = fasta.fetch(chrom, start_pos - 1, end_pos).upper()
    
    for i, ref_base in enumerate(ref_seq):
        pos = start_pos + i
        # Skip if reference base is not a standard base
        if ref_base not in bases:
            continue
        
        # Generate all possible alternate alleles (excluding the reference)
        alt_bases = bases - {ref_base}
        for alt in alt_bases:
            variants.append((chrom, pos, ref_base, alt))
    
    return variants

def count_variants(df: pd.DataFrame, 
                  fasta: pysam.FastaFile,
                  variant_col: str,
                  chrom: str,
                  start_pos: int,
                  end_pos: int,
                  fxn_allow_2ndvar: List[ Callable[[str,int,str,str], bool]] = [lambda chrom,pos,ref,alt: False] ) -> pd.DataFrame:
    """
    Count variants in both single-variant and any-context scenarios.
    
    Args:
        df: Input DataFrame with variant lists
        fasta: pysam.FastaFile object for the reference genome
        variant_col: Name of column containing variant lists
        chrom: Chromosome to analyze
        start_pos: Start position of interval
        end_pos: End position of interval
        fxn_allow_2ndvar: List of lambda functions(chrom,pos,ref,alt) to allow 2nd variant in any context
        
    Returns:
        DataFrame with variant counts
    """
    # Generate all possible variants in the interval
    all_variants = generate_single_base_variants(fasta, chrom, start_pos, end_pos)
    
    # Initialize count dictionaries
    single_counts = {(c, p, r, a): 0 for c, p, r, a in all_variants}
    any_counts = single_counts.copy()

    strand = None
    
    # Count variants in the data
    for i,r in df.iterrows():
        variant_list = r[variant_col]
        variants = get_variants_from_list(variant_list)

        # Skip empty variant lists
        if not variants:
            continue

        l_nonignored = [ v for v in variants if not any(
            [ fxn(v[0],v[1],v[2],v[3]) for fxn in fxn_allow_2ndvar])]

        for v in variants:
            if not strand:
                strand = v[4]
            else:
                if v[4] != strand:
                    raise ValueError(f"Variant {v} has different strand than previous variants")
            
            key = (v[0], v[1], v[2], v[3])
            if key in any_counts:
                any_counts[key] += 1
        
                if len(l_nonignored) == 1:
                    single_counts[key] += 1        
            
    # Create output DataFrame
    results = []
    for (c, p, r, a) in all_variants:
        results.append({
            'chrom': c,
            'position': p,
            'ref': r,
            'alt': a,
            'single_variant_count': single_counts[(c, p, r, a)],
            'any_context_count': any_counts[(c, p, r, a)]
        })

    results = pd.DataFrame(results)

    results['strand'] = strand
    
    return results

def plot_variant_counts_base_stacked(df: pd.DataFrame,
                                    chrom: str,
                                    start_pos: int,
                                    end_pos: int,
                                    gene_model_mapper: GeneModelMapper = None) -> alt.Chart:
    
    """Create a bar plot of variant counts by position using Altair."""
    # Melt the DataFrame to get counts in a single column
    melted = df.melt(
        id_vars=['chrom', 'position', 'ref', 'alt', 'strand'],
        value_vars=['single_variant_count', 'any_context_count'],
        var_name='count_type',
        value_name='count'
    )

    lc = []
    
    for count_type in melted['count_type'].unique():
        mcur = melted.loc[melted['count_type'] == count_type]
    
        c1 = alt.Chart(mcur).mark_bar(
            width=1,
            binSpacing=0,
        ).encode(
            x=alt.X('position:Q', title='Position').scale( reverse=melted.iloc[0]['strand'] == '-' ).scale(zero=False),
            y=alt.Y('count:Q', title='# barcodes'),
            color=alt.Color('alt:N', title='Alternate Base'),
            tooltip=['chrom', 'position', 'ref', 'alt', 'count_type', 'count', 'strand']
        ).properties(
            width=1200,
            height=250,
            title='All SNVs, #unique barcodes per position'
        )
    
        if gene_model_mapper:
            models_of_interest = gene_model_mapper.gene_model.loc[
                (gene_model_mapper.gene_model['genome_chrom'] == chrom) &
                (gene_model_mapper.gene_model['genome_start'] <= end_pos) &
                (gene_model_mapper.gene_model['genome_end'] >= start_pos)
            ]

            exonbox = alt.Chart(models_of_interest).mark_rect(
                fill='gray'
            ).encode(
                x=alt.X('genome_start:Q').scale(zero=False),
                x2=alt.X2('genome_end:Q'),
                y=alt.YValue(0),
                y2=alt.Y2Value(50),
            )

            c1 =  c1 + exonbox
        
        lc.append(c1)

    chartf = alt.vconcat(*lc)
        
    return chartf

def plot_variant_counts_base_row(df: pd.DataFrame,
                                    chrom: str,
                                    start_pos: int,
                                    end_pos: int,
                                    gene_model_mapper: GeneModelMapper = None) -> alt.Chart:
    
    """Create a bar plot of variant counts by position using Altair."""
    # Melt the DataFrame to get counts in a single column
    melted = df.melt(
        id_vars=['chrom', 'position', 'ref', 'alt', 'strand'],
        value_vars=['single_variant_count', 'any_context_count'],
        var_name='count_type',
        value_name='count'
    )

    lc = []
    
    for altbase in sorted(melted['alt'].unique()):
        mcur = melted.loc[melted['alt'] == altbase]
    
        c1 = alt.Chart(mcur).mark_bar(
            width=1,
            binSpacing=0,
        ).encode(
            x=alt.X('position:Q', title='Position').scale( reverse=melted.iloc[0]['strand'] == '-' ).scale(zero=False),
            y=alt.Y('count:Q', title='# barcodes'),
            color=alt.Color('count_type:N', title='Count Type'),
            tooltip=['chrom', 'position', 'ref', 'alt', 'count_type', 'count', 'strand']
        ).properties(
            width=1200,
            height=250,
            title=f'#barcodes with SNV mut={altbase}, by position'
        )
    
        if gene_model_mapper:
            models_of_interest = gene_model_mapper.gene_model.loc[
                (gene_model_mapper.gene_model['genome_chrom'] == chrom) &
                (gene_model_mapper.gene_model['genome_start'] <= end_pos) &
                (gene_model_mapper.gene_model['genome_end'] >= start_pos)
            ]

            exonbox = alt.Chart(models_of_interest).mark_rect(
                fill='gray'
            ).encode(
                x=alt.X('genome_start:Q').scale(zero=False),
                x2=alt.X2('genome_end:Q'),
                y=alt.YValue(0),
                y2=alt.Y2Value(50),
            )

            c1 =  c1 + exonbox
        
        lc.append(c1)

    chartf = alt.vconcat(*lc)
        
    return chartf

def plot_mutfreq_waterfall(df: pd.DataFrame,
                                    chrom: str,
                                    start_pos: int,
                                    end_pos: int,
                            **kwargs):

    """Create a bar plot of variant counts by position using Altair."""
    # Melt the DataFrame to get counts in a single column
    f, ax = plt.subplots(1,1,**kwargs)

    df = df.loc[ (df['chrom'] == chrom) & (df['position'] >= start_pos) & (df['position'] <= end_pos) ]

    v = np.array(df['single_variant_count'])
    v = v[ np.argsort(-v) ]
    plt.scatter( y=v,
                 x=np.arange( v.shape[0] ),
                 label='SNV, alone',
                 s=4,
    )

    v = np.array(df['any_context_count'])
    v = v[ np.argsort(-v) ]
    plt.scatter( y=v,
                 x=np.arange( v.shape[0] ),
                 label='SNV, allowing for other variants',
                 s=4,
    )

    plt.ylabel('#barcodes/variant')
    plt.xlabel('variant rank')
    plt.legend()

    return plt.gcf()   
    


def main():
    parser = argparse.ArgumentParser(description='Analyze variant coverage in a genomic interval')
    parser.add_argument('--ref_fasta', type=str, required=True,
                      help='Reference genome FASTA file path')
    parser.add_argument('--input_file', type=str, required=True,
                      help='Input table file path')
    parser.add_argument('--variant_col', type=str, default="variant_list_genome",
                      help='Name of column containing variant lists')

    parser.add_argument('--chrom', type=str, required=True,
                      help='Chromosome to analyze')
    parser.add_argument('--cloned_start_pos', type=int, required=True,
                      help='Start position of interval')
    parser.add_argument('--cloned_end_pos', type=int, required=True,
                      help='End position of interval')

    parser.add_argument('--mut_targ_start_pos', type=int, required=True,
                      help='Start position of interval')
    parser.add_argument('--mut_targ_end_pos', type=int, required=True,
                      help='End position of interval')

    parser.add_argument('--gene_model_file', type=str, default=None,
                      help='Gene model file path')
    parser.add_argument('--output_table', type=str, default=None,
                      help='Output table file path')
    parser.add_argument('--output_plot_base', type=str, default=None,
                      help='Output plot file path')
    parser.add_argument('--lfxn_allow_2ndvar', type=str, default=None,
                      help='List of lambda functions(chrom,pos,ref,alt) to allow 2nd variant in any context')
    args = parser.parse_args()
    
    # Open reference genome
    fasta = pysam.FastaFile(args.ref_fasta)
    
    # Read input data
    df = pd.read_table(args.input_file)
    
    if args.lfxn_allow_2ndvar:
        l_fxn_allow_2ndvar = eval(args.lfxn_allow_2ndvar)
    else:
        l_fxn_allow_2ndvar = []

    if args.gene_model_file:
        gene_model_mapper = GeneModelMapper.from_file(args.gene_model_file)
    else:
        gene_model_mapper = None

    # Process variants
    result_df = count_variants(
        df,
        fasta,
        args.variant_col,
        args.chrom,
        args.cloned_start_pos,
        args.cloned_end_pos,
        l_fxn_allow_2ndvar,
    )
    
    # Save results
    if args.output_table:
        result_df.to_csv(args.output_table, sep='\t', index=False)
    
    # Create and save plot
    if args.output_plot_base:
        chart = plot_variant_counts_base_stacked(result_df, args.chrom, args.cloned_start_pos, args.cloned_end_pos, gene_model_mapper)
        chart.save(f'{args.output_plot_base}.single_any.html')
    
        chart = plot_variant_counts_base_row(result_df, args.chrom, args.cloned_start_pos, args.cloned_end_pos, gene_model_mapper)
        chart.save(f'{args.output_plot_base}.bybase.html')

        wfall = plot_mutfreq_waterfall(result_df, args.chrom, args.mut_targ_start_pos, args.mut_targ_end_pos, figsize=(8,4))
        wfall.savefig(f'{args.output_plot_base}.waterfall.png')

    # Clean up
    fasta.close()

if __name__ == '__main__':
    main() 