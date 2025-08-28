import cmd
import os.path as op
import argparse
import pandas as pd
import numpy as np
from typing import Tuple, Dict

import pysam
import Bio.Seq

from splshared.ssshared import *

class GeneModelMapper:
    
    @staticmethod
    def from_file(gene_model_file):
        gene_model = pd.read_csv(gene_model_file, sep='\t')
        return GeneModelMapper(gene_model)

    def __init__(self, tbl=None):
        if tbl is None:
            self.tbl = gene_model_tbl_checker.new_odict_of_lists()
        else:
            self.tbl = tbl
            gene_model_tbl_checker.check(tbl)       

        
    # to implement - convert frmo c.XXX to genomic coordinates 

class VectorExonTable:
    @staticmethod
    def from_file(vector_exon_file):
        vector_exon_tbl = pd.read_csv(vector_exon_file, sep='\t')
        return VectorExonTable(vector_exon_tbl)

    def __init__(self, tbl=None):
        if tbl is None:
            self.tbl = vector_exon_tbl_checker.new_odict_of_lists()
        else:
            self.tbl = tbl
            vector_exon_tbl_checker.check(tbl)

    def get_const_exon_positions(
        self, 
    ) -> Dict[str,Tuple[Tuple[int,int],Tuple[int,int]]]:

        mvec_constupdn = {}

        for vecname in self.tbl['vector_chrom'].unique():
            li_up = (self.tbl['vector_chrom'] == vecname) & (self.tbl['exon_name'] == 'constup') 
            li_dn = (self.tbl['vector_chrom'] == vecname) & (self.tbl['exon_name'] == 'constdn') 

            if li_up.sum()!=1: raise ValueError(f"Expected 1 constup exon for {vecname}, found {li_up.sum()}")
            if li_dn.sum()!=1: raise ValueError(f"Expected 1 constdn exon for {vecname}, found {li_dn.sum()}")

            mvec_constupdn[ vecname ] = ( 
                            ( self.tbl.loc[ li_up ].iloc[0][ 'vector_start'], self.tbl.loc[ li_up ].iloc[0][ 'vector_end'] ),
                            ( self.tbl.loc[ li_dn ].iloc[0][ 'vector_start'], self.tbl.loc[ li_dn ].iloc[0][ 'vector_end'] ),
            )

        return mvec_constupdn
    
    def check_known_exon(self, vecname, exon_start, exon_end):
        """
        Check if a given exon is known in the vector exon table.
        
        Args:
            vecname (str): Vector chromosome name
            exon_start (int): Start position of the exon (1-based)
            exon_end (int): End position of the exon (1-based)
            
        Returns:
            bool: True if exon is known, False otherwise
        """
        li_ex = (self.tbl['vector_chrom'] == vecname) & (self.tbl['vector_start'] == exon_start) & (self.tbl['vector_end'] == exon_end)
        if li_ex.sum() == 0:
            return None
        elif li_ex.sum() == 1:
            return self.tbl.loc[ li_ex ].iloc[0]
        else:
            raise ValueError(f"Multiple exons found for {vecname}:{exon_start}-{exon_end}")

class GenomePlasmidMapper:
    """
    A class to handle coordinate translations between genomic coordinates and cloned fragment coordinates.
    
    This class loads a coordinate mapping table and provides methods to convert positions between
    different coordinate systems (genomic, vector/cloned fragment, cDNA/HGVS).
    """
    
    @staticmethod
    def from_file(mapping_file):
        """
        Load mapping table from a TSV file and create a new GenomePlasmidMapper instance.
        
        Args:
            mapping_file (str): Path to TSV file containing coordinate mappings
            
        Returns:
            GenomePlasmidMapper: New mapper instance initialized with the loaded table
            
        Raises:
            FileNotFoundError: If mapping_file does not exist
            ValueError: If mapping table format is invalid
        """
        mapping_table = pd.read_csv(mapping_file, sep='\t')
        return GenomePlasmidMapper(mapping_table)

    def __init__(self, mapping_table):
        """
        Initialize the CoordinateMapper with a mapping table.
        
        Args:
            mapping_table (pd.DataFrame, optional): DataFrame containing coordinate mappings with columns:
                - chrom: Chromosome name
                - genome_start: Start position in genomic coordinates (1-based)
                - genome_end: End position in genomic coordinates
                - vector_start: Start position in vector coordinates (1-based)
                - vector_end: End position in vector coordinates
                - cdna_start: Start position in cDNA coordinates (1-based)
                - cdna_end: End position in cDNA coordinates
                - genome_strand: Strand ('+' or '-')
                - exon_number: Exon number
        """

        if mapping_table is None:
            raise ValueError("mapping_table is required")
        
        self.mapping_table = mapping_table
        cloned_fragments_tbl_checker.check(self.mapping_table)
        self._validate_mapping_table()
        
    def _validate_mapping_table(self):
        valid_strands = {'+', '-'}
        invalid_strands = set(self.mapping_table['genome_strand']) - valid_strands
        if invalid_strands:
            raise ValueError(f"Invalid strand values: {invalid_strands}. Must be '+' or '-'")

    def check_sequences(self, 
                        ref_fasta: pysam.FastaFile,
                        plas_fasta: pysam.FastaFile, ):
        """
        check that the sequences in the mapping table are consistent with the sequences in the reference fasta files
        """

        for _, row in self.mapping_table.iterrows():
            ref_seq = plas_fasta.fetch(row['vector_chrom'],row['vector_start']-1, row['vector_end'])
            if row['genome_strand'] == '-' :
                ref_seq = str(Bio.Seq.Seq(ref_seq).reverse_complement())
            plas_seq = ref_fasta.fetch(row['genome_chrom'],row['genome_start']-1, row['genome_end'])
            ref_seq = ref_seq.upper()
            plas_seq = plas_seq.upper()
            if ref_seq != plas_seq:
                raise ValueError(f"Sequence mismatch at {row['vector_chrom']}:{row['vector_start']}-{row['vector_end']}")

    def genomic_to_vector(self, chrom, position, strand=None):
        """
        Convert genomic coordinates to vector coordinates.
        
        Args:
            chrom (str): Chromosome name
            position (int): Genomic position (1-based)
            strand (str, optional): Strand ('+' or '-'). If not provided, will use mapping table.
            
        Returns:
            int: Vector coordinate (1-based) or None if position is outside mapped regions
        """
        if self.mapping_table.empty:
            raise ValueError("No mapping table loaded")
            
        # Filter for chromosome
        region = self.mapping_table[self.mapping_table['genome_chrom'] == chrom]
        if region.empty:
            return None
            
        # Find relevant mapping region
        region = region[
            (region['genome_start'] <= position) & 
            (region['genome_end'] >= position)
        ]
        
        if region.empty:
            return None
            
        if len(region) > 1:
            raise ValueError(f"{chrom}:{position} maps to multiple regions")
            
        region = region.iloc[0]
        
        # If strand is provided, verify it matches
        if strand and strand != region['genome_strand']:
            raise ValueError(f"Provided strand {strand} doesn't match mapping {region['genome_strand']}")
            
        # Calculate relative position and convert
        if region['genome_strand'] == '+':
            relative_pos = position - region['genome_start']
            vector_pos = region['vector_start'] + relative_pos
        else:
            relative_pos = region['genome_end'] - position
            vector_pos = region['vector_start'] + relative_pos
            
        return ( region['vector_chrom'], vector_pos )

    def vector_to_genomic(self, vector_chrom, vector_position):
        """
        Convert vector coordinates to genomic coordinates.
        
        Args:
            vector_chrom (str): Vector chromosome name
            vector_position (int): Vector position (1-based)
            
        Returns:
            tuple: (chrom, position, strand) or None if position is outside mapped regions
        """
        if self.mapping_table.empty:
            raise ValueError("No mapping table loaded")
            
        # Find relevant mapping region
        region = self.mapping_table[
            (self.mapping_table['vector_chrom'] == vector_chrom) & 
            (self.mapping_table['vector_start'] <= vector_position) & 
            (self.mapping_table['vector_end'] >= vector_position)
        ]
        
        if region.empty:
            # raise ValueError(f"Vector {vector_chrom}:{vector_position} not found in this map")
            return None
            
        if len(region) > 1:
            raise ValueError(f"Vector {vector_chrom}:{vector_position} maps to multiple regions")
            
        region = region.iloc[0]
        
        # Calculate relative position and convert
        relative_pos = vector_position - region['vector_start']
        
        if region['genome_strand'] == '+':
            genomic_pos = region['genome_start'] + relative_pos
        else:
            genomic_pos = region['genome_end'] - relative_pos
            
        return (region['genome_chrom'], genomic_pos, region['genome_strand'])

