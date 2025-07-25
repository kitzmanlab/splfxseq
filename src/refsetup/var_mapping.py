import pandas as pd
import numpy as np

class CoordinateMapper:
    """
    A class to handle coordinate translations between genomic coordinates and cloned fragment coordinates.
    
    This class loads a coordinate mapping table and provides methods to convert positions between
    different coordinate systems (genomic, vector/cloned fragment, cDNA/HGVS).
    """
    
    def __init__(self, mapping_table=None):
        """
        Initialize the CoordinateMapper with a mapping table.
        
        Args:
            mapping_table (pd.DataFrame, optional): DataFrame containing coordinate mappings with columns:
                - chrom: Chromosome name
                - genomic_start: Start position in genomic coordinates (1-based)
                - genomic_end: End position in genomic coordinates
                - vector_start: Start position in vector coordinates (1-based)
                - vector_end: End position in vector coordinates
                - cdna_start: Start position in cDNA coordinates (1-based)
                - cdna_end: End position in cDNA coordinates
                - strand: Strand ('+' or '-')
                - exon_number: Exon number
        """
        self.mapping_table = mapping_table if mapping_table is not None else pd.DataFrame()
        self._validate_mapping_table()
        
    def _validate_mapping_table(self):
        """Validate the mapping table has required columns and correct data types."""
        required_cols = [
            'chrom', 'genomic_start', 'genomic_end', 'vector_start', 
            'vector_end', 'cdna_start', 'cdna_end', 'strand', 'exon_number'
        ]
        
        if not self.mapping_table.empty:
            missing_cols = [col for col in required_cols if col not in self.mapping_table.columns]
            if missing_cols:
                raise ValueError(f"Missing required columns: {missing_cols}")
            
            # Ensure numeric columns are numeric
            numeric_cols = ['genomic_start', 'genomic_end', 'vector_start', 
                          'vector_end', 'cdna_start', 'cdna_end', 'exon_number']
            for col in numeric_cols:
                self.mapping_table[col] = pd.to_numeric(self.mapping_table[col])
                
            # Validate strand values
            valid_strands = {'+', '-'}
            invalid_strands = set(self.mapping_table['strand']) - valid_strands
            if invalid_strands:
                raise ValueError(f"Invalid strand values: {invalid_strands}. Must be '+' or '-'")

    def load_mapping_table(self, mapping_table):
        """Load a new mapping table."""
        self.mapping_table = mapping_table
        self._validate_mapping_table()

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
        region = self.mapping_table[self.mapping_table['chrom'] == chrom]
        if region.empty:
            return None
            
        # Find relevant mapping region
        region = region[
            (region['genomic_start'] <= position) & 
            (region['genomic_end'] >= position)
        ]
        
        if region.empty:
            return None
            
        if len(region) > 1:
            raise ValueError(f"Position {position} maps to multiple regions")
            
        region = region.iloc[0]
        
        # If strand is provided, verify it matches
        if strand and strand != region['strand']:
            raise ValueError(f"Provided strand {strand} doesn't match mapping {region['strand']}")
            
        # Calculate relative position and convert
        if region['strand'] == '+':
            relative_pos = position - region['genomic_start']
            vector_pos = region['vector_start'] + relative_pos
        else:
            relative_pos = region['genomic_end'] - position
            vector_pos = region['vector_start'] + relative_pos
            
        return vector_pos

    def vector_to_genomic(self, position):
        """
        Convert vector coordinates to genomic coordinates.
        
        Args:
            position (int): Vector position (1-based)
            
        Returns:
            tuple: (chrom, position, strand) or None if position is outside mapped regions
        """
        if self.mapping_table.empty:
            raise ValueError("No mapping table loaded")
            
        # Find relevant mapping region
        region = self.mapping_table[
            (self.mapping_table['vector_start'] <= position) & 
            (self.mapping_table['vector_end'] >= position)
        ]
        
        if region.empty:
            return None
            
        if len(region) > 1:
            raise ValueError(f"Position {position} maps to multiple regions")
            
        region = region.iloc[0]
        
        # Calculate relative position and convert
        relative_pos = position - region['vector_start']
        
        if region['strand'] == '+':
            genomic_pos = region['genomic_start'] + relative_pos
        else:
            genomic_pos = region['genomic_end'] - relative_pos
            
        return (region['chrom'], genomic_pos, region['strand'])

    def vector_to_cdna(self, position):
        """
        Convert vector coordinates to cDNA coordinates.
        
        Args:
            position (int): Vector position (1-based)
            
        Returns:
            tuple: (cdna_position, exon_number) or None if position is outside mapped regions
        """
        if self.mapping_table.empty:
            raise ValueError("No mapping table loaded")
            
        # Find relevant mapping region
        region = self.mapping_table[
            (self.mapping_table['vector_start'] <= position) & 
            (self.mapping_table['vector_end'] >= position)
        ]
        
        if region.empty:
            return None
            
        if len(region) > 1:
            raise ValueError(f"Position {position} maps to multiple regions")
            
        region = region.iloc[0]
        
        # Calculate relative position and convert
        relative_pos = position - region['vector_start']
        cdna_pos = region['cdna_start'] + relative_pos
            
        return (cdna_pos, region['exon_number'])

    def cdna_to_vector(self, position):
        """
        Convert cDNA coordinates to vector coordinates.
        
        Args:
            position (int): cDNA position (1-based)
            
        Returns:
            tuple: (vector_position, exon_number) or None if position is outside mapped regions
        """
        if self.mapping_table.empty:
            raise ValueError("No mapping table loaded")
            
        # Find relevant mapping region
        region = self.mapping_table[
            (self.mapping_table['cdna_start'] <= position) & 
            (self.mapping_table['cdna_end'] >= position)
        ]
        
        if region.empty:
            return None
            
        if len(region) > 1:
            raise ValueError(f"Position {position} maps to multiple regions")
            
        region = region.iloc[0]
        
        # Calculate relative position and convert
        relative_pos = position - region['cdna_start']
        vector_pos = region['vector_start'] + relative_pos
            
        return (vector_pos, region['exon_number'])

    def get_exon_boundaries(self, coordinate_system='genomic'):
        """
        Get all exon boundaries in the specified coordinate system.
        
        Args:
            coordinate_system (str): One of 'genomic', 'vector', or 'cdna'
            
        Returns:
            list of tuples: [(start, end, exon_number), ...] sorted by position
        """
        if self.mapping_table.empty:
            raise ValueError("No mapping table loaded")
            
        if coordinate_system == 'genomic':
            start_col, end_col = 'genomic_start', 'genomic_end'
        elif coordinate_system == 'vector':
            start_col, end_col = 'vector_start', 'vector_end'
        elif coordinate_system == 'cdna':
            start_col, end_col = 'cdna_start', 'cdna_end'
        else:
            raise ValueError("coordinate_system must be one of: genomic, vector, cdna")
            
        boundaries = list(zip(
            self.mapping_table[start_col],
            self.mapping_table[end_col],
            self.mapping_table['exon_number']
        ))
        
        return sorted(boundaries)

    def is_exonic(self, position, coordinate_system='genomic'):
        """
        Check if a position is within an exon.
        
        Args:
            position (int): Position to check
            coordinate_system (str): One of 'genomic', 'vector', or 'cdna'
            
        Returns:
            tuple: (bool, exon_number) or (False, None) if intronic/outside
        """
        boundaries = self.get_exon_boundaries(coordinate_system)
        
        for start, end, exon_num in boundaries:
            if start <= position <= end:
                return (True, exon_num)
                
        return (False, None)