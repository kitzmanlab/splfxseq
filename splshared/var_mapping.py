import os.path as op
import argparse
import pandas as pd
import numpy as np


from splshared.ssshared import *

class GeneModelMapper:
    
    @staticmethod
    def from_file(gene_model_file):
        gene_model = pd.read_csv(gene_model_file, sep='\t')
        return GeneModelMapper(gene_model)

    def __init__(self, gene_model):
        gene_model_tbl_checker.check(gene_model)
        self.gene_model = gene_model

    # to implement - convert frmo c.XXX to genomic coordinates 

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


# def main():

#     opts = argparse.ArgumentParser()

#     opts.add_argument('--in_align_tbl', dest='in_align_tbl',
#         help='table of subassembly reads aligned to reference', required=True)

#     # opts.add_argument('--ref_fasta', dest='ref_fasta',
#     #     required=True)

#     # opts.add_argument('--col_libname', dest='col_libname',
#     #     help='column with library/sample name', required=True)

#     # opts.add_argument('--col_bam', dest='col_bam',
#     #     help='column with bam path', required=True)

#     # opts.add_argument('--region', dest='region',
#     #     help='chrom:start-stop to report over', required=True)

#     # opts.add_argument('--out_tbl', dest='out_tbl',
#     #     help='output table here', required=True)

#     o = opts.parse_args()

#     chrom, corng = o.region.split(':')[0], (int(o.region.split(':')[1].split('-')[0]),int(o.region.split(':')[1].split('-')[1]))


#     ref_fasta = pysam.Fastafile(o.ref_fasta)

#     samptbl = pd.read_csv(o.in_align_tbl,sep='\t')

#     hdrout = [
#         'chrom','pos', 'ref',
#         'sample',
#         'reads_all','reads_fwd','reads_rev',
#         'matches','matches_fwd','matches_rev',
#         'mismatches','mismatches_fwd','mismatches_rev',
#         'deletions','deletions_fwd','deletions_rev',
#         'insertions','insertions_fwd','insertions_rev',
#         'A','C','G','T',
#         'A_fwd','A_rev',
#         'C_fwd','C_rev',
#         'G_fwd','G_rev',
#         'T_fwd','T_rev',
#         'N_fwd','N_rev' ]

#     hdr_copydirectly = [
#         'chrom', 'pos', 
#         'reads_all','reads_fwd','reads_rev',
#         'matches','matches_fwd','matches_rev',
#         'mismatches','mismatches_fwd','mismatches_rev',
#         'deletions','deletions_fwd','deletions_rev',
#         'insertions','insertions_fwd','insertions_rev',
#         'A','C','G','T',
#         'A_fwd','A_rev',
#         'C_fwd','C_rev',
#         'G_fwd','G_rev',
#         'T_fwd','T_rev',
#         'N_fwd','N_rev' ]

#     tbl_out = { c:[] for c in hdrout }


#     for _, r in samptbl.iterrows():

#         bam = pysam.AlignmentFile( r[o.col_bam] )
#         libname = r[o.col_libname]

#         cvgcur=pss.load_variation_strand(
#               bam,
#               fafile=ref_fasta,
#               chrom=chrom,
#               start=corng[0]-1,
#               end=corng[1],
#               pad=True,
#               truncate=True,
#               max_depth=100000)
            
#         tbl_out['sample'].extend( [libname]*len(cvgcur['chrom']) )

#         reflist = cvgcur['ref']
#         for i in range(len(reflist)):
#             reflist[i] = reflist[i].decode()

#         for col in hdr_copydirectly:
#             tbl_out[col].extend( cvgcur[col].tolist() )

#         tbl_out['ref'].extend( reflist )

#     tbl_out = pd.DataFrame( tbl_out )

#     # make it 1 based
#     tbl_out['pos']+=1
#     # fix encoding
#     for col in ['chrom','ref']:
#         tbl_out[col]=tbl_out[col].str.decode('ascii')

#     tbl_out.to_csv(o.out_tbl,sep='\t',index=False)




# if __name__ == '__main__':
#     main()

