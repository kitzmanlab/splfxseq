import pysam
import os.path as op
import argparse
import pandas as pd
from splshared.ssshared import revComp, gene_model_tbl_checker
from splshared.var_mapping import GenomePlasmidMapper, GeneModelMapper, VectorExonTable

up_ex_std = "CAAGCTGCACGCTCTAGAGTCGACCCAGCA"
dn_ex_std = "ACCTGGAGATCTCCCGAGGGGACCCGACAG"

def main():

    opts = argparse.ArgumentParser( description='manage mapping between genomic and vector coordinates' )

    opts.add_argument('--genomic_exon_tbl', required=True, help='table of exons with genomic coordinates and strand', dest='genomic_exon_tbl')
    opts.add_argument('--cloned_frag_tbl', required=True, help='table with cloned fragments\' genomic and vector coordinates', dest='cloned_frag_tbl')
    opts.add_argument('--fasta_plasmid', required=True, help='plasmid fasta file', dest='fasta_plasmid')
    opts.add_argument('--fasta_genome', required=True, help='genome fasta file', dest='fasta_genome')
    opts.add_argument('--out_vec_exon_tbl', required=True, help='output table of vector exon locations', dest='out_vec_exon_tbl')

    opts.add_argument('--custom_up_ex_seq', default=None, help='custom upstream exon sequence', dest='custom_up_ex_seq')
    opts.add_argument('--custom_dn_ex_seq', default=None, help='custom downstream exon sequence', dest='custom_dn_ex_seq')

    o = opts.parse_args()

    # create the mapper
    mapper = GenomePlasmidMapper.from_file(o.cloned_frag_tbl)
    plas_fasta=pysam.FastaFile(o.fasta_plasmid)
    mapper.check_sequences(ref_fasta=pysam.FastaFile(o.fasta_genome), plas_fasta=plas_fasta)
    
    # 
    genomic_exon_tbl = GeneModelMapper.from_file(o.genomic_exon_tbl)
    vector_exon_tbl = VectorExonTable()

    for i, r in genomic_exon_tbl.tbl.iterrows():
        res = mapper.genomic_to_vector(r['genome_chrom'], r['genome_start'], r['genome_strand'])
        assert res is not None, f"no vector coordinates found for {r['genome_chrom']}:{r['genome_start']}-{r['genome_end']}"
        vec_chrom, vec_start = res
        res = mapper.genomic_to_vector(r['genome_chrom'], r['genome_end'], r['genome_strand'])
        assert res is not None, f"no vector coordinates found for {r['genome_chrom']}:{r['genome_start']}-{r['genome_end']}"
        _, vec_end = res
        vector_exon_tbl.tbl['vector_chrom'].append(vec_chrom)
        vector_exon_tbl.tbl['vector_start'].append(vec_start)
        vector_exon_tbl.tbl['vector_end'].append(vec_end)

        for cn in genomic_exon_tbl.tbl.columns:
            if cn in vector_exon_tbl.tbl:
                vector_exon_tbl.tbl[cn].append(r[cn])
        
    for plas_chrom in set( vector_exon_tbl.tbl['vector_chrom'] ):
        plas_seq = plas_fasta.fetch(plas_chrom)

        # look for constant up / dn exon sequences in vector
        up_ex_pos = plas_seq.find(up_ex_std)
        dn_ex_pos = plas_seq.find(dn_ex_std)
        assert up_ex_pos != -1, f"upstream exon sequence {up_ex_std} not found in plasmid"
        assert dn_ex_pos != -1, f"downstream exon sequence {dn_ex_std} not found in plasmid"
        
        up_ex = o.custom_up_ex_seq if o.custom_up_ex_seq else up_ex_std
        dn_ex = o.custom_dn_ex_seq if o.custom_dn_ex_seq else dn_ex_std

        up_ex_pos = plas_seq.find(up_ex)
        dn_ex_pos = plas_seq.find(dn_ex)

        assert up_ex_pos != -1, f"upstream exon sequence {up_ex} not found in plasmid"
        assert dn_ex_pos != -1, f"downstream exon sequence {dn_ex} not found in plasmid"

        vector_exon_tbl.tbl['vector_chrom'].append(plas_chrom)
        vector_exon_tbl.tbl['vector_start'].append(up_ex_pos+1)
        vector_exon_tbl.tbl['vector_end'].append(up_ex_pos+1+len(up_ex)-1)
        vector_exon_tbl.tbl['exon_name'].append('constup')
        
        vector_exon_tbl.tbl['vector_chrom'].append(plas_chrom)
        vector_exon_tbl.tbl['vector_start'].append(dn_ex_pos+1)
        vector_exon_tbl.tbl['vector_end'].append(dn_ex_pos+1+len(dn_ex)-1)
        vector_exon_tbl.tbl['exon_name'].append('constdn')

        # set entries for the genomic columns to null since there is no genomic correlate of the const exons
        for cn in vector_exon_tbl.tbl: 
            if cn not in 'vector_chrom,vector_start,vector_end,exon_name'.split(','):
                vector_exon_tbl.tbl[cn].append(None)
                vector_exon_tbl.tbl[cn].append(None)

    vector_exon_tbl.tbl = pd.DataFrame(vector_exon_tbl.tbl)

    vector_exon_tbl.tbl['genome_start']=[ str(int(x)) if not pd.isnull(x) else "" for x in vector_exon_tbl.tbl['genome_start']  ]
    vector_exon_tbl.tbl['genome_end']=[ str(int(x)) if not pd.isnull(x) else "" for x in vector_exon_tbl.tbl['genome_end']  ]

    vector_exon_tbl.tbl.to_csv(o.out_vec_exon_tbl, sep='\t', index=False)

if __name__ == '__main__':
    main()

