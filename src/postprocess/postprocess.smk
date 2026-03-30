######################################################
# Configuration values you must specify:
#
#
# > lib_table: varfx report, output from var_effect_indiv_samps.smk  
#
#   one row per sequencing library to be processed

#   required columns
#   - libname   (this nust be uniq; if you are running the same lib across multiple batches you need to 
#   - var_rpt   (path to var_rpt table for the lib )
#   - var_rpt_snvonly (path to SNV only var_rpt for lib)
#   - var_bc_rpt (path to per var bcrpt for lib)
#   - exon      (exon/ref to which the lib is aligned and for which you want aggregated stats to be generated)

# > map_table: a table for mapping between genome, vector and cdna coordinates, one row per ref/exon
# > outdir
# 
#######################################################
#  optional config values:
#
#   - num_bc_per_var    (minimum number of singleton only bc required to support each variant, per sample, defaults to 5)
#   - cumd_cutoff       (cutoff for cumulative distribution i.e. first x amount of reads for a given sample, defaults to 0.90)
#   - intronic_bp       (number of basepair from exon boundary beyond which we define as "null" in test for stat sig, defaults to 10)
#   - splai_table         (path to a splai table to merge with summary table, must contain a chrom, pos, ref, alt columns)
#   - clinvar_table       (path to a clinvar table to merge with summary table, must contain chr, pos, ref, alt columns)
#   - dms_table          (path to a dms table to merge with summary table, must contain Variant and 'LOF score' columns, it will then get protein variants for each bp change)
#   - gnomad_table         (path to a gnomad table to merge with summary table, must contain CHROM, POS, REF, ALT columns)

#######################################################

import os.path as op
import os
import pandas as pd
import pysam 
import Bio.Seq
import post_process_stats

import altair as alt 

########

# localrules: outtbl

assert 'lib_table' in config, 'must specify sample table'


PREFIX = config['prefix'] if 'prefix' in config  else ''


NUM_BC_PER_VAR = config['num_bc_per_var'] if 'num_bc_per_var' in config else 5
CUMD_CUTOFF = config['cumd_cutoff'] if 'cumd_cutoff' in config else 0.90
INTRONIC_BP = config['intronic_bp'] if 'intronic_bp' in config else 10

tbl_args = ''
if config.get('splai_table'):
    tbl_args += f'--splai_tbl {config["splai_table"]} '  
if config.get('clinvar_table'):   
    tbl_args += f'--clinvar_tbl {config["clinvar_table"]} ' 
if config.get('dms_table'):       
    tbl_args += f'--dms_tbl {config["dms_table"]} ' 
if config.get('gnomad_table'):     
    tbl_args +=  f'--gnomad_tbl {config["gnomad_table"]} '   

########
#
# load and check sample table.

l_reqd_cols = [ 'libname', 'var_rpt', 'var_rpt_snvonly', 'var_bc_rpt', 'exon' ]
tbl_lib = pd.read_table( config['lib_table'] )
assert all( [ col in tbl_lib.columns for col in l_reqd_cols ] ), 'lib table must have columns: '+','.join(l_reqd_cols)
assert len(set(tbl_lib['libname'])) == tbl_lib.shape[0], 'all libname entries must be unique'
  
lLibs = tbl_lib['libname'].unique()
lExons = tbl_lib['exon'].unique()

########
#
# load and check map table

lreqd_cols_map_tbl = ["genome_chrom", "vector_chrom", "gene_model_id", "exon", "vector_genome_start", "vector_genome_end", "vector_start", "vector_end", "ex_genome_start", "ex_genome_end", "cdna_start", "genome_strand"]
tbl_map = pd.read_table( config['map_table'] )
assert all( [ col in tbl_map.columns for col in lreqd_cols_map_tbl ] ), 'map table must have columns: '+','.join(lreqd_cols_map_tbl)
assert set(lExons).issubset(set(tbl_map['exon'].unique())), "Not all exons from lib table are present in mapping table" 

# expected output files

assert 'outdir' in config, 'must specify output directory'
OUT_DIR = config['outdir']
# outputs

l_out_combined_long_rpt = expand('{}/perexon_rpts/long/{}{{exon}}.long.txt'.format(OUT_DIR,PREFIX), exon=lExons)
l_out_combined_wide_rpt = expand('{}/perexon_rpts/wide/{}{{exon}}.wide.txt'.format(OUT_DIR,PREFIX), exon=lExons)
l_out_summary  = expand('{}/perexon_rpts/summary/{}{{exon}}.summary.txt'.format(OUT_DIR,PREFIX), exon=lExons)

outRpt = OUT_DIR+'/'+PREFIX+'summary_rpt.txt'
lOutFiles = [outRpt] + l_out_combined_long_rpt + l_out_combined_wide_rpt + l_out_summary

########
# ensure rerun if the file paths change
def exon_sample_inputs(wc):
    rows = tbl_lib.loc[tbl_lib["exon"] == wc.exon]
    return rows["var_rpt_snvonly"].tolist() + rows["var_bc_rpt"].tolist()

########

rule all:
    input:
        lOutFiles

rule per_exon_process:
    input:
        varfx_tbl=config["lib_table"],
        map_tbl=config["map_table"],
        sample_files=exon_sample_inputs,
    output:
        long_rpt = op.join(OUT_DIR, 'perexon_rpts/long/'+PREFIX+'{exon}.long.txt'),
        wide_rpt = op.join(OUT_DIR, 'perexon_rpts/wide/'+PREFIX+'{exon}.wide.txt'),
        summary_rpt = op.join(OUT_DIR, 'perexon_rpts/summary/'+PREFIX+'{exon}.summary.txt'),
    threads: 1
    params: tbl_args=tbl_args
    resources:
        mem_per_cpu="10gb", 
        cpus="1", 
        time="0:30:00",
    run:


        shell("""
            post_process_stats \
            --varfx_tbl {input.varfx_tbl} \
            --map_tbl {input.map_tbl} \
            {params.tbl_args} \
            --exon {wildcards.exon} \
            --num_bc_per_var {NUM_BC_PER_VAR} \
            --cumd_cutoff {CUMD_CUTOFF} \
            --intronic_bp {INTRONIC_BP} \
            --out_long {output.long_rpt} \
            --out_wide {output.wide_rpt} \
            --out_summary {output.summary_rpt}
        """)


rule out_rpt:
    input:
        long_rpt = expand(rules.per_exon_process.output.long_rpt, exon=lExons),
        wide_rpt = expand(rules.per_exon_process.output.wide_rpt, exon=lExons),
        summary_rpt = expand(rules.per_exon_process.output.summary_rpt, exon=lExons),
    output:
        out_rpt = op.join(OUT_DIR,PREFIX+'summary_rpt.txt'),
    run:
        out_rpt = {k:[] for k in 'exon,long_rpt,wide_rpt,summary_rpt'.split(',')}
        out_rpt['exon'] = lExons
        out_rpt['long_rpt'] = list(input.long_rpt)
        out_rpt['wide_rpt'] = list(input.wide_rpt)
        out_rpt['summary_rpt'] = list(input.summary_rpt)
        out_rpt = pd.DataFrame(out_rpt)
        out_rpt.to_csv(output.out_rpt, sep='\t', index=False)
