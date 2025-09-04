######################################################
# Configuration values you must specify:
#
#
# > lib_table  
#   one row per sequencing library to be processed
#   required columns
#   - libname   (this nust be uniq; if you are running the same lib across multiple batches you need to 
#   - bam  (path to barcode-sorted bam file)
#   - ref_seq_name (name of the reference sequence)


# vector_exon_tbl 
#
# TODO DOC
#
# > outdir
# 
#######################################################
#  optional config values:
#
#   prefix 
#
#   combine_rna_samp_options, additional options to pass to combine_rna_samp.py
#
#   process_rna_aligns_opts, additional options to pass to process_rna_aligns.py

#######################################################
# For example, to run:
#  

import os.path as op
import os
import pandas as pd
import pysam 
import Bio.Seq

import altair as alt 

########

localrules: outtbl

assert 'lib_table' in config, 'must specify sample table'

PREFIX = config['prefix'] if 'prefix' in config  else ''
OUT_DIR = config['outdir']

VECTOR_EXON_TBL = config['vector_exon_tbl']

if 'combine_rna_samp_options' in config:
    COMBINE_RNA_SAMP_OPTIONS = config['combine_rna_samp_options']
else:
    COMBINE_RNA_SAMP_OPTIONS = "  --otherisos_perbc_min_read_count 1 --otherisos_perbc_min_psi 0.025 --otherisos_minbc_withinallsamp 5 --otherisos_minbc_sumacrosssamp 5 "

if 'process_rna_aligns_opts' in config:
    PROCESS_RNA_ALIGNS_OPTS = config['process_rna_aligns_opts']
else:
    PROCESS_RNA_ALIGNS_OPTS = " --spl-tol 3 --indel-tol 5 --min-matches-for 20 --min-matches-rev 20 --max-soft-clip-for 10 --max-soft-clip-rev 10 "

########
#
# load and check sample table.

l_reqd_cols = [ 'libname', 'bam', 'ref_seq_name' ]
tbl_lib = pd.read_table( config['lib_table'] )
assert all( [ col in tbl_lib.columns for col in l_reqd_cols ] ), 'lib table must have columns: '+','.join(l_reqd_cols)
assert len(set(tbl_lib['libname'])) == tbl_lib.shape[0], 'all libname entries must be unique'
    
tbl_lib = tbl_lib.set_index( 'libname',drop=False )

lLibs = tbl_lib['libname'].unique()

l_refs = tbl_lib['ref_seq_name'].unique()

########
# expected output files

assert 'outdir' in config, 'must specify output directory'

# final output bam

l_out_bcxiso = expand('{}/perlib/{}{{libname}}.reads_by_bc_x_iso.txt.gz'.format(OUT_DIR,PREFIX), libname=lLibs)
l_out_overall_status_rpt = expand('{}/perlib/{}{{libname}}.reads_by_status.txt'.format(OUT_DIR,PREFIX), libname=lLibs)
l_out_overall_iso_rpt = expand('{}/perlib/{}{{libname}}.reads_by_isoform.txt'.format(OUT_DIR,PREFIX), libname=lLibs)

outRpt = OUT_DIR+'/'+PREFIX+'count_report.txt'

lOutFiles = [outRpt] + l_out_bcxiso + l_out_overall_status_rpt + l_out_overall_iso_rpt

outSampXIsoPlot = op.join( OUT_DIR, 'joined_samp_x_isogrp.html' )
outSampXIsoTbl = op.join( OUT_DIR, 'joined_samp_x_isogrp.txt' )

lOutFiles += [outSampXIsoPlot, outSampXIsoTbl]

# by ref:
l_temp_byref = expand('{}/perref/{{ref_seq_name}}/samples_to_agg.txt'.format(OUT_DIR,PREFIX), ref_seq_name=l_refs)


lOutFiles += l_temp_byref

########

rule all:
    input:
        lOutFiles

rule count_isoform_events:
    input:
        bam = lambda wc: tbl_lib.loc[ wc.libname ][ 'bam' ],
        vector_exon_tbl = config['vector_exon_tbl'],
    output:
        bcxiso = op.join( OUT_DIR, 'perlib/'+PREFIX+'{libname}.reads_by_bc_x_iso.txt.gz' ),
        overall_status_rpt = op.join( OUT_DIR, 'perlib/'+PREFIX+'{libname}.reads_by_status.txt' ),
        overall_iso_rpt =  op.join( OUT_DIR, 'perlib/'+PREFIX+'{libname}.reads_by_isoform.txt' ),
    threads: 4
    params:
        vector_exon_tbl = config['vector_exon_tbl'],
        process_rna_aligns_opts = PROCESS_RNA_ALIGNS_OPTS,
    resources:
        mem_per_cpu="4gb", 
        cpus="4", 
        time="1:00:00"
    run:
        shell("""
            process_rna_aligns \
            --bam {input.bam} \
            --vector_exon_tbl {input.vector_exon_tbl} \
            --fn_rpt_base {OUT_DIR}/perlib/{PREFIX}{wildcards.libname} \
            --rpt-extra-kv 'libname:{wildcards.libname}' \
            {params.process_rna_aligns_opts}
        """)

def get_libnames_for_ref(wc):
    llibs = list(tbl_lib.loc[tbl_lib['ref_seq_name']==wc.ref_seq_name, 'libname'])
    return llibs
    # return expand(op.join( OUT_DIR, 'perlib/'+PREFIX+'{libname}.reads_by_bc_x_iso.txt.gz' ), libname=llibs)

rule per_ref_combine_counts:
    input:
        # bcxiso = lambda wc: get_libnames_for_ref(wc)
        bcxiso = lambda wc: expand(rules.count_isoform_events.output.bcxiso, libname=get_libnames_for_ref(wc))
    output:
        per_ref_tbl = op.join( OUT_DIR, 'perref/{ref_seq_name}/samples_to_agg.txt' ),
        isogrp_tbl = op.join( OUT_DIR, 'perref/{ref_seq_name}/isogroups.txt' ),
    run:    
        print(wildcards.ref_seq_name)
        print(get_libnames_for_ref(wildcards))
        print(input.bcxiso)
        print(output.per_ref_tbl)

        tbl_per_ref = { cn:[] for cn in 'libname,bc_status_table'.split(',') }
        tbl_per_ref['bc_status_table'] = input.bcxiso
        tbl_per_ref['bc_status_table_withisogrp'] = [ fn.replace('.reads_by_bc_x_iso.txt.gz', '.reads_by_bc_x_isogrp.txt.gz') for fn in input.bcxiso ]
        tbl_per_ref['per_samp_isogrp_rpt'] = [ op.join(OUT_DIR, 'perlib/'+PREFIX+f'{libname}.reads_by_isogrp.txt') for libname in get_libnames_for_ref(wildcards) ]
        
        tbl_per_ref['libname'] = list(tbl_lib.loc[tbl_lib['ref_seq_name']==wildcards.ref_seq_name, 'libname'])
        tbl_per_ref = pd.DataFrame(tbl_per_ref)
        tbl_per_ref.to_csv(output.per_ref_tbl,sep='\t',index=False)

        shell("""
        combine_rna_counts \
            --samplesheet {output.per_ref_tbl} \
            --out_isogrps {output.isogrp_tbl} \
            --seq_name {wildcards.ref_seq_name} \
            --vector_exon_tbl {VECTOR_EXON_TBL} \
            {COMBINE_RNA_SAMP_OPTIONS}
            """)

rule outtbl:
    input:
        per_ref_tbl = expand(rules.per_ref_combine_counts.output.per_ref_tbl, ref_seq_name=l_refs),
        isogrp_tbl = expand(rules.per_ref_combine_counts.output.isogrp_tbl, ref_seq_name=l_refs),
    output:
        table_out=outRpt,
        samp_x_iso_plot=outSampXIsoPlot,
        samp_x_iso_tbl=outSampXIsoTbl,
    run:
        # concatenate the sample tables grouped by each ref
        tbl_out = []

        for ref_seq_name, fn_isogrptbl, fn in zip(l_refs, input.isogrp_tbl, input.per_ref_tbl):
            t = pd.read_table(fn)
            t['ref_seq_name'] = ref_seq_name
            t['isogrp_tbl'] = fn_isogrptbl
            tbl_out.append(t)
        
        tbl_out = pd.concat(tbl_out)
        tbl_out.to_csv(output.table_out,sep='\t',index=False)

        # join sample x isogrp tables, make plot
        shell('plot_sample_x_isogrp --sample_tbl {output.table_out} --out_plot {output.samp_x_iso_plot} --out_joined_tbl {output.samp_x_iso_tbl}')