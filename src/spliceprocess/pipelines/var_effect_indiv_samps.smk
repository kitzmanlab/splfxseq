######################################################
# Configuration values you must specify:
#
#
# > lib_table: library counts report, output from count_isoform_events_singleref.smk  
#
#   one row per sequencing library to be processed

#   required columns
#   - libname   (this nust be uniq; if you are running the same lib across multiple batches you need to 
#   - isogrp_tbl (path to isoform group table for this reference)
#   - bc_status_table_withisogrp  (path to isoform group counts table per barcode)
#   - ref_seq_name (name of the reference sequence)
#
#   - pairing_tbl_name (name of the pairing table) <- THIS ONE MUST BE ADDED AND MUST MATCH A ROW IN pairing_tbls

# > pairing_tables:  table of pairing tables
#   required columns
#   - pairing_tbl_name (name of the pairing table)
#   - pairing_tbl_path (path to the pairing table)
#
# > outdir
# 
#######################################################
#  optional config values:
#
#   prefix 
#
#   suppress_perbc_rpt (default False)

#######################################################

import os.path as op
import os
import pandas as pd
import pysam 
import Bio.Seq

import altair as alt 

########

# localrules: outtbl

assert 'lib_table' in config, 'must specify sample table'
assert 'pairing_tables' in config, 'must specify pairing table table'

PREFIX = config['prefix'] if 'prefix' in config  else ''
OUT_DIR = config['outdir']

MAKE_PERBC_RPT = ('suppress_perbc_rpt' not in config)

########
#
# load and check sample table.

l_reqd_cols = [ 'libname', 'bc_status_table_withisogrp', 'ref_seq_name', 'pairing_tbl_name' ]
tbl_lib = pd.read_table( config['lib_table'] )
assert all( [ col in tbl_lib.columns for col in l_reqd_cols ] ), 'lib table must have columns: '+','.join(l_reqd_cols)
assert len(set(tbl_lib['libname'])) == tbl_lib.shape[0], 'all libname entries must be unique'
  
tbl_lib = tbl_lib.set_index( 'libname',drop=False )

lLibs = tbl_lib['libname'].unique()

########
#
# load and check pairing table

lreqd_cols_pairing_tbls = [ 'pairing_tbl_name', 'pairing_tbl_path' ]
tbl_pairing_tbls = pd.read_table( config['pairing_tables'] )
assert all( [ col in tbl_pairing_tbls.columns for col in lreqd_cols_pairing_tbls ] ), 'pairing table must have columns: '+','.join(lreqd_cols_pairing_tbls)
assert len(set(tbl_pairing_tbls['pairing_tbl_name'])) == tbl_pairing_tbls.shape[0], 'all pairing_tbl_name entries must be unique'
tbl_pairing_tbls = tbl_pairing_tbls.set_index( 'pairing_tbl_name' )
l_pairing_tbl_names = list(tbl_pairing_tbls.index.unique())

####
# merge libs w/ pairing tables

tbl_lib = pd.merge( 
    tbl_lib,
    tbl_pairing_tbls,
    left_on='pairing_tbl_name',
    right_on='pairing_tbl_name',
    how='left',
)

libs_pairing_tbl_notfound = list(tbl_lib.loc[  
    tbl_lib['pairing_tbl_path'].isnull()
].libname)

if len(libs_pairing_tbl_notfound) > 0:
    raise ValueError('pairing tables not found for these entries ' + ','.join(libs_pairing_tbl_notfound))

tbl_lib = tbl_lib.set_index('libname',drop=False)

########
# expected output files

assert 'outdir' in config, 'must specify output directory'

# outputs

l_out_var_rpt = expand('{}/pervar_rpts/{}{{libname}}.var_rpt.txt'.format(OUT_DIR,PREFIX), libname=lLibs)
l_out_var_rpt_snvonly = expand('{}/pervar_rpts_snvonly/{}{{libname}}.var_rpt_snvonly.txt'.format(OUT_DIR,PREFIX), libname=lLibs)
l_out_snv_plot  = expand('{}/pervar_plots_snvonly/{}{{libname}}.var_rpt_snvonly.html'.format(OUT_DIR,PREFIX), libname=lLibs)
l_out_per_bc_rpt  = expand('{}/pervar_bcrpts/{}{{libname}}.bybc.txt.gz'.format(OUT_DIR,PREFIX), libname=lLibs) if MAKE_PERBC_RPT else []

outRpt = OUT_DIR+'/'+PREFIX+'varfx_rpt.txt'
lOutFiles = [outRpt] + l_out_var_rpt + l_out_var_rpt_snvonly + l_out_snv_plot + l_out_per_bc_rpt

########

rule all:
    input:
        lOutFiles

rule per_var_process:
    input:
        isogrp_tbl = lambda wc: tbl_lib.loc[ wc.libname ][ 'isogrp_tbl' ],
        bc_x_iso_status_tbl = lambda wc: tbl_lib.loc[ wc.libname ][ 'bc_status_table_withisogrp' ],
        bc_pairing_tbl = lambda wc: tbl_lib.loc[ wc.libname ][ 'pairing_tbl_path' ],
    params:
        seq_name = lambda wc: tbl_lib.loc[ wc.libname ][ 'ref_seq_name' ],
    output:
        var_rpt = op.join(OUT_DIR, 'pervar_rpts/'+PREFIX+'{libname}.var_rpt.txt'),
        per_bc_rpt = op.join(OUT_DIR, 'pervar_bcrpts/'+PREFIX+'{libname}.bybc.txt.gz') if MAKE_PERBC_RPT else None,
    threads: 1
    resources:
        mem_per_cpu="10gb", 
        cpus="1", 
        time="1:00:00"
    run:
        if MAKE_PERBC_RPT:
            shell("""
                rna_counts_add_var_info \
                --isogrp_tbl {input.isogrp_tbl} \
                --bc_x_iso_status_tbl {input.bc_x_iso_status_tbl} \
                --seq_name {params.seq_name} \
                --bc_pairing_tbl {input.bc_pairing_tbl} \
                --libname {wildcards.libname} \
                --out_var_rpt {output.var_rpt} \
                --out_perbc {output.per_bc_rpt}
            """)
        else:
            shell("""
                rna_counts_add_var_info \
                --isogrp_tbl {input.isogrp_tbl} \
                --bc_x_iso_status_tbl {input.bc_x_iso_status_tbl} \
                --seq_name {params.seq_name} \
                --bc_pairing_tbl {input.bc_pairing_tbl} \
                --libname {wildcards.libname} \
                --out_var_rpt {output.var_rpt}
            """)

rule filt_to_snv:
    input:
        var_rpt = op.join(OUT_DIR, 'pervar_rpts/'+PREFIX+'{libname}.var_rpt.txt')
    output:
        var_rpt_snvonly = op.join(OUT_DIR, 'pervar_rpts_snvonly/'+PREFIX+'{libname}.var_rpt_snvonly.txt'),
    run:
        shell("""
            varfx_singlesamp filt_to_snv --var_tbl_in {input.var_rpt} --var_tbl_out {output.var_rpt_snvonly}
        """)

rule mk_varfx_plot:
    input:
        var_rpt_snvonly = op.join(OUT_DIR, 'pervar_rpts_snvonly/'+PREFIX+'{libname}.var_rpt_snvonly.txt'),
        isogrp_tbl = lambda wc: tbl_lib.loc[ wc.libname ][ 'isogrp_tbl' ],
    params:
        seq_name = lambda wc: tbl_lib.loc[ wc.libname ][ 'ref_seq_name' ],
    output:
        varfx_plot = op.join(OUT_DIR, 'pervar_plots_snvonly/'+PREFIX+'{libname}.var_rpt_snvonly.html')
    run:
        shell("""
            varfx_singlesamp singlerep_plot --var_tbl {input.var_rpt_snvonly} --isogrp_tbl {input.isogrp_tbl} --ref_seq_name {params.seq_name} --out_plot {output.varfx_plot}
        """)

rule out_rpt:
    input:
        per_var_rpt = expand(rules.per_var_process.output.var_rpt, libname=lLibs),
        per_var_snv_rpt = expand(rules.filt_to_snv.output.var_rpt_snvonly, libname=lLibs),
    output:
        var_rpt_out = op.join(OUT_DIR,PREFIX+'varfx_rpt.txt')
    run:
        out_rpt = {k:[] for k in 'libname,var_rpt,var_rpt_snvonly'.split(',')}
        out_rpt['libname'] = lLibs
        out_rpt['var_rpt'] = list(input.per_var_rpt)
        out_rpt['var_rpt_snvonly'] = list(input.per_var_snv_rpt)
        out_rpt = pd.DataFrame(out_rpt)
        out_rpt.to_csv(output.var_rpt_out, sep='\t', index=False)
