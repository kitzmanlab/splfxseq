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

l_reqd_cols = [ 'libname', 'isogrp_tbl', 'bc_status_table_withisogrp', 'ref_seq_name', 'pairing_tbl_name' ]
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
l_out_var_rpt = expand(
    op.join(OUT_DIR, 'perwt_rpts', PREFIX + '{libname}.wt_rpt.txt'),
    libname=lLibs
)
outRpt = OUT_DIR+'/'+PREFIX+'wtfx_rpt.txt'
lOutFiles = [outRpt] + l_out_var_rpt 

########

rule all:
    input:
        lOutFiles

rule per_wt_process:
    input:
        isogrp_tbl = lambda wc: tbl_lib.loc[ wc.libname ][ 'isogrp_tbl' ],
        bc_x_iso_status_tbl = lambda wc: tbl_lib.loc[ wc.libname ][ 'bc_status_table_withisogrp' ],
        bc_pairing_tbl = lambda wc: tbl_lib.loc[ wc.libname ][ 'pairing_tbl_path' ],
    params:
        seq_name = lambda wc: tbl_lib.loc[ wc.libname ][ 'ref_seq_name' ],
    output:
        wt_rpt = op.join(OUT_DIR, 'perwt_rpts/'+PREFIX+'{libname}.wt_rpt.txt'),
        per_bc_rpt = op.join(OUT_DIR, 'perwt_bcrpts/'+PREFIX+'{libname}.bybc.txt.gz') if MAKE_PERBC_RPT else [],
    threads: 1
    resources:
        mem_per_cpu="10gb", 
        cpus="1", 
        runtime="1h"
    run:
        if MAKE_PERBC_RPT:
            shell("""
                rna_counts_add_wt_info \
                --isogrp_tbl {input.isogrp_tbl} \
                --bc_x_iso_status_tbl {input.bc_x_iso_status_tbl} \
                --seq_name {params.seq_name} \
                --bc_pairing_tbl {input.bc_pairing_tbl} \
                --libname {wildcards.libname} \
                --out_wt_rpt {output.wt_rpt} \
                --out_perbc {output.per_bc_rpt}
            """)
        else:
            shell("""
                rna_counts_add_wt_info \
                --isogrp_tbl {input.isogrp_tbl} \
                --bc_x_iso_status_tbl {input.bc_x_iso_status_tbl} \
                --seq_name {params.seq_name} \
                --bc_pairing_tbl {input.bc_pairing_tbl} \
                --libname {wildcards.libname} \
                --out_wt_rpt {output.wt_rpt}
            """)

rule out_rpt:
    input:
        per_wt_rpt = expand(rules.per_wt_process.output.wt_rpt, libname=lLibs),
        per_bc_rpt = expand(rules.per_wt_process.output.per_bc_rpt,libname=lLibs)
    output:
        wt_rpt_out = op.join(OUT_DIR,PREFIX+'wtfx_rpt.txt')
    run:
        out_rpt = {k:[] for k in 'libname,wt_rpt,wt_bc_rpt,exon'.split(',')}
        out_rpt['libname'] = lLibs
        out_rpt['wt_rpt'] = list(input.per_wt_rpt)
        out_rpt["wt_bc_rpt"] = list(input.per_bc_rpt) if input.per_bc_rpt else [None] * len(lLibs)
        out_rpt['exon'] = tbl_lib.loc[lLibs, 'pairing_tbl_name'].tolist()
        out_rpt = pd.DataFrame(out_rpt)
        out_rpt.to_csv(output.wt_rpt_out, sep='\t', index=False)
