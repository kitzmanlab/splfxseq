######################################################
#
# Assess barcode overlap between MPSA RNAseq sequencing library and barcode pairing table
#
#
# Configuration values you must specify:
#
#
# > lib_table  
#   one row per sequencing library to be processed
#   required columns
#   - libname   (this must be unique)
#   - pairing_tbl_name (name of pairing table to use)
#   - [colname_histo] (the name column name of the barcode histogram; default 'bc_histo', specified by option colname_histo)
#
# > pairing_tbls
#   one row per pairing table to be processed
#   required columns
#   - pairing_tbl_name (name of pairing table to use)
#   - pairing_tbl_path (path to pairing table)
#
# > outdir
# 
#
#######################################################
#  optional config values:
#
# - colname_histo (the name column name of the barcode histogram to look for in the library table; default 'bc_histo')
# - prefix (prefix for output files; default '')
# - lmin_barcode (minimum barcode counts, comma separated list of values; default '1,5,10')
# - lqile (quantile values, comma separated list of values; default '0.75,0.90')

import os.path as op
import os
import pandas as pd
import pysam 
import Bio.Seq

import altair as alt 

########

localrules: outtbl

assert 'lib_table' in config, 'must specify sample table'
assert 'pairing_tbls' in config, 'must specify pairing table table'

###

assert 'outdir' in config, 'must specify output directory'  
OUT_DIR = config['outdir'].strip()
COLNAME_HISTO = config['colname_histo'].strip() if 'colname_histo' in config else 'bc_histo'
PREFIX = config['prefix'].strip() if 'prefix' in config else ''
LMIN_BARCODE = config['lmin_barcode'].strip() if 'lmin_barcode' in config else '1,5,10'
LQILE = config['lqile'].strip() if 'lqile' in config else '0.75,0.90'

## load pairing table

l_reqd_cols_pairing_tbls = [ 'pairing_tbl_name', 'pairing_tbl_path' ]
tbl_pairing_tbls = pd.read_table( config['pairing_tbls'] )
assert all( [ col in tbl_pairing_tbls.columns for col in l_reqd_cols_pairing_tbls ] ), 'pairing table must have columns: '+','.join(l_reqd_cols_pairing_tbls)
assert len(set(tbl_pairing_tbls['pairing_tbl_name'])) == tbl_pairing_tbls.shape[0], 'all pairing_tbl_name entries must be unique'

tbl_pairing_tbls = tbl_pairing_tbls.set_index( 'pairing_tbl_name' )

l_pairing_tbl_names = list(tbl_pairing_tbls.index.unique())

########
# load and check sample table.

l_reqd_cols = [ 'libname', COLNAME_HISTO, 'pairing_tbl_name' ]
tbl_lib = pd.read_table( config['lib_table'] )
assert all( [ col in tbl_lib.columns for col in l_reqd_cols ] ), 'lib table must have columns: '+','.join(l_reqd_cols)
assert len(set(tbl_lib['libname'])) == tbl_lib.shape[0], 'all libname entries must be unique'

tbl_libpairing = pd.merge( tbl_lib, tbl_pairing_tbls, on='pairing_tbl_name', how='left' )

if tbl_libpairing['pairing_tbl_path'].isnull().any():
    raise ValueError( 'some libraries have pairing table entires that are not found' + ','.join(set( tbl_lib.loc[tbl_libpairing['pairing_tbl_path'].isnull()]['libname'])))

tbl_lib = tbl_libpairing

tbl_lib = tbl_lib.set_index( 'libname',drop=False )

lLibs = list(tbl_lib['libname'].unique())

########
# expected output files

l_knownbc_overlap_summary = expand('{}/indiv_sample_rpts/{}{{libname}}_indivsamp_overlap_knownbc.tsv'.format(OUT_DIR,PREFIX), libname=lLibs)
l_knownbc_overlap_qile_summary = expand('{}/indiv_sample_rpts/{}{{libname}}_indivsamp_overlap_knownbc_qiles.tsv'.format(OUT_DIR,PREFIX), libname=lLibs)

outRpt = OUT_DIR+'/'+PREFIX+'per_sample_overlap_rpt.tsv'
outRpt_qiles = OUT_DIR+'/'+PREFIX+'per_sample_overlap_rpt_qiles.tsv'

lOutFiles = [outRpt, outRpt_qiles] + l_knownbc_overlap_summary + l_knownbc_overlap_qile_summary

########

rule all:
    input:
        lOutFiles


# def lut(t,i,col):
#     print(t.loc[i])
#     return t.loc[i][col]

rule check_overlap_rnabc_knownbc:  
    input:
        bc_histo = lambda wc: tbl_lib.loc[ wc.libname ][ COLNAME_HISTO ],
        pairing_tbl = lambda wc: tbl_lib.loc[ wc.libname ][ 'pairing_tbl_path' ],
        # bc_histo = lut(tbl_lib,lambda wc: wc.libname,COLNAME_HISTO),
        # pairing_tbl = lut(tbl_lib,lambda wc: wc.libname,'pairing_tbl_path'),
    output:
        overlap_summary = op.join(OUT_DIR,'indiv_sample_rpts/'+PREFIX+'{libname}_indivsamp_overlap_knownbc.tsv'),
        overlap_qile_summary = op.join(OUT_DIR,'indiv_sample_rpts/'+PREFIX+'{libname}_indivsamp_overlap_knownbc_qiles.tsv'),
    params:
        libname = lambda wc: wc.libname,
    run:
        shell(f"""
        check_overlap_rnabc_knownbc singlesamp --samp_name {{params.libname}} \
            --pairing_tbl {{input.pairing_tbl}} \
            --bc_histo {{input.bc_histo}} \
            --out_base {OUT_DIR}/indiv_sample_rpts/{PREFIX}{wildcards.libname} \
            --lmin_barcode {LMIN_BARCODE} \
            --lqile {LQILE}"""
        )


rule outtbl:
    input:
        overlap_summary = expand(rules.check_overlap_rnabc_knownbc.output.overlap_summary, libname=lLibs),
        overlap_qile_summary = expand(rules.check_overlap_rnabc_knownbc.output.overlap_qile_summary, libname=lLibs),
    output:
        table_out=outRpt,
        table_qile_out=outRpt_qiles,
    run:

        loverlap_summary = [ pd.read_table(fn) for fn in input.overlap_summary ]
        loverlap_qile_summary = [ pd.read_table(fn) for fn in input.overlap_qile_summary ]

        loverlap_summary = pd.concat(loverlap_summary,axis=0)
        loverlap_qile_summary = pd.concat(loverlap_qile_summary,axis=0)

        loverlap_summary['pairing_tbl_name'] = pd.merge(loverlap_summary,tbl_lib,left_on='libname',right_index=True,how='left')['pairing_tbl_name']
        loverlap_qile_summary['pairing_tbl_name'] = pd.merge(loverlap_qile_summary,tbl_lib,left_on='libname',right_index=True,how='left')['pairing_tbl_name']

        loverlap_summary.to_csv( output.table_out, sep='\t', index=False )
        loverlap_qile_summary.to_csv( output.table_qile_out, sep='\t', index=False )

        shell(f"""
        check_overlap_rnabc_knownbc summarize_samps --overlap_table {output.table_out} \
            --overlap_qile_table {output.table_qile_out} \
            --out_base {OUT_DIR}/{PREFIX}"""
        )
