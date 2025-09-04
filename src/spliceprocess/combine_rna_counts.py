from collections import defaultdict
from typing import Dict,List,Set
import pandas as pd
from splshared.var_mapping import VectorExonTable


def gather_counts_across_samples(
    samptbl: pd.DataFrame,
    ref_seq_name: str,
    otherisos_perbc_min_read_count: int,
    otherisos_perbc_min_psi: float,
):

    # list of files to combine
    # load files in chunk mode
    # apply min read count and min psi criteria
    # per file, return iso -> [read count] list , for qualifying barcodes

    mlib_iso_readcts = {}

    for libname, bcstatusfn in zip( samptbl['libname'], samptbl['bc_status_table'] ):
        
        print('processing',libname)

        miso_bcreadcts = defaultdict( list )
        for bcstat_chunk in pd.read_csv( bcstatusfn, sep='\t', chunksize=100000, compression='gzip' ):

            bcstat_chunk['psi'] = bcstat_chunk['ok_readcount'] / bcstat_chunk['totalrd_ok']

            li_in_ref = bcstat_chunk['isoform'].str.split(':', expand=True)[0] == ref_seq_name
            li_ct_ok = bcstat_chunk['ok_readcount'] >= otherisos_perbc_min_read_count
            li_psi_ok = bcstat_chunk['psi'] >= otherisos_perbc_min_psi
            bcstat_chunk_filter = bcstat_chunk[ li_in_ref & li_ct_ok & li_psi_ok ] 

            for iso, isotbl in bcstat_chunk_filter.groupby( 'isoform' ):
                lreadcts = isotbl['ok_readcount'].tolist()
                miso_bcreadcts[iso] += lreadcts
        
        mlib_iso_readcts[libname] = miso_bcreadcts

    return mlib_iso_readcts



def categorize_isoforms(
    mlib_iso_readcts: Dict[str, Dict[str, List[int]]],
    vec_exons_known: VectorExonTable,
    ref_seq_name: str,
    otherisos_minbc_withinallsamp: int,
    otherisos_minbc_sumacrosssamp: int,
    known_exons_strict: bool = False,
):

    """
    go through isoforms.  
    - can we name using one of the known exons? if so, keep    
    - if not, consider whether it can be in the OTHER category
    """

    s_all_isos = set()
    for libname in mlib_iso_readcts:
        s_all_isos.update( mlib_iso_readcts[libname].keys() )

    miso_known = defaultdict(list)  # known exon name(s) comma sepd -> list of isoform names
    siso_otheraccepted = set()  # 
    siso_otherother = set()  # 

    for iso in s_all_isos:
        seqname = ref_seq_name
        assert seqname == iso.split(':')[0]
        lcorng = [ (int(se.split('_')[0]), int(se.split('_')[1])) for se in iso.split(':')[1].split(',') if len(se)>0 ]

        # keep track of exons - known (by name), unknown, and both together in sorted order
        lknown, lunk, lku = [], [], []
        for corng in lcorng:
            exon_info = vec_exons_known.check_known_exon( seqname, corng[0], corng[1] )
            if exon_info is not None:
                lknown.append( exon_info['exon_name'] )
                lku.append( exon_info['exon_name'] )
            else:
                lunk.append( '%d_%d' % (corng[0], corng[1]) )
                lku.append( '%d_%d' % (corng[0], corng[1]) )

        if len(lku)==0:
            # this is skip
            miso_known['SKIP'].append(iso)
            continue

        if len(lknown)>0:
            if len(lunk)==0 or (not known_exons_strict):
                miso_known[','.join(lknown) ].append( iso )
                continue
            
        # if we get here, check to see if this meets criteria for "OTHER" isoform group
        mlib_nbc = { lib:len(mlib_iso_readcts[lib][iso]) for lib in mlib_iso_readcts }
        passes_every_indivsamp = all( [mlib_nbc[lib] >= otherisos_minbc_withinallsamp for lib in mlib_nbc ] )
        passes_sum_across_samp = sum( mlib_nbc.values() ) >= otherisos_minbc_sumacrosssamp

        if passes_every_indivsamp and passes_sum_across_samp:
            siso_otheraccepted.add(iso)
        else:
            siso_otherother.add(iso)

    return miso_known, siso_otheraccepted, siso_otherother


def agg_persamp_isogrp_counts(
    libname: str,
    ref_seq_name: str,
    fn_in_bcstatus: str,
    fn_out_bcstatus: str,
    misogrp_to_isos: Dict[str, List[str]],
    siso_otheraccepted: Set[str],
):
    # read in input bcstatus table
    # aggregate by isogrp
    # annotate bcstatus table with isogrp name and write back out

    miso_to_grp = {}
    for ig in misogrp_to_isos:
        for iso in misogrp_to_isos[ig]:
            miso_to_grp[iso] = ig

    out_rpt = { k:[] for k in ['libname',
                              'overall_nrd', 'overall_nbc',
                               'overall_nrd_ok','overall_nbc_ok', 
                               'frac_rd_ok', 'frac_bc_ok',
                               'overall_nrd_counted','overall_nbc_counted',
                               'frac_okrd_counted', 'frac_okbc_counted',
                               'seq_name',
                               'isogrp_name',
                               'isogrp_nrd_ok', 'isogrp_nbc_ok',
                               'isogrp_nrd_counted', 'isogrp_nbc_counted',
                               'isogrp_psi_rd',
                               'isogrp_psi_bc'] }


    nread, nbc = 0, 0
    nread_ok, nbc_ok = 0, 0
    nread_counted, nbc_counted = 0, 0
    misogrp_nrd_ok, misogrp_nbc_ok = defaultdict(int), defaultdict(int)
    misogrp_nrd_counted, misogrp_nbc_counted = defaultdict(int), defaultdict(int)

    isfirst = True

    for bcstat_chunk in pd.read_csv( fn_in_bcstatus, sep='\t', chunksize=100000, compression='gzip' ):

        lout = []

        for iso, isotbl in bcstat_chunk.groupby( 'isoform' ):
            
            cur_nread = sum( [ isotbl[f'totalrd_{cn}'].sum() for cn in ['ok','bad_ends','bad_starts','secondary','unpaired','unmapped','soft_clipped'] ] )
            cur_nread_ok = isotbl['totalrd_ok'].sum()
            cur_nbc_ok = ( isotbl['ok_readcount']>0 ).sum()

            nread += cur_nread
            nbc += isotbl.shape[0]
            nread_ok += cur_nread_ok
            nbc_ok += cur_nbc_ok

            if iso in miso_to_grp:
                isogrp = miso_to_grp[iso]
            elif iso in siso_otheraccepted:
                isogrp = 'OTHER'
            else:
                isogrp = None

            if isogrp is not None:
                cur_nread_counted = cur_nread_ok
                cur_nbc_counted = isotbl.shape[0]
            else:
                cur_nread_counted, cur_nbc_counted = 0, 0
            
            nread_counted += cur_nread_counted
            nbc_counted += cur_nbc_counted

            misogrp_nrd_ok[isogrp] += cur_nread_ok
            misogrp_nbc_ok[isogrp] += cur_nbc_ok
            misogrp_nrd_counted[isogrp] += cur_nread_counted
            misogrp_nbc_counted[isogrp] += cur_nbc_counted

            if isogrp is not None:
                isotbl2=isotbl.copy()
                isotbl2['isogrp_name'] = isogrp
                lout.append( isotbl2 )

        if len(lout)>0:
            lout=pd.concat( lout )
            lout=lout.sort_values(by='bc')
            if isfirst:
                lout.to_csv( fn_out_bcstatus, sep='\t', index=False, compression='gzip' )
                isfirst = False
            else:
                lout.to_csv( fn_out_bcstatus, sep='\t', index=False, mode='a', header=False, compression='gzip' )

            
    for isogrp in list(misogrp_to_isos.keys())+['OTHER']:
        out_rpt['libname'].append( libname )
        out_rpt['overall_nrd'].append( nread )
        out_rpt['overall_nbc'].append( nbc )
        out_rpt['overall_nrd_ok'].append( nread_ok )
        out_rpt['overall_nbc_ok'].append( nbc_ok )
        out_rpt['frac_rd_ok'].append( nread_ok / nread if nread > 0 else 0 )
        out_rpt['frac_bc_ok'].append( nbc_ok / nbc if nbc > 0 else 0 )
        out_rpt['overall_nrd_counted'].append( nread_counted )
        out_rpt['overall_nbc_counted'].append( nbc_counted )
        out_rpt['frac_okrd_counted'].append( nread_counted / nread_ok if nread_ok > 0 else 0 )
        out_rpt['frac_okbc_counted'].append( nbc_counted / nbc_ok if nbc_ok > 0 else 0 )
        out_rpt['seq_name'].append( ref_seq_name )
        out_rpt['isogrp_name'].append( isogrp )
        out_rpt['isogrp_nrd_ok'].append( misogrp_nrd_ok[isogrp] )
        out_rpt['isogrp_nbc_ok'].append( misogrp_nbc_ok[isogrp] )
        out_rpt['isogrp_nrd_counted'].append( misogrp_nrd_counted[isogrp] )
        out_rpt['isogrp_nbc_counted'].append( misogrp_nbc_counted[isogrp] )
        out_rpt['isogrp_psi_rd'].append( misogrp_nrd_counted[isogrp] / sum(misogrp_nrd_counted.values()) if sum(misogrp_nrd_counted.values()) > 0 else 0 )
        out_rpt['isogrp_psi_bc'].append( misogrp_nbc_counted[isogrp] / sum(misogrp_nbc_counted.values()) if sum(misogrp_nbc_counted.values()) > 0 else 0 )

    out_rpt = pd.DataFrame( out_rpt )
    return out_rpt


def main():
    import argparse

    parser = argparse.ArgumentParser(description='combine processed rna count tables across samples/replicates, single reference target sequence only')
    
    parser.add_argument('--samplesheet', help='sample sheet with library names and paths to processed rna count tables; expected columns: libname, bc_status_table, per_samp_isogrp_stats\n\
        with bc_status_table giving INPUT path to per-barcode isoform counts table\n\
        with bc_status_table_withisogrp giving OUTPUT path to per-sample isoform group stats table\n\
        and per_samp_isogrp_rpt giving OUTPUT path to per-sample isoform group stats table', 
        dest='samplesheet' )

    parser.add_argument('--out_isogrps', help='output file name for isoform group table', dest='out_isogrps' )

    parser.add_argument('--seq_name', help='name of target reference sequence', dest='seq_name' )

    parser.add_argument('--vector_exon_tbl', help='vector exon table', dest='vector_exon_tbl' )

    parser.add_argument('--otherisos_perbc_min_read_count', default=1, type=int, help='min reads for bc to contribute to OTHER isoform', dest='otherisos_perbc_min_read_count' )
    parser.add_argument('--otherisos_perbc_min_psi', default=0.025, type=float, help='min within-bc psi for bc to contribute to OTHER isoform', dest='otherisos_perbc_min_psi' )
    parser.add_argument('--otherisos_minbc_withinallsamp', default=1, type=int, help='to be counted as OTHER isoform group, min #bcs within EVERY sample', dest='otherisos_minbc_withinallsamp' )
    parser.add_argument('--otherisos_minbc_sumacrosssamp', default=1, type=int, help='to be counted as OTHER isoform group, min #bcs, summed across ALL samples', dest='otherisos_minbc_sumacrosssamp' )
    parser.add_argument('--lenient_named_isoforms', default=True, action='store_false', help='lenient mode for OTHER isoform: require all exons to be known to be named isoform', dest='strict' )

    args = parser.parse_args()

    samptbl = pd.read_table( args.samplesheet )

    vec_exons = VectorExonTable.from_file( args.vector_exon_tbl )

    lib_iso_cts = gather_counts_across_samples(
        samptbl,
        args.seq_name,
        args.otherisos_perbc_min_read_count,
        args.otherisos_perbc_min_psi
    )

    miso_known, siso_otheraccepted, siso_otherother = categorize_isoforms(
        lib_iso_cts,
        vec_exons,
        args.seq_name,
        args.otherisos_minbc_withinallsamp,
        args.otherisos_minbc_sumacrosssamp,
        args.strict
    )

    isotbl = {k:[] for k in ['seq_name', 'isogrp_name', 'isoform']}
    for isogrp in set(miso_known):
        for iso in miso_known[isogrp]:
            isotbl['seq_name'].append( args.seq_name )
            isotbl['isogrp_name'].append( isogrp )
            isotbl['isoform'].append( iso )

    if 'SKIP' not in miso_known:
        isotbl['seq_name'].append( args.seq_name )
        isotbl['isogrp_name'].append( 'SKIP' )
        isotbl['isoform'].append( f'{args.seq_name}:' )

    for iso in siso_otheraccepted:
        isotbl['seq_name'].append( args.seq_name )
        isotbl['isogrp_name'].append( 'OTHER' )
        isotbl['isoform'].append( iso )

    isotbl = pd.DataFrame( isotbl )
    isotbl.to_csv( args.out_isogrps, sep='\t', index=False )

    for _, r in samptbl.iterrows():
        per_samp_rpt = agg_persamp_isogrp_counts(
            r['libname'],
            args.seq_name,
            r['bc_status_table'],
            r['bc_status_table_withisogrp'],
            miso_known,
            siso_otheraccepted
        )
        per_samp_rpt=per_samp_rpt.sort_values( by=['isogrp_psi_bc'], ascending=False )
        per_samp_rpt.to_csv( r['per_samp_isogrp_rpt'], sep='\t', index=False )

if __name__ == '__main__':
    main()