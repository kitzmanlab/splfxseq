import sys
import pysam    
import itertools
from collections import defaultdict, Counter
from typing import Tuple
import pandas as pd
from splshared.var_mapping import VectorExonTable

def clean_jns_pe( read1: pysam.AlignedSegment,
                  read2: pysam.AlignedSegment,
                  corng_upstream_ex: Tuple[int,int],
                  corng_dnstream_ex: Tuple[int,int],
                  spl_tol: int,
                  indel_tol: int,
                  min_matches_for: int,
                  min_matches_rev: int,
                  max_soft_clip_for: int,
                  max_soft_clip_rev: int ) -> Tuple[str,str]:

    if read1.get_blocks()[ 0 ][ 0 ] > read2.get_blocks()[ 0 ][ 0 ]:
        temp_r1 = read1
        read1 = read2
        read2 = temp_r1

    if read1.reference_name != read2.reference_name:
        return ('unmapped',None)

    refname = read1.reference_name

    r12_jns = [ read1.get_blocks(), read2.get_blocks() ]
    r12_cigs = [ read1.cigartuples, read2.cigartuples ]

    #if one of the reads is unmapped skip the pair
    if not all( r12_jns ):
        return ('unmapped',None)

    r12_jns = [ list( r12_jns[ 0 ] ), list( r12_jns[ 1 ] ) ]

    us_cnst = ds_cnst = False

    #screens out reads with too many bases soft_clipped
    if sum( jn[ 1 ] - jn[ 0 ] + 1 for jn in r12_jns[ 0 ] ) < min_matches_for \
    or sum( jn[ 1 ] - jn[ 0 ] + 1 for jn in r12_jns[ 1 ] ) < min_matches_rev:
        return ('soft_clipped',None)

    if sum( [t[1] for t in r12_cigs[0] if t[0] == 4] ) > max_soft_clip_for \
    or sum( [t[1] for t in r12_cigs[1] if t[0] == 4] ) > max_soft_clip_rev:
        return ('soft_clipped',None)

    #this section will collapse any small indels
    for i,r_cig in enumerate( r12_cigs ):

        #there are no indels - all is well
        if not any( c_o == 1 or c_o == 2  for c_o,c_l in r_cig ):
            continue

        ri_jn = []
        j = 0
        indel = False

        for c_o,c_l in r_cig:

            if indel:
                indel = False
                continue

            #if its a match add to append jn
            #don't add any very short alignments since those are unreliable
            if c_o == 0 and c_l >= spl_tol:
                ri_jn.append( r12_jns[ i ][ j ] )
                j += 1
            #if its an indel within tolerance, add it to previous junction
            elif ri_jn and ( c_o == 1 or c_o == 2 ) and c_l <= indel_tol:
                indel = True
                ri_jn[ -1 ] = ( ri_jn[ -1 ][ 0 ], r12_jns[ i ][ j ][ 1 ] )
                j += 1

        r12_jns[ i ] = ri_jn

    #this section combines the reads
    #if the reads don't overlap combine them
    #so ( 1267, 1312 ), ( 1480, 1512 ) becomes ( 1267, 1512 )
    if r12_jns[ 0 ][ -1 ][ 1 ] <= r12_jns[ 1 ][ 0 ][ 0 ]:
        r12_jns[ 0 ][ -1 ] = ( r12_jns[ 0 ][ -1 ][ 0 ], r12_jns[ 1 ][ 0 ][ 1 ] )
        r12_jns = tuple( sorted( r12_jns[ 0 ] + r12_jns[ 1 ][ 1: ] ) )

    else:
        r12_jns = tuple( sorted( r12_jns[ 0 ] + r12_jns[ 1 ] ) )

    r12_jn_comb = []
    for _jn in r12_jns:

        #checks if read splices properly from upstream constant donor
        if _jn[ 1 ] >= corng_upstream_ex[ 1 ] - spl_tol and _jn[ 1 ] <= corng_upstream_ex[ 1 ] + spl_tol:
            us_cnst = True
            continue

        #checks if read splices properly into the downstream constant acceptor
        if _jn[ 0 ] >= corng_dnstream_ex[ 0 ] - spl_tol and _jn[ 0 ] <= corng_dnstream_ex[ 0 ] + spl_tol:
            ds_cnst = True
            continue

        #changes to 1-based inclusive numbering
        jn = ( _jn[ 0 ] + 1, _jn[ 1 ] )

        #if the jn is outside of the constant exons
        if jn[ 0 ] >= corng_upstream_ex[ 1 ] and jn[ 1 ] < corng_dnstream_ex[ 0 ] + 1:

            if not r12_jn_comb or jn[ 0 ] > r12_jn_comb[ -1 ][ 1 ]:
                r12_jn_comb.append( jn )
            #if the forward and reverse reads cross - combine them
            #so ( 1267, 1312 ), ( 1297, 1512 ) becomes ( 1267, 1512 )
            else:
                r12_jn_comb[ -1 ] = ( r12_jn_comb[ -1 ][ 0 ], jn[ 1 ] )

    r12_jn_comb = tuple( r12_jn_comb )

    r12_jn_comb = f'{refname}:' + ','.join([f'{jn[0]}_{jn[1]}' for jn in r12_jn_comb])

    if not ds_cnst:
        return ('bad_ends',r12_jn_comb)
    elif not us_cnst:
        return ('bad_starts',r12_jn_comb)
    else:      
        return ('ok',r12_jn_comb)
    

def process_PE_rna_bam(
    bam: pysam.AlignmentFile,
    vec_exons: VectorExonTable,
    spl_tol: int,
    indel_tol: int,
    min_matches_for: int,
    min_matches_rev: int,
    max_soft_clip_for: int,
    max_soft_clip_rev: int ,
    fn_bam_out_reject_bad_ends:str=None,
    fn_bam_out_reject_bad_starts:str=None,
    fn_bam_out_reject_secondary:str=None,
    fn_bam_out_reject_softclip:str=None,
    fn_bam_out_reject_unmapped:str=None,
    ):

    """
    Process a MPSA RNAseq bam file, which is expected to be sorted by read name and contain paired-end data with one alignment per read.
    Supplementary alignments should be excluded from the input. 
    Output is a dictionary 
        status_count -> per-read status, count of readpairs
        bc_status -> per-barcode, count of readpairs for each status
        bc_iso -> map from barcode to isoform to count, for "ok" status read pairs only
        iso_count -> counter of isoform read count, for "ok" status read pairs only
    """

    rejreason_bam = {}

    if fn_bam_out_reject_bad_ends:
        rejreason_bam['bad_ends'] = pysam.AlignmentFile( fn_bam_out_reject_bad_ends, "wb", template = bam )
    if fn_bam_out_reject_bad_starts:
        rejreason_bam['bad_starts'] = pysam.AlignmentFile( fn_bam_out_reject_bad_starts, "wb", template = bam )
    if fn_bam_out_reject_secondary:
        rejreason_bam['secondary'] = pysam.AlignmentFile( fn_bam_out_reject_secondary, "wb", template = bam )
    if fn_bam_out_reject_softclip:
        rejreason_bam['soft_clipped'] = pysam.AlignmentFile( fn_bam_out_reject_softclip, "wb", template = bam )
    if fn_bam_out_reject_unmapped:
        rejreason_bam['unmapped'] = pysam.AlignmentFile( fn_bam_out_reject_unmapped, "wb", template = bam )

    mvec_constupdn = vec_exons.get_const_exon_positions()

    mbc_status = defaultdict( lambda:Counter( {'ok':0, 'bad_ends':0, 'bad_starts':0, 'secondary':0, 'unpaired':0, 'unmapped':0, 'soft_clipped':0} ) )
    mbc_iso = defaultdict( lambda:Counter() )
    miso_rdcnt = Counter()
    mstatus_rdcnt = Counter( {'ok':0, 'bad_ends':0, 'bad_starts':0, 'secondary':0, 'unpaired':0, 'unmapped':0, 'soft_clipped':0} )

    Npairs_ctr = 0

    for readname, reads in itertools.groupby( bam, key=lambda x: x.query_name ):  

        Npairs_ctr += 1
        if Npairs_ctr % 1000000 == 0:
            sys.stderr.write( f'{Npairs_ctr} pairs...\n' )
            sys.stderr.flush()

        reads = list( reads )

        lr1 = [ r for r in reads if r.is_read1 ]
        lr2 = [ r for r in reads if r.is_read2 ]
    
        if len(lr1) < 1:
            raise ValueError( f'read {readname} has no read1 alignment' )
        if len(lr2) < 1:
            raise ValueError( f'read {readname} has no read2 alignment' )
        
        if len(lr1) > 1:
            raise ValueError( f'read {readname} has too many read1 alignments' )
        if len(lr2) > 1:
            raise ValueError( f'read {readname} has too many read2 alignments' )

        read1 = lr1[0]
        read2 = lr2[0]

        bc = read1.get_tag( 'BC' )
        if bc!=read2.get_tag( 'BC' ):
            raise ValueError( f'read {readname} has different barcodes in read1, read2' )

        rej = None
        if read1.is_secondary or read2.is_secondary:
            rej = 'secondary'
        elif not read1.is_paired or not read2.is_paired:
            rej = 'unpaired'
        elif read1.is_unmapped or read2.is_unmapped or read1.mate_is_unmapped:
            rej = 'unmapped'

        chromname = read1.reference_name
        if read1.reference_name != read2.reference_name:
            raise ValueError( f'read {readname} has different reference names in read1, read2' )
        corng_upstream_ex, corng_dnstream_ex = mvec_constupdn[chromname]

        if rej:
            if rej in rejreason_bam:
                rejreason_bam[ rej ].write( read1 )
                rejreason_bam[ rej ].write( read2 )

            mbc_status[ bc ][ rej ] += 1
            mstatus_rdcnt[ rej ] += 1
            continue

        else:
            jns = clean_jns_pe( read1,
                                read2,
                                corng_upstream_ex,
                                corng_dnstream_ex,
                                spl_tol,
                                indel_tol,
                                min_matches_for,
                                min_matches_rev,
                                max_soft_clip_for,
                                max_soft_clip_rev )

            rej = jns[0]
            if jns[0] != 'ok':
                if rej in rejreason_bam:
                    rejreason_bam[ rej ].write( read1 )
                    rejreason_bam[ rej ].write( read2 )
            else:
                mbc_iso[ bc ][ jns[1] ] += 1
                miso_rdcnt[ jns[1] ] += 1

            mbc_status[ bc ][ rej ] += 1
            mstatus_rdcnt[ rej ] += 1

        rej = None
        
    return  {'bc_status':mbc_status, 'bc_iso':mbc_iso, 'iso_count':miso_rdcnt, 'status_count':mstatus_rdcnt}


def process_PE_rna_bam_write_rpt(
    bam: pysam.AlignmentFile,
    fn_rpt_base: str, 
    vec_exons: VectorExonTable,
    spl_tol: int,
    indel_tol: int,
    min_matches_for: int,
    min_matches_rev: int,
    max_soft_clip_for: int,
    max_soft_clip_rev: int ,
    fn_bam_out_reject_bad_ends:str=None,
    fn_bam_out_reject_bad_starts:str=None,
    fn_bam_out_reject_secondary:str=None,
    fn_bam_out_reject_softclip:str=None,
    fn_bam_out_reject_unmapped:str=None,
    rpt_extra_kv:dict=None,
    ):

    res = process_PE_rna_bam(
        bam,
        vec_exons,
        spl_tol,
        indel_tol,
        min_matches_for,
        min_matches_rev,
        max_soft_clip_for,
        max_soft_clip_rev,
        fn_bam_out_reject_bad_ends,
        fn_bam_out_reject_bad_starts,
        fn_bam_out_reject_secondary,
        fn_bam_out_reject_softclip,
        fn_bam_out_reject_unmapped,
    )

    # write report of read counts by status
    mstatus_rdcnt = pd.Series( res['status_count'] )
    if rpt_extra_kv:
        for k in rpt_extra_kv:
            mstatus_rdcnt[k]=rpt_extra_kv[k]
    mstatus_rdcnt.to_frame().T.to_csv( fn_rpt_base + '.reads_by_status.txt', sep='\t', index=False )

    # write report of read counts by isoform
    miso_rdcnt = pd.Series( res['iso_count'] ).to_frame()
    miso_rdcnt.index.name = 'isoform'
    miso_rdcnt.columns = ['ok_readcount']
    if rpt_extra_kv:
        for k in rpt_extra_kv:
            miso_rdcnt[k]=rpt_extra_kv[k]
    miso_rdcnt = miso_rdcnt.sort_values( by='ok_readcount', ascending=False )
    miso_rdcnt.to_csv( fn_rpt_base + '.reads_by_isoform.txt', sep='\t', index=True )

    # write (bc,iso) long table
    mbc_iso = res['bc_iso'] 
    bcxiso = pd.DataFrame(
    itertools.chain.from_iterable(
        [ [ (bc,iso,mbc_iso[bc][iso]) for iso in mbc_iso[bc] ] for bc in mbc_iso ]
    ),
    columns=['bc','isoform','ok_readcount'])

    mbc_status = pd.DataFrame( res['bc_status'] ).T
    mbc_status.columns = [f'totalrd_{c}' for c in mbc_status.columns]
    mbc_status = mbc_status.fillna(0).astype(int)
    mbc_status.index.name = 'bc'

    bcxiso_wstatus = pd.merge( bcxiso, mbc_status, left_on='bc', right_index=True, how='left' )

    bcxiso_wstatus.to_csv( fn_rpt_base + '.reads_by_bc_x_iso.txt.gz', sep='\t', index=False, compression='gzip' )

def main():
    import argparse

    parser = argparse.ArgumentParser(description='Process paired-end RNA alignments and write reports')
    
    parser.add_argument('--bam', help='Input BAM file')
    parser.add_argument('--vector_exon_tbl', required=True, help='Vector exon table')
    parser.add_argument('--fn_rpt_base', required=True, help='Base filename for output reports')
    parser.add_argument('--spl-tol', type=int, default=3, help='Splice junction tolerance')
    parser.add_argument('--indel-tol', type=int, default=5, help='Indel tolerance')
    parser.add_argument('--min-matches-for', type=int, default=20, help='Minimum matches for forward read')
    parser.add_argument('--min-matches-rev', type=int, default=20, help='Minimum matches for reverse read')
    parser.add_argument('--max-soft-clip-for', type=int, default=10, help='Maximum soft clipping for forward read')
    parser.add_argument('--max-soft-clip-rev', type=int, default=10, help='Maximum soft clipping for reverse read')
    parser.add_argument('--reject-bad-ends', default=None, help='BAM file for rejected reads with bad ends')
    parser.add_argument('--reject-bad-starts', default=None, help='BAM file for rejected reads with bad starts')
    parser.add_argument('--reject-secondary', default=None, help='BAM file for rejected secondary alignments')
    parser.add_argument('--reject-softclip', default=None, help='BAM file for rejected soft clipped reads')
    parser.add_argument('--reject-unmapped', default=None, help='BAM file for rejected unmapped reads')
    parser.add_argument('--rpt-extra-kv', default='', help='Extra key-value pairs for reports')

    args = parser.parse_args()

    vec_exons = VectorExonTable.from_file( args.vector_exon_tbl )

    # Convert rpt_extra_kv list to dict if provided
    rpt_extra_kv = dict( [ kv.split(':') for kv in args.rpt_extra_kv.split(',') ] ) if args.rpt_extra_kv and len(args.rpt_extra_kv) > 0 else None

    process_PE_rna_bam_write_rpt(
        pysam.AlignmentFile( args.bam ),
        args.fn_rpt_base,
        vec_exons,
        args.spl_tol,
        args.indel_tol,
        args.min_matches_for,
        args.min_matches_rev,
        args.max_soft_clip_for,
        args.max_soft_clip_rev,
        args.reject_bad_ends,
        args.reject_bad_starts,
        args.reject_secondary,
        args.reject_softclip,
        args.reject_unmapped,
        rpt_extra_kv
    )



if __name__ == '__main__':
    main()