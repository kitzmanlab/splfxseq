import pysam    
from typing import Tuple

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
    