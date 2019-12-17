import pysam
import os
import pandas as pd
import numpy as np
import networkx as nx
from collections import Counter
import matplotlib.pyplot as plt

import itertools

import pybedtools as pbt


def make_junction_graph(exon_bed):

    """Create a graph of all possible traversals of specified exons

    Args:
        exon_bed (PyBedTool): bed file with exon records.  The first record is assumed to be the constant first exon,
            and the last record is assumed to be the constant final exon

    Returns:
        dict (str -> tuple): dictionary of named path (ie isoform) traversing exons.  Keys are path names,
            values are tuple (series) of tuples corresponding to exons on that path.  The exon coordinates returned are
            1-based, inclusive.

    """

    gr = nx.DiGraph()

    nex = len(exon_bed)

    gr.add_node( exon_bed[0].start+1 )
    gr.add_node( exon_bed[0].end )

    # add donors and acceptors
    for iex in range(1, len(exon_bed)):
        gr.add_node( exon_bed[iex].start+1 ) # +1 to make this 1-based,inclusive
        gr.add_node( exon_bed[iex].end )     # +0 to make this 1-based,inclusive

    # loop through donors in order
    for iex in range(len(exon_bed)):
        gr.add_edge( exon_bed[iex].start+1, exon_bed[iex].end )
        # at each donor, add a path to any acceptor which is >= donor position
        for jex in range(len(exon_bed)):
            if exon_bed[jex].start+1 > exon_bed[iex].end:
                gr.add_edge( exon_bed[iex].end, exon_bed[jex].start+1 )

    gr.add_node( exon_bed[nex-1].start+1 )
    gr.add_node( exon_bed[nex-1].end )

    lpaths = list( [path for path in
                      nx.all_simple_paths( gr, exon_bed[0].start+1, exon_bed[nex-1].end )] )

    # convert paths to series of (exonstart,exonend) for legibility

    lpaths = [ tuple( [ (path[i*2],path[i*2+1]) for i in range(int(len(path)/2)) ] )
               for path in lpaths ]

    lpathnames = ['iso{:02d}'.format( i ) for i in range(len(lpaths)) ]

    return dict(zip(lpathnames, lpaths))


def trunc_isoforms_by_readlayout_SE(
    pathdict,
    read_start,
    read_length ):

    """Create "compatibility" group of isoforms cropped by a given single end read of fixed length and
    fixed starting position

    Args:
        pathdict (dict<str> -> tuple): isoform name to path dictionary
        read_start (int): fixed starting position of a read; 1-based, inclusive
        read_length (int): fixed read length

    Returns:
        - dict (str -> tuple): dictionary of isoform truncated by read position and length. Keys = compatbility group name, values = exon (or partial exon) starts and ends, in 1-based inclusive coordiantes

        - dict (str -> list of str) - within each named compatbility group (keys) --> which named isoforms from pathdict are equivalent (values)

    """

    # compat_grp --> ( visible path, [ isonames from paths input ] )

    m_iso_vispath = {}

    for pathname in pathdict:
        path = pathdict[pathname]
        found_ex = False
        for iex in range(int(len(path))):
            coord_ex = ( path[iex][0], path[iex][1] )
            if read_start >= coord_ex[0] and read_start <= coord_ex[1]:
                found_ex=True
                break

        if not found_ex:
            # this read start position does not exist on this path.
            m_iso_vispath[ pathname ] = None
            continue

        # truncate path to what is 'visible' given the proposed layout
        path_vis = []
        read_bp_remain = read_length

        in_first_ex = True
        while read_bp_remain > 0 and iex < int(len(path)):
            coord_ex = ( path[iex][0], path[iex][1] )

            if in_first_ex :
                coord_start = read_start
                in_first_ex = False
            else:
                coord_start = coord_ex[0]

            coord_end = min( coord_ex[1], coord_start+read_bp_remain-1 )

            path_vis.append( (coord_start,coord_end) )

            read_bp_remain -= (coord_end - coord_start + 1)

            iex += 1

        path_vis = tuple(path_vis) # change to hashable type

        m_iso_vispath[ pathname ] = path_vis

    m_vispath_name = {}
    m_name_vispath = {}

    n_vispath = 0
    for vispath in set(m_iso_vispath.values()):
        vispathname = 'isogrp{:02d}'.format(n_vispath)
        m_name_vispath[vispathname] = vispath
        m_vispath_name[vispath] = vispathname
        n_vispath+=1

    m_vispath_liso = {}
    for pathname in pathdict:
        vispath = m_iso_vispath[ pathname ]
        vispathname = m_vispath_name[vispath]

        if vispathname not in m_vispath_liso:
            m_vispath_liso[ vispathname ] = [ pathname ]
        else:
            m_vispath_liso[ vispathname] .append( pathname )

    return m_name_vispath, m_vispath_liso


def compute_isoform_counts(
    bam,
    isogrpdict,
    read_start_coord,
    tol_first_last=0,
    min_reads=10,
    count_otherisos=False
):
    """
    Create per-barcode counts of reads matching each isoform. Resulting dataset contains a count of total reads, unmapped
    reads, reads with a bad start (greater than 5 (default) basepairs away from most common start), how many reads
    match each isoform, and the psis for each isoform (matching reads/usable reads).

    Args:
        bam (pysam.AlignmentFile):  a barcode tag-grouped bam file of spliced alignments to reporter sequence
        isogrpdict (dict): isoform compatbility group name --> list of exon start, end coodinates (1-based, inclusive)
        read_start_coord (int): expected fixed read start position
        tol_first_last (int): allow up to this many bases tolerance at the first and last exons
        min_reads (int): require at least this many reads to be present per barcode group
        count_otherisos (bool): should reads not matching any of the known isoforms be counted in the total?

    Returns:
        Pandas data frame with per-barcode read counts
    """
    rowList = []

    ctr_bcs,ctr_reads=0,0

    for tag, _reads in itertools.groupby( bam, lambda _r:_r.get_tag( 'RX' )):

        reads=list(_reads)

        if len(reads)<min_reads:
            continue

        n_umapped=0
        n_badstart=0
        n_nomatch = 0
        ln_matches=[0] * len( isogrpdict )

        for read in reads:
            if read.is_unmapped:
                n_umapped+=1
            elif abs(read.reference_start+1-read_start_coord)>tol_first_last:
                n_badstart+=1
            else:
                cur_matches = check_junctions2( read, isogrpdict, tol_first_last )

                if sum(cur_matches)==0:
                    n_nomatch+=1
                else:
                    assert sum(cur_matches) <= 1, ( str(read),str(cur_matches) )

                    # for now, not dealing with read groups matching multiple isoform groups - this should not happen

                    for i in range(len(cur_matches)):
                        ln_matches[i]+=cur_matches[i]

        total = n_umapped + n_badstart + n_nomatch + sum(ln_matches)

        ctr_bcs+=1
        ctr_reads+=len(reads)

        # for debugging only go through a few
        # if ctr_bcs==10:break

        if ctr_bcs % 1000 == 0 :
            print('processed {} bcs, {} reads'.format( ctr_bcs, ctr_reads ))

        rowList.append( [tag,total,n_umapped,n_badstart,n_nomatch]+ln_matches )

    # from IPython.core.debugger import set_trace
    # set_trace()

    psidf=pd.DataFrame(rowList, columns=['barcode','num_reads','unmapped_reads','bad_starts','other_isoform']+list(isogrpdict.keys()))

    if not count_otherisos:
        psidf['usable_reads']=psidf.num_reads-(psidf.unmapped_reads+psidf.bad_starts+psidf.other_isoform)

        for iso in isogrpdict:
            psidf[iso+'_psi'] = psidf[iso] / psidf.usable_reads
    else:
        psidf['usable_reads']=psidf.num_reads-(psidf.unmapped_reads+psidf.bad_starts)

        for iso in isogrpdict:
            psidf[iso+'_psi'] = psidf[iso] / psidf.usable_reads

        psidf['other_isoform_psi'] = psidf['other_isoform'] / psidf.usable_reads

    psidf = psidf.set_index('barcode')

    return psidf


# check coverage of exon and of junctions
def check_junctions2( read, isogrpdict, tol_first_last=0 ):

    """Check an individual aligned read vs a list of isoforms

    Args:
        read (pysam.AlignedSegment): a single aligned read
        isogrpdict (dict): isoform compatbility group name --> list of exon start, end coodinates (1-based, inclusive)
        tol_first_last (int): allow up to this many bases tolerance at the first and last exons

    Returns:
        list of bools indicating whether this alignment is consistent with each of the provided isoform compat. groups.
    """

    # get blocks of reference coverage from the read.
    l_refcoord_blks = read.get_blocks()

    # coalesce any adjacent blocks; fix from [00) coordinate system to [11]
    l_refcoord_blks_joined = []
    for i in range(len(l_refcoord_blks)):
        l_refcoord_blks_joined.append( (l_refcoord_blks[i][0]+1,l_refcoord_blks[i][1]) )
        if len(l_refcoord_blks_joined) >= 2:
            if l_refcoord_blks_joined[-1][0] == l_refcoord_blks_joined[-2][1]+1:
                l_refcoord_blks_joined[-2:]=[ (l_refcoord_blks_joined[-2][0], l_refcoord_blks_joined[-1][1]) ]

    l_refcoord_blks_joined=tuple(l_refcoord_blks_joined)

    lmatches=[]

    # now compare to isogrps
    for isogrpname in isogrpdict:
        isogrp = isogrpdict[isogrpname]

        match = False
        lblks = len( l_refcoord_blks_joined )
        if lblks == len( isogrp ):
            # from IPython.core.debugger import set_trace
            # set_trace()

            # do the aligned blocks of this read strictly equal the isoform coord blocks?
            if isogrp == l_refcoord_blks_joined:
                # yes, great, it's a match
                match=True
            elif tol_first_last>0:
                # no, but check for mismatch at the very first or last base within the specified tol.
                match=all([c_read[0] in range(c_ref[0]-tol_first_last,c_ref[0]+tol_first_last)
                           and c_read[1] in range(c_ref[1]-tol_first_last,c_ref[1]+tol_first_last)
                           for c_ref,c_read in zip(isogrp,l_refcoord_blks_joined)])

                #commented out by CS - wasn't working properly

                """if lblks>1:
                    if (isogrp[1:-1] == l_refcoord_blks_joined[1:--1]) and \
                       l_refcoord_blks_joined[0][1]==isogrp[0][1] and \
                       l_refcoord_blks_joined[0][0]>=isogrp[0][0]-tol_first_last and \
                       l_refcoord_blks_joined[0][0]<=isogrp[0][0]+tol_first_last and \
                       l_refcoord_blks_joined[-1][1]>=isogrp[1][1]-tol_first_last and \
                       l_refcoord_blks_joined[-1][1]<=isogrp[1][1]+tol_first_last :

                        match=True
                else:
                    if l_refcoord_blks_joined[0][1]==isogrp[0][1] and \
                       l_refcoord_blks_joined[0][0]>=isogrp[0][0]-tol_first_last and \
                       l_refcoord_blks_joined[0][0]<=isogrp[0][0]+tol_first_last and \
                       l_refcoord_blks_joined[-1][1]>=isogrp[1][1]-tol_first_last and \
                       l_refcoord_blks_joined[-1][1]<=isogrp[1][1]+tol_first_last :
                        match=True"""

        lmatches.append( match )

    return lmatches

def filter_on_barc_len( jxnbybctbl, max_len=35 ):
    li = jxnbybctbl.index.str.len() <= max_len
    subtbl = jxnbybctbl.loc[li].copy()
    return subtbl


def combine_isogrps(
    new_grpnames_to_old_grpnames,
    jxnbybctbl
 ):
    """ combine barcode groups

    Arguments:
        new_grpnames_to_old_grpnames {[type]} -- [description]
    """

    newtbl = pd.DataFrame()
    for c in ['num_reads','unmapped_reads','bad_starts','other_isoform','usable_reads']:
        newtbl[c] = jxnbybctbl[c]

    for newgrp in new_grpnames_to_old_grpnames:
        oldgrp = new_grpnames_to_old_grpnames[ newgrp ]

        if type(oldgrp)==str:
            oldgrp=[oldgrp]
        else:
            assert type(oldgrp)==list

        if len(oldgrp)==1:
            newtbl[newgrp] = jxnbybctbl[oldgrp]
        else:
            newtbl[newgrp] = jxnbybctbl[oldgrp].sum(axis=1)

    for newgrp in new_grpnames_to_old_grpnames:
        newtbl[newgrp+'_psi'] = newtbl[newgrp] / newtbl['usable_reads']

    return newtbl
