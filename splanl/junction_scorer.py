import pysam
import os
import pandas as pd
import numpy as np
import networkx as nx
from collections import Counter
import matplotlib.pyplot as plt
import copy
import itertools
import pybedtools as pbt

def adjust_junctions(pysam_align,
                    refseq,
                    cnst_exons,
                    outfile):

    for bc, _reads in itertools.groupby( pysam_align, lambda _r: _r.get_tag( 'RX' ) ):

        reads = list(_reads)

        for read in reads:

            iso = read.get_blocks()

            #don't check junctions for unmapped reads or if there's only one block
            if read.is_unmapped or len( iso )==1:
                outfile.write( read )
                continue

            #check if the end of any junction is more downstream than the expected end of the downstream exon
            #check if the end of any junction is at the expected end of the downstream exon
            if any( jn[1] < cnst_exons[0][1] for jn in iso ) and not( any( jn[1] == cnst_exons[0][1] for jn in iso ) ):
                #grabs the last junction which is less than the expected end of the constant exon
                bad_jn = [ jn for jn in iso if jn[1] < cnst_exons[0][1] ][-1]
                #grabs junction following the screwed up junction
                next_jn = iso[ iso.index(bad_jn)+1 ]

                #if the next exonic region is smaller than the difference from the bad to expected junction
                #write the read as is and move to next read
                if (next_jn[1] - next_jn[0]) < (cnst_exons[0][1] - bad_jn[1]):
                    outfile.write( read )
                    continue

                #if its not also check if the suffix of whats missing from the constant exon matches the prefix of the next junction
                if refseq[ bad_jn[1]:cnst_exons[0][1] ] == refseq[ next_jn[0]:next_jn[0] + ( cnst_exons[0][1] - bad_jn[1] ) ]:
                    cigarlist = [ list(tup) for tup in read.cigartuples ]
                    #gets indices that are a match and need editing
                    match_idx = [ idx for idx, tup in enumerate( read.cigartuples ) if tup[0] == 0 ]
                    bad_jn_idx = match_idx[ iso.index(bad_jn) ]

                    assert len(match_idx) > iso.index(bad_jn), bc+'\nNo matched segments after bad junction'
                    next_jn_idx = match_idx[ iso.index(bad_jn) + 1 ]

                    #adds bases back to the bad junction to force it to end at the expected constant exon end
                    cigarlist[ bad_jn_idx ][1]+=cnst_exons[0][1] - bad_jn[1]
                    #subtract bases from the next matched segment
                    cigarlist[ next_jn_idx ][1]-=cnst_exons[0][1] - bad_jn[1]
                    read.cigartuples = tuple( tuple(l) for l in cigarlist )

            #check if the start of any junction is farther upstream than the expected start of the constant upstream exon
            #check if the start of any junction is at the expected start of the constant upstream exon
            elif any( jn[0] > cnst_exons[1][0] for jn in iso ) and not( any( jn[0] == cnst_exons[1][0] for jn in iso ) ):
                #grabs the first junction which is less than the expected end of the constant exon
                bad_jn = [ jn for jn in iso if jn[0] > cnst_exons[1][0] ][0]
                #grabs junction preceeding the screwed up junction
                prev_jn = iso[ iso.index(bad_jn)-1 ]

                #if the previous exonic region is smaller than the difference from the bad to expected junction
                #write the read as is and move to next read
                if (prev_jn[1] - prev_jn[0]) < (bad_jn[0] - cnst_exons[1][0]):
                    outfile.write( read )
                    continue

                #finally check if the prefix lost on the constant upstream matches the suffix on the previous junction
                if refseq[ cnst_exons[1][0]:bad_jn[0] ] == refseq[ prev_jn[1]-( bad_jn[0] - cnst_exons[1][0] ):prev_jn[1] ]:
                    cigarlist = [ list(tup) for tup in read.cigartuples ]
                    #gets indices that are a match and need editing
                    match_idx = [ idx for idx, tup in enumerate( read.cigartuples ) if tup[0] == 0 ]
                    bad_jn_idx = match_idx[ iso.index(bad_jn) ]

                    assert iso.index(bad_jn) > 0, bc+'\nNo matched segments before bad junction'

                    prev_jn_idx = match_idx[ iso.index(bad_jn) - 1 ]

                    #adds bases back to the bad junction to force it to start at the expected constant exon start
                    cigarlist[ bad_jn_idx ][1]+=bad_jn[0] - cnst_exons[1][0]
                    #removes bases from the previous matched segment
                    cigarlist[ prev_jn_idx ][1]-=bad_jn[0] - cnst_exons[1][0]
                    read.cigartuples = tuple( tuple(l) for l in cigarlist )

            #sanity check that we didn't alter the length of the read
            assert ( read.infer_query_length() == read.query_length ), bc+'\nWrong read length - invalid cigar'

            outfile.write( read )

    #make sure to flush the buffer!
    outfile.close()


def clean_jns(jns,
            cnst_exons):
    """Removes any junctions within the constant exons,
    joins adjacent junctions,
    and removes the last junction if its <= 3 bp long

    Args:
        jns (list of tuples): list of tuples from a PySam read_blocks() function call
        cnst_exons (list of tuples): coords of known constant exons

    Returns:
        jns_joined (list of tuples): same list of tuples as input with
        junctions within the constant exons removed
    """
    jns_joined = []

    prev_jn = None

    for jn in jns:
        if jn[0] >= cnst_exons[0][1] and jn[1] < cnst_exons[1][0]:
            jns_joined.append(jn)
            cur_jn = jn
            if prev_jn:
                if cur_jn[0] == prev_jn[1]+1:
                    jns_joined[-2:]=[ (prev_jn[0], cur_jn[1]) ]
            prev_jn=cur_jn

    #removes the last junction if its less than 3 bp long
    #hard to align a so few reads accurately so ignore them
    if len(jns_joined)>0 and ( jns_joined[-1][1] - jns_joined[-1][0] ) <= 3:
        jns_joined = jns_joined[:-1]

    return(jns_joined)

def make_full_bed(jn_df,
                    cnst_exons,
                    chrom_name):

    in_df = jn_df.copy()

    in_df = in_df.set_index('junction').sort_index()

    blnk_bed = { 'chrom': [chrom_name]*( in_df.shape[0]+2 ),
                 'start': [cnst_exons[0][0]],
                 'end': [cnst_exons[0][1]],
                 'name': ['upstream_cnst']}

    for i,jn in enumerate(in_df.index):
        blnk_bed['start'].append(jn[0])
        blnk_bed['end'].append(jn[1])
        blnk_bed['name'].append('junction_'+str(i))

    blnk_bed['start'].append(cnst_exons[1][0])
    blnk_bed['end'].append(cnst_exons[1][1])
    blnk_bed['name'].append('downstream_cnst')

    bed_df = pd.DataFrame(blnk_bed)

    return(bed_df)

def remove_cnst_ex(blk,
                    cnst_exons):

    cnst_min = min( min( cnst_exons ) )
    cnst_max = max( max( cnst_exons ) )

    for jn in copy.deepcopy(blk):
        #if the junction starts before the start of the constant exon remove it
        if jn[0] <= cnst_min:
            blk.remove(jn)
            continue
        #if the junction ends after the end of the constant exon remove it
        if jn[1] >= cnst_max:
            blk.remove(jn)
            continue
        for ex in cnst_exons:
            if jn[0] in range(ex[0],ex[1]) or jn[1] in range(ex[0],ex[1]):
                blk.remove(jn)
                break

    return(blk)


def get_all_isoforms(pysam_align,
                    cnst_exons):
    """Gets all isoforms within a pysam alignment file

    Args:
        pysam_align (PySam AlignmentFile): a pysam alignment file with all of the reads for the sample

    Returns:
        iso_df (pandas dataframe): Dataframe containing all isoforms seen within the data along with the
        number of reads associated with that isoform
    """

    #creates a counter object of all of the exon junctions
    all_isos_cnt = Counter( [ tuple( clean_jns( read.get_blocks(), cnst_exons ) ) for read in pysam_align ] )

    iso_df = pd.DataFrame.from_dict(all_isos_cnt, orient='index').reset_index()
    iso_df = iso_df.rename(columns={'index':'isoform', 0:'read_count'}).sort_values(by='read_count', ascending=False)

    return(iso_df)

def number_isoforms(iso_df):

    out_tbl = iso_df.copy()

    out_tbl['iso_num'] = ['iso'+str(i).zfill(4) for i in range(iso_df.shape[0])]

    out_tbl = out_tbl.set_index('iso_num')

    return out_tbl

def check_tolerance(var_pos,
                    unique_jns,
                    tol):

    within_tol = False

    for jn in unique_jns:
        if var_pos >= ( jn - tol ) and var_pos <= ( jn + tol ):
            within_tol = True
            break

    return within_tol

def cluster_vars(var_bc_counts,
                min_cluster_size,
                tol,
                clustered_bc_ave):

    clusters = {}
    i = 0

    for variants in var_bc_counts:
        for var in variants.split(','):
            var_pos =  int( var.split(':')[1] )
            if not clusters:
                clusters[ 'cluster'+str(i) ] = [ (var_pos,), var_bc_counts[ variants ] ]
                i+=1
                continue
            for cluster, info in clusters.items():
                var_clustered = False
                #if the variant is nearby first variant in cluster add to cluster
                if check_tolerance( var_pos, [ info[0][0] ], tol):
                    clusters[ cluster ][0]+= ( var_pos, )
                    clusters[ cluster ][1] += var_bc_counts[ variants ]
                    var_clustered = True
                    #if we've reached the minimum cluster count and the average bc count is at threshold exit
                    if len( clusters[ cluster ][0] ) == min_cluster_size and \
                    ( clusters[ cluster ][1] / len( clusters[ cluster ][0] ) ) >=  clustered_bc_ave:
                        return( True )
                    #either way we clustered the variant so break outta the loop
                    break
                #otherwise make a new cluster for it
            if not var_clustered:
                clusters[ 'cluster'+str(i) ] = [ (var_pos,), var_bc_counts[ variants ] ]
                i+=1

    return( False )


def summarize_isos_by_var_bc(pysam_align,
                             cnst_exons,
                                satbl,
                                iso_df,
                                unique_jns,
                                print_count=5,
                                min_maxbc_count=100,
                                tol=10,
                                min_cluster_size=3,
                                clustered_bc_ave=3,
                                min_bc_max_reads=(2,2)):
    """Gets counts for number of variants and barcodes per isoform. Prints variants with the most number of barcodes.

    Args:
        pysam_align (PySam AlignmentFile): a tag sorted pysam alignment file with all of the reads for the sample
        cnst_exons (list of tuples): coords of known constant exons
        satbl (pandas dataframe): subassembly dataframe which maps barcodes to variants
                                    Indexed on readgroupid
        iso_df (pandas dataframe): dataframe of isoforms, isoform groups, and read counts created from get_all_isoforms
                                    and create_isogrps
                                    Indexed on isogrp
        print_count (int): number of desired variants to print out per isoform - program will print the variants
                            with the highest number of barcodes

    Returns:
        iso_df (pandas dataframe): Original iso_df with the addition of number of variants per isoform and number
                                    of barcodes per isoform
    """
    out_tbl = iso_df.copy()
    satbl_c = satbl.copy()

    satbl_c = satbl_c.dropna(subset=['variant_list'])

    isogrp_to_bc = {}

    bc_cnt, read_cnt = 0,0

    for bc, _reads in itertools.groupby( pysam_align, lambda _r: _r.get_tag( 'RX' ) ):

        if bc not in satbl_c.index:
            continue

        for iso, reads in itertools.groupby( _reads, lambda _r: tuple( clean_jns( _r.get_blocks(), cnst_exons ))):

            r=list(reads)

            isogrp = out_tbl.loc[out_tbl['isoform']==iso].index

            assert len(isogrp)==1, 'Isoforms are not unique'

            isogrp=isogrp[0]

            if isogrp not in isogrp_to_bc:
                isogrp_to_bc[ isogrp ]={}

            if bc not in isogrp_to_bc[ isogrp ]:
                isogrp_to_bc[ isogrp ][ bc ] = len(r)
            else:
                isogrp_to_bc[ isogrp ][ bc ] += len(r)

            read_cnt+=len(r)

        bc_cnt+=1

        if bc_cnt%10000==0:
            print('Barcodes processed:',bc_cnt,'Reads processed:',read_cnt)

    #takes 16 minutes to get here
    isogrp_var_count = { isogrp: [ len( isogrp_to_bc[ isogrp ] ),
                                    #should count number of variants
                                    len( set( satbl_c.loc[ bc ].variant_list
                                    for bc in isogrp_to_bc[ isogrp ] ) ),
                                    max( isogrp_to_bc[ isogrp ].values() )
                                  ]
                        for isogrp in isogrp_to_bc
                        }

    chucked_bcs, chucked_reads = 0,0
    missed_bcs, missed_reads = 0,0

    for isogrp in isogrp_to_bc:

        bcs = [ bc for bc in isogrp_to_bc[ isogrp ] ]

        var_to_bc = {}
        for bc in bcs:
            var = satbl_c.loc[ bc ].variant_list
            if var in var_to_bc:
                var_to_bc[ var ].add(bc)
            else:
                var_to_bc[ var ] = { bc }

        var_bc_counts = { var: len(var_to_bc[ var ]) for var in var_to_bc }

        var_bc_counts_sort = { var: count for var, count in sorted(var_bc_counts.items(), key=lambda item: -item[1])}

        #get the variant with the highest number of barcodes
        top_var = list( itertools.islice( var_bc_counts_sort.items(), 1 ) ) [0]

        #get max bcs from one variant for each isogrp
        isogrp_var_count[ isogrp ].append( top_var[1] )

        print( 'For isoform:', out_tbl.loc[ isogrp ].isoform )
        print( 'The variants with the top', print_count, 'number of barcodes are:')
        print( list(itertools.islice( var_bc_counts_sort.items(), print_count ) ) )

        #get a list of all variants
        var_list = [ v for var in var_bc_counts_sort for v in var.split(',') ]
        var_list_unique = list( set( var_list ) )

        #this gets the first position where the isoform differs from expected junctions
        jn_diff = list( { jn for jns in out_tbl.loc[isogrp].isoform for jn in jns }
                        - set( unique_jns ) )
        if len( jn_diff )>1:
            jn_diff = [ min( jn_diff ) ]

        #lets keep the best isoforms
        #first check if the max bc per isoform is greater than 100
        if top_var[1] >= min_maxbc_count:
            isogrp_var_count[ isogrp ].append( 1 )
        #requires the top variant to have at least x barcodes and the top barcode to have at least y reads
        elif top_var[1] <= min_bc_max_reads[0] and isogrp_var_count[ isogrp ][2] <= min_bc_max_reads[1]:
            chucked_bcs += isogrp_var_count[ isogrp ][0]
            chucked_reads += out_tbl.loc[ isogrp ].read_count
            isogrp_var_count[ isogrp ].append( 0 )
            continue
        #check if the var with the max bc is within tolerance of the where the isoform differs
        elif any( [ check_tolerance( int( var.split(':')[1] ),
                                 jn_diff,
                                 tol )
                    for var in top_var[0].split(',') ] ):
            isogrp_var_count[ isogrp ].append( 2 )
        #check if at least min_cluster_size of variants are within a tol of each other
        elif cluster_vars(  var_bc_counts_sort,
                            min_cluster_size,
                            tol,
                            clustered_bc_ave ):
            isogrp_var_count[ isogrp ].append( 3 )
        #check if a double variant is also listed as a single variant
        elif len( var_list ) != len( var_list_unique ):
            isogrp_var_count[ isogrp ].append( 4 )
        else:
            isogrp_var_count[ isogrp ].append( 0 )
            missed_bcs += isogrp_var_count[ isogrp ][0]
            missed_reads += out_tbl.loc[ isogrp ].read_count

    print('%i (%.2f%%) barcodes failed the min_bc_max_reads filter' \
            % ( chucked_bcs, 100*(chucked_bcs/sum(i[0] for i in isogrp_var_count.values() ) ) ) )
    print('%i (%.2f%%) reads failed the min_bc_max_reads filter' \
            % ( chucked_reads, 100*(chucked_reads/out_tbl.read_count.sum() ) ) )
    print('%i (%.2f%%) barcodes did not fulfill any filter' \
            % ( missed_bcs, 100*(missed_bcs/sum(i[0] for i in isogrp_var_count.values() ) ) ) )
    print('%i (%.2f%%) reads did not fulfill any filter'\
            % ( missed_reads, 100*(missed_reads/out_tbl.read_count.sum() ) ) )

    iso_cnts = pd.DataFrame.from_dict(isogrp_var_count, orient='index',
                                    columns=['num_bcs','num_vars','max_reads_per_bc','max_bc_per_var','filter'])

    out_tbl = pd.merge(out_tbl, iso_cnts, left_index=True, right_index=True).sort_values(by=['read_count'],ascending=False)

    return(out_tbl)

def combine_isoforms(iso_df,
                    cryptic_exons):

    out_tbl = iso_df.copy()

    out_tbl.drop( columns = ['filter'], inplace=True )

    isoform_comb = []

    for iso in out_tbl.isoform:

        cur_iso = []
        prev_jn = None

        for jn in iso:

            #don't add junctions overlapping the cryptic_exons
            if not ( any( [ jn[0] in range( ce[0],ce[1] ) or jn[1] in range( ce[0],ce[1] )
                            for ce in cryptic_exons ] ) ):
                cur_iso.append( jn )
                cur_jn = jn
                if prev_jn:
                    #combine adjacent junctions
                    if cur_jn[0] == prev_jn[1]:
                        cur_iso[-2:]=[ (prev_jn[0], cur_jn[1]) ]
                prev_jn=cur_jn

        cur_iso = tuple( cur_iso )
        isoform_comb.append( cur_iso )

    out_tbl[ 'comb_isoform' ] = isoform_comb

    #sums across matching isoforms
    #print(out_tbl.groupby(['comb_isoform'],as_index=False).max())
    sums_df = out_tbl.groupby(['comb_isoform'],as_index=False).sum().drop(columns=['max_reads_per_bc','max_bc_per_var'])
    #max across matching isoforms
    max_df = out_tbl.groupby(['comb_isoform'],as_index=False).max().drop(columns=['read_count','num_bcs','num_vars','isoform'])
    #counts number of isoforms per isogrp
    counts_df = out_tbl.groupby(['comb_isoform'],as_index=False).count()[['comb_isoform','isoform']].rename(columns={"isoform": "isoform_counts"})
    #concatenates isoform numbers to track from previous tables
    iso_nums_df = out_tbl.reset_index().groupby(['comb_isoform'])['index'].apply(','.join).reset_index()
    #concatenates all isoforms to track from previous table
    isoform_df = out_tbl.reset_index().groupby(['comb_isoform'])['isoform'].apply(tuple).reset_index()

    out_tbl = pd.merge(sums_df,max_df,on='comb_isoform')
    out_tbl = pd.merge(out_tbl,counts_df,on='comb_isoform')
    out_tbl = pd.merge(out_tbl,isoform_df,on='comb_isoform')
    out_tbl = pd.merge(out_tbl,iso_nums_df,on='comb_isoform').rename(columns={"index": "comb_iso_nums"})
    out_tbl.sort_values(by=['read_count'], ascending=False, inplace=True)

    out_tbl['index'] = ['isogrp'+str(i).zfill(3) for i in range(out_tbl.shape[0])]
    out_tbl.set_index('index',inplace=True)

    return (out_tbl)

def check_individual_variants(pysam_align,
                                bc_list,
                                suspicious_isos=[]):

    iso_dict = {}
    full_reads = {}

    for bc, _reads in itertools.groupby( pysam_align, lambda _r: _r.get_tag( 'RX' ) ):

        if bc not in bc_list:
            continue

        reads = list(_reads)

        for r in reads:
            iso = tuple( r.get_blocks() )
            if iso not in iso_dict:
                iso_dict[ iso ] = {bc: 1}
            elif bc not in iso_dict[ iso ]:
                iso_dict[ iso ][ bc ] = 1
            else:
                iso_dict[ iso ][ bc ] += 1
            if iso in suspicious_isos:
                if iso not in full_reads:
                    full_reads[ iso ] = [ r ]
                else:
                    full_reads[ iso ].append(r)

    return([iso_dict, full_reads])

def add_junctions(iso,
                    con_exons,
                    read_len,
                    start_site):

    #add in first constant exon
    iso = ( ( start_site, con_exons[0][1] ), ) + iso

    used_bases = sum( jn[1]-jn[0] for jn in iso )

    if used_bases > read_len:
        print(iso)

    #assert used_bases <= read_len, str(iso)+'\nMore bases in isoform than read length'

    #if we still have bases left, add final constant exon
    if used_bases < read_len:
        iso_end = con_exons[1][0] + (read_len - used_bases)
        iso = iso + ( ( con_exons[1][0], iso_end ), )

    #make 1 - based inclusive
    isos_out = tuple( [ ( jn[0]+1, jn[1] ) for jn in iso ] )

    return( isos_out )


def create_iso_dict(iso_df,
                    con_exons,
                    read_len,
                    start_site):

    iso_df_c = iso_df.copy()

    isogrpdict = { iso: add_junctions( iso_df_c.loc[ iso ].isoform,
                                        con_exons,
                                        read_len,
                                        start_site )
                    for iso in iso_df_c.index}

    return( isogrpdict )

def create_iso_dict_no_cnst(iso_df):

    iso_df_c = iso_df.copy()

    #make 1-based inclusive
    isogrpdict = { iso: tuple( [ ( jn[0]+1, jn[1] )
                    for jn in iso_df_c.loc[ iso ].isoform ] )
                    for iso in iso_df_c.index}

    return( isogrpdict )


def filter_junctions(jn_df,
                    thresh):
    """Removes junctions with read counts lower than the threshold.

    Args:
        jn_df (Pandas data frame): Pandas data frame containing junctions and their read col_read_counts
        thresh (int): Number of reads required to remain in the data frame

    Returns:
        jn_df_filt (Pandas data frame): Pandas data frame with only junctions that meet the read count threshold
    """
    return( jn_df[ jn_df[ 'read_count' ]>=thresh ] )

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
    #print(read)
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
    #print('data',l_refcoord_blks_joined)
    lmatches=[]

    # now compare to isogrps
    for isogrpname in isogrpdict:
        isogrp = isogrpdict[isogrpname]
        #print(isogrpname,isogrp)

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
                match=all( [ ( c_read[0] in range( c_ref[0]-tol_first_last,c_ref[0]+tol_first_last )
                           and c_read[1] in range( c_ref[1]-tol_first_last,c_ref[1]+tol_first_last ) )
                           if ( c_ref==isogrp[0] or c_ref==isogrp[-1] )
                           else ( c_read == c_ref )
                           for c_ref,c_read in zip(isogrp,l_refcoord_blks_joined) ] )
                #print(match)

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
