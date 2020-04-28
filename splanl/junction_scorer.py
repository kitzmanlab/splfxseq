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
from ast import literal_eval

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

    for _jn in jns:
        #changes to 1-based inclusive numbering
        jn = (_jn[0]+1, _jn[1])
        if jn[0] >= cnst_exons[0][1] and jn[1] < ( cnst_exons[1][0]+1 ):
            jns_joined.append(jn)
            cur_jn = jn
            #joins adjacent junctions
            if prev_jn:
                if cur_jn[0] == prev_jn[1]+1:
                    jns_joined[-2:]=[ (prev_jn[0], cur_jn[1]) ]
            prev_jn=cur_jn

    #removes the last junction if its less than 3 bp long
    #hard to align a so few reads accurately so ignore them
    if len(jns_joined)>0 and ( jns_joined[-1][1] - jns_joined[-1][0] ) <= 3:
        jns_joined = jns_joined[:-1]

    return(jns_joined)

def get_all_isoforms(align_by_samp_dict,
                    cnst_exons):
    """Gets all isoforms within a pysam alignment file

    Args:
        align_dict (dictionary): dictionary with sample names as keys and pysam alignment files as values

    Returns:
        out_tbls (dictionary): dictionary with sample names as keys and pandas dataframes with isoform tuples as the
        index and number of reads representing each isoform as the only column
    """
    out_tbls = {}

    for samp in align_by_samp_dict:

        print(samp)

        #creates a counter object of all of the exon junctions
        all_isos_cnt = Counter( [ tuple( clean_jns( read.get_blocks(), cnst_exons ) ) for read in align_by_samp_dict[ samp ] ] )

        iso_df = pd.DataFrame.from_dict( all_isos_cnt, orient='index' ).reset_index()
        iso_df = iso_df.rename(columns={'index':'isoform', 0:'read_count'}).sort_values(by='read_count', ascending=False)
        iso_df = iso_df.set_index('isoform')

        out_tbls[ samp ] = iso_df

    return( out_tbls )

def number_and_merge_isoforms(isodf_by_samp_dict):

    out_tbl = { samp+'_read_count': [] for samp in isodf_by_samp_dict }

    out_tbl[ 'isoform' ] = list( set ( [ iso for samp in isodf_by_samp_dict
                                              for iso in isodf_by_samp_dict[ samp ].index ] ) )

    #number all isoforms so they are of the form iso0001
    out_tbl[ 'isonum' ] = [ 'iso'+str(i).zfill( len( str( len( out_tbl[ 'isoform' ] ) ) ) )
                            for i in range( len( out_tbl[ 'isoform' ] ) ) ]

    for samp in isodf_by_samp_dict:

        print(samp)

        #set all the read counts as zero
        out_tbl[ samp+'_read_count' ] = [ 0 for i in range( len( out_tbl[ 'isoform' ] ) ) ]

        for idx, iso in enumerate( isodf_by_samp_dict[ samp ].index ):

            out_tbl[ samp+'_read_count' ][ out_tbl[ 'isoform' ].index( iso ) ] = isodf_by_samp_dict[ samp ].iloc[ idx ].read_count

    out_tbl = pd.DataFrame( out_tbl )

    #reorder columns so isoform is first column
    ordered_col = ['isoform'] + [ col for col in out_tbl.columns if col != 'isoform' ]
    out_tbl = out_tbl[ ordered_col ]

    out_tbl = out_tbl.set_index('isonum')

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


def summarize_isos_by_var_bc(align_by_samp_dict,
                             cnst_exons,
                                satbl,
                                iso_df,
                                unique_jns,
                                canonical_isos,
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

    for samp in align_by_samp_dict:

        print(samp)

        isogrp_to_bc = {}

        bc_cnt, read_cnt = 0,0

        for bc, _reads in itertools.groupby( align_by_samp_dict[ samp ], lambda _r: _r.get_tag( 'RX' ) ):

            if bc not in satbl_c.index:
                continue

            for iso, reads in itertools.groupby( _reads, lambda _r: tuple( clean_jns( _r.get_blocks(), cnst_exons ))):

                r=list(reads)

                isogrp = out_tbl.loc[out_tbl['isoform']==iso].index

                assert len(isogrp)==1, 'The isoform matches '+str( len( isogrp ) )+' isogroups instead of the expected 1'

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
            #first check if the isoform exactly matches an isoform we're expecting to see
            if any( [ out_tbl.loc[ isogrp ].isoform == can_iso for can_iso in canonical_isos ] ):
                isogrp_var_count[ isogrp ].append( 1 )
            #then check if the max bc per isoform is greater than 100
            elif top_var[1] >= min_maxbc_count:
                isogrp_var_count[ isogrp ].append( 2 )
            #requires the top variant to have at least x barcodes and the top barcode to have at least y reads
            elif top_var[1] <= min_bc_max_reads[0] and isogrp_var_count[ isogrp ][2] <= min_bc_max_reads[1]:
                chucked_bcs += isogrp_var_count[ isogrp ][0]
                chucked_reads += out_tbl.loc[ isogrp ][ samp+'_read_count' ]
                isogrp_var_count[ isogrp ].append( 0 )
                continue
            #check if the var with the max bc is within tolerance of the where the isoform differs
            elif any( [ check_tolerance( int( var.split(':')[1] ),
                                 jn_diff,
                                 tol )
                    for var in top_var[0].split(',') ] ):
                isogrp_var_count[ isogrp ].append( 3 )
            #check if at least min_cluster_size of variants are within a tol of each other
            elif cluster_vars(  var_bc_counts_sort,
                            min_cluster_size,
                            tol,
                            clustered_bc_ave ):
                isogrp_var_count[ isogrp ].append( 4 )
            #check if a double variant is also listed as a single variant
            elif len( var_list ) != len( var_list_unique ):
                isogrp_var_count[ isogrp ].append( 5 )
            else:
                isogrp_var_count[ isogrp ].append( 0 )
                missed_bcs += isogrp_var_count[ isogrp ][0]
                missed_reads += out_tbl.loc[ isogrp ][ samp+'_read_count' ]

        print('%i (%.2f%%) barcodes failed the min_bc_max_reads filter' \
            % ( chucked_bcs, 100*(chucked_bcs/sum(i[0] for i in isogrp_var_count.values() ) ) ) )
        print('%i (%.2f%%) reads failed the min_bc_max_reads filter' \
            % ( chucked_reads, 100*(chucked_reads/out_tbl[ samp+'_read_count' ].sum() ) ) )
        print('%i (%.2f%%) barcodes did not fulfill any filter' \
            % ( missed_bcs, 100*(missed_bcs/sum(i[0] for i in isogrp_var_count.values() ) ) ) )
        print('%i (%.2f%%) reads did not fulfill any filter'\
            % ( missed_reads, 100*(missed_reads/out_tbl[ samp+'_read_count' ].sum() ) ) )

        cols = [ samp + suffix for suffix in ['_num_bcs','_num_vars','_max_reads_per_bc','_max_bc_per_var','_filter'] ]
        iso_cnts = pd.DataFrame.from_dict(isogrp_var_count, orient='index',
                                    columns=cols)

        out_tbl = pd.merge(out_tbl, iso_cnts, left_index=True, right_index=True, how="outer")

        out_tbl.index.name = 'isonum'

    #if iso_counts doesn't contain all the isoforms, the values will be null for those rows
    #change to zeros
    out_tbl = out_tbl.fillna(0)

    #get totals across samples
    for suffix in [ '_read_count', '_num_bcs', '_num_vars' ]:

        #sum of each row across the columns
        col = [ col for col in out_tbl.columns if suffix in col ]
        out_tbl[ 'total'+suffix ] = out_tbl[ col ].sum( axis=1 )

    #counts number of samples passing filter for each isoform
    #specifically, counting non zero filter values
    filter_col = [ col for col in out_tbl.columns if 'filter' in col ]
    out_tbl[ 'total_passfilt' ] = out_tbl[ filter_col ].astype( bool ).sum( axis=1 )

    return(out_tbl)

def combine_isoforms(iso_df,
                    cryptic_exons):

    out_tbl = iso_df.copy()

    filter_cols = [ col for col in out_tbl.columns if 'filter' in col ]
    out_tbl.drop( columns = filter_cols, inplace=True )

    isoform_comb = []

    for iso in out_tbl.isoform:

        #this will make sure any strings are converted to tuples
        iso = literal_eval(iso)

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
    max_cols = [ col for col in out_tbl.columns if 'max' in col ]
    sums_df = out_tbl.groupby(['comb_isoform'],as_index=False).sum().drop(columns=max_cols)
    #max across matching isoforms
    other_cols = [ col for col in out_tbl.columns for name in ['read_count','num_bcs','num_vars','isoform'] if name in col and col!= 'comb_isoform' ]
    max_df = out_tbl.groupby(['comb_isoform'],as_index=False).max().drop(columns=other_cols)
    #counts number of isoforms per isogrp
    counts_df = out_tbl.groupby(['comb_isoform'],as_index=False).count()[['comb_isoform','isoform']].rename(columns={"isoform": "isoform_counts"})
    #concatenates isoform numbers to track from previous tables
    #iso_nums_df = out_tbl.reset_index().groupby(['comb_isoform'])['index'].apply(','.join).reset_index()
    iso_nums_df = out_tbl.reset_index().groupby(['comb_isoform'])['isonum'].apply(','.join).reset_index()
    #concatenates all isoforms to track from previous table
    isoform_df = out_tbl.reset_index().groupby(['comb_isoform'])['isoform'].apply(tuple).reset_index()

    print('sums',sums_df)
    print('max',max_df)
    out_tbl = pd.merge(sums_df,max_df,on='comb_isoform')
    print('one',out_tbl.head())
    out_tbl = pd.merge(out_tbl,counts_df,on='comb_isoform')
    print('two',out_tbl.head())
    out_tbl = pd.merge(out_tbl,isoform_df,on='comb_isoform')
    print('three',out_tbl.head())
    out_tbl = pd.merge(out_tbl,iso_nums_df,on='comb_isoform').rename(columns={"index": "comb_iso_nums"})
    print('four',out_tbl.head())
    out_tbl = out_tbl.sort_values(by='total_read_count', ascending=False)
    print('five',out_tbl.head())


    out_tbl['isogrp'] = ['isogrp'+str(i).zfill( len( str( out_tbl.shape[0] ) ) ) for i in range( out_tbl.shape[0] ) ]
    out_tbl.set_index('isogrp',inplace=True)

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

def create_iso_dict_no_cnst(iso_df):

    iso_df_c = iso_df.copy()

    isogrpdict = { iso: iso_df_c.loc[ iso ].isoform
                    for iso in iso_df_c.index }

    return( isogrpdict )

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
    cnst_exons,
    tol_first_last=0,
    min_reads=1,
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
            #if the end of the upstream constant exon isn't within the tolerance count it as a bad start
            elif not( any( jn[1] <= cnst_exons[ 0 ][ 1 ] + tol_first_last
                            or jn[1] > cnst_exons[ 0 ][ 1 ] - tol_first_last
                            for jn in read.get_blocks() ) ):
                n_badstart+=1
            else:
                cur_matches = check_junctions2( read, isogrpdict, cnst_exons, tol_first_last )

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
def check_junctions2( read, isogrpdict, cnst_exons, tol_first_last=0 ):

    """Check an individual aligned read vs a list of isoforms

    Args:
        read (pysam.AlignedSegment): a single aligned read
        isogrpdict (dict): isoform compatbility group name --> list of exon start, end coodinates (1-based, inclusive)
        tol_first_last (int): allow up to this many bases tolerance at the first and last exons

    Returns:
        list of bools indicating whether this alignment is consistent with each of the provided isoform compat. groups.
    """
    # get blocks of reference coverage from the read.
    l_refcoord_blks_cleaned = clean_jns( read.get_blocks(), cnst_exons )

    l_refcoord_blks_cleaned=tuple(l_refcoord_blks_cleaned)

    lmatches=[]

    # now compare to isogrps
    for isogrpname in isogrpdict:
        isogrp = isogrpdict[isogrpname]
        #print(isogrpname,isogrp)

        match = False
        lblks = len( l_refcoord_blks_cleaned )
        if lblks == len( isogrp ) and isogrp == l_refcoord_blks_cleaned:
            match=True

        lmatches.append( match )

    return lmatches

def filter_on_barc_len( jxnbybctbl, max_len=35 ):
    li = jxnbybctbl.index.str.len() <= max_len
    subtbl = jxnbybctbl.loc[li].copy()
    return subtbl

def create_named_isogrps(iso_grp_dict,
                        iso_names_dict,
                        remove_exon_coords,
                        upstream_ex_len,
                        read_len,
                        tol = 0):

    named_iso_grps = { isoname: [] for isoname in iso_names_dict }
    #set up an other category for isoforms that don't match
    named_iso_grps['OTHER']=[]

    for isonum,iso in iso_grp_dict.items():

        iso_cleaned = []
        for jn in iso:
            #check if the junction overlaps the exons to remove
            if not( any ( jn[0] in range( r_exon[0], r_exon[1] ) or jn[1] in range( r_exon[0], r_exon[1] )
                    for r_exon in remove_exon_coords ) ):
                iso_cleaned.append(jn)

        #more likely to be off by a base if we are near the read length
        bases_used = sum( jn[1] - jn[0] +1 for jn in iso_cleaned ) + upstream_ex_len

        #this will fail if there's two possible isoform matches
        #would that ever happen?
        match=False
        for isoname, iso_to_match in iso_names_dict.items():
            #check if it perfectly matches the named isoform
            if iso_cleaned == iso_to_match:
                match = isoname
                break

            #if they're different lengths don't bother going through all the tolerance checks
            elif len( iso_to_match ) != len( iso_cleaned ):
                continue

            #if we have are using a tolerance and we are near the read length
            elif tol > 0 and bases_used in range( read_len - tol, read_len +tol ):
                #check if the final junction is within the tolerance (allows for deletions in the upstream exon)

                #gets index of all mismatched junctions
                mismatched_jn = [ iso_cleaned.index( jn_c )
                                for jn_m, jn_c in zip( iso_to_match, iso_cleaned )
                                if jn_m != jn_c ]

                #we only want to allow mismatches at the end of the last junction
                #so first check the only mismatch is the last junction
                #then check the start of the junction matches
                #then check if the final junction is within tolerance
                if mismatched_jn[0] == len( iso_cleaned ) - 1 \
                and iso_cleaned[-1][0] == iso_to_match[-1][0] \
                and iso_cleaned[-1][1] in range(iso_to_match[-1][1] - tol, iso_to_match[-1][1] + tol):
                    match = isoname
                    break

        if match:
            named_iso_grps[ match ].append( isonum )
        else:
            named_iso_grps[ 'OTHER' ].append( isonum )

    return(named_iso_grps)

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
