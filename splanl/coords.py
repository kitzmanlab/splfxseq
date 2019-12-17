

def pos_to_hgvspos(
    lvecpos,
    vec_corange_cloned,
    vec_corange_exons,
    cdna_corange_exons ):

    """
    Convert from vector position to HGVS cDNA location

    Arguments:
        lvecpos {pandas dataframe column} -- pandas dataframe column with the vector positions
        vec_corange_cloned {tuple of ints} -- range of vector coordinates which contain the cloned insert
        vec_corange_exons {list of tuples of ints} -- range of vector coordinates for each exon
        cdna_corange_exons {list of tuples of ints} -- range of cdna coordinates for each exon

    Returns:
        list -- list of cdna coordinates
    """


    assert len(vec_corange_exons)==len(cdna_corange_exons), 'must be same # of exons in both parameters'

    for (corng_vec,corng_ex) in zip(vec_corange_exons,cdna_corange_exons):
        assert corng_ex[1]-corng_ex[0]+1 == corng_vec[1]-corng_vec[0]+1 , 'exons must be same lengths'
        assert corng_vec[0] >= vec_corange_cloned[0] and corng_vec[1] <= vec_corange_cloned[1], 'exon must be entirely within cloned region'

    loutpos=[]

    for vecpos in lvecpos:
        if vecpos < vec_corange_cloned[0] or vecpos > vec_corange_cloned[1]:
            coordstr=None
        else:
            # is there only one exon?
            if len( vec_corange_exons ) == 1:
                # is the location before the exon, after, or inside?
                if vecpos < vec_corange_exons[0][0]:
                    coordstr = 'c.{:d}-{:d}'.format(
                        cdna_corange_exons[0][0],
                        vec_corange_exons[0][0]-vecpos )
                elif vecpos > vec_corange_exons[0][1]:
                    coordstr = 'c.{:d}+{:d}'.format(
                        cdna_corange_exons[0][1],
                        vecpos-vec_corange_exons[0][1] )
                else:
                    coordstr = 'c.{:d}'.format( cdna_corange_exons[0][0] + vecpos - vec_corange_exons[0][0] )
            else:
                # not handling multiple exons yet; if in an intron, we'd need to figure out which is nearer, then create the coordinate relative to that
                assert 1==0, 'not implemented yet'

        loutpos.append(coordstr)

    return loutpos
