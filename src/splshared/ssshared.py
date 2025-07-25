import numpy as np
import pandas as pd
from collections import OrderedDict

def zopen(fn,mode,encoding='utf-8'):
    # return gzip.GzipFile(fn,mode=mode) if fn.endswith('gz') else open(fn,mode=mode)
    return gzip.open(fn,mode=mode,encoding=encoding) if fn.endswith('gz') else open(fn,mode=mode,encoding=encoding)

rcMapRNA = {'A':'U', 'C': 'G', 'G': 'C', 'U': 'A', 'T': 'A', 'N': 'N', '-': '-', ' ': ' ',
            'a':'u', 'c': 'g', 'g': 'c', 'u': 'a', 't': 'a', 'n': 'n'}
rcMapDNA = {'A':'T', 'C': 'G', 'G': 'C', 'U': 'A', 'T': 'A', 'N': 'N', '-': '-', ' ': ' ',
            'a':'t', 'c': 'g', 'g': 'c', 'u': 'a', 't': 'a', 'n': 'n'}

rcMapDNAAmbigs = { 'A':'T', 'C':'G', 'G':'C', 'T':'A',
                   'M':'K', 'R':'Y', 'W':'W', 'S':'S',
                   'Y':'R', 'K':'M', 'V':'B', 'H':'D',
                   'D':'H', 'B':'V', 'N':'N' }

def revComp(seq, useRNA=False):
    return ''.join( [useRNA and rcMapRNA[b] or rcMapDNA[b] for b in reversed(seq)])

def isWithin(cooCoord, cooRange, fIsInclusive=True):
    if fIsInclusive: return (cooCoord <= max(cooRange) and cooCoord >= min(cooRange))
    else:            return (cooCoord <  max(cooRange) and cooCoord >  min(cooRange))

def chrStartStop(s):
    return s.split(':')[0],int(s.split(':')[1].split('-')[0]),int(s.split(':')[1].split('-')[1])
    
def chrTupStartStop(s):
    return s.split(':')[0],(int(s.split(':')[1].split('-')[0]),int(s.split(':')[1].split('-')[1]))

def coordConversion( typeFrom='[01]', typeTo='[00]' ):
    intvCorrectionFrom=[0,0]
    intvCorrectionFrom[0] = ( typeFrom[0]=='(' and 1 or
                          typeFrom[0]=='[' and 0 or
                          0 ) +\
                        ( typeFrom[1]=='0' and 0 or
                          typeFrom[1]=='1' and -1 or
                          0 ) 
    intvCorrectionFrom[1] = ( typeFrom[3]==')' and 0 or
                          typeFrom[3]==']' and 1 or
                          0 ) +\
                        ( typeFrom[1]=='0' and 0 or
                          typeFrom[1]=='1' and -1 or
                          0 ) 

    intvCorrectionTo=[0,0]
    intvCorrectionTo[0] = ( typeTo[0]=='(' and 1 or
                          typeTo[0]=='[' and 0 or
                          0 ) +\
                        ( typeTo[1]=='0' and 0 or
                          typeTo[1]=='1' and -1 or
                          0 ) 
    intvCorrectionTo[1] = ( typeTo[3]==')' and 0 or
                          typeTo[3]==']' and 1 or
                          0 ) +\
                        ( typeTo[1]=='0' and 0 or
                          typeTo[1]=='1' and -1 or
                          0 ) 

    return ( intvCorrectionFrom[0]-intvCorrectionTo[0], 
             intvCorrectionFrom[1]-intvCorrectionTo[1] )


def chrStartStopArg_11incl_to_00incl(s):
    try:
        chrom,start,stop=chrStartStop(s)
        return (chrom,
        	    coordConversion('[11]','[00]')[0]+start,
        	    coordConversion('[11]','[00]')[0]+stop )
    except:
        raise argparse.ArgumentTypeError('range must be chrom:start-stop')



def mkodict( lcols ):
    return OrderedDict( [(c,[]) for c in lcols] )



# col1:type,col2:type,...
# type - None means don't check; else, use a numpy dtype string.

class DframeFormatChecker(object):
    
    """check that the names and types of a dataframe columns match the expected format"""
    def __init__(self, formatstr=None, addl_cols_ok=True):

        self.mcolname_type = OrderedDict( [  ( col.split(':')[0],None) if ((':' not in col) or (col.split(':')[1]=='None')) else
                                             ( col.split(':')[0],np.dtype(col.split(':')[1])) 
                                             for col in formatstr.split(',')] )
        self.addl_cols_ok = addl_cols_ok

    def check( self, dframe, convert=True ):

        df_colnames = set([ str(cn) for cn in dframe.columns] )
        check_colnames = set( [ cdef for cdef in self.mcolname_type ] )

        extra_colnames = df_colnames.difference( check_colnames )
        if not self.addl_cols_ok and len(extra_colnames) > 0:
            raise ValueError('disallowed extra columns: {}'.format( ','.join(list(extra_colnames)) ))

        missing_colnames = check_colnames.difference( df_colnames )

        if len(missing_colnames) > 0:
            raise ValueError('missing required columns: {}'.format( ','.join(list(missing_colnames)) ))

        bad_dts = False
        wrong_dt = []

        if dframe.shape[0]>0:
            for cn in check_colnames.intersection(df_colnames):
                if self.mcolname_type[ cn ] is not None:
                    if dframe[cn].dtype != self.mcolname_type[ cn ]:
                        if convert:
                            dframe[cn] = dframe[cn].astype(self.mcolname_type[ cn ])
                        else:
                            wrong_dt.append( '{}: found {}, expected {}'.format(cn, dframe[cn].dtype, self.mcolname_type[ cn ] ))
                            bad_dts = True

            if bad_dts:
                raise ValueError('bad column types: {}'.format( '; '.join(list(wrong_dt)) ))

        return

    def new_odict_of_lists( self ):
        return mkodict( self.mcolname_type.keys() )



"""
cloned_fragments_tbl:  table defining genomic source of cloned fragments - to translate from genomic to construct coordinates.  All coordinates are 1 based inclusive. Construct fragment is assumed on positive strand.
    - genome_chrom          genomic chromosome name of the fragment
    - genome_start         genome start coordinate of the fragment
    - genome_end           genome end coordinate of the fragment
    - genome_strand        genome strand of the fragment
    - vector_chrom      construct chromosome name of the fragment
    - vector_start       construct start coordinate of the fragment
    - vector_end         construct end coordinate of the fragment
"""
cloned_fragments_tbl_checker = DframeFormatChecker( 'genome_chrom:str,genome_start:int,genome_end:int,genome_strand:None,vector_chrom:None,vector_start:int,vector_end:int', True )

"""
gene model table:  table with annotations of cloned exons within genome coodinates
    - gene_model_id        gene model id, a unique identifier for the gene model
    - exon_name            exon name, a unique identifier for the exon
    - genome_chrom          genomic chromosome name of the fragment
    - genome_start         genome start coordinate of the fragment
    - genome_end           genome end coordinate of the fragment
    - genome_strand        genome strand of the fragment
    - cdna_start           cdna start coordinate of the exon, hgvs formatted, e.g., 1234
"""
gene_model_tbl_checker = DframeFormatChecker( 'gene_model_id:str,exon_name:str,genome_chrom:str,genome_start:int,genome_end:int,genome_strand:None,cdna_start:int', True )

