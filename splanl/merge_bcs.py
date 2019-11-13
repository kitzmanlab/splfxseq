
# coding: utf-8

# In[18]:


import os
import pandas as pd
import numpy as np


# In[135]:


def MergeBCs(varFile,filesToMerge,readCountFilter,numMuts=None):
    """Merges barcodes into a dataset with each row representing one variant. Barcodes with number of reads less
    than the read count filter are not included in final counts. You can specify the number of mutations in each
    variant - 0 invokes a wild type dataset.
    
    NOTE: Function has not been tested (and I believe will fail) if numMuts>1.
    
    Args:
        varFile (str): '/path/and/name/to/file.txt'
        filesToMerge (list of str): ['/path/and/name/to/file1.txt','/path/and/name/to/file2.txt'] 
        readCountFilter (int): Number of reads required to be included in counts
        numMuts (int or None): Specifies the number of mutations per barcode for final file.

    Returns:
        Nothing.
        
    Saves:
        text file with each row representing a variant in the same directory as the filesToMerge.
        
    Usage:
        MergeBCs('/path/and/name/to/varfile.txt',['/path/and/name/to/file1.txt','/path/and/name/to/file2.txt'],
        20,0)
    """
    mergeddf=MergeReps(filesToMerge)
    vardf,varBCDict=ProcessVarFile(varFile,numMuts)
    dfList=[]
    for var in varBCDict:
        if numMuts!=0:
            pos=var.split(':')[1]
            ref=var.split(':')[2]
            alt=var.split(':')[3]
        countDict=ComputeBCStats(varBCDict[var],mergeddf,readCountFilter)
        if numMuts!=0:
            rowList=ColumnAppend([],[var,numMuts,pos,ref,alt,[[varBCDict[var]]],2*len(varBCDict[var]),
                                 [val for val in countDict.values()]])
        else:
            rowList=ColumnAppend([],[var,numMuts,[[varBCDict[var]]],2*len(varBCDict[var]),
                                 [val for val in countDict.values()]])
        dfList.append(rowList)
    if numMuts!=0:
        varpsidf=pd.DataFrame(dfList, columns=['variant','num_variants','position','ref_allele','alt_allele','BCs',
                                          'total_BCs']+[key for key in countDict.keys()])
    else:
        varpsidf=pd.DataFrame(dfList, columns=['variant','num_variants','BCs','total_BCs']+
                              [key for key in countDict.keys()])
    varDirectory,file=os.path.split(varFile)
    BCDirectory,junk=os.path.split(filesToMerge[0])
    if numMuts!=0:
        outfile='/PSIdata_byVar_'+file.split('.')[0]+'_filter_'+str(readCountFilter)+'.txt'
    else:
        outfile='/PSIdata_byVar_WT_'+file.split('.')[0]+'_filter_'+str(readCountFilter)+'.txt'
    print(BCDirectory+outfile)
    if not os.path.exists(BCDirectory+outfile):
        varpsidf.to_csv(BCDirectory+outfile, index=None, mode='a', sep='\t')
    else:
        varpsidf.to_csv(BCDirectory+outfile, index=None, mode='w', sep='\t') 


# In[156]:


def ComputeBCStats(BCList,BCdf,readCountFilter):
    """Removes barcodes with less reads than the read filter. Computes number of reads matching each isoform for a
    single variant with all barcodes representing that variant in the BCList. Also computes other useful information
    like total number of reads, total number of barcodes, barcodes passing filter, reads unmapped, reads with bad starts.

    Args:
        BCList (list of str): ['BC1','BC2','BC3']
        BCdf (pandas dataframe): Dataframe by barcode with all reps merged and columns from each rep ending in '_repN'
        readCountFilter (int): Number of reads required to be included in counts

    Returns:
        countDict (dictionary - str:int or list of ints): dictionary with keys representing column names and values
        as the counts across barcodes that pass the read filter.
        
    Usage:
        countDict=ComputeBCStats(['BC1','BC2','BC3'],BCdf,20)
    """
    colList=list(BCdf.columns)
    countDict={}
    countDict.update({'BCused':0, 'BCpresent':0})
    countDict.update({col.rsplit('_',1)[0]:[] for col in colList if col[-1]=='1' and 'psi' in col})
    countDict.update({col.rsplit('_',1)[0]:0 for col in colList if (col[-1]=='1' and 'psi' not in col)})
    repList=list(set([col.rsplit('_',1)[1] for col in colList if len(col.rsplit('_',1))>1]))
    for rep in repList:
        for BC in BCList:
            BCRow=BCdf.loc[(BCdf['barcode']==BC) & (BCdf['usable_reads_'+rep]>0)]
            if BCRow.empty:
                continue
            else:
                countDict['BCpresent']+=1
                countDict['num_reads']+=BCRow['num_reads_'+rep].values[0]
                countDict['unmapped_reads']+=BCRow['unmapped_reads_'+rep].values[0]
                countDict['bad_starts']+=BCRow['bad_starts_'+rep].values[0]
                BCRow=BCRow.loc[BCRow['usable_reads_'+rep]>readCountFilter]
                if BCRow.empty:
                    continue
                else:
                    countDict['BCused']+=1
                    countDict['usable_reads']+=BCRow['usable_reads_'+rep].values[0]
                    for key in countDict.keys():
                        if 'isoform' in key and 'psi' not in key:
                            countDict[key]+=BCRow[key+'_'+rep].values[0]
                        elif 'psi' in key:
                            countDict[key].append(BCRow[key+'_'+rep].values[0])
    for key in list(countDict.keys()):
        if 'psi' in key:
            if len(countDict[key])>1:
                countDict[key+'_median']=np.median(np.array(countDict[key]))
            elif len(countDict[key])==1:
                countDict[key+'_median']=countDict[key][0]
            else:
                countDict[key+'_median']=np.nan
            del countDict[key]
        elif 'isoform' in key:
            if countDict['usable_reads']>0:
                countDict[key+'_psi_mean']=countDict[key]/countDict['usable_reads']
            else:
                countDict[key+'_psi_mean']=np.nan
    return countDict


# In[121]:


def ProcessVarFile(varFile,numMuts):
    """Takes a file and opens it into a pandas dataframe. File must contain columns: 'n_variants_passing', 
    'readgroupid', and 'variant_list' or it will throw an error. If the number of mutants is nonzero, then the 
    pandas dataframe is filtered to only contain variants with that number of mutations. If number of mutants is 
    left unspecified, the pandas dataframe is filtered to contain only variants with a nonzero number of mutations.
    If number of mutants is 0, then the pandas dataframe is filtered to contain only variants listed as 
    'no_variants_input'. Then a dictionary of [var]=[BC1,BC2,BC3] is created to map all the barcodes to their 
    variant. If number of mutants is 0, then the dictionary is of length 1 and the var='WT'

    Args:
        varFile (str): '/path/and/name/to/file.txt'
        numMuts (int or None): Specifies the number of mutations per barcode for final file.

    Returns:
        vardf (pandas dataframe): dataframe filtered by number of mutations allowed per barcode
        varBCDict (dictionary - str:list of str): variant from 'variant list' column and a list of barcodes
        corresponding to that variant
        
    Usage:
        vardf,varBCDict=ProcessVarFile('/path/and/name/to/file.txt',0)
    """
    vardf=pd.read_csv(varFile, sep='\t')
    if 'n_variants_passing' not in list(vardf.columns):
        raise ValueError('n_variants_passing is not a column name in '+varFile)
    elif 'readgroupid' not in list(vardf.columns):
        raise ValueError('readgroupid is not a column name in '+varFile)
    elif 'variant_list' not in list(vardf.columns):
        raise ValueError('variant_list is not a column name in '+varFile)
    varBCDict={}
    if numMuts and numMuts>0:
        vardf=vardf[vardf['n_variants_passing']==numMuts]
    elif numMuts is None:
        vardf=vardf[vardf['n_variants_passing']>0]
    elif numMuts==0:
        vardf=vardf[vardf['status']=='no_variants_input']
    for row in vardf.itertuples():
        if numMuts==0:
            var='WT'
        else:
            var=row.variant_list
        if var not in varBCDict:
            varBCDict[var]=[row.readgroupid]
        else:
            varBCDict[var].append(row.readgroupid)
    return(vardf,varBCDict)


# In[26]:


def MergeReps(filesToMerge):
    """Merges files by barcode using the 'rep' name as the suffix. Thus, all columns from the rep1 dataset will 
    end in '_rep1' within the merged dataset. Returns a value error if 'barcode' is not a column in all the files 
    to merge. Replaces any NaN's from dividing by zero as 0's.
    
    NOTE: This theoretically works to merge n different datasets but has only been tested using two datasets.

    Args:
        filesToMerge (list of str): ['/path/and/name/to/file1.txt','/path/and/name/to/file2.txt'] 

    Returns:
        mergeddf (pandas dataframe): dataframe of each original dataframe merged by the barcode
        
    Usage:
        mergeddf=MergeReps(['/path/and/name/to/file1.txt','/path/and/name/to/file2.txt'])
    """
    if len(filesToMerge)<2:
        raise ValueError('Need to include at least two files to merge together')
    dfDict={}
    for idx,file in enumerate(filesToMerge):
        #extracts repN from file name in the format '/path/to/file_repN.something.bam'
        rep=file.rsplit('_',2)[1]
        dfDict[rep+'df']=pd.read_csv(file, sep='\t')
        if 'barcode' not in list(dfDict[rep+'df'].columns):
            raise ValueError('Barcode is not a column name in '+file)
    #merges the first two dataframes with '_repN' as the suffix
    mergeddf=pd.merge(dfDict[list(dfDict.keys())[0]],dfDict[list(dfDict.keys())[1]], on='barcode', how='outer',
                      suffixes=('_'+list(dfDict.keys())[0][:-2],'_'+list(dfDict.keys())[1][:-2]))
    if len(filesToMerge)>2:
        for idx in range(2,len(filesToMerge)):
            mergeddf=pd.merge(mergeddf,dfDict[list(dfDict.keys())[idx]], on='barcode', how='outer',
                              suffixes=('','_'+list(dfDict.keys())[idx][:-2])) 
    #turns any missing values (result of dividing by 0 when 0 mapped reads) into zero values
    mergeddf.fillna(0,inplace=True)
    return(mergeddf)


# In[16]:


def ColumnAppend(rowList,colList):
    """Appends items to an existing list. Not a difficult task but cleans up main section of code.

    Args:
        rowList (list): List to append to
        colList (list): Items to append - if the items are a list, they are incorporated into existing list. 

    Returns:
        rowList (list): original list with additional items added
        
    Usage:
        rowList=ColumnAppend(rowList,colList)
        ColumnAppend([2,4],[1,4,6])
        [2,4,1,4,6]
        ColumnAppend([2,4],[[1,1,0],4,2])
        [2,4,1,1,0,4,2]
    """
    for col in colList:
        if isinstance(col,list):
            rowList+=col
        else:
            rowList.append(col)
    return(rowList)


# In[157]:


get_ipython().run_cell_magic('time', '', "MergeBCs('/nfs/kitzman2/smithcat/proj/campersplice/subassembly/JKP555.haps.final.txt',\n         ['/nfs/kitzman2/smithcat/proj/campersplice/rna-seq/exon10/PSIdata_byBC_JKP555_cDNA_fixed_rep1_exon10.txt',\n      '/nfs/kitzman2/smithcat/proj/campersplice/rna-seq/exon10/PSIdata_byBC_JKP555_cDNA_fixed_rep2_exon10.txt'],\n         20,1)")


# In[163]:


get_ipython().run_cell_magic('time', '', "MergeBCs('/nfs/kitzman2/smithcat/proj/campersplice/subassembly/JKP556.haps.final.txt',\n         ['/nfs/kitzman2/smithcat/proj/campersplice/rna-seq/exon11/PSIdata_byBC_JKP556_cDNA_fixed_rep1_exon11.txt',\n      '/nfs/kitzman2/smithcat/proj/campersplice/rna-seq/exon11/PSIdata_byBC_JKP556_cDNA_fixed_rep2_exon11.txt'],\n         20,1)")


# In[159]:


get_ipython().run_cell_magic('time', '', "MergeBCs('/nfs/kitzman2/smithcat/proj/campersplice/subassembly/JKP555.haps.final.txt',\n         ['/nfs/kitzman2/smithcat/proj/campersplice/rna-seq/exon10/PSIdata_byBC_JKP555_cDNA_fixed_rep1_exon10.txt',\n      '/nfs/kitzman2/smithcat/proj/campersplice/rna-seq/exon10/PSIdata_byBC_JKP555_cDNA_fixed_rep2_exon10.txt'],\n         20,0)")


# In[ ]:


get_ipython().run_cell_magic('time', '', "MergeBCs('/nfs/kitzman2/smithcat/proj/campersplice/subassembly/JKP556.haps.final.txt',\n         ['/nfs/kitzman2/smithcat/proj/campersplice/rna-seq/exon11/PSIdata_byBC_JKP556_cDNA_fixed_rep1_exon11.txt',\n      '/nfs/kitzman2/smithcat/proj/campersplice/rna-seq/exon11/PSIdata_byBC_JKP556_cDNA_fixed_rep2_exon11.txt'],\n         20,0)")


# In[161]:


get_ipython().run_cell_magic('time', '', "MergeBCs('/nfs/kitzman2/smithcat/proj/campersplice/subassembly/JKP555.haps.final.txt',\n         ['/nfs/kitzman2/smithcat/proj/campersplice/rna-seq/exon10/PSIdata_byBC_JKP555_cDNA_fixed_rep1_exon10.txt',\n      '/nfs/kitzman2/smithcat/proj/campersplice/rna-seq/exon10/PSIdata_byBC_JKP555_cDNA_fixed_rep2_exon10.txt'],\n         0)")


# In[ ]:


get_ipython().run_cell_magic('time', '', "MergeBCs('/nfs/kitzman2/smithcat/proj/campersplice/subassembly/JKP556.haps.final.txt',\n         ['/nfs/kitzman2/smithcat/proj/campersplice/rna-seq/exon11/PSIdata_byBC_JKP556_cDNA_fixed_rep1_exon11.txt',\n      '/nfs/kitzman2/smithcat/proj/campersplice/rna-seq/exon11/PSIdata_byBC_JKP556_cDNA_fixed_rep2_exon11.txt'],\n         0)")

