#!/usr/bin/env python3
"""
Description: Prepare feature matrix from Costanzo 2016 GI data
- Remove duplicate gene pairs
- Add pathway information to create label (yes/no co-functional)

"""
__author__ = "Kenia Segura Abá"


import datatable as dt
import pandas as pd
import numpy as np
from collections import Counter

def add_pwys(
    df1, df2, new_col, how="left", left_on="Query Gene", right_on="Accession-1", ):
    df1 = df1.merge(df2, how=how, left_on=left_on, right_on=right_on)
    df1.rename(columns={"Pathways of gene":new_col}, inplace=True)
    df1.drop("Accession-1", axis=1, inplace=True)
    return df1

if __name__=="__main__":
    ## Combine datasets
    # read in data    
    costanzo = "/mnt/home/seguraab/Shiu_Lab/Co-function/Data/Costanzo_2016/S1"
    EE = dt.fread(f"{costanzo}/SGA_ExE.txt") # essential x essential gene pairs
    EN = dt.fread(f"{costanzo}/SGA_ExN_NxE.txt") # essential x nonessential gene pairs
    NN = dt.fread(f"{costanzo}/SGA_NxN.txt") # nonessential x nonessential gene pairs

    # add interaction type column
    EE["Type"] = "ExE"
    EN["Type"] = "ExN"
    NN["Type"] = "NxN"

    # combine and save
    data = dt.rbind(EE, EN, NN)
    data.to_csv(f"{costanzo}/SGA_all.csv")
    data = data.to_pandas()
    # data.shape # 19313654 rows

    ## Add pathway information
    metacyc = "/mnt/home/seguraab/Shiu_Lab/Co-function/Data/MetaCyc"
    pwys = pd.read_csv(f"{metacyc}/All-genes-pathways-S288c.txt", sep="\t")
    pwys = pwys[["Accession-1", "Pathways of gene"]]
    # pwys.shape
    # len(pwys['Accession-1'].unique()) # no. of unique genes
    # pwys.dropna().shape # genes with pathway info
    
    # for query genes
    query = data["Query Strain ID"].str.split("_", n=2, expand=True)
    data.insert(0, "Query Gene", query[0])
    data = add_pwys(data, pwys, new_col="Query Pathway")

    # for array genes
    array = data["Array Strain ID"].str.split("_", n=2, expand=True)
    data.insert(1, "Array Gene", array[0])
    data = add_pwys(data, pwys, left_on="Array Gene", new_col="Array Pathway")

    ## Aggregate fitness data for duplicate pairs (A-B, B-A)
    data.drop(["Accession-1_x", "Accession-1_y"], axis=1, inplace=True)
    data.shape # (19313654, 16)
    data.to_csv(f"{costanzo}/SGA_all_pwys.csv", index=False) # save unfiltered dataset
    # data.groupby(['Query Gene', 'Array Gene']).size() # counts of duplicate gene pairs
    group = data[["Query Gene", "Array Gene"]].apply(frozenset, axis=1) # gene pairs
    data_med = data.iloc[:,[0,1,9,10,11]].groupby(group).aggregate("median") # median fitness values
    data_med.reset_index(inplace=True) # reset index
    data_med["index"] = data_med["index"].apply(lambda x: tuple(x)) # coerce frozenset to tuple
    data_med.index = pd.MultiIndex.from_tuples(data_med["index"], names=["Query Gene", "Array Gene"]) # create multi-index
    data_med.drop("index", axis=1, inplace=True) # drop "index" column containing tuples
    data_med.reset_index(inplace=True) # remove multiindex
    data_med.shape # (12936843, 5) no. of gene pairs with SMF and DMF data

    ## Delete remaining duplicate rows by gene pairs
    # Check how many duplicate pairs exist
    data_med.duplicated().sum() # 0
    
    ## Add pathway info again
    data_med = add_pwys(data_med, pwys, new_col="Query Pathway")
    data_med = add_pwys(data_med, pwys, left_on="Array Gene", new_col="Array Pathway")
    data_med = data_med.dropna() # drop gene pairs with missing pathway info
    data_med.shape # (258400, 7)

    ## Assign label to gene pairs
    # 1 (co-functional): query & array pathways are the same
    # 0 (not co-functional): query & array pathways differ
    # split Query Pathway and check if at least one pathway is in Array Pathway
    data_med['Index'] = np.linspace(0, data_med.shape[0], data_med.shape[0], endpoint=False, dtype="int32")
    data_med["Co-Function"] = data_med['Index'].\
        apply(lambda x: 1 if str(bool(\
            [pwy for pwy in data_med.iloc[x,5].split(" // ") \
            if (pwy in data_med.iloc[x,6].split(" // "))]
        ))=='True' else 0)
    data_med["Co-Function"].values.sum() # 5026 co-functional gene pairs, 253374 non-cofunctional gene pairs
    data_med.set_index("Index", drop=True, inplace=True)

    ## Ensure non-cofunctional gene pairs are not in the same gene family
    neg_set = data_med.loc[data_med["Co-Function"]==0,:]

    # gene family data

    # append gene family data
    # to query 
    # to array

    # check the gene pairs are not in the same gene family
    # drop gene pairs that are in the same gene family

    ## save final unblanced dataset
    # add positive class