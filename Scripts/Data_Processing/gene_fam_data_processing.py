#!/usr/bin/env python3
"""
Description: Prepare feature matrix from Costanzo 2016 GI data
- Remove duplicate gene pairs
- Add pathway information to create label (yes/no co-functional)
- Filter gene pairs by gene family to get final label

"""
__author__ = "Kenia Segura Abá"

import os
import datatable as dt
import pandas as pd
import numpy as np
from collections import Counter
os.chdir("/mnt/home/seguraab/Shiu_Lab/Co-function/")

if __name__=="__main__":
    # Read in feature table with co-functional label based on shared pathway annotations
    df = dt.fread("Data/Costanzo_2016/S1/SGA_pwys_fitness.csv")
    df = df.to_pandas()

    ## Ensure non-cofunctional gene pairs are not in the same gene family and assign new label
    fam = pd.read_csv("Data/PANTHER/pantherGeneListProtein_forencoding.txt", sep="\t")
    fam.shape # (5047, 13)
    len(fam.systematic_name.unique()) # 5011 unique genes
    len(fam.PANTHER_fam_subfam_num.unique()) # 3919 unique family IDs
    len(fam.PANTHER_fam_subfam.unique()) # 3895 unique family descriptions
    df["shared_fam"] = 0 # whether the genes in a pair are in the same gene family
    df["new_Co-Function"] = 0 # new co-function label based on pathway and gene family
    to_check = [] # indices to check when genes map to multiple gene families
    # [88197, 88202, 88965, 88990, 89024, 89119, 89128, 89141, 89153, 89166, 89238, 89290, 89311, 101492, 235411, 235425, 235439]
    for i in range(df.shape[0]):
        current_label = df["Co-Function"].values[i]
        if current_label==True: # update label for co-functional gene pairs according to pathways
            gene1 = df["Query Gene"].values[i]
            gene2 = df["Array Gene"].values[i]
            fam1 = fam[fam["systematic_name"]==gene1]["PANTHER_fam_subfam_num"].to_list()
            fam2 = fam[fam["systematic_name"]==gene2]["PANTHER_fam_subfam_num"].to_list()
            if (len(fam1)==1 and len(fam2)==1 and fam1==fam2): # in the same family
                df["new_Co-Function"].values[i] = 0 # change label to not co-functional
            elif (len(fam1)==1 and len(fam2)==1 and fam1!=fam2):
                df["new_Co-Function"].values[i] = 1 # keep the original label
            if (len(fam1)>1 or len(fam2) >1): # when either gene maps to multiple families
                to_check.append(i)
    df["new_Co-Function"].sum() # 4132 co-functional gene-pairs based on both pathways and gene families
    df["Co-Function"].sum() # 5026 co-functional gene-pairs based only on pathways
    df.to_csv("Data/Costanzo_2016/S1/SGA_pwys_gene_fam_fitness.csv", chunksize=1000) # save final unblanced dataset
    # there are 258400 gene pairs total and they should be unique
    # 4132 of these are co-functional. This dataset is extremely unbalanced.
    # All the gene pairs should have pathway and gene family annotations, but 
    # I didn't explicitly remove those with empty gene family annotation, 
    # probably need to add this additional condition and re-run the above code

    # ??? should I drop gene pairs that are in the same gene family, either way they're still co-functional

    # Create a balanced dataset
    # What approaches should I do?
    # add positive class