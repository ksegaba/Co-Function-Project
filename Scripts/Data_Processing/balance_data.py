"""
Create a balanced feature table to build models with
"""
__author__ = "Kenia Segura Abá"

import os
import datatable as dt
import pandas as pd
import matplotlib.pyplot as plt

os.chdir("/mnt/home/seguraab/Shiu_Lab/Co-function")

if __name__ == "__main__":
    # Read in pre-processed & imbalanced co-function and fitness feature table
    df = dt.fread("Data/Costanzo_2016/S1/SGA_pwys_gene_fam_fitness.csv")
    df = df.to_pandas()
    df = df[["Query Gene", "Array Gene", "Query single mutant fitness (SMF)", 
            "Array SMF", "Double mutant fitness", "new_Co-Function"]]
    
    # Check class imbalance
    ax = df.groupby(["new_Co-Function"]).count().plot.bar()
    plt.savefig("Figures/class_imbalanced.pdf")
    plt.close()

    # Check fitness distributions in each class
    ax = df.hist(by=df["new_Co-Function"], bins=20, alpha=0.3, width=0.04, figsize=(11,5))
    plt.savefig("Figures/class_imbalanced_fitness.pdf")
    plt.close()