#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 15 13:53:59 2021

@author: erikzhi
"""

import os
import pandas as pd
from functools import reduce
import numpy as np

def get_TPM_from_rsem_counts(counts_file, SRR_name, target_genes):
    """

    Returns:
      -counts - pd dataframe with TPM (transcripts per million) counts
    -------
    Requires:
      -counts_file - the output RSEM file with counts per gene (.genes.results)
      -SRR_name - name of the experiment
      -target_genes - list of genes of interest 

    """
    counts = pd.\
    read_csv(counts_file,\
    delimiter="\t")
    counts['gene_label'] =\
        counts['gene_id'].\
        str.split('=',3,True)[1].\
        str.split('_',2,True)[0]
    counts =\
        counts.rename(columns = {'TPM':SRR_name})
    counts = counts[counts['gene_label'].\
                    isin(target_genes)]
    counts = pd.DataFrame(counts, columns =\
    ['gene_label', SRR_name])
    
    return counts

def merge_rsem_TPMs(meta_df, matrix_dir, target_genes):
    """

    Returns:
      -reduced_TPMs - concatenated pd daraframe with TPMs
    -------
    Requires:
      -meta_df - pd datarame with metadata coumns: GSM, SRR, gene_label
          where:
              -GSM - GSM identities of GEO experiments
              -SRR - SRR file name of the corresponding experiments
              -gene_labels - experiment labels, each label corresponds
                             to experiment where the gene was perturbed
      -matrix_dir - path to dir where RSEM output files are located

    """
    TPMs = []

    for SRR in meta_df['SRR']:
        counts_file=os.path.join(matrix_dir,\
                                 SRR+".genes.results")
        print(counts_file)
        TPMs.append\
        (get_TPM_from_rsem_counts(counts_file, SRR, target_genes)) 
    
    reduced_TPMs =\
        reduce(lambda x, y:\
        pd.merge(x, y, on = 'gene_label'),\
            TPMs)

    return reduced_TPMs


def calculate_FC_from_TPMs(TPM_matrix, ctrl_vector):
    """
    
    Returns:
      -FC_matrix - log2FC matrix
    -------
    Requires:
      -TPM_matrix - pd with TPM values per gene per experiment
      -ctrl_vector - pd column with TPM values per gene

    """

    eps = 1e-7

    #rep1
    ctrl_vector_vals = ctrl_vector.values.tolist()
    ctrl_vector_vals = [item for sublist in ctrl_vector_vals for item in sublist]
    FC_matrix = TPM_matrix.div(ctrl_vector_vals, axis=0)
    
    FC_matrix = (FC_matrix + eps).apply(np.log2)
    FC_matrix.sort_index(inplace=True)
    
    FC_matrix =FC_matrix.drop(columns=\
    [col for col in FC_matrix if col not in FC_matrix.index.tolist()])
        
    return FC_matrix

def save_FC_to_GS(list_with_matrices):
    """
    
    Returns:
      -Y (expression), P (perturbation) matrices in GeneSpider format
     -------
    Requires:
      -list_with_matrices - list of log2FC matrices where
                            each matrix represents one replicate 
    """
    
    Y_list = pd.concat(list_with_matrices,axis=1)
    
    P_list = pd.DataFrame(np.tile(np.eye(len(list_with_matrices[0].columns)),\
                                  (1,len(Y_list))))
    
    return Y_list, P_list
