#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 16 15:13:17 2021
@author: erik
"""
import csv
import os
import sys
import pandas as pd
import argparse
from parse_rsem_output import *

def process_rsem_counts_to_GS_matrix(your_label,\
                                     your_rsem_dir, your_meta_file1,\
                                     your_meta_file2, your_out_dir):
    '''
    To wrap-up everything
    your_label: experiment label for output matrix
    your_rsem_dir: dir with rsem counts
    meta_file1: data/processed/metaGSM.txt
    meta_file2: output of first rule
    your_out_dir: dir with tables
    '''

    #handle metadata and target genes
    meta1 = pd.\
    read_csv(\
    os.path.join(os.getcwd(),\
    your_meta_file1),\
        header=1, delimiter="\t")
    meta1.columns = ['GSM','gene_label']
    meta1.GSM = meta1.GSM.str.strip()
    
    meta2 = pd.\
    read_csv(\
    os.path.join(os.getcwd(),\
    your_meta_file2),\
        header=0, delimiter=" ")
    meta2.columns = ['GSM','SRR']
    meta2.GSM = meta2.GSM.str.strip()
    meta_info = pd.merge(meta1, meta2,
                   on='GSM',
                   how='left')
    meta_ctrl =\
    meta_info[meta_info['gene_label']\
    .str.contains('WT', na=False)]
    meta_pert =\
    meta_info[~meta_info.\
    GSM.isin(meta_ctrl.GSM)]
    target_genes = meta_pert['gene_label'].tolist()

    #handle control and perturbation experiments 
    all_ctrl_df = merge_rsem_TPMs(meta_ctrl,\
        os.path.join(os.getcwd(), your_rsem_dir), target_genes)
    
    all_pert_df = merge_rsem_TPMs(meta_pert,\
        os.path.join(os.getcwd(), your_rsem_dir), target_genes)
    all_pert_df.index = all_pert_df['gene_label']
    all_pert_df=\
    all_pert_df.drop(columns=['gene_label'])

    all_pert_df.columns=\
    meta_pert.gene_label.tolist()
    
    #ctrls
    ctrl1 = meta_ctrl['SRR'].tolist()[0]
    ctrl2 = meta_ctrl['SRR'].tolist()[1]
    
    #save expression matrices for each replicate
    Y_matrix1 = calculate_FC_from_TPMs(all_pert_df, all_ctrl_df[[ctrl1]])
    Y_matrix2 = calculate_FC_from_TPMs(all_pert_df, all_ctrl_df[[ctrl2]])
    Y_list = [Y_matrix1, Y_matrix2]

    #export matrices as csv-files in GeneSpider format
    Y_file, P_file = save_FC_to_GS(Y_list)
    Y_file.to_csv\
            (str(os.path.join(os.getcwd(),str(your_out_dir),\
                          str(your_label)+"_y.csv")),\
             index=True, header=True, sep = '\t')
    P_file.to_csv\
            (str(os.path.join(os.getcwd(),str(your_out_dir),\
                          str(your_label)+"_p.csv")),\
             index=True, header=True, sep = '\t')
    return True



def main():
    '''Function to convert transcript counts of multiple experiments (RSEM)
    to GeneSpider expression matrix format
    Requires:
    Returns: Y,P-matrices
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument("-l", "--label", required=True, help="experiment label")
    parser.add_argument("-r", "--rsem_dir", required=True, help="rsem counts directory")
    parser.add_argument("-m1", "--meta_GSM_file", required=True, help="meta file 1")
    parser.add_argument("-m2", "--meta_SRR_file", required=True, help="meta file 2")
    parser.add_argument("-o", "--out_dir", required=True, help="output directory")
    args = parser.parse_args()

    process_rsem_counts_to_GS_matrix(args.label, args.rsem_dir, args.meta_GSM_file, args.meta_SRR_file, args.out_dir)

if __name__ == "__main__":
    main()
