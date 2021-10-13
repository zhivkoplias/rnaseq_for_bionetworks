# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.

import matlab.engine
eng = matlab.engine.start_matlab()
tf = eng.isprime(37)
print(tf)

"""
#import functions
import os
import numpy as np
import pandas as pd
from PandasGCTXParserL1000 import PandasGCTXParserL1000
gparser = PandasGCTXParserL1000()

#params
output_dir = '/home/erikz/sonnhammer/work-in-progress/GCTX_counts_L1000/matrices/'
data_dir = '/home/erikz/sonnhammer/work-in-progress/GCTX_counts_L1000/data/'
list_of_cell_lines = ['A375', 'A549', 'HA1E', 'HCC515', 'HEPG2',\
                      'HT29', 'MCF7', 'PC3']

#parse level5
for cell_line in list_of_cell_lines:
    #read data
    exp_data_lvl5, ctrl_data_lvl5 = gparser.read_gctx_data(cell_line,\
    L1000_gctx_file =os.path.join(data_dir,\
                    'GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx'),\
    inst_info_file = os.path.join(data_dir, 'GSE92742_Broad_LINCS_sig_info.txt'),
    gene_info_file = os.path.join(data_dir, 'GSE92742_Broad_LINCS_gene_info.txt'), level=5)
    
    #select columns with 3 reps
    exp_data_lvl5_matrix = gparser.gctoo2matrices_lvl5(exp_data_lvl5)
    
    #save as y-matrix
    exp_data_lvl5_matrix.to_csv\
            (output_dir+cell_line+'_lvl5_y.csv',\
             index=True, header=True, sep = '\t')

#parse level3 data
list_of_cell_lines = ['HA1E', 'HCC515', 'HEPG2',\
                      'HT29', 'A549', 'MCF7','PC3', 'A375']

shRNA_num = 2
for cell_line in list_of_cell_lines:
    print(cell_line)
    #read data
    exp_data_lvl3, ctrl_data_lvl3 = gparser.read_gctx_data(cell_line,\
    L1000_gctx_file =os.path.join(data_dir,\
                    'GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx'),\
    inst_info_file = os.path.join(data_dir, 'GSE92742_Broad_LINCS_inst_info.txt'),
    gene_info_file = os.path.join(data_dir, 'GSE92742_Broad_LINCS_gene_info.txt'), level=3)
    
    #subset
    exp_data_lvl3_subset = gparser.filter_gctx_data(exp_data_lvl3, ["X1","X2","X3"], 3)
    
    #select 'best shRNA experiment' (highest correlation between replicates)
    exp_data_lvl3_subset_bestshRNA_rep1, exp_data_lvl3_subset_bestshRNA_rep2,\
        exp_data_lvl3_subset_bestshRNA_rep3 = \
        gparser.select_one_perturbator(exp_data_lvl3_subset, ["X1","X2","X3"])
    
    #compute average across different shRNAs // comment out if not needed
    exp_data_lvl3_subset_bulk_rep1, exp_data_lvl3_subset_bulk_rep2,\
        exp_data_lvl3_subset_bulk_rep3 = \
        gparser.merge_all_perturbators(exp_data_lvl3_subset, ["X1","X2","X3"], shRNA_num)
    
    #exp_data_lvl3_subset_bulk_rep1 = exp_data_lvl3_subset_bestshRNA_rep1
    #exp_data_lvl3_subset_bulk_rep2 = exp_data_lvl3_subset_bestshRNA_rep2
    #exp_data_lvl3_subset_bulk_rep3 = exp_data_lvl3_subset_bestshRNA_rep3
    
    
    #prepare replicates // 
    #rep_counts is the threshold for the lowest number of shRNA per gene per experiment
    exp_data_lvl3_subset_bulk_rep1_matrix =\
        gparser.gctoo2matrices_lvl5(exp_data_lvl3_subset_bulk_rep1)
    exp_data_lvl3_subset_bulk_rep2_matrix =\
        gparser.gctoo2matrices_lvl5(exp_data_lvl3_subset_bulk_rep2)
    exp_data_lvl3_subset_bulk_rep3_matrix =\
        gparser.gctoo2matrices_lvl5(exp_data_lvl3_subset_bulk_rep3)
    
    #genes that are the same betwen experiments
    common_genes = np.intersect1d(np.intersect1d(\
                   exp_data_lvl3_subset_bulk_rep1_matrix.index,\
                   exp_data_lvl3_subset_bulk_rep2_matrix.index),\
                   exp_data_lvl3_subset_bulk_rep3_matrix.index).tolist()
    
    #filter out experiments that are not present in all three reps
    exp_data_lvl3_subset_bulk_rep1_matrix =\
        gparser.filter_overlapping_experiments(exp_data_lvl3_subset_bulk_rep1_matrix,\
                                               common_genes)
    exp_data_lvl3_subset_bulk_rep2_matrix =\
        gparser.filter_overlapping_experiments(exp_data_lvl3_subset_bulk_rep2_matrix,\
                                               common_genes)
    exp_data_lvl3_subset_bulk_rep3_matrix =\
        gparser.filter_overlapping_experiments(exp_data_lvl3_subset_bulk_rep3_matrix,\
                                               common_genes)
    
    #prepare control
    ctrl_data_lvl3.data_df.index = ctrl_data_lvl3.\
                row_metadata_df['pr_gene_symbol']
    ctrl_data_lvl3.data_df.sort_index(inplace=True)
    ctrl_data_lvl3 = ctrl_data_lvl3.data_df
    ctrl_data_lvl3 = ctrl_data_lvl3[ctrl_data_lvl3.index.\
                            isin(common_genes)]
    ctrl_data_lvl3['mean'] = ctrl_data_lvl3.mean(axis=1)
    
    #calculate FC
    exp_data_lvl3_subset_bulk_rep1_matrix_FC =\
        gparser.calculate_FC(exp_data_lvl3_subset_bulk_rep1_matrix,\
                             ctrl_data_lvl3['mean'])
    exp_data_lvl3_subset_bulk_rep2_matrix_FC =\
        gparser.calculate_FC(exp_data_lvl3_subset_bulk_rep2_matrix,\
                             ctrl_data_lvl3['mean'])
    exp_data_lvl3_subset_bulk_rep3_matrix_FC =\
        gparser.calculate_FC(exp_data_lvl3_subset_bulk_rep3_matrix,\
                             ctrl_data_lvl3['mean'])
    
    #save FC matrices
    test_all = pd.concat([exp_data_lvl3_subset_bulk_rep1_matrix_FC,\
                          exp_data_lvl3_subset_bulk_rep2_matrix_FC,\
                          exp_data_lvl3_subset_bulk_rep3_matrix_FC], axis=1)
    
    #save experimental data
    #test_all = pd.concat([exp_data_lvl3_subset_bulk_rep1_matrix,\
    #                      exp_data_lvl3_subset_bulk_rep2_matrix,\
    #                      exp_data_lvl3_subset_bulk_rep3_matrix], axis=1)
    
    test_all.to_csv\
            (output_dir+cell_line+'_lvl3_exp.csv',\
             index=True, header=True, sep = '\t')


            
            
            
            
#parse level4 data
for cell_line in list_of_cell_lines:
    #read data
    HEPG2_data4, HEPG2_data_ctrl4 = gparser.read_gctx_data(cell_line,\
    L1000_gctx_file =\
    "/scratch/erikzhi/L1000_data/level4/GSE92742_Broad_LINCS_Level4_ZSPCINF_mlr12k_n1319138x12328.gctx",\
    inst_info_file = \
    "/scratch/erikzhi/L1000_data/LIMS/GSE92742_Broad_LINCS_inst_info.txt")
    #select columns with 3 reps
    HEPG_data_subset = gparser.filter_gctx_data(HEPG2_data4, ["X1","X2","X3"], 3)
    #subset plates
    #one shRNA
    #HEPG_data_subset_p1_one, HEPG_data_subset_p2_one, HEPG_data_subset_p3_one = \
    #    gparser.select_one_perturbator(HEPG_data_subset, ["X1","X2","X3"])
    #bulk shRNA
    HEPG_data_subset_p1_bulk, HEPG_data_subset_p2_bulk, HEPG_data_subset_p3_bulk = \
        gparser.merge_all_perturbators(HEPG_data_subset, ["X1","X2","X3"])
    
    #split
    test1 = gparser.gctoo2matrices_lvl5(HEPG_data_subset_p1_bulk, rep_counts=1)
    test2 = gparser.gctoo2matrices_lvl5(HEPG_data_subset_p2_bulk, rep_counts=1)
    test3 = gparser.gctoo2matrices_lvl5(HEPG_data_subset_p3_bulk, rep_counts=1)
    
    #prepare control
    HEPG2_data_ctrl4.data_df.index = HEPG2_data_ctrl4.\
                row_metadata_df['pr_gene_symbol']
    HEPG2_data_ctrl4.data_df.sort_index(inplace=True)
    HEPG2_data_ctrl4 = HEPG2_data_ctrl4.data_df
    HEPG2_data_ctrl4 = HEPG2_data_ctrl4[HEPG2_data_ctrl4.index.\
                                    isin(test1.index.tolist())]
    HEPG2_data_ctrl4['mean'] = HEPG2_data_ctrl4.mean(axis=1)
    
    #calculate FC
    test11 = gparser.calculate_FC(test1, HEPG2_data_ctrl4['mean'])
    test22 = gparser.calculate_FC(test2, HEPG2_data_ctrl4['mean'])
    test33 = gparser.calculate_FC(test3, HEPG2_data_ctrl4['mean'])
    
    #save matrix  
    test_all = pd.concat([test1, test2, test3], axis=1)
    test_all.to_csv\
            (output_dir+cell_line+'_y_s.csv',\
             index=True, header=True, sep = '\t')