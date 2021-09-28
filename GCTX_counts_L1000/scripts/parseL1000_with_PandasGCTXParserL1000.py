# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.

import matlab.engine
eng = matlab.engine.start_matlab()
tf = eng.isprime(37)
print(tf)

"""

#parse level5 data
from PandasGCTXParserL1000 import PandasGCTXParserL1000
gparser = PandasGCTXParserL1000()

#params
output_dir = '/scratch/erikzhi/L1000/L1000_data/matrices/level5/'

list_of_cell_lines = ['A375', 'A549', 'HA1E', 'HCC515', 'HEPG2',\
                      'HT29', 'MCF7', 'PC3']

#parse
for cell_line in list_of_cell_lines:
    
    #read data
    HEPG2_data5, HEPG2_data_ctrl5 = gparser.read_gctx_data(cell_line,\
    L1000_gctx_file =\
    "/scratch/erikzhi/L1000/L1000_data/level5/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx",\
    inst_info_file = \
    "/scratch/erikzhi/L1000/L1000_data/LIMS/GSE92742_Broad_LINCS_sig_info.txt")
    
    #select columns with 3 reps
    test = gparser.gctoo2matrices_lvl5(HEPG2_data5)
    
    #save as y matrices
    test.to_csv\
            (output_dir+cell_line+'_lvl5_y.csv',\
             index=True, header=True, sep = '\t')

            
#parse level4 data
import pandas as pd            
from PandasGCTXParserL1000 import PandasGCTXParserL1000
gparser = PandasGCTXParserL1000()

#params
output_dir = '/scratch/erikzhi/L1000/L1000_data/matrices/level4/'
list_of_cell_lines = ['HEPG2']
list_of_cell_lines = ['A375', 'A549', 'HA1E', 'HCC515', 'HEPG2',\
                      'HT29', 'MCF7', 'PC3']

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


#parse level3 data
import pandas as pd    
from PandasGCTXParserL1000 import PandasGCTXParserL1000
gparser = PandasGCTXParserL1000()

list_of_cell_lines = ['HA1E', 'HCC515', 'HEPG2',\
                      'HT29']
#list_of_cell_lines = ['A549']
#problems: 'MCF7','A549', 'PC3', 'A375'
#list_of_cell_lines = ['HEPG2']
output_dir = '/scratch/erikzhi/L1000/L1000_data/matrices/level3/'

for cell_line in list_of_cell_lines:
    print(cell_line)
    #read data
    HEPG2_data, HEPG2_data_ctrl = gparser.read_gctx_data(cell_line)
    #subset
    HEPG_data_subset = gparser.filter_gctx_data(HEPG2_data, ["X1","X2","X3"], 3)
    #save bulk
    HEPG_data_subset_p1_one, HEPG_data_subset_p2_one, HEPG_data_subset_p3_one = \
        gparser.select_one_perturbator(HEPG_data_subset, ["X1","X2","X3"])
    HEPG_data_subset_p1_bulk, HEPG_data_subset_p2_bulk, HEPG_data_subset_p3_bulk = \
        gparser.merge_all_perturbators(HEPG_data_subset, ["X1","X2","X3"])
    
    #prepare replicates
    test1 = gparser.gctoo2matrices_lvl5(HEPG_data_subset_p1_one, rep_counts=1)
    test2 = gparser.gctoo2matrices_lvl5(HEPG_data_subset_p2_one, rep_counts=1)
    test3 = gparser.gctoo2matrices_lvl5(HEPG_data_subset_p3_one, rep_counts=1)
    
    #prepare control
    HEPG2_data_ctrl.data_df.index = HEPG2_data_ctrl.\
                row_metadata_df['pr_gene_symbol']
    HEPG2_data_ctrl.data_df.sort_index(inplace=True)
    HEPG2_data_ctrl = HEPG2_data_ctrl.data_df
    HEPG2_data_ctrl = HEPG2_data_ctrl[HEPG2_data_ctrl.index.\
                                    isin(test1.index.tolist())]
    HEPG2_data_ctrl['mean'] = HEPG2_data_ctrl.mean(axis=1)
    
    #calculate FC
    test1 = gparser.calculate_FC(test1, HEPG2_data_ctrl['mean'])
    test2 = gparser.calculate_FC(test2, HEPG2_data_ctrl['mean'])
    test3 = gparser.calculate_FC(test3, HEPG2_data_ctrl['mean'])
    
    #save FC matrices
    test_all = pd.concat([test1, test2, test3], axis=1)
    test_all.to_csv\
            (output_dir+cell_line+'_y_s.csv',\
             index=True, header=True, sep = '\t')
