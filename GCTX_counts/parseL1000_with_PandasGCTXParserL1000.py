# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.

import matlab.engine
eng = matlab.engine.start_matlab()
tf = eng.isprime(37)
print(tf)

"""

from PandasGCTXParserL1000 import PandasGCTXParserL1000
gparser = PandasGCTXParserL1000()

list_of_cell_lines = ['A375', 'A549', 'HA1E', 'HCC515', 'HEPG2',\
                      'HT29', 'MCF7', 'PC3']
    
list_of_cell_lines = ['HEPG2', 'HA1E']
   
#list_of_cell_lines = ['A375']

#list_of_cell_lines = ['A549', 'HA1E', 'HCC515', 'HEPG2',\
#                      'HT29', 'MCF7', 'PC3']

for cell_line in list_of_cell_lines:

    #read data
    HEPG2_data, HEPG2_data_ctrl = gparser.read_gctx_data(cell_line)
    #subset
    HEPG_data_subset = gparser.filter_gctx_data(HEPG2_data, ["X1","X2","X3"], 3)
    
    #save bulk
    HEPG_data_subset_p1_bulk, HEPG_data_subset_p2_bulk, HEPG_data_subset_p3_bulk = \
        gparser.merge_all_perturbators(HEPG_data_subset, ["X1","X2","X3"])
    #save best shRNA
    HEPG_data_subset_p1_one, HEPG_data_subset_p2_one, HEPG_data_subset_p3_one = \
        gparser.select_one_perturbator(HEPG_data_subset, ["X1","X2","X3"])
    
    #plot PCA selected
    gparser.plot_PCA(HEPG_data_subset_p1_one,\
                     HEPG_data_subset_p2_one,\
                         HEPG_data_subset_p3_one,\
                             HEPG2_data_ctrl,\
                             cell_line, '/scratch/erikzhi/plots/L1000/selected')
    #plot PCA bulk
    gparser.plot_PCA(HEPG_data_subset_p1_bulk,\
                     HEPG_data_subset_p2_bulk,\
                         HEPG_data_subset_p3_bulk,\
                             HEPG2_data_ctrl,\
                             cell_line, '/scratch/erikzhi/plots/L1000/bulk')
        
    #save data
    list_of_reps_s = [HEPG_data_subset_p1_one, HEPG_data_subset_p2_one,\
                                 HEPG_data_subset_p3_one]
    list_of_reps_b = [HEPG_data_subset_p1_bulk, HEPG_data_subset_p2_bulk,\
                                 HEPG_data_subset_p3_bulk]
    gparser.save2genespyder(list_of_reps_s, HEPG2_data_ctrl, cell_line,\
                            '/scratch/erikzhi/L1000_data/matrices/selected/',
                            fold_change=False)
    gparser.save2genespyder(list_of_reps_b, HEPG2_data_ctrl, cell_line,\
                            '/scratch/erikzhi/L1000_data/matrices/bulk/',
                            fold_change=False)
