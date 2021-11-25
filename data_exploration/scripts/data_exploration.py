#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 23 12:19:45 2021

@author: erik
"""
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 23 10:31:40 2021

@author: erikzhi
"""
#load libs and func
import pandas as pd
import os
import numpy as np
import data_exploration_plots as f_explore


#define dir
data_dir =\
    "/home/erik/sweden/sonnhammer/work-in-progress/data_exploration/data"
output_dir = "/home/erik/sweden/sonnhammer/work-in-progress/data_exploration/plots"

##load ecoli data (parker)
expression_data = pd.\
    read_csv(\
    os.path.join(data_dir,\
    "Parker_y.csv"),\
        header=0, delimiter=",")
        
cols = expression_data.select_dtypes(exclude=['float']).columns
expression_data[cols] = expression_data[cols].apply(pd.to_numeric, downcast='float', errors='coerce')
cell_line = 'ecoli'

#create plots
f_explore.plot_all_explorative_plots(expression_data, output_dir, cell_line, 2)
#----------------


#load yeast data (kemmeren)
expression_data = pd.\
    read_csv(\
    os.path.join(data_dir,\
    "Kemmeren_y.csv"),\
        header=0, delimiter="\t")
        
cols = expression_data.select_dtypes(exclude=['float']).columns
expression_data[cols] = expression_data[cols].apply(pd.to_numeric, downcast='float', errors='coerce')
cell_line = 'yeast'

#create plots
f_explore.plot_all_explorative_plots(expression_data, output_dir, cell_line, 2)
#----------------


#handle ENCODE data
expression_data = pd.\
    read_csv(\
    os.path.join(data_dir,\
    "ENCODE_hepg2_y.csv"),\
        header=0, delimiter="\t")
        
expression_data['GENE_SYMBOL'] = expression_data.index
first_column = expression_data.pop('GENE_SYMBOL')
expression_data.insert(0, 'GENE_SYMBOL', first_column)
expression_data.reset_index(drop=True, inplace=True)

cols = expression_data.select_dtypes(exclude=['float']).columns
expression_data[cols] = expression_data[cols].apply(pd.to_numeric, downcast='float', errors='coerce')
cell_line = 'ENCODE'

#create plots
f_explore.plot_all_explorative_plots(expression_data, output_dir, cell_line, 2)


#Load L1000 level 5
expression_data = pd.\
    read_csv(\
    os.path.join(data_dir,"L1000_A375_lvl5_y_ready.csv"),\
        header=0, delimiter="\t")
        
expression_data.rename(columns={'pr_gene_symbol':'GENE_SYMBOL'}, inplace=True)
cols = expression_data.select_dtypes(exclude=['float']).columns
expression_data[cols] = expression_data[cols].apply(pd.to_numeric, downcast='float', errors='coerce')
cell_line = 'L1000'


#create plots
f_explore.plot_all_explorative_plots(expression_data, output_dir, cell_line, 3)
