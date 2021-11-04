#!/usr/bin/env python3
# -*- coding: utf-8 -*-

class PandasGCTXParserL1000:
    """
    A set of wrappers to parse L1000 data with pandasGEXpress package/
    https://clue.io/cmapPy/pandasGEXpress.html
    /
    LINCS user guide: https://docs.google.com/document/d/1q2gciWRhVCAAnlvF2iRLuJ7whrGP6QjpsCMq1yWz7dU/
    Created on Tue Dec  1 11:37:19 2020

    @author: Erik Zhivkoplias
    """
    
    def __init__(self):
        
        """import libraries
        """
        
        print('loaded')
        return

    def read_gctx_data(self, cell_line,\
                       L1000_gctx_file, gene_info_file, inst_info_file, level, hrs="96"):
        """
    
        Parameters
        ----------
        L1000_gctx_file : name of gctx file, str.
        gene_info_file : name of level3 annotation file, str
        inst_info_file : name of level3 annotation file, str
        cell_line : name of the cell line, str.
    
        !see list of files here: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92742    
    
        Returns
        -------
        gctoo cell-line specific instances (exp and ctrl).
    
        """
        self.cell_line = cell_line
        
        #import libs
        import pandas as pd
        from cmapPy.pandasGEXpress.parse import parse
        
        #read meta info
        inst_info = pd.read_csv\
        (inst_info_file,\
         sep="\t", dtype=str)
        
        gene_info = pd.read_csv\
            (gene_info_file,\
             sep="\t", dtype=str)
        #ids
        if level==5:
            gctx_ids = 'sig_id'
        else:
            gctx_ids = 'inst_id'
        #lm genes
        landmark_gene = gene_info[gene_info["pr_is_lm"] == "1"]  
        landmark_gene_row_ids = landmark_gene["pr_gene_id"]
        landmark_gene_names = landmark_gene["pr_gene_symbol"]
       
        #parse meta info
        #experiments
        trt_sh = inst_info[(inst_info["pert_type"] == "trt_sh") &\
                           (inst_info["pert_time"] == hrs) &\
                           (inst_info["cell_id"] == cell_line) &\
                            (inst_info["pert_iname"].isin(landmark_gene_names))]
               
        
        trt_sh_ids = inst_info[gctx_ids][(inst_info["pert_type"] == "trt_sh") &\
                                           (inst_info["pert_time"] == hrs) &\
                                           (inst_info["cell_id"] == cell_line) &\
                            (inst_info["pert_iname"].isin(landmark_gene_names))]
        
        #print(trt_sh_ids)
        #ctrl
        ctrl = inst_info[(inst_info["pert_type"] == "ctl_vector") &\
                                                  (inst_info['pert_iname'] == "EMPTY_VECTOR") &\
                                                      (inst_info["pert_time"] == hrs) &\
                                                  (inst_info["cell_id"] == cell_line)]
        ctrl_ids = inst_info[gctx_ids][(inst_info["pert_type"] == "ctl_vector") &\
                                                  (inst_info['pert_iname'] == "EMPTY_VECTOR") &\
                                                      (inst_info["pert_time"] == hrs) &\
                                                  (inst_info["cell_id"] == cell_line)]
        #print(ctrl_ids)
        
        
                
        #subset gctx files
        level3_data_experiments = parse\
                    (L1000_gctx_file,\
                                  rid = landmark_gene_row_ids, cid = trt_sh_ids)          
        level3_data_ctrl = parse\
                    (L1000_gctx_file,\
                                  rid = landmark_gene_row_ids, cid = ctrl_ids)
                        
                
        #annotate instances
        trt_sh.set_index(gctx_ids, inplace=True) 
        ctrl.set_index(gctx_ids, inplace=True)
        landmark_gene.set_index("pr_gene_id", inplace=True)
        
        level3_data_experiments.col_metadata_df = trt_sh
        level3_data_experiments.row_metadata_df = landmark_gene
        level3_data_ctrl.col_metadata_df = ctrl
        level3_data_ctrl.row_metadata_df = landmark_gene
    
                
        return level3_data_experiments, level3_data_ctrl
    
    def filter_gctx_data(self, gctoo_instance, list_of_plates,\
                         num_of_plates):
        """

        Parameters
        ----------
        gctoo_instance : TYPE
            DESCRIPTION.
        list_of_plates : TYPE
            DESCRIPTION.
        num_of_plates : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        
        #import libs
        import pandas as pd
        import cmapPy.pandasGEXpress.subset_gctoo as sg
        
        
        #subset experiments performed on the plates with valid names
        col_meta_data = gctoo_instance.col_metadata_df
        col_meta_data['plate_num'] =\
            col_meta_data['rna_plate'].str.split('_',3,True)[3]  
        all_plate_ids =\
        col_meta_data.index[col_meta_data.plate_num.isin(list(list_of_plates))]
        
        gctoo_instance = sg.subset_gctoo(gctoo_instance,\
                                          cid=list(all_plate_ids))
        col_meta_data = gctoo_instance.col_metadata_df
   


        #select pertubated genes with at least "num_of_plates" plates
        mask_col_data = col_meta_data[['rna_well', 'pert_id']]
        mask_col_data['key_to_filter'] = mask_col_data.values.tolist()
        #create dict with cids
        dict_with_cids = pd.Series(tuple(mask_col_data.key_to_filter.values), index=mask_col_data.index).to_dict()     
        del mask_col_data['key_to_filter'] 

        #create list of keys
        mask = mask_col_data.groupby(mask_col_data.columns.tolist()).size().reset_index().\
            rename(columns={0:'counts'})
        mask_pert_id_filter = mask[mask.counts>=num_of_plates]
        mask_pert_id_filter = mask_pert_id_filter[["rna_well", 'pert_id']]

        mask_pert_id_filter['key_to_filter'] = mask_pert_id_filter.values.tolist()
        list_of_keys = mask_pert_id_filter['key_to_filter'].values.tolist()

        #create list of cids
        list_of_cids = []
        
        #add cid to a list if key is present in dictionary with cids
        for key in list_of_keys:
            if [k for k,v in dict_with_cids.items() if v == key]:
                list_of_cids.append([k for k,v in dict_with_cids.items() if v == key])

        list_of_cids = [item for sublist in list_of_cids for item in sublist]

        #subset gctoo instance with list of cids
        gctoo_instance_subset = sg.subset_gctoo\
            (gctoo_instance, cid=list_of_cids)
            
        return gctoo_instance_subset
    
    
    def merge_tech_duplicates(self, gctoo_instance, column_name, min_shRNAs_num=1):
        """
        

        Parameters
        ----------
        gctoo_instance : TYPE
            DESCRIPTION.
        column_name : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        #import libs
        import pandas as pd
        import numpy as np
        import cmapPy.pandasGEXpress.GCToo as GCToo
        import cmapPy.pandasGEXpress.subset_gctoo as sg
        import cmapPy.pandasGEXpress.concat as cg 
        
        
        #select pert_iname duplicated on each plate
        col_meta_data = gctoo_instance.col_metadata_df
        mask_pert_ids_non_dup = \
                gctoo_instance.col_metadata_df.groupby\
                    ([column_name]).size().loc[lambda x: x == 1].index        
        mask_pert_ids_dup = \
                gctoo_instance.col_metadata_df.groupby\
                ([column_name]).size().loc[lambda x: x > min_shRNAs_num].index
                
        #subset without duplicates
        data_rep1_nondup_ids = col_meta_data.index\
                    [col_meta_data[column_name].isin(list(mask_pert_ids_non_dup))]
        data_rep1_nondup = sg.subset_gctoo\
                    (gctoo_instance, cid=list(data_rep1_nondup_ids))
                
        #subset duplicates
        data_rep1_dup_ids = col_meta_data.index\
                [col_meta_data[column_name].isin(list(mask_pert_ids_dup))]
        data_rep1_dup = sg.subset_gctoo\
                    (gctoo_instance, cid=list(data_rep1_dup_ids))
        
        #merged gctoo instances
        new_gctoo_instances = []
        
        #calculate means for each duplicated perturbator
        for pert_dup in mask_pert_ids_dup:
            #select cis
            pert_dup_ids = list(col_meta_data.index\
                            [col_meta_data[column_name] == pert_dup])
            #calculate mean expression value for them
            merged_pert_dup_array = data_rep1_dup.\
                    data_df[pert_dup_ids].mean(axis=1).values
            
            cid_index = np.where(data_rep1_dup.col_metadata_df.index.values\
                                 == pert_dup_ids[0])
            #create cid pd
            cid_db = data_rep1_dup.col_metadata_df.iloc[cid_index]
            
            #create pd for expression values
            merged_pert_dup_db =\
                    pd.DataFrame(merged_pert_dup_array,
                        index = data_rep1_nondup.row_metadata_df.index.values,
                        columns=[pert_dup_ids[0]])
            
            #create new gctoo instances
            merged_pert_dup_gctoo =\
                    GCToo.GCToo(data_df=pd.DataFrame(merged_pert_dup_db), 
                        row_metadata_df=data_rep1_nondup.row_metadata_df.copy(), 
                        col_metadata_df=cid_db, 
                        make_multiindex=True)
            
            #add to the list
            new_gctoo_instances.append(merged_pert_dup_gctoo)
        if min_shRNAs_num == 1:
            new_gctoo_instances.append(data_rep1_nondup)
        
        data_rep_cleaned = cg.hstack(new_gctoo_instances)
            
        return data_rep_cleaned
    
    def merge_all_perturbators(self, gctoo_instance, list_of_plates, min_shRNAs_num):
        """Merging shRNA experiments (with the same pert_iname)
            Input: gctoo_instance
            Returns
        -------
        data_rep1, data_rep2, data_rep3: 3 gctoo instances, plate-selected 
        
        
        """
        #import libs
        import cmapPy.pandasGEXpress.subset_gctoo as sg
        
        #gctoo_instance = HEPG2_data_3_reps
        col_meta_data = gctoo_instance.col_metadata_df
        
        
        #select reps
        rep1_ids = col_meta_data.index[col_meta_data["plate_num"] ==\
                                       list_of_plates[0]]
        data_rep1 = sg.subset_gctoo\
                (gctoo_instance, cid=list(rep1_ids))
                    
        rep2_ids = col_meta_data.index[col_meta_data["plate_num"] ==\
                                       list_of_plates[1]]
        data_rep2 = sg.subset_gctoo\
                (gctoo_instance, cid=list(rep2_ids))
        
        rep3_ids = col_meta_data.index[col_meta_data["plate_num"] ==\
                                       list_of_plates[2]]
        data_rep3 = sg.subset_gctoo\
                (gctoo_instance, cid=list(rep3_ids))
        
        data_rep1_all = self.merge_tech_duplicates(data_rep1, "pert_iname", min_shRNAs_num)
        data_rep2_all = self.merge_tech_duplicates(data_rep2, "pert_iname", min_shRNAs_num)
        data_rep3_all = self.merge_tech_duplicates(data_rep3, "pert_iname", min_shRNAs_num)
        
        return data_rep1_all, data_rep2_all, data_rep3_all
    
    def select_one_perturbator(self, gctoo_instance, list_of_plates):
        
        """

        Parameters
        ----------
        gctoo_instance: gctoo instance
        col_meta_data: meta info 
    
        Returns
        -------
        data_rep1, data_rep2, data_rep3: 3 gctoo instances, plate-selected
    
        """
        #import libs
        import cmapPy.pandasGEXpress.subset_gctoo as sg
        import cmapPy.pandasGEXpress.concat as cg
        
        col_meta_data = gctoo_instance.col_metadata_df
        pert_gctoo_names = []
        
        for pert_name in col_meta_data['pert_iname'].unique():
            
            #subset shRNA
            DDR_ids = col_meta_data.index[col_meta_data["pert_iname"] == pert_name]
            data_DDR = sg.subset_gctoo(gctoo_instance, cid=list(DDR_ids))
            col_meta_data_DDR = data_DDR.col_metadata_df
            shRNA_types = col_meta_data_DDR['pert_id'].unique()
        
            dict_shRNA_type_var = {}
            for shRNA_type in shRNA_types:
                #print(shRNA_type)
                #store cids and subset data
                col_meta_data_DDR_oneshRNA_ids = \
                col_meta_data_DDR.index[col_meta_data_DDR['pert_id'] == shRNA_type]
                data_DDR_oneshRNA = sg.subset_gctoo\
                (gctoo_instance, cid=list(col_meta_data_DDR_oneshRNA_ids))
                #data_DDR_oneshRNA.data_df
                
            
                #calculate variance between shRNAs replicates (sum all genes)
                shRNA_type_var = sum(data_DDR_oneshRNA.data_df.var(axis=1))/len(data_DDR_oneshRNA.data_df)
                #add shRNA type and variance to dictionary
                dict_shRNA_type_var[shRNA_type] = shRNA_type_var
        
            #select shRNA with the min variance
            dict_shRNA_type_var_sorted = {k: v for k, v in sorted(dict_shRNA_type_var.items(),\
                                                     key=lambda item: item[1])}
            dict_shRNA_type_var_keys = list(dict_shRNA_type_var_sorted.keys())
            #print(dict_shRNA_type_var_sorted)
            #choose the one with the least variance across plates
            min_var_shRNA = dict_shRNA_type_var_keys[0]
        
            shRNA_id = col_meta_data_DDR.index\
            [col_meta_data_DDR['pert_id'] == min_var_shRNA]
            
            data_DDR_oneshRNA = sg.subset_gctoo\
                (gctoo_instance, cid=list(shRNA_id))
            pert_gctoo_names.append(data_DDR_oneshRNA)
        
        merged_instance = cg.hstack(pert_gctoo_names)
        col_meta_data = merged_instance.col_metadata_df
        
        #select reps
        rep1_ids = col_meta_data.index[col_meta_data["plate_num"] ==\
                                       list_of_plates[0]]
        data_rep1 = sg.subset_gctoo\
                (merged_instance, cid=list(rep1_ids))
        #data_rep1.col_metadata_df.sort_values('pert_id').head(40)
                
        rep2_ids = col_meta_data.index[col_meta_data["plate_num"] ==\
                                       list_of_plates[1]]
        data_rep2 = sg.subset_gctoo\
                (merged_instance, cid=list(rep2_ids))
        
        rep3_ids = col_meta_data.index[col_meta_data["plate_num"] ==\
                                       list_of_plates[2]]
        data_rep3 = sg.subset_gctoo\
                (merged_instance, cid=list(rep3_ids))
                
        data_rep1_all = self.merge_tech_duplicates(data_rep1, "pert_id")
        data_rep2_all = self.merge_tech_duplicates(data_rep2, "pert_id")
        data_rep3_all = self.merge_tech_duplicates(data_rep3, "pert_id")
        
        return data_rep1_all, data_rep2_all, data_rep3_all
    
    def plot_PCA(self, gctoo_instance_1, gctoo_instance_2, gctoo_instance_3,\
                 gctoo_instance_ctrl, cell_line, output_dir):
        """
        

        Parameters
        ----------
        gctoo_instance_1 : TYPE
            DESCRIPTION.
        gctoo_instance_2 : TYPE
            DESCRIPTION.
        gctoo_instance_3 : TYPE
            DESCRIPTION.
        cell_line : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        
        #import libs
        import pandas as pd
        import numpy as np
        import plotly.graph_objs as go
        from plotly.graph_objs import Scene, XAxis, YAxis, ZAxis
        from sklearn.decomposition import PCA
        from sklearn.preprocessing import StandardScaler
        
        gctoo_instance_1 = gctoo_instance_1.data_df.add_prefix('rep1_')
        gctoo_instance_2 = gctoo_instance_2.data_df.add_prefix('rep2_')
        gctoo_instance_3 = gctoo_instance_3.data_df.add_prefix('rep3_')
        gctoo_instance_ctrl = gctoo_instance_ctrl.data_df.add_prefix('ctrl')
            
        labels = {'experiments':list(gctoo_instance_1.columns)+\
                      list(gctoo_instance_2.columns)+list(gctoo_instance_3.columns)+\
                          list(gctoo_instance_ctrl.columns),
                      'plate':(['plate 1'] * len(gctoo_instance_1.columns))+\
                          (['plate 2'] * len(gctoo_instance_2.columns))+\
                              (['plate 3'] * len(gctoo_instance_3.columns))+\
                              (['control'] * len(gctoo_instance_ctrl.columns))} 
                
        my_labels = pd.DataFrame(labels)
            
        my_data = pd.concat([gctoo_instance_1, gctoo_instance_2,\
                             gctoo_instance_3, gctoo_instance_ctrl],\
                                axis=1)
                
          
        #convert data for PCA
        scaler = StandardScaler()
        X_train = my_data.T
        X_train_scl = scaler.fit_transform(X_train)
            
        #explained var
        components = 10
        pca = PCA(n_components=components)
        Y = pca.fit(X_train_scl)
        var_exp = Y.explained_variance_ratio_
        cum_var_exp = np.cumsum(var_exp)
            
        # Plot the explained variance
        x = ["PC%s" %i for i in range(1,components)]
        trace1 = go.Bar(
                x=x,
                y=list(var_exp),
                name="Explained Variance")
            
        trace2 = go.Scatter(
                x=x,
                y=cum_var_exp,
                name="Cumulative Variance")
            
        layout = go.Layout(
                title='Explained variance '+str(cell_line)+"",
                xaxis=dict(title='Principle Components', tickmode='linear'))
            
        data = [trace1, trace2]
        fig = go.Figure(data=data, layout=layout)
        fig.write_html(output_dir+str(cell_line)+'_explained.html',\
                           auto_open=True)
            
        # Project first three components
        Y_train_pca = pca.fit_transform(X_train_scl)
            
        traces = []
        for name in ['plate 1', 'plate 2', 'plate 3', 'control']:
            trace = go.Scatter3d(
                    x=Y_train_pca[my_labels.plate==name,0],
                    y=Y_train_pca[my_labels.plate==name,1],
                    z=Y_train_pca[my_labels.plate==name,2],
                    mode='markers',
                    name=name,
                    marker=go.Marker(size=10, line=go.Line(width=1),opacity=1))
                
            traces.append(trace)
            
        layout = go.Layout(
            showlegend=True,
            scene=Scene(
                    xaxis=XAxis(title='PC1: '\
                                + str(round(pca.explained_variance_ratio_[0]*100, 3))\
                                    + '%'),
                    yaxis=YAxis(title='PC2: '\
                                + str(round(pca.explained_variance_ratio_[1]*100, 3))\
                                    + '%'),
                    zaxis=ZAxis(title='PC3: '\
                                + str(round(pca.explained_variance_ratio_[2]*100, 3))\
                                    + '%')
                    ),
                title="Projection of First Three Principle Components "+str(cell_line)+""
                )
            
        data = traces
        fig = go.Figure(data=data, layout=layout)
        fig.write_html(output_dir+str(cell_line)+'_PCA.html',\
                           auto_open=True)
            
        return True


    def gctoo2matrices(self, gctoo_instance_cleaned_experiments,
                   gctoo_instance_control, fold_change=True):
        """
        
        function to convert pre-processed experimental data to matrices
    
        Parameters
        ----------
        gctoo_instance_cleaned_experiments : cell line perturbed with shRNA,
                                            one rep subsetm gctoo instance
        gctoo_instance_control : cell line perturbed with empty vectors,
                                gctoo instance
    
        Returns
        -------
        pd dataframe with normalized values (Y matrix)
        pd dataframe with perturbation design (P matrix)
    
        """
        #import libs
        import pandas as pd
        import numpy as np
        
        #annotate and sort ctrl
        gctoo_instance_control.data_df.index = gctoo_instance_control.row_metadata_df['pr_gene_symbol']
        ctrl_data =gctoo_instance_control.data_df
        ctrl_data.sort_index(inplace=True)
        ctrl_data['mean'] = ctrl_data.mean(axis=1)
        
        #annotate shRNA experiments
        gctoo_instance_cleaned_experiments.data_df.columns =\
        gctoo_instance_cleaned_experiments.col_metadata_df['pert_iname']             
        gctoo_instance_cleaned_experiments.data_df.index =\
        gctoo_instance_cleaned_experiments.row_metadata_df['pr_gene_symbol']
        pert_data = gctoo_instance_cleaned_experiments.data_df
        
        #sort shRNA experiments
        pert_data = pert_data.\
                       reindex(sorted(pert_data.columns),axis=1)
        pert_data.sort_index(inplace=True)
        
        
        #create expression matrix
        eps = 1e-7
        
        if fold_change==False:
            expression_data_matrix_log = pert_data.copy()
            
        else:
            expression_data_matrix_log = pert_data.div(ctrl_data['mean'].tolist(),\
                                           axis=0)
                
        expression_data_matrix_log = expression_data_matrix_log + eps
        expression_data_matrix_log = expression_data_matrix_log.apply(np.log2)
        
        #sort expression matrix             
        expression_data_matrix_log = expression_data_matrix_log.drop(\
        columns=[col for col in expression_data_matrix_log if col\
               not in expression_data_matrix_log.index.tolist()])
        expression_data_matrix_log = expression_data_matrix_log[expression_data_matrix_log.index.\
                                    isin(expression_data_matrix_log.columns.tolist())]

        
        expression_data_matrix_log = expression_data_matrix_log.\
                       reindex(sorted(expression_data_matrix_log.columns),axis=1)
        expression_data_matrix_log.sort_index(inplace=True)
        
        
        #plt.hist(my_data1_log.to_numpy().flatten(), bins=200)
        
        #z-score normalization
        #my_data1_log_z = my_data1_log
        
        #data_matrix_log_minmax = data_matrix_log.copy()
        
        #annotate experiments
        #data_matrix_log_minmax.data_df.index = data_matrix_log_minmax.row_metadata_df['pr_gene_symbol']
        #data_matrix_log_minmax.data_df.columns = data_matrix_log_minmax.col_metadata_df['pert_iname']
        #data_matrix_log_minmax.data_df.sort_index(inplace=True)
        #data_matrix_log_minmax.data_df = data_matrix_log_minmax.data_df.\
        #        reindex(sorted(data_matrix_log_minmax.data_df.columns),axis=1)
        #data_matrix_log_minmax = data_matrix_log_minmax[data_matrix_log_minmax.\
         #                                                               index.isin(data_matrix_log_minmax.columns.tolist())]
        #data_matrix_log_minmax.sort_index(inplace=True)
        #data_matrix_log_minmax = data_matrix_log_minmax.reindex(sorted(data_matrix_log_minmax.columns),axis=1)
        
      #  for row in data_matrix_log_minmax.iterrows():
      #      data_matrix_log_minmax[row] = (data_matrix_log_minmax[row] -\
      #                           data_matrix_log_minmax[row].mean())\
       #         /data_matrix_log_minmax[row].std(ddof=0)
        #plt.hist(my_data1_log_z.to_numpy().flatten(), bins=200)
        
        
        #min-max normalization // well not exactly, between -2 and 1
        #minmax*(0.5 - -2) + -2 
        #data_matrix_log_minmax = data_matrix_log.copy()
        
        #for col in data_matrix_log_minmax.columns:
        #    data_matrix_log_minmax[col] = (data_matrix_log_minmax[col].values\
        #                                - data_matrix_log_minmax[col].min())/\
       #                          (data_matrix_log_minmax[col].max()\
        #                          - data_matrix_log_minmax[col].min())
        #    data_matrix_log_minmax[col] = (data_matrix_log_minmax[col]*(2.5))-2
                
        #plt.hist(my_data1_log_minmax.to_numpy().flatten(), bins=200)
        
        #create pert matrix
        pert_data_matrix_log = expression_data_matrix_log.copy()
        #print(pert_data_matrix_log_minmax)
        
        for column in pert_data_matrix_log.columns:
            boolean_true_ind = pert_data_matrix_log[column].index == column
            real_pert = [i for i, val in enumerate(boolean_true_ind) if val]
            pert_data_matrix_log[column].iloc[real_pert[0]] = 1
        
        pert_data_matrix_log[pert_data_matrix_log != 1] = 0
        pert_data_matrix_log[pert_data_matrix_log == 1] = -1
        
        #return matrices
        return expression_data_matrix_log, pert_data_matrix_log


    def save2genespyder(self, gctoo_instance_cleaned_experiments_list,
                   gctoo_instance_control, cell_line, output_dir,\
                       inst_info_file = \
   "/scratch/erikzhi/L1000_data/LIMS/GSE92742_Broad_LINCS_inst_info.txt",
                       gene_info_file =\
   "/scratch/erikzhi/L1000_data/LIMS/GSE92742_Broad_LINCS_gene_info.txt",
   fold_change=True):
        """
        Save Y and P matrices as csv file

        Parameters
        ----------
        list_of_matrices : TYPE
            DESCRIPTION.
        output_dir : TYPE
            DESCRIPTION.

        Returns
        -------
        bool
            DESCRIPTION.

        """
        #import libs
        import numpy as np
        import pandas as pd
        
        #define lm genes
        #gene_info = pd.read_csv\
        #    (gene_info_file,\
        #     sep="\t", dtype=str)
        #landmark_gene = gene_info[gene_info["pr_is_lm"] == "1"] 
        
        #save matrices
        
        y1, p1 = self.gctoo2matrices(gctoo_instance_cleaned_experiments_list[0],
                                gctoo_instance_control)
        y2, p2 = self.gctoo2matrices(gctoo_instance_cleaned_experiments_list[1],
                                gctoo_instance_control)
        y3, p3 = self.gctoo2matrices(gctoo_instance_cleaned_experiments_list[2],
                                gctoo_instance_control)
        
        #columns to drop
        cols_to_drop = np.concatenate\
            ((np.setdiff1d(np.union1d(y3.columns, y1.columns), np.intersect1d(y3.columns, y1.columns)),\
              np.setdiff1d(np.union1d(y3.columns, y2.columns), np.intersect1d(y3.columns, y2.columns))))
        
        #final sort
        for matrix in [y1,p1,y2,p2,y3,p3]:
            matrix.drop(cols_to_drop, inplace=True, axis=1, errors='ignore')
            matrix.sort_index(axis=1, inplace=True)
            
        y_all = pd.concat([y1,y2,y3],axis=1)
        p_all = pd.concat([p1,p2,p3],axis=1)
        
        y_all.to_csv\
            (output_dir+cell_line+'_y.csv',\
             index=True, header=True, sep = '\t')
        
        p_all.to_csv\
            (output_dir+cell_line+'_p.csv',\
             index=True, header=True, sep = '\t')
                
        return True
    
    
    def gctoo2matrices_lvl5(self, gctoo_instance_lvl5, rep_counts=1):
        """
        
        function to convert pre-processed experimental data (level 5) to matrices
    
        Parameters
        ----------
        gctoo_instance_lvl5 : cell line perturbed with shRNA,
                            FC values
        gctoo_instance_control : cell line perturbed with empty vectors,
                                gctoo instance
    
        Returns
        -------
        pd dataframe with expression values (Y matrix) in GS format
    
        """
        import pandas as pd
        import cmapPy.pandasGEXpress.subset_gctoo as sg
        
        mask_col_data = gctoo_instance_lvl5.col_metadata_df['pert_iname']
        
        mask_col_data_ktf = mask_col_data.values.tolist()
        #create dict with cids
        dict_with_cids = pd.Series(tuple(mask_col_data_ktf),\
                                   index=mask_col_data.index).to_dict()     
    
        #create list of keys
        mask = mask_col_data.groupby(mask_col_data.tolist()).size().reset_index().\
            rename(columns={'pert_iname':'counts'})
        mask_pert_id_filter = mask[mask.counts==rep_counts]
    
        list_of_keys = mask_pert_id_filter['index'].values.tolist()
    
        #create list of cids
        list_of_cids = []
            
        #add cid to a list if key is present in dictionary with cids
        for key in list_of_keys:
            if [k for k,v in dict_with_cids.items() if v == key]:
                list_of_cids.append([k for k,v in dict_with_cids.items() if v == key])
    
        list_of_cids = [item for sublist in list_of_cids for item in sublist]
    
        #subset gctoo instance with list of cids
        gctoo_instance_lvl5 = sg.subset_gctoo\
            (gctoo_instance_lvl5, cid=list_of_cids)
            
        #subset pr genes
        pr_genes = gctoo_instance_lvl5.row_metadata_df['pr_gene_symbol']
        dict_with_cids = pd.Series(tuple(pr_genes.tolist()), index=pr_genes.index).to_dict() 
        
        col_names = gctoo_instance_lvl5.col_metadata_df['pert_iname'].values.tolist()
        
        #create list of cids
        list_of_rids = []
            
        #add cid to a list if key is present in dictionary with cids
        for key in col_names:
            if [k for k,v in dict_with_cids.items() if v == key]:
                list_of_rids.append([k for k,v in dict_with_cids.items() if v == key])
    
        list_of_rids = list(set([item for sublist in list_of_rids for item in sublist]))
        
        gctoo_instance_lvl5 = sg.subset_gctoo\
            (gctoo_instance_lvl5, rid=list_of_rids)
            
        
        #annotate
        gctoo_instance_lvl5.data_df.index = gctoo_instance_lvl5.\
            row_metadata_df['pr_gene_symbol']
        gctoo_instance_lvl5.data_df.sort_index(inplace=True)
        
        gctoo_instance_lvl5.data_df.columns = gctoo_instance_lvl5.col_metadata_df['pert_iname'].\
        values+'_'+gctoo_instance_lvl5.col_metadata_df['pert_id'].values+\
        gctoo_instance_lvl5.col_metadata_df['pert_dose'].values
        
        #drop dupl columns
        selected_expression = gctoo_instance_lvl5.data_df
        selected_expression.columns = pd.io.parsers.\
            ParserBase({'names':selected_expression.columns}).\
            _maybe_dedup_names(selected_expression.columns)
        #test = test.loc[:,~test.columns.duplicated()]
        
        #sort columns
        selected_expression = selected_expression.\
            reindex(sorted(selected_expression.columns), axis=1)
        
        #return matrices
        return selected_expression
    
    def calculate_FC(self, your_dataset, column_with_ctrl):
        import numpy as np
        #create expression matrix
        eps = 1e-7
        #print(your_dataset.index)
        #print(column_with_ctrl.index)
        column_with_ctrl = column_with_ctrl + eps
        
        your_dataset = your_dataset[your_dataset.index.isin(column_with_ctrl.index.tolist())]
        #your_dataset = your_dataset.drop(columns=[col for col in your_dataset if col not in column_with_ctrl.index.tolist()])
        #your_dataset.sort_index(inplace=True)
        column_with_ctrl = column_with_ctrl[column_with_ctrl.index.isin(your_dataset.index.tolist())]
        
        FC_dataset = your_dataset.div(column_with_ctrl.tolist(),axis=0)
                
        FC_dataset = FC_dataset + eps
        FC_dataset = FC_dataset.apply(np.log2)
    
        return FC_dataset
    
    def filter_overlapping_experiments(self, rep_matrix, common_labels):
        """
        filter out genes-esperiments pairs that are not present in all three reps
        """
        all_experiments = rep_matrix.columns.values.tolist()
        experiment_labels = [experiment.split('_',1)[0] for experiment in all_experiments]
        experiments_dictionary = dict(zip(all_experiments, experiment_labels))
        
        for key, value in experiments_dictionary.copy().items():
            if value not in common_labels:
                del experiments_dictionary[key]
        
        rep_matrix =\
        rep_matrix[rep_matrix.\
                            columns.intersection(list(experiments_dictionary.keys()))]
        rep_matrix =\
        rep_matrix[rep_matrix.index.\
                                isin(list(experiments_dictionary.values()))]
        
        return rep_matrix
        
