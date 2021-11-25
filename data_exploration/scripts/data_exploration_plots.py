"""
Scripts to generate a series of descriptive plots for perturbation-based expression data. 
Requires: sorted pandas df (GS format), number of replicates, folder to save plot
"""
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 23 12:19:45 2021

@author: erik
"""
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 23 13:23:05 2021

@author: erikzhi
"""


#1. clustergram
def plot_expression_clustergram(your_data):
    import numpy as np
    from heatmapcluster import heatmapcluster
    import matplotlib.pyplot as plt
    from varname import nameof
    
    data_name = nameof(your_data)
    
    fig = plt.figure(figsize=(20,15))
    
    font = {'family' : 'DejaVu Sans',
        'weight' : 'bold',
        'size'   : 22}

    plt.rc('font', **font)
    
    your_data_np = np.round(np.matrix(your_data.loc\
     [:, your_data.columns != 'GENE_SYMBOL']), 4)
    
    h = heatmapcluster(your_data_np, your_data.index, your_data.columns[1:],
                   num_row_clusters=3, num_col_clusters=3,
                   label_fontsize=8,
                   xlabel_rotation=-75,
                   cmap=plt.cm.coolwarm,
                   show_colorbar=True,
                   top_dendrogram=True)
    plt.savefig(str(data_name)+'_plot_expression_clustergram.svg')
    plt.show()
    
    plt.close()
    
    return

#2 pca1
def plot_pca_expression_pc1_vs_pc2(your_data):
    import numpy as np
    import matplotlib.pyplot as plt
    from sklearn.preprocessing import StandardScaler
    from sklearn.decomposition import PCA
    from sklearn.pipeline import Pipeline
    from varname import nameof
    
    data_name = nameof(your_data)
    
    your_data_np = np.round(np.matrix(your_data.loc\
     [:, your_data.columns != 'GENE_SYMBOL']), 4)
    
    pipeline = Pipeline([('scaling', StandardScaler()), ('pca', PCA(n_components=10))])
    pca_data = pipeline.fit_transform(your_data_np.T)
    explained_variance = np.var(pca_data, axis=0)/sum(np.var(pca_data, axis=0))
    
    #1250 for lasso 500 for lesshub
    
    #2. PCA
    test = np.array_split(pca_data,2)
    rep1 = test[0]
    rep2 = test[1]
    
    #pca1
    fig = plt.figure(figsize=(20,15))
    ax = fig.add_subplot(111)
    
    font = {'family' : 'DejaVu Sans',
        'weight' : 'bold',
        'size'   : 22}

    plt.rc('font', **font)
    
    ax.scatter(rep1[:, 0], rep1[:, 1], alpha=0.7, s=75, marker='o', edgecolor ="black",
                linewidths = 0.5, c='green')
    ax.scatter(rep2[:, 0], rep2[:, 1], alpha=0.7, s=75, marker='^', edgecolor ="black",
                linewidths = 0.5, c='blue')
    plt.gca().legend(('rep1','rep2'), loc='upper right', shadow=True, ncol=2)
    
    ax.set_xlabel(('PC1: '+ str(round(explained_variance[0], 2)*100)+'%'))
    ax.set_ylabel(('PC2: '+ str(round(explained_variance[1], 2)*100)+'%'))
    plt.title("PCA, -"+ str(len(rep1)) + " experiments")
    
    fig.savefig(str(data_name)+'_plot_pca_expression_pc1_vs_pc2.svg')
    plt.close()
    
    return

#2 pca2
def plot_pca_expression_pc2_vs_pc3(your_data):
    import numpy as np
    import matplotlib.pyplot as plt
    from sklearn.preprocessing import StandardScaler
    from sklearn.decomposition import PCA
    from sklearn.pipeline import Pipeline
    from varname import nameof
    
    data_name = nameof(your_data)
    
    your_data_np = np.round(np.matrix(your_data.loc\
     [:, your_data.columns != 'GENE_SYMBOL']), 4)
    
    pipeline = Pipeline([('scaling', StandardScaler()), ('pca', PCA(n_components=10))])
    pca_data = pipeline.fit_transform(your_data_np.T)
    explained_variance = np.var(pca_data, axis=0)/sum(np.var(pca_data, axis=0))
    
    #1250 for lasso 500 for lesshub
    
    #2. PCA
    test = np.array_split(pca_data,2)
    rep1 = test[0]
    rep2 = test[1]
    
    #pca1
    fig = plt.figure(figsize=(20,15))
    ax = fig.add_subplot(111)
    
    font = {'family' : 'DejaVu Sans',
        'weight' : 'bold',
        'size'   : 22}

    plt.rc('font', **font)
    
    ax.scatter(rep1[:, 1], rep1[:, 2], alpha=0.7, s=75, marker='o', edgecolor ="black",
                linewidths = 0.5, c='green')
    ax.scatter(rep2[:, 1], rep2[:, 2], alpha=0.7, s=75, marker='^', edgecolor ="black",
                linewidths = 0.5, c='blue')
    plt.gca().legend(('rep1','rep2'), loc='upper right', shadow=True, ncol=2)
    
    ax.set_xlabel(('PC2: '+ str(round(explained_variance[1], 2)*100)+'%'))
    ax.set_ylabel(('PC3: '+ str(round(explained_variance[2], 2)*100)+'%'))
    plt.title("PCA, -"+ str(len(rep1)) + " experiments")
    
    fig.savefig(str(data_name)+'_plot_pca_expression_pc2_vs_pc3.svg')
    plt.close()
    
    return

#4 singular and eigenvals distribution
def plot_singular_and_eigenvals_distr(your_data, reps=2):
    import numpy as np
    import scipy
    import matplotlib.pyplot as plt
    from varname import nameof
    
    data_name = nameof(your_data)
    
    
    your_data_np = np.round(np.matrix(your_data.loc\
     [:, your_data.columns != 'GENE_SYMBOL']), 4)
    
    cols = int(your_data_np.shape[1]/reps)
    rows = int(your_data_np.shape[0])
    ind=0
    
    #2 svd distribution
    svd_distr = scipy.linalg.svdvals(your_data_np)
    svd_distr.sort()

    reshaped_data = np.empty((rows,cols,reps))
    list_eigen = np.empty((cols, reps))
    for i in range(reps):
        reshaped_data[:,:,i] = your_data_np[:,ind:cols]
        eigvals_distr = scipy.linalg.eigvals(your_data_np[:,ind:cols])
        list_eigen[:,i] = eigvals_distr
        ind = ind+cols
        cols = ind+cols

    av_eigenvals = np.abs(np.mean(list_eigen, axis=1))
    av_eigenvals.sort()
    x = np.linspace(1,len(av_eigenvals),len(av_eigenvals))
    
    fig = plt.figure(figsize=(20,15))
    
    font = {'family' : 'DejaVu Sans',
        'weight' : 'bold',
        'size'   : 22}

    plt.rc('font', **font)
    
    plt.plot(x, av_eigenvals, '-o')  # arguments are passed to np.histogram
    plt.plot(x, svd_distr, '-ro')  # arguments are passed to np.histogram
    plt.legend(["eigenvalues", "singular values"], loc='upper center',\
           bbox_to_anchor=(0.5, 1.00), shadow=True, ncol=2)
    plt.title("Dimensionality")
    
    fig.savefig(str(data_name)+'_plot_singular_and_eigenvals_distr.svg')
    plt.close()
    
    return


#5 Standard Error (in expression between replicates) and Y-expression
def plot_median_expression_vs_sd_error(your_data, reps=2):
    
    import numpy as np
    import matplotlib.pyplot as plt
    from varname import nameof
    
    data_name = nameof(your_data)
    
    your_data_np = np.round(np.matrix(your_data.loc\
     [:, your_data.columns != 'GENE_SYMBOL']), 4)
    
    cols = int(your_data_np.shape[1]/reps)
    rows = int(your_data_np.shape[0])
    ind=0
    
    reshaped_data = np.empty((rows,cols,reps))
    #list_eigen = np.empty((cols, reps))
    
    #reshaped_data = np.array(np.hsplit(your_data_np, reps))
    if reps==2:
        for i in range(reps):
            reshaped_data[:,:,i] = your_data_np[:,ind:cols]
            #eigvals_distr = scipy.linalg.eigvals(your_data_np[:,ind:cols])
            #list_eigen[:,i] = eigvals_distr
            ind = ind+cols
            cols = ind+cols
    

    av_variance = np.var(reshaped_data, axis=2)
    av_mean = np.mean(reshaped_data, axis=2)
    
    fig = plt.figure(figsize=(20,15))
    ax = fig.add_subplot(111)
    
    font = {'family' : 'DejaVu Sans',
        'weight' : 'bold',
        'size'   : 22}

    plt.rc('font', **font)
    
    plt.hist(av_mean.flatten(), bins=cols, color='b', alpha=0.4)  # arguments are passed to np.histogram
    
    plt.axvline(x = np.median(av_mean.flatten()), color = 'black', label = 'average var: '+' '+\
                str(np.round(av_variance.mean(), decimals=4)))

    ax.axvspan(np.median(av_mean.flatten())-av_variance.mean(), np.median(av_mean.flatten())+av_variance.mean(),\
               alpha=0.5, color='red')
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.00), shadow=True, ncol=2)
    plt.title("logFC histogram and average variance")
    
    fig.savefig(str(data_name)+'_plot_median_expression_vs_av_variance.svg')
    plt.close()
    
    return


#6. F-matrix score distribution // how much was performed the target gene in all experiments
def plot_pert_score(your_data, reps=3):
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    from varname import nameof
    
    data_name = nameof(your_data)
    
    #get intended perturbation values
    if reps==2:
        step=your_data.shape[0]+1
        rep1 = your_data.iloc[:,1:step]
        rep2 = your_data.iloc[:,step:your_data.shape[1]]
        real1 = pd.Series(np.diag(rep1), index=[rep1.index, rep1.columns])
        real2 = pd.Series(np.diag(rep2), index=[rep2.index, rep2.columns])
        real_all = real1.append(real2)
    else:
        step=your_data.shape[0]+1
        rep1 = your_data.iloc[:,1:step]
        rep2 = your_data.iloc[:,step:step+step-1]
        rep3 = your_data.iloc[:,step+step-1:your_data.shape[1]]

        real1 = pd.Series(np.diag(rep1), index=[rep1.index, rep1.columns])
        real2 = pd.Series(np.diag(rep2), index=[rep2.index, rep2.columns])
        real3 = pd.Series(np.diag(rep3), index=[rep3.index, rep3.columns])
        real_all = real1.append(real2)
        real_all = real_all.append(real3)
    
    #get most negative values
    fmat_real = pd.DataFrame({'real_min': your_data.loc[:, your_data.columns != 'GENE_SYMBOL'].min()})
    fmat_real['int_min'] = real_all.values
    #calculate pert score
    fmat_real['pert score'] = fmat_real['int_min']/fmat_real['real_min']
    
    
    #plot it 
    fig = plt.figure(figsize=(20,15)) 
    font = {'family' : 'DejaVu Sans',
        'weight' : 'bold',
        'size'   : 22}

    plt.rc('font', **font)
    
    plt.hist(fmat_real['pert score'], bins=your_data.shape[0]*3, color='b', alpha=0.4)
    plt.title("Perturbation off-target score")
    fig.savefig(str(data_name)+'_plot_pert_score.svg')
    plt.close()
    
    return

def calc_pert_rank(your_dataset, reps=2):
    
    from scipy.stats import rankdata as rd
    import matplotlib.pyplot as plt
    from varname import nameof
    import numpy as np
    
    fig = plt.figure(figsize=(20,15))
    
    font = {'family' : 'DejaVu Sans',
        'weight' : 'bold',
        'size'   : 50}

    plt.rc('font', **font)
    
    if reps==3:
        pert_mat = np.concatenate((np.eye(your_dataset.shape[0]),
                                   np.eye(your_dataset.shape[0])),\
                                  axis=1)*(-1)   
    
        pert_mat = np.concatenate((pert_mat, np.eye(your_dataset.shape[0])),\
                                  axis=1)*(-1)   
    else:
        pert_mat = np.concatenate((np.eye(your_dataset.shape[0]),np.eye(your_dataset.shape[0])),axis=1)*(-1)
    
    rank_within_sought_r1 = 0
    rank_within_sought_r2 = 0
    rank_within_sought_r3 = 0
    rank_within_sought_r4 = 0
    rank_within_sought_r5 = 0
    rank_within_sought_r6 = 0
    
    print(pert_mat.shape)
    perts = np.nonzero(pert_mat)
    your_dataset = np.round(np.matrix(your_dataset.loc\
     [:, your_dataset.columns != 'GENE_SYMBOL']), 4)
    data_rank = rd(your_dataset,axis=0,method='min')
    
    for i in range(data_rank.shape[1]):
        if data_rank[perts[0][i],perts[1][i]] == 1:
            rank_within_sought_r1 +=1
        elif data_rank[perts[0][i],perts[1][i]] == 2:
            rank_within_sought_r2 +=1
        elif data_rank[perts[0][i],perts[1][i]] == 3:
            rank_within_sought_r3 +=1
        elif data_rank[perts[0][i],perts[1][i]] == 4:
            rank_within_sought_r4 +=1
        elif data_rank[perts[0][i],perts[1][i]] == 5:
            rank_within_sought_r5 +=1
        elif data_rank[perts[0][i],perts[1][i]] <= 6:
            rank_within_sought_r6 +=1
    
    rank_within_sought = [rank_within_sought_r1,rank_within_sought_r2,rank_within_sought_r3,\
                          rank_within_sought_r4,rank_within_sought_r5,rank_within_sought_r6]
    rank_within_sought = list(map(lambda x: (x/data_rank.shape[1])*100, rank_within_sought))
    
    
    names = ['1', '2', '3', "4", '5', '6']
    plt.bar(names, rank_within_sought, width=0.7, color=plt.get_cmap('tab20c').colors, edgecolor='k', 
        linewidth=2)
    plt.xlabel("ranks")
    plt.ylabel("pecrentage of experiments")
    
    data_name = nameof(your_dataset)
    fig.savefig(str(data_name)+'_pert_ranks.svg')
    plt.show()
    
    return

#function to run all plots
def plot_all_explorative_plots(your_dataset, output_dir, cell_line_name, num_of_reps):
    import os
    
    home_dir = '/home/erik/sweden/sonnhammer/scripts'
    directory = cell_line_name
    path = os.path.join(output_dir, directory)
    print(path)
    if not os.path.exists(path):
        os.mkdir(path)
    os.chdir(path)
    
    #call all functions
    plot_expression_clustergram(your_dataset)
    plot_pca_expression_pc1_vs_pc2(your_dataset)
    plot_pca_expression_pc2_vs_pc3(your_dataset)
    #plot_singular_and_eigenvals_distr(your_dataset)
    plot_median_expression_vs_sd_error(your_dataset, num_of_reps)
    plot_pert_score(your_dataset, num_of_reps)
    calc_pert_rank(your_dataset, num_of_reps)
    
    #return to script dir
    os.chdir(home_dir)
    return