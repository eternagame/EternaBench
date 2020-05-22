import numpy as np
import pandas as pd

import seaborn as sns
import matplotlib.pyplot as plt

from glob import glob
import sys, os
from copy import deepcopy

def plot_reac_heatplot(df, ind_range=None, **kwargs):
    '''Plot heatplot image of reactivities.
    Input: full_df style dataframe.'''

    if ind_range is None:
    	start,finish = 0,-1
    else:
    	start, finish=ind_range

    max_len = np.max([len(x) for x in df['reactivity'][start:finish]])
    arr= []

    for x in df['reactivity'][start:finish]:
        arr.append(np.concatenate([x,np.zeros(max_len-len(x))]))
    plt.imshow(np.array(arr).T, cmap='gist_heat_r',vmin=0,vmax=1.5, **kwargs)
    plt.xlabel('Construct')
    plt.ylabel('Sequence position')

def plot_estimate_heatplot(df, ind_range=None, package='vienna_2', **kwargs):
    '''Plot heatplot image of predicted p(unp) values.
    Input: full_df style dataframe.'''

    if ind_range is None:
    	start,finish = 0,-1
    else:
    	start, finish=ind_range
    	
    max_len = np.max([len(x) for x in df['reactivity'][start:finish]])
    arr= []

    for x in df['p_%s'% package][start:finish]:
        arr.append(np.concatenate([x,np.zeros(max_len-len(x))]))
    plt.imshow(np.array(arr).T, cmap='gist_heat_r',vmin=0,vmax=1, **kwargs)
    plt.xlabel('Construct')
    plt.ylabel('Sequence position')

def plot_EB(corr_df, package_list=None, split_by_nucleotides=False):
    '''Plot in eternabench manuscript style from correlation dataframe.'''
    
    packages = ['vienna_1', 'vienna_2', 'vienna_2_nodangles',
     'vienna_langdon_pars', 'vienna_rnasoft_pars',
    'vienna_2_60C',
    'nupack_95', 'nupack_95_nodangles', 'nupack_99', 'nupack_99_nodangles',
    'rnastructure', 'rnastructure_nocoax',
    'contrafold_1', 'contrafold_2', 'contrafold_2_nc',  
    'rnasoft_99','rnasoft_99_nodangles', 'rnasoft_07',
    'rnasoft_blstar', 'rnasoft_bl_nodangles', 'rnasoft_lam-cg', 'rnasoft_nom-cg',
    'eternafold_A','eternafold_B','eternafold_C','eternafold_D','eternafold_E',
    'eternafold_F', 'eternafold_G']

    package_titles = ['Vienna 1', 'Vienna 2', 'Vienna 2, no dangles',
     'Vienna 2, `Langdon 2018` params', 'Vienna 2, `RNASoft 2007` params',
      'Vienna 2, 60C',
    'NUPACK 1995', 'NUPACK 1995, no dangles', 'NUPACK 1999', 'NUPACK 1999, no dangles',
    'RNAstructure', 'RNAstructure, no coaxial stacking',
    'CONTRAfold 1', 'CONTRAfold 2', 'CONTRAfold 2, noncomplementary',  
    'RNAsoft 1999', 'RNAsoft 1999, no dangles', 'RNAsoft 2007',
    'RNAsoft BLstar', 'RNAsoft BL, no dangles', 'RNAsoft LAM-CG', 'RNAsoft NOM-CG','EternaFold_SRR',
    'EternaFold_SCRR','Eternafold_C','EternaFold_D','EternaFold_E','EternaFold_F','EternaFold_G']
    
    sort_order = ['vienna', 'nupack','rnastructure','contrafold','rnasoft','eternafold']
    
    corr_df['package_type'] = corr_df.apply(lambda row: row['package'].split('_')[0], axis=1)

    new_cdf = pd.DataFrame()
    
    for pkg in sort_order:
        tmp_df = corr_df.loc[corr_df['package_type']==pkg]
        new_cdf = new_cdf.append(tmp_df, ignore_index = True)
        
    corr_df = new_cdf
    
    nuc_palette = [[255/255, 211/255, 0],
                   [113/255, 188/255, 120/255],
                   [1, 102/255, 102/255],
                   [51/255, 153/255, 1],]

    if split_by_nucleotides and 'nucl' not in corr_df.keys():
        raise ValueError('DF must have nucleotide field to plot correlations across nucleotides.')

    if not package_list:
        package_list = corr_df['package'].unique()

    ind_val, unsorted_pkg_list = [], []
    for pkg in package_list:
        ind_val.append(packages.index(pkg))
        unsorted_pkg_list.append(pkg)

    sorted_pkg_list = [x for _,x in sorted(zip(ind_val, unsorted_pkg_list))]

    print(sorted_pkg_list)
    
    tmp_pal=sns.color_palette('Set1',5)
    standardized_palette = [tmp_pal[x] for x in [1,2,0,4,3]] + [[0,0,0]] 
    
    print(corr_df['package_type'].unique())
    if split_by_nucleotides:
        
        sns.pointplot(x='package', y='Corr', data=corr_df, order = sorted_pkg_list, hue='nucl',
                      join=False, ci='sd', marker='.', scale=0.5, dodge=0.5, palette=nuc_palette)

    else:
        sns.pointplot(x='package', y='Corr', data=corr_df, order = sorted_pkg_list,
                      hue='package_type', palette=standardized_palette,
                      join=False, ci='sd', marker='.', scale=0.5, dodge=False)

    plt.xticks([x for x in range(len(sorted_pkg_list))],
                                 [package_titles[packages.index(x)] for x in sorted_pkg_list], rotation=45, horizontalalignment='right')
    plt.grid()
    plt.legend(frameon=False)
    plt.ylabel('Correlation(Chem. Map. data, p(unpaired))')


    
