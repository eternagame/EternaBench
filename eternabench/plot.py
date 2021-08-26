import numpy as np
import pandas as pd

import seaborn as sns
import matplotlib.pyplot as plt
from copy import deepcopy
from glob import glob
import sys, os
from copy import deepcopy
from scipy.stats import pearsonr,spearmanr
import warnings
warnings.filterwarnings("ignore")


blue, orange, green, red, purple, brown, _, _, _, _=sns.color_palette()
palette = [brown, blue, red, green,orange, purple, [0,0,0]]

def get_packages():
    tmp = pd.read_csv(os.environ['ETERNABENCH_PATH']+'/eternabench/package_metadata.csv')
    tmp = tmp.set_index('package')
    return tmp

def get_external_Dataset_data():
    return pd.read_csv(os.environ['ETERNABENCH_PATH']+'/eternabench/external_dataset_metadata.csv')

def get_palette():
    return palette, ['Other','ViennaRNA','NUPACK','RNAstructure','CONTRAfold','RNAsoft', 'EternaFold']

def reactivity_heatmap(df, ind_range=None, **kwargs):
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
    plt.imshow(np.array(arr), cmap='gist_heat_r',vmin=0,vmax=1.5, **kwargs)
    plt.ylabel('Construct')
    plt.xlabel('Sequence position')

def punpaired_heatmap(df, ind_range=None, package='vienna_2', **kwargs):
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
    plt.imshow(np.array(arr), cmap='gist_heat_r',vmin=0,vmax=1, **kwargs)
    plt.xlabel('Construct')
    plt.ylabel('Sequence position')

def create_array(df_o, col1, col2, value, col1_subset=None, col2_subset=None):
    '''Create array of `value` from dataframe. 
    First index is col1, subset, second index is col2 '''
    
    if col1_subset is None:
        col1_subset = list(df_o[col1].unique())
    if col2_subset is None:
        col2_subset = list(df_o[col2].unique())
        

    arr = np.nan*np.ones([len(col1_subset),len(col2_subset)])
    for i, col1_val in enumerate(col1_subset):
        for j, col2_val in enumerate(col2_subset):
            tmp = df_o.loc[df_o[col1]==col1_val][df_o[col2] == col2_val]
            if len(tmp)>0:
                assert len(tmp)==1
                arr[i,j] = tmp[value].values[0]
                    
    return arr, col1_subset, col2_subset
        
def single_barplot(df, cat, val, err, titles, **kwargs):
    u = df[cat].unique()
    x = np.arange(len(u))
    offsets = (np.arange(len(u))-np.arange(len(u)).mean())/(len(u)+1.)
    plt.barh(x,df[val].values,  xerr=df[err].values, **kwargs)
    plt.ylim([-0.5,len(x)-0.5])

def ranked_heatmap(zscores, ranking, metric='pearson_zscore_by_Dataset_mean', package_order=None,
                   dataset_field='Dataset', 
                   dataset_order=None,
                   figsize=(7,5), width_ratios=[2,1], cbar_loc=[-0.1,0.05], fig=None, axes=None,
                   vmin=None,vmax=None, ext=False,dataset_labels=None):
    '''Plot heatmap of packages ranked over datasets, with best at top.'''

    if package_order is not None:
        zsc_copy=pd.DataFrame()
        r_copy = pd.DataFrame()
        for pkg in package_order:
            zsc_copy= zsc_copy.append(zscores.loc[zscores.package==pkg],ignore_index=True)
            r_copy= r_copy.append(ranking.loc[ranking.package==pkg],ignore_index=True)

        zscores = deepcopy(zsc_copy)
        ranking = deepcopy(r_copy)

    if dataset_order is not None:
        zsc_copy=pd.DataFrame()
        for ds in dataset_order:
            zsc_copy= zsc_copy.append(zscores.loc[zscores[dataset_field]==ds],ignore_index=True)

        zscores = deepcopy(zsc_copy)

    tmp, package_list, dataset_list = create_array(zscores, 'package', dataset_field, metric)

    if axes is None:
        gridkw = dict(width_ratios=width_ratios)
        fig, (ax1, ax2) = plt.subplots(1,2, gridspec_kw=gridkw, figsize=figsize)
    else:
        ax1, ax2 = axes
    
    package_data = get_packages()
    palette, hue_order = get_palette()

    package_titles=[]
    for pkg in package_list:
        try:
            title = package_data.loc[pkg]['title']
            package_titles.append(title)
        except:
            package_titles.append(pkg)

    color_range = np.max(np.abs(tmp))

    if vmin is None:
        vmin=-1*color_range
    if vmax is None:
        vmax = color_range

    im = ax1.imshow(tmp, cmap='seismic_r', vmin = vmin, vmax=vmax,origin='upper')
    ax1.set_yticks(range(len(package_titles)))
    ax1.set_yticklabels(package_titles)
    
    ax1.set_xticks(range(len(dataset_list)))

    if ext:
        ext_data = get_external_Dataset_data()
        ax1.set_xticklabels([ext_data.loc[ext_data['Dataset']==k].Title.values[0] for k in dataset_list],rotation=45,horizontalalignment='right')
    else:
        if dataset_labels is not None:
            ax1.set_xticklabels(dataset_labels, rotation=45, horizontalalignment='right')
        else:
            ax1.set_xticklabels(dataset_list,rotation=45,horizontalalignment='right')

    cbaxes = fig.add_axes([*cbar_loc, 0.15, 0.03]) 
    fig.colorbar(im, cax = cbaxes,orientation='horizontal', label='Z-score')

    ranking = ranking.merge(package_data, on='package')
    zscores = zscores.merge(package_data, on='package')

    ranking['category'] = [hue_order[i] for i in ranking['category']]
    zscores['category'] = [hue_order[i] for i in zscores['category']]

    sns.barplot(y='title', x=metric, hue='category', dodge=False, data=zscores,
                ax=ax2, palette=palette, hue_order=hue_order, linewidth=0)
    ax2.yaxis.set_ticks_position('right')
    ax2.set_ylabel('')
    ax2.axvline(0,color='k',linewidth=0.5,linestyle=':')
    ax2.set_xlabel('Avg. Z-score')
    ax2.set_xlim([vmin,vmax])
    ax2.legend([],[], frameon=False)


def ranked_heatmap_w_bar_overhead(zscores, ranking, metric='pearson_zscore_by_Dataset_mean', package_order=None,
                   dataset_field='Dataset', figsize=(7,5), width_ratios=[2,1], cbar_loc=[-0.1,0.05],
                   vmin=None,vmax=None, ext=False,barplot_ymax=0.75, barplot_ylabel='Mean Corr.',RMSE=False, rainbow=False):
    '''Plot heatmap of packages ranked over datasets, with best at top.'''

    if RMSE:
        package_order = list(reversed(package_order))
        zscores[metric] *= -1
        ranking[metric] *= -1
        
    if package_order is not None:
        zsc_copy=pd.DataFrame()
        r_copy = pd.DataFrame()
        for pkg in package_order:
            zsc_copy= zsc_copy.append(zscores.loc[zscores.package==pkg],ignore_index=True)
            r_copy= r_copy.append(ranking.loc[ranking.package==pkg],ignore_index=True)

        zscores = deepcopy(zsc_copy)
        ranking = deepcopy(r_copy)



    tmp, package_list, dataset_list = create_array(zscores, 'package', dataset_field, metric)

    gridkw = dict(width_ratios=width_ratios, height_ratios = [.5,3])
    fig, ax = plt.subplots(2,2, gridspec_kw=gridkw, figsize=figsize)

    package_data = get_packages()
    palette, hue_order = get_palette()

    package_titles=[]
    for pkg in package_list:
        title = package_data.loc[pkg]['title']
        if len(title)>0:
            package_titles.append(title)
        else:
            package_titles.append(pkg)

    color_range = np.max(np.abs(tmp))

    if vmin is None:
        vmin=-1*color_range
    if vmax is None:
        vmax = color_range

    ax0=ax[0][0]

    if rainbow:
        sns.barplot(x=dataset_field, y=metric.split('_')[0]+'_mean', data=zscores, palette='rainbow', ax=ax0)
    else:
        sns.barplot(x=dataset_field, y=metric.split('_')[0]+'_mean', data=zscores, color='grey', ax=ax0)

    ax0.set_xticks([])
    ax0.set_ylim([0,barplot_ymax])
    ax0.set_xlabel('')
    ax0.set_ylabel(barplot_ylabel)

    ax0.spines['right'].set_visible(False)
    ax0.spines['top'].set_visible(False)


    ax[0][1].axis('off')

    ax1=ax[1][0]

    im = ax1.imshow(tmp, cmap='seismic_r', vmin = vmin, vmax=vmax, aspect='auto', origin='upper left')

    ax1.set_yticks(range(len(package_titles)))
    ax1.set_yticklabels(package_titles)
    
    ax1.set_xticks(range(len(dataset_list)))

    if ext:
        ext_data = get_external_Dataset_data()
        ax1.set_xticklabels([ext_data.loc[ext_data['Dataset']==k].Title.values[0] for k in dataset_list],rotation=45,horizontalalignment='right')
    else:
        ax1.set_xticklabels(dataset_list,rotation=45,horizontalalignment='right')

    cbaxes = fig.add_axes([*cbar_loc, 0.15, 0.03]) 
    fig.colorbar(im, cax = cbaxes,orientation='horizontal', label='Z-score')

    ranking = ranking.merge(package_data, on='package')

    ranking['category'] = [hue_order[i] for i in ranking['category']]

    ax2=ax[1][1]
    # sns.barplot(y='title', x=metric, hue='category', dodge=False, data=ranking.iloc[::-1],
    #             ax=ax2, palette=palette, hue_order=hue_order, )

    color_order = [palette[hue_order.index(x)] for x in ranking ['category']]

    ax2.barh(np.arange(len(ranking)), ranking[metric].values,  xerr=ranking[metric.replace('mean','std')].values,color=color_order)
    ax2.set_ylim([-0.5,len(ranking)-0.5])

    ax2.yaxis.set_ticks_position('right')
    ax2.set_ylabel('')
    ax2.axvline(0,color='k',linewidth=0.5,linestyle=':')
    ax2.set_xlabel('Avg. Z-score')
    ax2.set_xlim([vmin-0.1,vmax+0.1])
    ax2.legend([],[], frameon=False)

def corrfunc(x,y, ax=None, **kws):

    r, pval = spearmanr(x, y, nan_policy='omit')

    ax = ax or plt.gca()    
    # m, b = np.poly1d(np.polyfit(x, y, 1))
    # xmin, xmax = ax.get_xlim()
    # plt.plot([xmin,xmax],[xmin*m+b,xmax*m+b],c='k',linestyle=':')
    # ax.set_xlim([xmin,xmax])
    ax.set_title(f'Spearman R = {r:.2f}',loc='right')

def jitterbox(**kwargs):
    'supply x, y, hue, data'
    sns.stripplot(**kwargs, dodge=True, alpha=0.5)
    ax = sns.boxplot(**kwargs, dodge=True, fliersize=0, zorder=10)
    
    for patch in ax.artists:
        r, g, b, a = patch.get_facecolor()
        patch.set_facecolor((r, g, b, 0))
        patch.set_zorder(10)

    handles, labels = ax.get_legend_handles_labels()
    l = plt.legend(handles[0:2], labels[0:2], bbox_to_anchor=(1, 1), frameon=False)
