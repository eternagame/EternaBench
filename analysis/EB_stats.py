import numpy as np
from scipy import stats
import pandas as pd
import matplotlib.pyplot as plt
from copy import copy
import textwrap

def bootstrap_inds(len_item):
    return np.random.choice(range(len_item), len_item)

def bootstrap_metric(df, x_data='logkd_nolig_scaled', y_data='log_pij', rmse=False, package_list=None,parametric=False, n_bootstraps=1):
    '''Input dataframe of experimental values and predicted balues. Bootstrap similarity metric.
    
    Returns:
    output_df: dictionary of package metrics over bootstrapped rounds.
    pval_dict: dictionary of probability that each package was the winner over all the bootstrapped rounds.
    pval_matrix: pairwise comparison matrix of p-values, computed by transforming correlation using fisher transform and 
        obtaining a z-score assuming normal distribution.
    empirical pval_matrix: pairwise comparison matrix of p-values, computed by counting number of times each package has a higher value 
        than the other.'''
    
    output_df = pd.DataFrame()
    pval_dict = {k:0 for k in package_list}
    dist_matrix = np.zeros([n_bootstraps, len(package_list), len(package_list)])
    empirical_pval_matrix = np.zeros([n_bootstraps, len(package_list), len(package_list)])
    
    xdata_vec = df[x_data].values
    
    for bs_ind in range(n_bootstraps):
        
        current_list = []
        for pkg in package_list:
            ydata_vec = df['%s_%s' % (y_data, pkg)].values

            bs_inds = bootstrap_inds(len(df))
            
            x = xdata_vec[bs_inds].astype(float)
            y = ydata_vec[bs_inds]
            
            if not rmse:
                val = np.corrcoef(x, y)[0][1]
            else:
                val = np.sqrt(np.mean(np.square(x-y)))

            current_list.append(val)
            if not rmse:
                output_df = output_df.append({'package': pkg, 'package_type': pkg.split('_')[0], 'C': val}, ignore_index=True)
            else:
                output_df = output_df.append({'package': pkg, 'package_type': pkg.split('_')[0], 'rmse': val}, ignore_index=True)

        if not rmse:
            winner = package_list[np.argmax(current_list)]
            x_C, y_C = np.meshgrid(np.tanh(current_list), np.tanh(current_list))
            
        else:
            winner = package_list[np.argmin(current_list)]
            x_C, y_C = np.meshgrid(current_list, current_list)
        

        dist_matrix[bs_ind] = (x_C - y_C)
        empirical_pval_matrix[bs_ind][x_C > y_C] = 1
        empirical_pval_matrix[bs_ind][x_C < y_C] = -1
        
        pval_dict[winner] +=1

    for k in package_list:
        pval_dict[k] /= n_bootstraps
        
    if parametric:
        z_score_matrix = np.mean(dist_matrix,axis=0)/np.std(dist_matrix,axis=0)
        pval_matrix = (stats.norm.cdf(z_score_matrix)-stats.norm.cdf(-1*z_score_matrix)) # 2-sided t test
    else:
        pval_matrix = np.mean(empirical_pval_matrix,axis=0)
            
    return output_df, pval_dict, pval_matrix

    
def plot_pairwise_matrices(dct, package_list, n_rows=4, n_cols=4, reorder=None,
                           empirical=False, titles=None, labels=True, dataset_titles=None, figsize=(10,10)):
    ctr=0
    n_packages = len(package_list)
    
    fig, ax = plt.subplots(nrows=n_rows, ncols=n_cols, subplot_kw={'aspect':'equal'},figsize=figsize)
    
    if dataset_titles is None:
        dataset_titles = [k for k in dct.keys()]
        
    if titles is None:
        titles=package_list
        
    if reorder is not None:
        titles= [titles[x] for x in reorder]

    matrices = [v for v in dct.values()]

    for i in range(n_rows):
        for j in range(n_cols):
            curr_ax = ax[i,j]

            if ctr < len(matrices):
                plot_arr = np.zeros([n_packages, n_packages])
                
                tmp = copy(matrices[ctr])
                if reorder is not None:
                    tmp = tmp[reorder]
                    tmp = tmp[:,reorder]
                    
                for ii in range(n_packages):
                    for jj in range(n_packages):
                        if ii>jj:
                            plot_arr[ii,jj] = tmp[ii,jj]

                heatmap = curr_ax.pcolor(plot_arr,cmap='seismic', vmin=-3, vmax=3)

                curr_ax.set_xticks(np.arange(n_packages)[:-1] + 0.5, minor=False)
                curr_ax.set_yticks(np.arange(n_packages)[1:]  + 0.5, minor=False)

                curr_ax.set_xticklabels(titles[:-1], minor=False, rotation=45, horizontalalignment='left')
                curr_ax.set_yticklabels(titles[1:], minor=False)
                curr_ax.xaxis.tick_top()
                   
                curr_ax.set_xlabel('\n'.join(textwrap.wrap(dataset_titles[ctr],20)))

                for typ in ['left','right','bottom','top']:
                    curr_ax.spines[typ].set_visible(False)

                ctr+=1
                
            else:
                fig.delaxes(curr_ax)           
    plt.tight_layout()
    plt.show()
    return fig
    
def bootstrap_all_datasets(df, x_data=None, y_data=None, package_list=None, n_bootstraps=100,rmse=False):
    
    dataset_list = df['dataset'].unique()
    n_datasets = len(dataset_list)
    n_packages = len(package_list)

    split_df = [df.loc[df['dataset']==x] for x in dataset_list]
    results = np.zeros([n_bootstraps, n_datasets, n_packages, n_packages])
    tot_list = []
    for bs_ind in range(n_bootstraps):
        # this gets indices to bootstrap the datasets.
        ds_bs_inds = bootstrap_inds(len(dataset_list))
        pval_matrix_list = []

        for j, ds_ind in enumerate(ds_bs_inds):
            # this bootstraps the individual dataset to compute the correlation coefficients.
            _, _, pval_matrix = bootstrap_metric(split_df[ds_ind], x_data = x_data, y_data = y_data, package_list=package_list, n_bootstraps=1, rmse=rmse)
            results[bs_ind, j] = pval_matrix
            
    return results


def sig_over_individual_datasets(df, x_data=None, y_data=None, package_list=None, n_bootstraps=100, rmse=False):
    
    dataset_list = df['dataset'].unique()
    n_datasets = len(dataset_list)
    n_packages = len(package_list)

    split_df = [df.loc[df['dataset']==x] for x in dataset_list]
    results = {k:[] for k in dataset_list}

    for j, ds_ind in enumerate(dataset_list):
        print(ds_ind)
        # this bootstraps the individual dataset to compute the correlation coefficients.
        _, _, pval_matrix = bootstrap_metric(split_df[j], x_data = x_data, y_data = y_data, package_list=package_list, n_bootstraps=n_bootstraps, rmse=rmse)
        results[ds_ind] = pval_matrix

    return results

def plot_net_overview(arr, package_list=None, titles=None, reorder=None, figsize=(4,4),fontsize=16, rmse=False):
    
    n_datasets=arr.shape[1]
    n_bs_inds=arr.shape[0]
    n_packages = arr.shape[2]
    tmp2 = copy(arr)
    tmp2[np.where(tmp2<0)] = 0
    tmp=np.sum(tmp2,axis=1)
    counter = np.zeros([n_bs_inds, n_packages, n_packages])
    win_inds = np.where(tmp > n_datasets/2)
    counter[win_inds] = 1
    mean_arr, std_arr = np.mean(counter,axis=0), np.std(counter,axis=0)

    fig, curr_ax = plt.subplots(subplot_kw={'aspect':'equal'},figsize=figsize)
    
    if titles is None:
        titles=package_list
        
    if reorder is not None:
        mean_arr = mean_arr[reorder]
        mean_arr = mean_arr[:,reorder]
        titles= [titles[x] for x in reorder]
                

    plot_arr = 0.5*np.ones([n_packages, n_packages])

    for ii in range(n_packages):
        for jj in range(n_packages):
            if ii>jj:
                plot_arr[ii,jj] = copy(mean_arr[ii,jj])

    if rmse:
        heatmap = curr_ax.pcolor(plot_arr,cmap='seismic_r', vmin=-1.5, vmax=2.5)
    else:
        heatmap = curr_ax.pcolor(plot_arr,cmap='seismic', vmin=-1.5, vmax=2.5)

    curr_ax.set_xticks(np.arange(n_packages)[:-1] + 0.5, minor=False)
    curr_ax.set_yticks(np.arange(n_packages)[1:] + 0.5, minor=False)

    curr_ax.set_xticklabels(titles[:-1], minor=False, rotation=45, horizontalalignment='left',fontsize=fontsize)
    curr_ax.set_yticklabels(titles[1:],  minor=False, fontsize=fontsize)
    curr_ax.xaxis.tick_top()
       
    for typ in ['left','right','bottom','top']:
        curr_ax.spines[typ].set_visible(False)

    for i in range(n_packages):
        for j in range(n_packages):
            
            if i<j:

                if rmse:
                    curr_ax.text(i+0.2,j+0.4,"%.3f" % (1-mean_arr[i,j]),fontsize=12)
                else:
                    curr_ax.text(i+0.2,j+0.4,"%.3f" % (mean_arr[i,j]),fontsize=12)
                    
                
    plt.tight_layout()
    plt.show()
    return fig