import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

def percentile(n):
    def percentile_(x):
        return np.percentile(x, n)
    percentile_.__name__ = 'percentile_%s' % n
    return percentile_

def bootstrap_inds(len_item):
    return np.random.choice(range(len_item), len_item)
    
def calculate_metric(df, x_data='reactivity', y_data='p', split_by_nucleotides=False, package_list=None, n_bootstraps=10, metric='pearson'):
    '''Input: concatenated dataframe or riboswitch dataframe.

    Calculate '''
    
    if metric not in ['spearman', 'pearson','RMSE']:
        raise RuntimeError('metric %s not understood.' % metric)

    corr_df=pd.DataFrame()
    xdata_vec = df[x_data].values

    for bs_ind in range(n_bootstraps):
        current_list=[]
        bs_inds = bootstrap_inds(len(xdata_vec))

        for k in package_list:
            ydata_vec = df['%s_%s' % (y_data, k)].values

            x = xdata_vec[bs_inds].astype(float)
            y = ydata_vec[bs_inds]

            #cat_p = df['p_%s' % k].values

            if metric=='spearman':
                C, _ = spearmanr(x, y)
            elif metric=='pearson':
                C, _ = pearsonr(x,y)
            elif metric=='RMSE':
                C = np.sqrt(np.mean(np.square(x-y)))

            current_list.append(C)
            corr_df = corr_df.append({metric: C, 'package': k,'bs_ind': bs_ind}, ignore_index=True)

    return corr_df
    
def calculate_Z_scores(df, package_list=None, dataset_list = None, metric='pearson', dataset_field='Dataset', ranking_category=None, sort=True, include_efold=True):
    '''Input:

    Dataframe with fields for `bootstrap`, `package`, correlation metric (like `pearson`), and `dataset_field`.
    sort (bool): sort the returned dataframes by ranking.

    Output
    stats_by_dataset: mean metric and z-score by dataset.
    ranking: mean metric and z-score over all datasets.

    '''

    if package_list is not None:
        df = df.loc[df.package.isin(package_list)]

    if not include_efold:
        df = df.loc[~df.package.str.contains('eternafold')]

    if dataset_list is not None:
        df = df.loc[df[dataset_field].isin(dataset_list)]

    df[metric+'_std_across_packages'] = df.groupby([dataset_field,'bs_ind'])[metric].transform('std')
    df[metric+'_mean_across_packages'] = df.groupby([dataset_field,'bs_ind'])[metric].transform('mean')
    df = df.loc[df[metric+'_std_across_packages']!=0]

    # zscore is calculated per dataset
    df[metric+'_zscore_by_%s' % dataset_field] = df.apply(lambda row: (row[metric] - row[metric+'_mean_across_packages'])/row[metric+'_std_across_packages'], axis=1)
    df = df.sample(frac=1)

    stats_by_dataset = df.groupby(['package',dataset_field])[[metric,metric+'_zscore_by_%s' % dataset_field]].agg([np.mean, np.std, percentile(2.5), percentile(97.5)])
    stats_by_dataset = stats_by_dataset.reset_index()

    stats_by_dataset.columns = ['_'.join(col) if len(col[-1])>0 else col[0] for col in stats_by_dataset.columns]

    ranking = df.groupby('package')[[metric, metric+'_zscore_by_%s' % dataset_field]].agg([np.mean, np.std, percentile(2.5), percentile(97.5)])

    ranking.columns = ['_'.join(col) if len(col[-1])>0 else col[0] for col in ranking.columns]
    ranking = ranking.reset_index()

    if sort: 
        ranking = ranking.sort_values( metric+'_zscore_by_%s_mean' % dataset_field)
        pkg_sort = list(reversed(list(ranking['package']))) # puts best packages at top of dataset / plot

        sorted_zscore_stats = pd.DataFrame()

        for pkg in pkg_sort:
            sorted_zscore_stats = sorted_zscore_stats.append(stats_by_dataset.loc[stats_by_dataset.package==pkg])

        return sorted_zscore_stats, ranking

    else:
        return stats_by_dataset, ranking
