import pandas as pd
import os, sys
from eternabench.stats import calculate_Z_scores


package_list=['vienna_2', 'vienna_2_60C', 'rnastructure', 'rnastructure_60C', 'rnasoft_blstar','contrafold_2','eternafold_B']

external_dataset_types = pd.read_csv(os.environ['ETERNABENCH_PATH']+'/eternabench/external_dataset_metadata.csv')

RNA_CLASSES = list(external_dataset_types.Class.unique())

EB_CM_bootstraps=pd.DataFrame()

for pkg in package_list:
    tmp = pd.read_json(os.environ['ETERNABENCH_PATH']+'/data/ChemMapping/bootstraps/CM_pearson_Dataset_%s_BOOTSTRAPS.json.zip' % pkg)
    EB_CM_bootstraps = EB_CM_bootstraps.append(tmp, ignore_index=True)

EB_CM_bootstraps = EB_CM_bootstraps.loc[EB_CM_bootstraps.Dataset=='RYOS_I']
EB_CM_bootstraps['Dataset'] = 'Leppek,2021 In-line-seq'

net_dataset_zscore_stats=pd.DataFrame()
net_ranking = pd.DataFrame()

for window_size in [300,600,900,1200]:

    df = pd.DataFrame()

    for pkg in package_list:
        #print('reading in file %s '% fil)
        tmp = pd.read_json("%s/data/Ext%d/CM_pearson_Dataset_%s_BOOTSTRAPS.json.zip" % (os.environ['ETERNABENCH_PATH'], window_size, pkg))
        df = df.append(tmp, ignore_index=True)

    df = df.append(EB_CM_bootstraps, ignore_index=True)
    df = df.merge(external_dataset_types, on='Dataset', how='left')

    for rna_class in RNA_CLASSES:

        print(window_size, rna_class)

        tmp = df.loc[df.Class==rna_class]

        if len(tmp)>0:

            zscore_stats, ranking = calculate_Z_scores(tmp, metric='pearson', dataset_field='Dataset',sort=False) #omitting package list to not throw error when no nupack for w=1200

            zscore_stats['window_size'] = window_size
            ranking['window_size']=window_size

            zscore_stats['Class'] = rna_class
            ranking['Class'] = rna_class

            net_dataset_zscore_stats= net_dataset_zscore_stats.append(zscore_stats, ignore_index=True)
            net_ranking= net_ranking.append(ranking, ignore_index=True)

net_dataset_zscore_stats.to_csv('EternaBench_external_zscores_by_dataset.csv',index=False)
net_ranking.to_csv('EternaBench_external_ranking.csv',index=False)


