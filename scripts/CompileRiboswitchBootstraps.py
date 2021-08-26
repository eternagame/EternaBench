from glob import glob
import pandas as pd
import os, sys
from eternabench.stats import calculate_Z_scores

for metric in ['pearson','spearman','RMSE']:

    bp_bootstraps = glob(os.environ['ETERNABENCH_PATH']+'/data/Riboswitch/RS_%s_*_bps_BOOTSTRAPS.json.zip' % metric)
    Z_bootstraps = glob(os.environ['ETERNABENCH_PATH']+'/data/Riboswitch/RS_%s_*_Z_BOOTSTRAPS.json.zip' % metric)

    #############################
    #### Load bps bootstraps ####
    #############################

    bootstraps = pd.DataFrame()

    print('%d bps files found' % len(bp_bootstraps))
    for x in sorted(bp_bootstraps):
        print(x)

    print('')
    print('%d Z files found' % len(Z_bootstraps))
    for x in sorted(Z_bootstraps):
        print(x)

    for fil in bp_bootstraps:
        tmp = pd.read_json(fil)
        bootstraps = bootstraps.append(tmp, ignore_index=True)
    #bootstraps.to_json(os.environ['ETERNABENCH_PATH']+'/scoring_data/Riboswitch_bps_%s_bootstraps.json.zip' % metric)

    ##############################################
    #### Calculate Zscores without eternafold ####
    ##############################################

    zscore_stats, ranking = calculate_Z_scores(bootstraps, metric=metric, dataset_field='Dataset', include_efold=False)

    zscore_stats.to_csv('RS_bps_%s_zscores.csv' % metric,index=False)
    ranking.to_csv('RS_bps_%s_rank.csv' % metric,index=False)

    #########################################################
    #### Calculate Zscores WITH eternafold over test set ####
    #########################################################

    efold_test_set_package_list = []

    with open(os.environ['ETERNABENCH_PATH']+'/scripts/package_benchmark_+eternafold.txt','r') as f:
        for lin in f.readlines():
            efold_test_set_package_list.append(lin.strip())

    print("Former length", len(bootstraps))
    bootstraps = bootstraps.loc[~bootstraps.Dataset.str.contains('Ribologic')]
    print("New length", len(bootstraps))
    zscore_stats, ranking = calculate_Z_scores(bootstraps, metric=metric, dataset_field='Dataset', package_list = efold_test_set_package_list)

    zscore_stats.to_csv('RS_bps_%s_zscores_Fig3_efold_testset.csv' % metric,index=False)
    ranking.to_csv('RS_bps_%s_rank_Fig3_efold_testset.csv' % metric,index=False)


    ##############################################
    #### Calculate Zscores with efold versions ###
    ##############################################

    #(note this is still without Ribologic)
    zscore_stats, ranking = calculate_Z_scores(bootstraps, metric=metric, dataset_field='Dataset', include_efold=True)

    zscore_stats.to_csv('RS_bps_%s_zscores_all_efold_versions.csv' % metric,index=False)
    ranking.to_csv('RS_bps_%s_rank_all_efold_versions.csv' % metric,index=False)


    ###########################
    #### Load Z bootstraps ####
    ###########################

    bootstraps = pd.DataFrame()

    for fil in Z_bootstraps:
        tmp = pd.read_json(fil)
        bootstraps = bootstraps.append(tmp, ignore_index=True)

    #bootstraps.to_json(os.environ['ETERNABENCH_PATH']+'/scoring_data/Riboswitch_Z_Bootstraps.json.zip')

    ##############################################
    #### Calculate Zscores without eternafold ####
    ##############################################

    for calc in ['logkd_nolig_scaled', 'logkd_lig_scaled', 'log_AR']:
        tmp = bootstraps.loc[bootstraps.Calculation==calc]
        zscore_stats, ranking = calculate_Z_scores(tmp, metric=metric, dataset_field='Dataset',include_efold=False)

        zscore_stats.to_csv('RS_Z_%s_%s_zscores.csv' % (metric, calc),index=False)
        ranking.to_csv('RS_Z_%s_%s_rank.csv' % (metric, calc),index=False)
