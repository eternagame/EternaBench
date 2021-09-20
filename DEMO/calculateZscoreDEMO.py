from eternabench.stats import calculate_Z_scores

import pandas as pd


package_list = [x.strip() for x in open('package_list.txt','r').readlines()]

CM_files = ['CM_pearson_Dataset_%s_BOOTSTRAPS.json.zip' % x for x in package_list]
RS_files = ['RS_pearson_Dataset_%s_bps_BOOTSTRAPS.json.zip' % x for x in package_list]

CM = pd.DataFrame()
for fil in CM_files:
	CM = CM.append(pd.read_json(fil), ignore_index=True)
print('Chem Mapping Rnd 1 scores')
cm_zscores, _ = calculate_Z_scores(CM)
print(cm_zscores[['package','pearson_mean','pearson_std','pearson_zscore_by_Dataset_mean']])

RS = pd.DataFrame()
for fil in RS_files:
	RS = RS.append(pd.read_json(fil), ignore_index=True)

print('')
print('Riboswitch "Ribologic FMN" scores')

rs_zscores, _ = calculate_Z_scores(RS)
print(rs_zscores[['package','pearson_mean','pearson_std','pearson_zscore_by_Dataset_mean']])
