import pandas as pd
from arnie.mea.mea_utils import score_two_dbns
from arnie.mfe import mfe

def calc_RNAstructure_pks(row, pseudo=False):
    print(row['name'], row['dataset'])
    if row['modifier']=='DMS':
        return mfe(row['sequence'], package='rnastructure', dms_signal = row['reactivity'], pseudo=pseudo)
    elif row['modifier']=='SHAPE':
        return mfe(row['sequence'], package='rnastructure', shape_signal = row['reactivity'], pseudo=pseudo)

def score_prediction(row, package):
	sen, ppv, mcc, fscore, N = score_two_dbns(row[package], row['GT_struct'])
	return mcc

df = pd.read_csv('shape_dataset_input.txt',delimiter='\t')

df['SHAPEknots'] = df.apply(lambda row: calc_RNAstructure_pks(row, pseudo=True), axis=1)
df['RNAstructure+SHAPE/DMS'] = df.apply(lambda row: calc_RNAstructure_pks(row, pseudo=False), axis=1)
df['RNAstructure'] = df.apply(lambda row: mfe(row['sequence'], package='rnastructure'), axis=1)

df['ViennaRNA (Washietl)'] = df.apply(lambda row: mfe(row['sequence'], package='vienna', probing_signal = row['reactivity']), axis=1)
df['ViennaRNA (Deigan)'] = df.apply(lambda row: mfe(row['sequence'], probing_signal=row['reactivity'], package='vienna',shapeMethod='D'), axis=1)
df['ViennaRNA (Zaringhalam)'] = df.apply(lambda row: mfe(row['sequence'], probing_signal=row['reactivity'], package='vienna',shapeMethod='Z'), axis=1)
df['ViennaRNA (RNAfold)'] = df.apply(lambda row: mfe(row['sequence'], package='vienna'), axis=1)

df['EternaFold'] = df.apply(lambda row: mfe(row['sequence'], package='eternafold'), axis=1)
df['EternaFold+SHAPE/DMS'] = df.apply(lambda row: mfe(row['sequence'], package='eternafold',probing_signal = row['reactivity'], probing_kws={'kappa': 0.1}), axis=1)

# df.to_csv('shape_dataset_with_predicted_structures.txt',sep='\t')

#### Score structures ####

pkg_list = ['SHAPEknots','RNAstructure+SHAPE/DMS','RNAstructure','ViennaRNA (Washietl)',
'ViennaRNA (Deigan)','ViennaRNA (Zaringhalam)','ViennaRNA (RNAfold)','EternaFold','EternaFold+SHAPE/DMS']

for pkg in pkg_list:
	print('scoring', pkg)
	df['MCC %s' % pkg] = df.apply(lambda row: score_prediction(row, pkg), axis=1)

df.to_csv('shape_dataset_with_structures_and_scores.txt',sep='\t')
