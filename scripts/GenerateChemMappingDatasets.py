from glob import glob
from eternabench.chemmapping_utils import *
from datetime import date

todaysdate = date.today().strftime("%d%b%Y")


rdat_files = {
'Round 00': 'ETERNA_R00_0000.rdat',
'Round 01': 'ETERNA_R69_0000.rdat',
'Round 02': 'ETERNA_R70_0000.rdat',
'Round 03': 'ETERNA_R71_0000.rdat',
'Round 04': 'ETERNA_R72_0000.rdat',
'Round 05': 'ETERNA_R73_0000.rdat',
'Round 06': 'ETERNA_R74_0000.rdat',
'Round 07': 'ETERNA_R75_0000.rdat',
'Round 08': 'ETERNA_R76_0000.rdat',
'Round 09': 'ETERNA_R77_0002.rdat',
'Round 10': 'ETERNA_R78_0001.rdat',
'Round 11': 'ETERNA_R79_0001.rdat',
'Round 12': 'ETERNA_R80_0001.rdat',
'Round 13': 'ETERNA_R81_0001.rdat',
'Round 14': 'ETERNA_R82_0001.rdat',
'Round 15': 'ETERNA_R83_0003.rdat',
'Round 16': 'ETERNA_R84_0000.rdat',
'Round 17': 'ETERNA_R85_0000.rdat',
'Round 18': 'ETERNA_R86_0000.rdat',
'Round 19': 'ETERNA_R87_0001.rdat',
'Round 20': 'ETERNA_R89_0000.rdat',
'Round 21': 'ETERNA_R91_0000.rdat',
'Round 22': 'ETERNA_R92_0000.rdat',
'Round 23': 'ETERNA_R94_0000.rdat',
'RYOS_I':'Ann_EteRNA_Lib062020_RYOS_Analysis3_NextSeq_EtOhPrecip.2.rdat'
}

full_df = pd.DataFrame()
filtered_df = pd.DataFrame()

def write_EternaScore(row):
	try:
		return float(row['EteRNA'].split(':')[-1])
	except:
		return row['EteRNA']

if not os.path.exists(os.environ['ETERNABENCH_PATH']+'/data/ChemMappingPreprocessing/raw_rdats/'):
	raise RuntimeError("Could not find /EternaBench/data/ChemMappingPreprocessing/raw_rdats/ . Navigate to that directory and unzip `raw_rdats.zip`")

for name, fil in rdat_files.items():
	df = write_rdat_as_dataframe(os.environ['ETERNABENCH_PATH']+'/data/ChemMappingPreprocessing/raw_rdats/%s' % fil, verbose=True)
	# remove FMN data
	df = filter_FMN_containing(df)

	# remove DMS data (only for R70)
	df = df.loc[df['modifier']!='DMS']
	print(name, len(df))
	print(df.head())

	df['Dataset'] = name

	df = df.dropna(subset=['sequence'])
	df['sequence'] = [x.replace('.','').replace('\n','') for x in df['sequence']]

	df = filter_dataframe_with_cdhit(df, fil.replace('.rdat',''))
	full_df = full_df.append(df,ignore_index=True)

full_df = full_df.loc[~full_df.Dataset.isna()] # Unsure why some get set to NaNs
full_df['EternaScore'] = full_df.apply(lambda row: write_EternaScore(row), axis=1)
full_df = full_df.drop(columns=['Part 2', 'EteRNA'])

filtered_df = full_df.loc[full_df['passed_CDHIT_filter']==True]
full_df.to_json(os.environ['ETERNABENCH_PATH']+'/data/EternaBench_ChemMapping_Full_%s.json.zip' % todaysdate)
filtered_df.to_json(os.environ['ETERNABENCH_PATH']+'/data/EternaBench_ChemMapping_Filtered_%s.json.zip' % todaysdate)
