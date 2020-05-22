from glob import glob
from RDAT_utils import *

rdat_files = [
'ETERNA_R69_0000.rdat',
'ETERNA_R70_0000.rdat',
'ETERNA_R71_0000.rdat',
'ETERNA_R72_0000.rdat',
'ETERNA_R73_0000.rdat',
'ETERNA_R74_0000.rdat',
'ETERNA_R75_0000.rdat',
'ETERNA_R76_0000.rdat',
'ETERNA_R77_0002.rdat',
'ETERNA_R78_0001.rdat',
'ETERNA_R79_0001.rdat',
'ETERNA_R80_0001.rdat',
'ETERNA_R81_0001.rdat',
'ETERNA_R82_0001.rdat',
'ETERNA_R83_0003.rdat',
'ETERNA_R84_0000.rdat',
'ETERNA_R85_0000.rdat',
'ETERNA_R86_0000.rdat',
'ETERNA_R87_0001.rdat',
'ETERNA_R89_0000.rdat',
'ETERNA_R91_0000.rdat',
'ETERNA_R92_0000.rdat',
'ETERNA_R94_0000.rdat',
]

full_df = pd.DataFrame()
filtered_df = pd.DataFrame()

for fil in rdat_files:
	df = write_rdat_as_dataframe(fil, verbose=True)
	# remove FMN data
	df = filter_FMN_containing(df)

	# remove DMS data (only for R70)
	df = df.loc[df['modifier']!='DMS']
	print(len(df))

	df['filename'] = [os.path.basename(x) for x in df['filename']]

	full_df = full_df.append(df,ignore_index=True)
	cdhit_df = filter_dataframe_with_cdhit(df, fil.replace('.rdat',''))
	filtered_df = filtered_df.append(cdhit_df, ignore_index=True)

full_df.to_json('Mar23_unfiltered_cloud_lab_rounds.json')
filtered_df.to_json('Mar23_filtered_cloud_lab_rounds.json')





