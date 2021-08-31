import sys, os, argparse
import eternabench.chemmapping_utils as utils
from eternabench.stats import calculate_metric
import pandas as pd
import numpy as np
from scipy.stats import pearsonr, spearmanr

def ScoreRiboswitches(data, xdata, ydata, agg_field = 'Dataset', n_bootstraps=10, package_list=None, metric='pearson'):
	'''Usage:
	Calculate metric for specified packages for data, split up by experimental datasets.

	This script concatenates data across all nucleotides. 

	Input:
	data: dataframe containing one row per construct
	agg_field: type to aggregate per-nucleotide correlations over. Can be 'filename' or 'project_name'.
	package_list: list of packages to calculate correlations for. If not provided, calculates over all packages (found by 'p_')
	n_bootstraps: number of bootstraps to use for correlation calculation.
	metric: type of metric metric to compute. default is 'pearson'. Options are 'pearson', 'spearman'.

	Dataframe required fields:
	'reac': vector of reactivity data
	'filename': specifies the experiment the construct is in
	'p_*': for each row, vector of probabilities same length as 'reactivity' field

	Output:
	Dataframe of correlation values 
	'''

	correlation_data = pd.DataFrame()

	if package_list is None:
		package_list=[k.replace(ydata+'_' ,'') for k in data.keys() if k.startswith(ydata+'_')]

	print('Package list:', package_list)

	# process data
	for kind in data[agg_field].unique():
		print(kind)

		tmp_data = data.loc[data[agg_field]==kind]
		
		corr_data = calculate_metric(tmp_data, x_data=xdata, y_data=ydata, n_bootstraps=n_bootstraps,
			package_list=package_list, metric=metric)

		corr_data[agg_field] = kind
		corr_data['Calculation'] = xdata
		correlation_data = correlation_data.append(corr_data,ignore_index=True)

	return correlation_data

if __name__=='__main__':

	p = argparse.ArgumentParser(description="""Score written predictions""")

	p.add_argument("infile", action="store", help="Input json")
	p.add_argument("--test", action='store_true', help='Tests first 3 constructs in each package.')
	p.add_argument("--metadata", action='store', help='Dataframe in .json format. If provided, will use metadata from this df\
		instead of the input dataframe (For instance, if there is an error, or if there are duplicates of sequences for\
		different conditions.')
	p.add_argument("--field_to_aggregate", action='store', default='Dataset', help="Field to aggregate results over. In chem mapping, 'filename' for experiments, or 'project_name' for projects.")
	p.add_argument("--n_bootstraps", action='store',type=int, default=10)
	p.add_argument("--metric", action='store',default='spearman', help='spearman, pearson, or rmse')
	p.add_argument("--method", action='store', default='bps', help="Type of riboswitch calculation to store (that is present in dataset). Must be `Z` or `bps`.")
	p.add_argument("-v", "--verbose", action="store_true", help="Verbose")
	p.add_argument("-o", action="store", dest='outfile', help='name of output json file (default is <input_name>_BOOTSTRAPS.json.zip')

	args = p.parse_args()

	basename = args.infile.split('/')[-1].split('.')[0]

	if args.outfile is None:
		outfile = basename
	else:
		outfile = args.outfile

	df = pd.read_json(sys.argv[1])

	if args.metadata:
		print('Current df length: %d' % len(df))
		metadata_df = pd.read_json(args.metadata)

		keys_to_add = [k for k in df.keys() if 'p_' in k]+['sequence']
		metadata_df = metadata_df.merge(df[keys_to_add], on='sequence')

		df = metadata_df
		print('Using metadata from %s, new df length = %d' % (args.metadata, len(df)))

	df = df.dropna(subset=[args.field_to_aggregate])

	print('read in', basename)

	if args.method=='bps':
		x_inputs=['logkd_nolig_scaled']
		y_inputs=['log_kfold_est_bp']

	elif args.method=='Z':
		x_inputs=['logkd_nolig_scaled', 'logkd_lig_scaled','log_AR']
		y_inputs=['log_kfold_est_nolig_Z', 'log_kfold_est_lig_Z', 'log_AR_est']

	out = pd.DataFrame()
	for x_input, y_input in list(zip(x_inputs, y_inputs)):
		correlation_data = ScoreRiboswitches(df, x_input, y_input, agg_field = args.field_to_aggregate,\
			n_bootstraps=args.n_bootstraps, package_list=None, metric=args.metric)

		out = out.append(correlation_data, ignore_index=True)

	out.to_json(outfile+'_BOOTSTRAPS.json.zip')
