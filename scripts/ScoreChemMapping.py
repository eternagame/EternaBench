import sys, os, argparse
import eternabench.chemmapping_utils as utils
from eternabench.stats import calculate_metric
import pandas as pd
import numpy as np
from scipy.stats import pearsonr, spearmanr

def ScoreChemMapping(data, x_data='reactivity', y_data='p', agg_field = 'filename', n_bootstraps=10,
	package_list=None, winsorize_cutoff=95, metric='pearson'):
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
	pval_data=pd.DataFrame()

	if package_list is None:
		package_list=[k.replace('p_','') for k in data.keys() if k.startswith('p_')]

	print('Package list:', package_list)

	# process data
	for kind in data[agg_field].unique():

		tmp_data = data.loc[data[agg_field]==kind]
		tmp_data = tmp_data.loc[~tmp_data[x_data].isna()]
		if len(tmp_data) == 0:
			break
		tmp_concat_data = utils.write_concatenated_dataframe(tmp_data, reactivity_field = x_data)
		# remove NaNs                                                                                                                            
		nan_subset = [x_data] + [x for x in tmp_concat_data.keys() if x.startswith('p_')]                                                 
		tmp_concat_data = tmp_concat_data.dropna(subset=nan_subset)                                                                              
		#print('HI2', kind, len(tmp_concat_data))                                                                                                
		tmp_concat_data = tmp_concat_data.fillna(value=np.nan) # IDK why I had to change Nones back to nans                                      
                               
		if len(tmp_concat_data) > 10: # to prevent weird things
			print('Analyzing %s' % kind)

			# filter reactivity outliers
			try:
				tmp_concat_data = utils.filter_data(tmp_concat_data, winsorize_cutoff=winsorize_cutoff)
			except:
				pass
				# filtering doesn't work on project Johan's sequences because all the reactivity values are 0
				# and can't take percentile


			# remove nucleotides preceded by 6 (or more) A's
			tmp_concat_data = tmp_concat_data.loc[tmp_concat_data['in_polyA']==0]
			tmp_concat_data[agg_field] = kind

			corr_data = calculate_metric(tmp_concat_data, x_data=x_data, y_data=y_data, n_bootstraps=n_bootstraps,
				package_list=package_list, metric=metric)

			corr_data[agg_field] = kind

			correlation_data = correlation_data.append(corr_data,ignore_index=True)

	return correlation_data

def percentile(n):
    def percentile_(x):
        return np.percentile(x, n)
    percentile_.__name__ = 'percentile_%s' % n
    return percentile_
    
if __name__=='__main__':

	p = argparse.ArgumentParser(description="""Score written predictions""")

	p.add_argument("infile", action="store", help="Input json")
	p.add_argument("--test", action='store_true', help='Tests first 3 constructs in each package.')
	p.add_argument("--subset", action='store_true', help='Only run for one representative set of options for each package.')
	p.add_argument("--metadata", action='store', help='Dataframe in .json format. If provided, will use metadata from this df\
		instead of the input dataframe (For instance, if there is an error, or if there are duplicates of sequences for\
		different conditions.')
	p.add_argument("--field_to_aggregate", action='store', default='Dataset', help="Field to aggregate results over. In chem mapping, 'filename' for experiments, or 'project_name' for projects.")
	p.add_argument("--n_bootstraps", action='store',type=int, default=10)
	p.add_argument("--metric", action='store',default='spearman', help='spearman, pearson, or rmse')
	p.add_argument("-v", "--verbose", action="store_true", help="Verbose")
	p.add_argument("-o", action="store", dest='outfile', help='name of output json file (default is <input_name>_BOOTSTRAPS.json.zip')

	p.add_argument('--reactivity_field', action='store',default='reactivity', help='DataFrame field that contains vector of reactivity values. Default is "reactivity".')
	
	args = p.parse_args()
	print(args)

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

	correlation_data = ScoreChemMapping(df, agg_field = args.field_to_aggregate, 
		n_bootstraps=args.n_bootstraps, package_list=None, metric=args.metric,
		x_data = args.reactivity_field, y_data='p')
	
	correlation_data.to_json(outfile+'_BOOTSTRAPS.json.zip')
