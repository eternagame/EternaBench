import pandas as pd
import os

data_path = os.environ['ETERNABENCH_PATH']

def load_CM_data():
	print("Loading EternaBench ChemMapping Filtered dataset (n=14339), last update Jul 10 2021")
	return pd.read_json(data_path+"/data/EternaBench_ChemMapping_Filtered_10Jul2021.json.zip")

def load_CM_unfiltered_data():
	print("Loading EternaBench ChemMapping FULL dataset (n=41876), last update Jul 10 2021")
	return pd.read_json(data_path+"/data/EternaBench_ChemMapping_Full_10Jul2021.json.zip")

def load_CM_example_calculations():
	print("Loading EternaBench ChemMapping Round 1 with example calculations")
	return pd.read_json(data_path+'/data/EternaBench_ChemMapping_Example_PunpVectors_Round1.json.zip')

def load_CM_example_correlations():
	print("Loading bootstrapped pearson correlations for EternaBench ChemMapping Round 1")
	return pd.read_csv(data_path+'/scoring_data/EternaBench_ChemMapping_Example_Correlations_Round1.csv')

def load_CM_scores():
	print("Loading scores DONE_17Apr2021ScoreDF.json.zip")
	df = pd.read_json(data_path+'/scoring_data/DONE_17Apr2021ScoreDF.json.zip')
	print("Fields: %s" %', '.join(list(df.keys())))
	return df

def load_CM_calculations():
	print('Warning! This is a big file')
	return pd.read_json(data_path+"/ChemMapping/data/DONE_17Apr2021.json.zip")
