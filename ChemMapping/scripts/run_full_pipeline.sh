# this is mainly for reference!  Recommend against running the full thing like this unless in a job on a cluster.

# writes all_CL_unfiltered.json and "/all_CL_filtered.json
#python write_preprocessed_cloud_lab_jsons.py

# Runs arnie, which uses pandas dataframe handling to write p(unpaired) vector estimates to the dataframe. saves as another dataframe.
# --parallel flag handles pandas parallelization using the package pandarallel.

python write_estimates_in_arnie.py ../data/preprocessing/all_cloud_labs_filtered.json -o ../data/full_dataframes/all_cloud_labs.json --parallel -v
python write_estimates_in_arnie.py ../data/preprocessing/CE_rounds.json -o ../data/full_dataframes/CE_rounds.json --parallel -v

# Concatenates reactivity data from individual constructs, nucleotide information, and p(unpaired) values, to a single column.
# Retains construct-level information (modifier, temperature, design name, etc.)
# also filters reactivity datapoints to throw out values below zero and outliers over specified percentile cutoff.

#python filter_concatenate_dataset.py ../data/full_dataframes/all_cloud_labs.json -o ../data/cleaned_dataframes/all_cloud_labs.json
#python filter_concatenate_dataset.py ../data/full_dataframes/CE_rounds.json -o ../data/cleaned_dataframes/CE_rounds.json

# Run analysis and produce plots.
#python write_cloud_lab_round_analysis.py
