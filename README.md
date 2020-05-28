# Organization:

This repository contains the EternaBench datasets as well as scripts to reproduce the analysis presented in Wayment-Steele & Das (2020).

# Data Origin

- Chemical Mapping RDAT files may be downloaded from www.rmdb.stanford.edu.

- Eterna riboswitch datasets are detailed in the supporting information of 

Andreasson, J. O., ... & Das, R., Greenleaf, W. J. (2019). Crowdsourced RNA design discovers diverse, reversible, efficient, self-contained molecular sensors. bioRxiv, 877183.

- Ribologic riboswitch dataset is detialed in the supporting information of

Wu, M. J., Andreasson, J. O., Kladwang, W., Greenleaf, W., & Das, R. (2019). Automated design of diverse stand-alone riboswitches. ACS synthetic biology, 8(8), 1838-1846.

# General overview of pipeline

The `write_estimates.py` scripts for each dataset call utils provided in `ChemMapping/scripts/RDAT_utils.py` and `Riboswitch/scripts/pandas_arnie_utils.py`, respectively. These use Pandarallel for parallelized computing on the pandas dataframe and arnie for stat mech calculations. If an engine other than arnie is desired, a new util can be written in the style of the functions in these.

# ChemMapping

Raw datasets are contained in `data/preprocessing`.  Datasets with p(unpaired) estimates are in `data/full_datasets`.

## To write `EternaBench-CM` dataset from raw `.rdat` files

`cd ChemMapping/data/preprocessing`
`unzip raw_rdats.zip`
`python ../../scripts/write_preprocessed_cloud_lab_jsons.py`

This creates the dataset `data/preprocessing/filtered_cloud_lab_rounds.json.zip`, filtered for sequence redundancy within datasets using `CD-HIT-EST`. The unfiltered dataset, `unfiltered_cloud_lab_rounds.zip`, is produced with the same script (comment out the call to `filter_dataframe_with_cdhit`).

## To write p(unp) estimates

To write with all packages for Round 69:

`cd ChemMapping/data/full_datasets`
`python ../../scripts/write_estimates_in_arnie.py ../preprocessing/Round69.json.zip -v --parallel`

To write with subset of packages for all data:

`cd ChemMapping/data/full_datasets`
`python ../../scripts/write_estimates_in_arnie.py ../preprocessing/filtered_cloud_lab_rounds.json.zip -v --subset --parallel`

## To perform correlation analysis

Analyses are described in jupyter notebook format in the `analysis` directory.

# Riboswitch

## To write `EternaBench-Switch` dataset

```
cd Riboswitch/data/datasets_without_predictions
python ../../scripts/write_ribologic_orig.py # writes ribologic_orig.json
python ../../scripts/write_eterna_filtered_json.py # writes processed_eterna_df.json
python ../../scripts/merge_jsons.py processed_eterna_df.json ribologic_orig.json all_switch_data.json
```

## To write k_MS2 estimates

To write estimates for all packages for Ribologic dataset using base pair estimation:

```
cd Riboswitch/data/datasets_with_predictions
python ../../scripts/write_estimates.py ../datasets_without_predictions/all_switch_data.json -v --parallel --ribologic_only --flanking -o all_packages_bp
```

To write estimates for all packages for Ribologic dataset using constrained partition function estimation:

```
cd Riboswitch/data/datasets_with_predictions
python ../../scripts/write_estimates.py ../datasets_without_predictions/all_switch_data.json -v --method Z --parallel --ribologic_only --flanking -o all_packages_bp
```

To write estimates for package subset for all datasets using base pair estimation:

```
cd Riboswitch/data/datasets_with_predictions
python ../../scripts/write_estimates.py ../datasets_without_predictions/all_switch_data.json -v --parallel --flanking -o all_packages_bp
```

To write estimates for package subset for all datasets using constrained partition function estimation:

```
cd Riboswitch/data/datasets_with_predictions
python ../../scripts/write_estimates.py ../datasets_without_predictions/all_switch_data.json -v --method Z --subset --parallel --flanking -o all_packages_bp
```

Note: input dataframe must already have sequence constraint strings for input, output, and input+output aptamers to use `Z` method.

`write_estimates.py` calls utils provided in `Riboswitch/scripts/pandas_arnie_utils.py`. These use Pandarallel for parallelized computing on the pandas dataframe and arnie for stat mech calculations. If an engine other than arnie is desired, a new util can be written in the style of the `pandas_arnie_utils`.

## To perform correlation and RMSE analysis

Analyses are described in jupyter notebook format in the `analysis` directory.

# ExternalDatasets

## ExternalChemMapping

This follows the same outline as the ChemMapping processing described above. Briefly,

```
python scripts/write_estimates_in_arnie.py data/full_df_window600.json
```

Analysis may be found in `analysis/external_chem_mapping.ipynb`.

## PUM-binding data from Jarmoskaite (2019)

Writing estimates and analysis may be found in `analysis/PUM.ipynb`.