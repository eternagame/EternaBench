# EternaBench

This repository contains the EternaBench datasets and accompanying code, which evaluate RNA structure prediction using diverse thermodynamic prediction tasks and high-throughput datasets (Wayment-Steele et al, 2021).


## Setup

Add to python path and point to datasets by adding to .bashrc:
```bash
export PYTHONPATH=/path/to/EternaBench
export ETERNABENCH_PATH=/path/to/EternaBench
```

## Use cases

### I want to tinker with the data in the paper figures, or use pre-calculated correlation or z-score data

Notebooks in `analysis` regenerate all the figures in the manuscript. Each figure cell indicates a path to a csv that contains the raw correlation and z-score data. Datasets for representative experiments and packages are also included.

### I want to benchmark my novel algorithm against the algorithms contained here

This code uses [Arnie](https://github.com/DasLab/arnie/) to wrap the algorithms tested in this work.

If you have an algorithm that you want to demo on these datasets, we recommend checking in a PR to Arnie to wrap your algorithm.

This will make benchmarking easier and will also make your algorithm immediately available for other Arnie-wielding RNA thermodynamics fans to use!

Instructions for linking base-pair probability calculations to Arnie are [here](https://github.com/eternagame/EternaBench/blob/master/docs/linkToArnie.md). Briefly, the algorithm just needs to provide a symmetric matrix of probabilities p(i:j) as a numpy array.

### I want to regenerate thermodynamic calculations and z-score calculations for an example chemical mapping and riboswitch dataset on a single core

1. Git clone [Arnie](https://github.com/DasLab/arnie/).

2. Follow the Arnie instructions [here](https://github.com/DasLab/arnie/blob/master/docs/setup_doc.md) to set up all the packages you want to rerun.

3. modify package_list.txt to iterate over the packages you wish to run.

4. Run the bash script: ./run_demo.sh

5. Successful completion will end in a call to calculate bootstrapped Pearson correlation coefficients and z-scores.

```
python calculateZscoreDEMO.py 
Chem Mapping Rnd 1 scores
        package  pearson_mean  pearson_std  pearson_zscore_by_Dataset_mean
1    eternafold      0.738693     0.001563                        0.855730
0  contrafold_2      0.718083     0.001729                        0.242594
2      vienna_2      0.672972     0.002110                       -1.098324

Riboswitch "Ribologic FMN" scores
        package  pearson_mean  pearson_std  pearson_zscore_by_Dataset_mean
1    eternafold      0.642524     0.013263                        1.004608
0  contrafold_2      0.502365     0.017499                       -0.017645
2      vienna_2      0.368747     0.020251                       -0.986963
```

Example outputs are in `DEMO/example_outputs_from_demo`.

### I want to regenerate thermodynamic calculations for all the datasets on a cluster

The Slurm scripts used to generate the data for this paper are contained in `/cluster_scripts`. You will need to modify the Slurm headers for your own environment.

```
cd ${ETERNABENCH_PATH}/cluster_scripts
./SubmitParallelChemMapping.sh
./SubmitParallelRiboswitch.sh
./SubmitParallelExternalData.sh
```

### I want to regenerate the filtered EternaBench datasets from the raw data

1. Git clone [RDatKit](https://github.com/ribokit/RDATKit.git) and follow instructions there to add it your python path.

2. Git clone [CD-HIT](https://github.com/weizhongli/cdhit) and export its path:

```
CDHIT_PATH='/path/to/cdhit'
```

3. Run the below python scripts.
```
python ${ETERNABENCH_PATH}/scripts/GenerateChemMappingDatasets.py
python ${ETERNABENCH_PATH}/scripts/GenerateRiboswitchDatasets.py
```

Takes about 12 minutes runtime to regenerate both. Example intermediate CDHIT outputs are provided in `cluster_scripts/CDHIT_example_output`.

# Organization 

`analysis`: python notebooks to reproduce paper figures.

`cluster_scripts`: scripts to run entire EternaBench benchmarking using SLURM cluster system.

`data`:

-	`DEMO_ChemMapping.json.zip`: Input data for "Cloud lab Round 1", example chemical mappingdataset discussed in main text.
-	`DEMO_Riboswitch.json.zip`: Input data for "Ribologic FMN" dataset, example riboswitch dataset discussed in main text. 
-  `EternaBench_*.json.zip`: Full and filtered EternaBench datasets without calculations.
-	`ChemMappingPreprocessing`: Initial datasets used to create chem mapping benchmark.
-	`RiboswitchPreprocessing`: Initial datasets used to create riboswitch benchmark.
-	`RiboswitchCalculations`: Example datasets with K_fold calculations (see notebooks in `analysis` for example calls to plot these).
-	`ChemMappingCalculations`: Example datasets with p(unpaired) calculations (see notebooks in `analysis` for example calls to plot these).
-	`ExternalData`: inputs and calculations for external collected datasets (see notebooks in `analysis` for example calls to plot these).

`DEMO`: One non-parallelized script to regenerate observable calculations for one representative dataset from Chem Mapping and Riboswitch, and calculate significance.

`docs`: Documentation.

`eternabench`: EternaBench API source.

`scoring_data`: CSVs containing all data and metrics used for evaluation in the paper. These data are the raw input to the figures plotted in the notebooks in `analysis`.

`scripts`: Scripts to calculate observables and bootstrap correlation significance over datasets.


# Data Origin


- *Chemical Mapping* RDAT files may be downloaded from www.rmdb.stanford.edu.

- *Riboswitch* datasets are detailed in the supporting information of 

Andreasson, J. O., ... & Das, R., Greenleaf, W. J. (2019). Crowdsourced RNA design discovers diverse, reversible, efficient, self-contained molecular sensors. bioRxiv, 877183.

- Ribologic riboswitch dataset is detailed in the supporting information of

Wu, M. J., Andreasson, J. O., Kladwang, W., Greenleaf, W., & Das, R. (2019). Automated design of diverse stand-alone riboswitches. ACS synthetic biology, 8(8), 1838-1846.


- Links to download the *Chemical mapping datasets from other papers*, which permit their re-use in the EternaBench dataset form, is included at [docs/Eternabench_external_dataset_license_info.csv](docs/Eternabench_external_dataset_license_info.csv).

# Contact

rhiju@stanford.edu

hannah.wayment.steele@gmail.com
