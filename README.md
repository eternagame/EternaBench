# EternaBench

This repository contains the EternaBench datasets as well as scripts to reproduce the analysis presented in Wayment-Steele et al. (2021).


## Setup

Add to python path and point to datasets by adding to .bashrc:
```bash
export PYTHONPATH=/path/to/EternaBench
export ETERNABENCH_PATH=/path/to/EternaBench
```

## Use cases (from least to most resource-intensive)

(An ask.) A clear use-case for this repository is to benchmark a novel algorithm against the algorithms contained here. This code uses [Arnie](https://github.com/DasLab/arnie/) to wrap the algorithms tested in this work.

If you have an algorithm for which you wish to demonstrate its superior performance on these datasets, we ask you might consider checking in a PR to Arnie to wrap it. This will not only make your benchmarking easier, it will also make your algorithm immediately available for other RNA thermodynamics aficionados to use with these data and the Arnie ecosystem!

Instructions for linking base-pair probability calculations to Arnie are [here](). Basically, the algorithm only needs to provide a symmetric matrix of probabilities p(i:j) as a numpy array.

### I want to play with the data behind the figures in the paper, or use pre-calculated correlation coefficient or z-score data

Notebooks in `analysis` regenerate all the figures in the manuscript. Each figure cell indicates a path to a csv that contains the raw correlation and z-score data. Datasets for representative experiments and packages are also included.

### I want to regenerate thermodynamic calculations and z-score calculations for a representative chemical mapping and/or riboswitch dataset on a single core

1. Git clone [Arnie](https://github.com/DasLab/arnie/).

2. Follow the Arnie instructions [here](https://github.com/DasLab/arnie/blob/master/docs/setup_doc.md) to set up all the packages you want to rerun.

3. modify runDemo.sh to iterate over the packages you wish to run.

4.

```
cd ${ETERNABENCH_PATH}/DEMO
chmod +x run_demo.sh
./run_demo.sh
```

### I want to regenerate thermodynamic calculations for all the datasets on a cluster

The slurm scripts used to generate the data for this paper are contained in `/cluster_scripts`.

```
./SubmitParallelChemMapping.sh
./SubmitParallelRiboswitch.sh
./SubmitParallelExternalData.sh
```

### I want to regenerate the filtered EternaBench datasets from the raw data

1. Git clone [RDatKit](https://github.com/ribokit/RDATKit.git) and follow instructions there to your python path.

2. Git clone [CD-HIT](https://github.com/weizhongli/cdhit) and export its path:

```
CDHIT_PATH='/path/to/cdhit'
```

3. Run the below python scripts.
```
python ${ETERNABENCH_PATH}/scripts/GenerateChemMappingDatasets.py
python ${ETERNABENCH_PATH}/scripts/GenerateRiboswitchDatasets.py
```

# Organization 

`analysis`: python notebooks to reproduce paper figures.

`cluster_scripts`: scripts to run entire EternaBench benchmarking using SLURM cluster system.

`data`: EternaBench datasets.

- 	`EternaBench_*.json.zip`: Full and filtered EternaBench datasets without calculations.
-	`ChemMappingPreprocessing`: raw datasets used to create chem mapping benchmark.
-	`RiboswitchPreprocessing`: raw datasets used to create riboswitch benchmark.
-	`RiboswitchCalculations`: Example datasets with K_fold calculations.
-	`ChemMappingCalculations`: Example datasets with p(unpaired) calculations.
-	`ExternalData`: inputs and calculations for external collected datasets.

`DEMO`: single script to regenerate observable calculations for one representative dataset from Chem Mapping and Riboswitch, and calculate significance.

`docs`: Documentation.

`eternabench`: EternaBench API source.

`scoring_data`: CSVs containing all data and metrics used for evaluation in the paper. These data are the raw input to the figures plotted in the notebooks in `analysis`.

`scripts`: scripts to calculate observables and bootstrap correlation significance over datasets. Full documentation is in `docs/RunBenchmarkREADME.md`.




# Data Origin

- Chemical Mapping RDAT files may be downloaded from www.rmdb.stanford.edu.

- Eterna riboswitch datasets are detailed in the supporting information of 

Andreasson, J. O., ... & Das, R., Greenleaf, W. J. (2019). Crowdsourced RNA design discovers diverse, reversible, efficient, self-contained molecular sensors. bioRxiv, 877183.

- Ribologic riboswitch dataset is detailed in the supporting information of

Wu, M. J., Andreasson, J. O., Kladwang, W., Greenleaf, W., & Das, R. (2019). Automated design of diverse stand-alone riboswitches. ACS synthetic biology, 8(8), 1838-1846.

