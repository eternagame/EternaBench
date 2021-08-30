# EternaBench

This repository contains the EternaBench datasets as well as scripts to reproduce the analysis presented in Wayment-Steele et al. (2021).


## Setup

Add to python path and point to datasets by adding to .bashrc:
```bash
export PYTHONPATH=/path/to/EternaBench
export ETERNABENCH_PATH=/path/to/EternaBench
```

## Use cases (from least to most resource-intensive)

(NB) A clear use-case for this repository is to benchmark a novel algorithm against the algorithms contained here. This code uses [Arnie](https://github.com/DasLab/arnie/) to wrap the algorithms tested in this work.

If you have an algorithm that you wish to demonstrate its superior performance on these datasets on, I recommend checking in a PR to Arnie to wrap it, this will make benchmarking it easier, and also make your algorithm immediately available for others to use with these data!

### I want to use pre-calculated correlation coefficient or z-score data from the paper, or regenerate figures in the preprint

Notebooks in `analysis` regenerate all the figures in the manuscript. Each figure cell indicates a path to a csv that contains the correlation mean, standard deviation and z-score statistics for each dataset and each package.

### I want to regenerate thermodynamic calculations and z-score calculations for a representative chemical mapping and/or riboswitch dataset on a single core

1. Git clone [Arnie](https://github.com/DasLab/arnie/).

2. Follow the Arnie instructions [here](https://github.com/DasLab/arnie/blob/master/docs/setup_doc.md) to set up all the packages you want to rerun.

3. modify runDemo.sh to iterate over the packages you wish to run.

4.

```
cd /path/to/eternabench
chmod +x runDemo.sh
./runEternaBench.sh
```

### I want to regenerate thermodynamic calculations for all the datasets on a cluster

The slurm scripts used to generate the data for this paper are contained in `/cluster_scripts`.

Once paths are set, the jobs can be started via
```
./runChemMapping.sh
./runRiboswitch.sh
./runExternal.sh
```

### I want to regenerate the filtered EternaBench datasets from the raw data

1. Git clone [CD-HIT](https://github.com/weizhongli/cdhit) and export its path:

```
CDHIT_PATH='/path/to/cdhit'
```

2. Run the below python scripts.
```
cd /path/to/EternaBench/data/chemmapping_preprocessing
python ../../scripts/GenerateChemMappingDatasets.py
cd /path/to/EternaBench/data/riboswitch_preprocessing
python ../../scripts/GenerateRiboswitchDatasets.py
```

# Organization 

`data`: datasets

`scripts`: scripts to regenerate benchmark. Full documentation is in `docs/RunBenchmarkREADME.md`.

`analysis`: python notebooks to reproduce the figures in (Wayment-Steele, 2021).

`eternabench`: EternaBench API source.


# Data Origin

- Chemical Mapping RDAT files may be downloaded from www.rmdb.stanford.edu.

- Eterna riboswitch datasets are detailed in the supporting information of 

Andreasson, J. O., ... & Das, R., Greenleaf, W. J. (2019). Crowdsourced RNA design discovers diverse, reversible, efficient, self-contained molecular sensors. bioRxiv, 877183.

- Ribologic riboswitch dataset is detailed in the supporting information of

Wu, M. J., Andreasson, J. O., Kladwang, W., Greenleaf, W., & Das, R. (2019). Automated design of diverse stand-alone riboswitches. ACS synthetic biology, 8(8), 1838-1846.

