## Chemical Mapping

### 1. Create input files

Runtime: 4 minutes on Sherlock Login Node

```python
cd EternaBench/ChemMapping/scripts
python GenerateChemMappingDatasets.py
```

Start of output, indicating Round 0 got processed correctly:

```
Output of MAPseeker v1.2
from data: 122112_RD_EteRNAPlayerProjects_1M7test
{'modifier': ['1M7'], 'MAPseq': ['tag:RTB000'], 'chemical': ['MgCl2:10mM', 'HEPES:50mM(pH8.0)'], 'temperature': ['24C'], 'processing': ['overmodificationCorrectionExact', 'ligationBiasCorrection:0.500', 'backgroundSubtraction', 'normalization:boxplot']}
Round 0 2058
                    EteRNA        ID        MAPseq  ... signal_to_noise                                          structure temperature
0  score:EteRNA_score:86.7  769718-1  [tag:RTB000]  ...    medium:1.437  .............................((((.((((.((((......       [24C]
1  score:EteRNA_score:68.0  769337-1  [tag:RTB000]  ...    medium:3.227  .............................((((.((((.((((......       [24C]
2  score:EteRNA_score:70.7  769330-1  [tag:RTB000]  ...    medium:2.590  .............................((((.((((.((((......       [24C]
3  score:EteRNA_score:73.3  769323-1  [tag:RTB000]  ...    medium:1.012  .............................((((.((((.((((......       [24C]
4  score:EteRNA_score:89.3  767147-1  [tag:RTB000]  ...    medium:2.130  .............................((((.((((.((((......       [24C]

[5 rows x 17 columns]
Running CD-HIT-EST on ETERNA_R00_0000...
Getting sequences from each cluster
```

Produces two files in `EternaBench/data`:

`EternaBench_ChemMapping_Full_<todaysdate>.json.zip`: Contains all constructs (N=41876, as of Jul 8 2021). 21 MB.

`EternaBench_ChemMapping_Filtered_<todaysdate>.json.zip`: Contains constructs filtered using CD-HIT-EST <80% sequence similarity cutoff (N=14339, Jul 8 2021). 7.3 MB.

To update what dataset the EternaBench API loads, change this in `eternabench/load_data.py`

```python
def load_CM_data():
	return pd.read_json(data_path+"/ChemMapping/scripts/Jul082021_filtered_cloud_lab_rounds.json")

def load_CM_example_calculations():
	return pd.read_json(data_path+'ChemMapping/data/Round69_10May2021.json.zip')

def load_CM_calculations():
	print('Warning! This is a big file')
	return pd.read_json(data_path+"/ChemMapping/data/DONE_17Apr2021.json.zip")
```

In filtered_df, 101 rows without project names.

### 2. generate p(unpaired) vector calculations

(This is best done as a job on a cluster. Sherlock, 30 GB memory, using Pandarallel parallelization: XX hours)

A. Command to test the script is running correctly:

```
cd EternaBench/scripts
python CalculatePunpVectors.py ../data/EternaBench_ChemMapping_Filtered_08Jul2021.json.zip --verbose --test -o example_output/TEST
```

This should print the following log:

```
#### ETERNABENCH: Calculate P(Unpaired) Vectors #####
Last Update July 2021, H Wayment-Steele
Running calculations using Arnie located at ['/home/users/hannahw1/arnie']
Arnie paths:
rnastructure: /home/groups/rhiju/hannahw1/secstruct_software/RNAstructure/exe
rnasoft: /home/groups/rhiju/hannahw1/secstruct_software/MultiRNAFold
contrafold_1: /home/groups/rhiju/hannahw1/secstruct_software/contrafold_v1_10/src
contrafold_2: /home/groups/rhiju/hannahw1/secstruct_software/contrafold_v2_02/src
contrafold_se: /home/groups/rhiju/hannahw1/secstruct_software/contrafold-se/src
linearpartition: /home/groups/rhiju/hannahw1/secstruct_software/LinearPartition/bin
linearfold: /home/groups/rhiju/hannahw1/secstruct_software/LinearFold/bin
vienna_2: /home/groups/rhiju/hannahw1/secstruct_software/ViennaRNA-2.4.10/src/bin
vienna_1: /home/groups/rhiju/hannahw1/secstruct_software/ViennaRNA-1.8.5/bin
nupack: /home/groups/rhiju/hannahw1/secstruct_software/nupack3.2.2/build/bin
eternafoldparams: /home/groups/rhiju/hannahw1/secstruct_software/EternaFoldParams/params.v1
TMP: /tmp
writing p_unp for vienna_2
writing p_unp for vienna_2_nodangles
writing p_unp for vienna_2_60C
writing p_unp for vienna_1
writing p_unp for nupack_99
writing p_unp for nupack_99_60C
writing p_unp for nupack_99_nodangles
writing p_unp for nupack_95
writing p_unp for nupack_95_nodangles
writing p_unp for rnastructure
writing p_unp for rnastructure_nocoax
writing p_unp for rnastructure_60C
writing p_unp for contrafold_1
writing p_unp for contrafold_2
writing p_unp for contrafold_2_nc
writing p_unp for vienna_langdon_pars
writing p_unp for vienna_rnasoft_pars
writing p_unp for rnasoft_99
writing p_unp for rnasoft_07
writing p_unp for rnasoft_blstar
writing p_unp for rnasoft_99_nodangles
writing p_unp for rnasoft_bl_nodangles
writing p_unp for rnasoft_lam-cg
writing p_unp for rnasoft_nom-cg
writing p_unp for eternafold_A_cfold
writing p_unp for eternafold_B_cfold
writing p_unp for eternafold_C_cfold
writing p_unp for eternafold_D_cfold
writing p_unp for eternafold_E_cfold
writing p_unp for eternafold_F_cfold
writing p_unp for eternafold_G_cfold
writing p_unp for eternafold_A
writing p_unp for eternafold_B
writing p_unp for eternafold_C
writing p_unp for eternafold_D
writing p_unp for eternafold_E
writing p_unp for eternafold_F
writing p_unp for eternafold_G
Writing output, chunked by Dataset, to TEST_Dataset_chunks_08Jul2021
```

This should produce the file `TEST.json.zip`, which contains p(unp) calculations for the first three constructs in the dataset.

When this script is run on all the cloud lab rounds, the resulting calculations file will be quite large. The dataset is also returned in chunks by Cloud lab rounds in the dir `TEST_chunks_<todaysdate>`. For this test case, the dir will only have one output file (`Round_0.json.zip`).

B. An example submission script to run on the full dataset is provided at `EternaBench/scripts/slurm_run_CM_benchmark.sh`. Relevant lines below:

```
source ~/.bash_profile

OUTNAME=${ETERNABENCH_PATH}/data/ChemMapping_PunpVectors_08Jul2021

cd ${ETERNABENCH_PATH}/scripts
python CalculatePunpVectors.py ${ETERNABENCH_PATH}/data/EternaBench_ChemMapping_Filtered_08Jul2021.json.zip --parallel --verbose -o $OUTNAME
```

### 3. Calculate correlation between p(unpaired) calculations and reactivity values.

## Riboswitches

### 1. Create input files

Runtime:

```
cd Riboswitch/data/datasets_without_predictions
python ../../scripts/write_ribologic_orig.py # writes ribologic_orig.json
python ../../scripts/write_eterna_filtered_json.py # writes processed_eterna_df.json
python ../../scripts/merge_jsons.py processed_eterna_df.json ribologic_orig.json all_switch_data.json
```

### 2. Run K_fold predictions

Runtime:

```
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
cd Riboswitch
python scripts/write_estimates.py data/preprocessing/all_switch_data.json -v --parallel --flanking -o all_packages_bp
```

To write estimates for package subset for all datasets using constrained partition function estimation:

```
cd Riboswitch/data/datasets_with_predictions
python ../../scripts/write_estimates.py ../datasets_without_predictions/all_switch_data.json -v --method Z --subset --parallel --flanking -o all_packages_bp
```

Note: input dataframe must already have sequence constraint strings for input, output, and input+output aptamers to use `Z` method.

`write_estimates.py` calls utils provided in `Riboswitch/scripts/pandas_arnie_utils.py`. These use Pandarallel for parallelized computing on the pandas dataframe and arnie for stat mech calculations. If an engine other than arnie is desired, a new util can be written in the style of the `pandas_arnie_utils`.
```

### 3. Score K_fold predictions