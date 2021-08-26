#!/bin/bash
#SBATCH --job-nam=compile.sh
#SBATCH -o /home/groups/rhiju/hannahw1/EternaBench/scoring_data/compile.sh.out
#SBATCH -p biochem,owners,normal
#SBATCH -n 1
#SBATCH --mail-user=hannahw1@stanford.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --mem-per-cpu=30G
#SBATCH --cpus-per-task=6
#SBATCH --time=10:00:00


source ~/.bash_profile

cd /home/groups/rhiju/hannahw1/EternaBench/scoring_data

python ../scripts/CompileBootstrappedResults.py ../data/Ext300/CM_pearson_Dataset_ -o Ext300 --dataset_field Dataset --metric pearson --calculate_Z_scores ../scripts/package_subset.txt
python ../scripts/CompileBootstrappedResults.py ../data/Ext600/CM_pearson_Dataset_ -o Ext600 --dataset_field Dataset --metric pearson --calculate_Z_scores ../scripts/package_subset.txt
python ../scripts/CompileBootstrappedResults.py ../data/Ext900/CM_pearson_Dataset_ -o Ext900 --dataset_field Dataset --metric pearson --calculate_Z_scores ../scripts/package_subset.txt
python ../scripts/CompileBootstrappedResults.py ../data/Ext1200/CM_pearson_Dataset_ -o Ext1200 --dataset_field Dataset --metric pearson --calculate_Z_scores ../scripts/package_subset.txt

python ../scripts/CompileBootstrappedResults.py ../data/ChemMapping/bootstraps/CM_pearson_Dataset -o EB_-efold  --dataset_field Dataset --metric pearson --calculate_Z_scores ../scripts/package_benchmark_-eternafold.txt
python ../scripts/CompileBootstrappedResults.py ../data/ChemMapping/bootstraps/CM_pearson_Dataset -o EB_+efold  --dataset_field Dataset --metric pearson --calculate_Z_scores ../scripts/package_benchmark_+eternafold.txt
python ../scripts/CompileBootstrappedResults.py ../data/ChemMapping/bootstraps/CM_pearson_Dataset -o EB_all_eternafold_versions  --dataset_field Dataset --metric pearson --calculate_Z_scores ../scripts/package_list.txt 
python ../scripts/CompileBootstrappedResults.py ../data/ChemMapping/bootstraps/CM_spearman_Dataset -o EB  --dataset_field Dataset --metric spearman --calculate_Z_scores ../scripts/package_benchmark_-eternafold.txt 
python ../scripts/CompileBootstrappedResults.py ../data/ChemMapping/bootstraps/CM_pearson_project -o EB_projects_-efold  --dataset_field project_name --metric pearson --calculate_Z_scores ../scripts/package_benchmark_-eternafold.txt 
python ../scripts/CompileBootstrappedResults.py ../data/ChemMapping/bootstraps/CM_pearson_project -o EB_projects_+efold  --dataset_field project_name --metric pearson --calculate_Z_scores ../scripts/package_benchmark_-eternafold.txt 

#riboswitches
