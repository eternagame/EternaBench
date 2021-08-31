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

cd ${ETERNABENCH_PATH}

# python scripts/CompileBootstrappedResults.py data/ChemMapping/bootstraps/CM_pearson_Dataset -o scoring_data/EB_-efold  --dataset_field Dataset --metric pearson\
#  --calculate_Z_scores scripts/package_benchmark_-eternafold.txt --dataset_list scripts/EternaBench_ChemMapping_datasets.txt

# python scripts/CompileBootstrappedResults.py data/ChemMapping/bootstraps/CM_pearson_Dataset -o scoring_data/EB_+efold_test  --dataset_field Dataset --metric pearson\
#  --calculate_Z_scores scripts/package_benchmark_+eternafold.txt --dataset_list scripts/EternaFold_test_datasets.txt

# python scripts/CompileBootstrappedResults.py data/ChemMapping/bootstraps/CM_pearson_Dataset -o scoring_data/EB_all_eternafold_versions  --dataset_field Dataset --metric pearson\
#  --calculate_Z_scores scripts/package_list.txt  --dataset_list scripts/EternaFold_test_datasets.txt

# python scripts/CompileBootstrappedResults.py data/ChemMapping/bootstraps/CM_spearman_Dataset -o scoring_data/EB  --dataset_field Dataset --metric spearman\
#  --calculate_Z_scores scripts/package_benchmark_-eternafold.txt  --dataset_list scripts/EternaBench_ChemMapping_datasets.txt

python scripts/CompileBootstrappedResults.py data/ChemMapping/bootstraps/CM_pearson_project -o scoring_data/EB_projects_-efold  --dataset_field project_name --metric pearson\
 --calculate_Z_scores scripts/package_benchmark_-eternafold.txt --dataset_list scripts/EternaBench_ChemMapping_datasets.txt

#python scripts/CompileBootstrappedResults.py data/ChemMapping/bootstraps/CM_pearson_project -o scoring_data/EB_projects_+efold  --dataset_field project_name --metric pearson --calculate_Z_scores scripts/package_benchmark_+eternafold.txt  --dataset_list scripts/EternaFold_test_datasets.txt

#riboswitches
