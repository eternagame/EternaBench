#!/bin/bash
#SBATCH --job-nam=Score
#SBATCH -o /home/groups/rhiju/hannahw1/EternaBench/scripts/example_output/Score.out
#SBATCH -p biochem,normal
#SBATCH -N 2
#SBATCH --mem-per-cpu 30GB
#SBATCH --mail-user=hannahw1@stanford.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --time=40:00:00


source ~/.bash_profile
grep -m 1 'cpu cores' /proc/cpuinfo

INPUT=${ETERNABENCH_PATH}/data/DONE_17Apr2021.json.zip
OUTPUT=${ETERNABENCH_PATH}/scoring_data/EternaBench_ChemMapping_by_project_pearson
METADATA=${ETERNABENCH_PATH}/data/EternaBench_ChemMapping_Filtered_10Jul2021.json.zip
cd ${ETERNABENCH_PATH}/scripts

python ScoreChemMapping.py $INPUT --n_bootstraps 1000 --similarity pearson -o ${OUTNAME}_score_by_project --field_to_aggregate project_name --metadata $METADATA
