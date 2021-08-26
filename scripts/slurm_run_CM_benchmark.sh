#!/bin/bash
#SBATCH --job-nam=CalculatePunpVectors
#SBATCH -o /home/groups/rhiju/hannahw1/EternaBench/scripts/example_output/CalculatePunpVectors.out
#SBATCH -p biochem,owners,normal
#SBATCH -n 3
#SBATCH --mem-per-cpu=30G
#SBATCH --mail-user=hannahw1@stanford.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --time=40:00:00


source ~/.bash_profile

OUTNAME=${ETERNABENCH_PATH}/data/ChemMapping_PunpVectors_08Jul2021

cd ${ETERNABENCH_PATH}/scripts
python CalculatePunpVectors.py ${ETERNABENCH_PATH}/data/EternaBench_ChemMapping_Filtered_08Jul2021.json.zip --parallel --verbose -o $OUTNAME
# python ScoreChemMapping.py $OUTNAME --n_bootstraps 1000 --similarity pearson -o ${OUTNAME}_score_by_expt
# python ScoreChemMapping.py $OUTNAME --n_bootstraps 1000 --similarity pearson --field_to_aggregate project_name -o ${OUTNAME}_score_by_project
