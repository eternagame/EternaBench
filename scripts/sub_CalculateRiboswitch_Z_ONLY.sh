#!/bin/bash
#SBATCH --job-nam=CalculateRiboswitch
#SBATCH -o /home/groups/rhiju/hannahw1/EternaBench/scripts/example_output/%j.out
#SBATCH -e /home/groups/rhiju/hannahw1/EternaBench/scripts/example_output/%j.err
#SBATCH -p biochem,normal
#SBATCH -N 2
#SBATCH --mem-per-cpu 30GB
#SBATCH --mail-user=hannahw1@stanford.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --time=40:00:00

source ~/.bash_profile
grep -m 1 'cpu cores' /proc/cpuinfo

run () {
	#prints the current job
    echo $1
    #$1 actually runs the job
    $1
}
cd ${ETERNABENCH_PATH}/data/Riboswitch

if [ ! -f RS_${2}_Z.json.zip ]; then
		if [ "$2" == "vienna_1" ]; then
			run "python ${ETERNABENCH_PATH}/scripts/CalculateRiboswitchKfold.py ${1} --package ${2} --verbose -o RS --flanking --method Z"
		else
			run "python ${ETERNABENCH_PATH}/scripts/CalculateRiboswitchKfold.py ${1} --package ${2} --parallel --verbose -o RS --flanking --method Z"
		fi
fi

if [ ! -f RS_pearson_Dataset_${2}_Z_BOOTSTRAPS.json.zip ]; then
run "python ${ETERNABENCH_PATH}/scripts/ScoreRiboswitches.py RS_${2}_Z.json.zip --method Z --metric pearson --n_bootstraps 1000 -o RS_pearson_Dataset_${2}_Z"
fi
if [ ! -f RS_RMSE_Dataset_${2}_Z_BOOTSTRAPS.json.zip ]; then
run "python ${ETERNABENCH_PATH}/scripts/ScoreRiboswitches.py RS_${2}_Z.json.zip --method Z --metric RMSE --n_bootstraps 1000 -o RS_RMSE_Dataset_${2}_Z"
fi
if [ ! -f RS_spearman_Dataset_${2}_Z_BOOTSTRAPS.json.zip ]; then
run "python ${ETERNABENCH_PATH}/scripts/ScoreRiboswitches.py RS_${2}_Z.json.zip --method Z --metric spearman --n_bootstraps 1000 -o RS_spearman_Dataset_${2}_Z"
fi

