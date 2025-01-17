#!/bin/bash
#SBATCH --job-nam=CalculatePunpVectors
#SBATCH -o /home/groups/rhiju/hannahw1/EternaBench/scripts/example_output/%j.out
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
if [ ! -f CM_${2}.json.zip ]; then
   run "python ${ETERNABENCH_PATH}/scripts/CalculatePunpVectors.py ${1} --package ${2} --parallel --verbose -o CM"
fi

for typ in "reactivity" "deg_50C" "deg_Mg_pH10" "deg_Mg_50C" "deg_pH10"; do

if [ ! -f CM_pearson_Dataset_${2}_${typ}_BOOTSTRAPS.json.zip ]; then
run "python ${ETERNABENCH_PATH}/scripts/ScoreChemMapping.py CM_${2}.json.zip --metric pearson --n_bootstraps 1000 -o CM_pearson_Dataset_${2}_${typ} --reactivity_field ${typ}"
fi
done

