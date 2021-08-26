#!/bin/bash
#SBATCH --job-nam=K_bps
#SBATCH -o /home/groups/rhiju/hannahw1/EternaBench/scripts/example_output/sub_Ribo_bps.out
#SBATCH -p biochem,biochem
#SBATCH -n 1
#SBATCH --mail-user=hannahw1@stanford.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --mem-per-cpu=30G
#SBATCH --cpus-per-task=6
#SBATCH --time=40:00:00

source ~/.bash_profile
grep -m 1 'cpu cores' /proc/cpuinfo

run () {
	#prints the current job
    echo $1
    #$1 actually runs the job
    $1
}
cd $ETERNABENCH_PATH/data

python ${ETERNABENCH_PATH}/scripts/CalculateRiboswitchKfold.py EternaBench_Riboswitch_Filtered_12Jul2021.json.zip --parallel --flanking -o EternaBench_Riboswitch_bps -v