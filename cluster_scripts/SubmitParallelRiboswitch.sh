#!/bin/bash

INPUT=${ETERNABENCH_PATH}/data/EternaBench_Riboswitch_Filtered_07Aug2021.json.zip

while read package; do
    echo ${package}
    sbatch ${ETERNABENCH_PATH}/cluster_scripts/sub_CalculateRiboswitch_BPS_ONLY.sh ${INPUT} ${package}
done<${ETERNABENCH_PATH}/cluster_scripts/package_list.txt

while read package; do
    echo ${package}
    sbatch ${ETERNABENCH_PATH}/cluster_scripts/sub_CalculateRiboswitch_Z_ONLY.sh ${INPUT} ${package}
done<${ETERNABENCH_PATH}/cluster_scripts/packages_for_Z.txt
