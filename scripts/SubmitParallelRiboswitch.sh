#!/bin/bash

INPUT=${ETERNABENCH_PATH}/data/EternaBench_Riboswitch_Filtered_07Aug2021.json.zip

while read package; do
    echo ${package}
    sbatch ${ETERNABENCH_PATH}/scripts/sub_CalculateRiboswitch_BPS_ONLY.sh ${INPUT} ${package}
done<${ETERNABENCH_PATH}/scripts/package_list.txt
#done<${ETERNABENCH_PATH}/scripts/rnasoft_pkgs.txt
while read package; do
    echo ${package}
    sbatch ${ETERNABENCH_PATH}/scripts/sub_CalculateRiboswitch_Z_ONLY.sh ${INPUT} ${package}
done<${ETERNABENCH_PATH}/scripts/packages_for_Z.txt
