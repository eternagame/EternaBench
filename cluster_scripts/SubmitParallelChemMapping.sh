#!/bin/bash

INPUT=${ETERNABENCH_PATH}/data/EternaBench_ChemMapping_Filtered_10Jul2021.json.zip

while read package; do
    echo ${package}
    sbatch ${ETERNABENCH_PATH}/cluster_scripts/sub_CalculatePunpVectors.sh ${INPUT} ${package}
done<${ETERNABENCH_PATH}/cluster_scripts/package_list.txt
