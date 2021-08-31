!/bin/bash

for i in 300 600 900 1200; do 
	INPUT=${ETERNABENCH_PATH}/data/ExternalData/ExternalData_window${i}_uniq.json.zip

	while read package; do
	    echo ${package}
	    sbatch ${ETERNABENCH_PATH}/cluster_scripts/sub_Ext_benchmarks.sh ${INPUT} ${package}
	done<${ETERNABENCH_PATH}/cluster_scripts/package_subset.txt
done

INPUT=${ETERNABENCH_PATH}/data/ExternalData/RYOS_FULL_11Aug2021.json.zip

while read package; do
    echo ${package}
    sbatch ${ETERNABENCH_PATH}/cluster_scripts/sub_RYOS.sh ${INPUT} ${package}
done<${ETERNABENCH_PATH}/cluster_scripts/package_subset.txt
