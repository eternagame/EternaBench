#!/bin/bash

CM_INPUT=${ETERNABENCH_PATH}/data/DEMO_ChemMapping.json.zip
RS_INPUT=${ETERNABENCH_PATH}/data/DEMO_Riboswitch.json.zip

while read pkg; do

	# Chem Mapping datasets; condensed from cluster_scripts/sub_CalculatePunpVectors.sh
	if [ ! -f CM_${pkg}.json.zip ]; then
	python ${ETERNABENCH_PATH}/scripts/CalculatePunpVectors.py $CM_INPUT --package $pkg --verbose -o CM
	fi
	if [ ! -f CM_pearson_Dataset_${pkg}_BOOTSTRAPS.json.zip ]; then
	python ${ETERNABENCH_PATH}/scripts/ScoreChemMapping.py CM_${pkg}.json.zip --metric pearson --n_bootstraps 1000 -o CM_pearson_Dataset_${pkg}
	fi

	# Riboswitch datasets; condensed from cluster_scripts/sub_CalculateRiboswitch_BPS_ONLY.sh
	if [ ! -f RS_${pkg}_bps.json.zip ]; then
	python ${ETERNABENCH_PATH}/scripts/CalculateRiboswitchKfold.py $RS_INPUT --package $pkg --verbose -o RS --flanking
	fi

	if [ ! -f RS_pearson_Dataset_${pkg}_bps_BOOTSTRAPS.json.zip ]; then
	python ${ETERNABENCH_PATH}/scripts/ScoreRiboswitches.py RS_${pkg}_bps.json.zip --metric pearson --n_bootstraps 1000 -o RS_pearson_Dataset_${pkg}_bps
	fi

done<package_list.txt

python calculateZscoreDEMO.py
