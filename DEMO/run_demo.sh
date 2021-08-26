#!/bin/bash

CM_INPUT=${ETERNABENCH_PATH}/data/DEMO_ChemMapping.json.zip
RS_INPUT=${ETERNABENCH_PATH}/data/DEMO_Riboswitch.json.zip

for pkg in vienna_2 contrafold_2 eternafold_B; do

	# Chem Mapping datasets; condensed from scripts/sub_CalculatePunpVectors.sh
	python ${ETERNABENCH_PATH}/scripts/CalculatePunpVectors.py $CM_INPUT --package $pkg --verbose -o CM
	python ${ETERNABENCH_PATH}/scripts/ScoreChemMapping.py CM_${pkg}.json.zip --metric pearson --n_bootstraps 1000 -o CM_pearson_Dataset_${pkg}

	# Riboswitch datasets; condensed from scripts/sub_CalculateRiboswitch_BPS_ONLY.sh
	python ${ETERNABENCH_PATH}/scripts/CalculateRiboswitchKfold.py $RS_INPUT --package $pkg --verbose -o RS --flanking
	python ${ETERNABENCH_PATH}/scripts/ScoreRiboswitches.py RS_${pkg}_bps.json.zip --metric pearson --n_bootstraps 1000 -o RS_pearson_Dataset_${pkg}_bps
done

