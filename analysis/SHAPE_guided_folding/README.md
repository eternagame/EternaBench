Scripts and analyses related to SHAPE-guided folding.

Datasets

`SHAPE_guided_all_constructs_6Jun2022.txt`: Contains input sequences, predicted structures, and MCC scores for package options.

Analysis

`analyze_zscores_plot.ipynb`: Plot, perform bootstrapping Z-score analysis

Scripts/inputs to reproduce

`shape_dataset_input.txt`: Clean spreadsheet of input sequences, ground truth structures, reactivity, and associated metadata

`run_SHAPE_guided_folding.py`: Perform SHAPE-guided folding and score MCC to ground truth for dataset in Arnie for RNAstructure, ViennaRNA, and EternaFold package options.

`write_DCP_dataset.py`: example code for extracting reactivity from RDATs.

Raw data

`raw_RDATs`: contains all RDATs used in analysis for Chen 2017, Kappel 2020, "DeepChemicalProfiling" datasets.

`ShapeKnots_DATA`: contains SHAPEknots dataset downloaded from https://weekslab.com/publications/, with a few errors manually corrected.
