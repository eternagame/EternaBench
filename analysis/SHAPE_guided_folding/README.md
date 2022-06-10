Scripts and analyses related to SHAPE-guided folding.

### Datasets

`SHAPE_guided_all_constructs_6Jun2022.txt`: Contains input sequences, predicted structures, and MCC scores for package options.

### Run SHAPE-directed folding in Arnie

`shape_dataset_input.txt`: Clean spreadsheet of input sequences, ground truth structures, reactivity, and associated metadata

`run_SHAPE_guided_folding.py`: Perform SHAPE-guided folding and score MCC to ground truth for dataset in Arnie for RNAstructure, ViennaRNA, and EternaFold package options.

Usage: 
```
python run_SHAPE_guided_folding.py
```

### Analysis

`analyze_zscores_plot.ipynb`: Plot, perform bootstrapping Z-score analysis

### Raw data

`raw_RDATs`: contains all RDATs used in analysis for datasets from Cheng et al. (2017) PNAS, Kappel et al. (2020) Nat. Methods, "DeepChemicalProfiling" (first published here).

`ShapeKnots_DATA`: contains SHAPEknots dataset downloaded from https://weekslab.com/publications/, with a few errors manually corrected.

`write_DCP_dataset.py`: example code for extracting reactivity from RDATs.
