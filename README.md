# Organization:

This repository contains the EternaBench datasets as well as scripts to reproduce the analysis presented in Wayment-Steele et al. (2021).

`data`: datasets

`scripts`: scripts to regenerate benchmark. Full documentation is in `docs/RunBenchmarkREADME.md`.

`analysis`: python notebooks to reproduce the figures in (Wayment-Steele, 2021).

`eternabench`: EternaBench API source.

# Quick Use

Add to python path and point to datasets by adding to .bashrc:
```bash
export PYTHONPATH=/path/to/EternaBench
export ETERNABENCH_PATH=/path/to/EternaBench
```

Load a dataset using

```python
import eternabench as eb
data = eb.load_CM_data()
```

# Data Origin

- Chemical Mapping RDAT files may be downloaded from www.rmdb.stanford.edu.

- Eterna riboswitch datasets are detailed in the supporting information of 

Andreasson, J. O., ... & Das, R., Greenleaf, W. J. (2019). Crowdsourced RNA design discovers diverse, reversible, efficient, self-contained molecular sensors. bioRxiv, 877183.

- Ribologic riboswitch dataset is detailed in the supporting information of

Wu, M. J., Andreasson, J. O., Kladwang, W., Greenleaf, W., & Das, R. (2019). Automated design of diverse stand-alone riboswitches. ACS synthetic biology, 8(8), 1838-1846.

