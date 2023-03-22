[![python-version](https://img.shields.io/badge/python-3.7+-blue.svg)](https://www.python.org/downloads/release/python-360/)

# TIMON 
Motif Co-Occurrence Network analysis (TIMON). is a tool to analyze motif co-occurrences in a set of peaks (ATAC, enhancer etc.). In addition, when compare the input data with background distal elements (summing 15 ENCODE cell types), Moe can identify significantly co-occuring motif pairs from the input peak set.

For questions on installation or usage, please open an issue, submit a pull request, or contact Rick Z. Li (zhl022@eng.ucsd.edu).

### Installation
```bash
git clone https://github.com/rzzli/TIMON.git
cd TIMON
pip install -e .
```
#### optional: install in a new conda environment
```bash
conda create -n TIMON_conda python=3.7.7
conda activate TIMON_conda
git clone https://github.com/rzzli/TIMON.git
cd TIMON
pip install -e .
```
add the new conda environment to jupyter notebook
```bash
conda activate TIMON_conda
python -m ipykernel install --user --name TIMON_conda
```

### Uninstall TIMON/conda environment
uninstall TIMON package
```bash
pip uninstall TIMON
```
remove the conda environment
```bash
conda env remove -n TIMON_conda
```

remove the conda environment on jupyter notebook
```bash
jupyter kernelspec uninstall TIMON_conda
```

### import TIMON 
```python
from TIMON.motifFreq import *
```

### Tutorials

This section is actively under development.

- [Fetal and Postnatal microglia motif co-occurrence analysis.](tutorials/microglia_cooccurrence.ipynb)

### License
[This project is licensed under MIT](https://github.com/rzzli/TIMON/blob/main/LICENSE)

### Contribution
TIMON was developed primarily by Rick Z Li, with contributions and suggestions from Claudia Han and supervision from Christopher Glass.
