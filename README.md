[![python-version](https://img.shields.io/badge/python-3.7+-blue.svg)](https://www.python.org/downloads/release/python-360/)
[![DOI](https://zenodo.org/badge/DOI/10.1016/j.immuni.2023.07.016)](https://doi.org/10.1016/j.immuni.2023.07.016))


# TIMON
Transcription Factor Interaction Inference from Motif Co-occurrence Networks (TIMON) is a tool to analyze transcription factor (TF) motif co-occurrences from epigenetic datasets (open chromatin, transcription factor binding sites, enhancers, etc.). The first step of TIMON is to identify non-overlapping motifs to identify the best matching TFs in the sequences of interest.  Next, TIMON builds co-occurrences matrices between TFs. Through comparing the TF co-occurrences in the cell type of interest with background cell types (50+ cell types from ENCODE), TIMON identifies significantly enriched TF pair in the cell type of interest and thus construct the TF co-occurrence networks. 

For questions on installation or usage, please open an issue, submit a pull request, or contact Rick Z. Li (zhl022@ucsd.edu)
<p align="center">
<img src="https://github.com/rzzli/TIMON/blob/main/image/TIMON.jpg" width="900" height="512">
</p>

### Installation
```bash
git clone https://github.com/rzzli/TIMON.git
cd TIMON
pip install -e .
```
#### optional: install in a new conda environment
```bash
conda create -n timon_conda python=3.7.7
conda activate timon_conda
git clone https://github.com/rzzli/TIMON.git
cd TIMON
pip install -e .
```
add the new conda environment to jupyter 
```bash
conda activate monet_conda
python -m ipykernel install --user --name monet_conda
```
The tutorial requires PyGraphviz package, to install:
```bash
conda install --channel conda-forge pygraphviz
```

### Uninstall TIMON/conda environment
uninstall TIMON package
```bash
pip uninstall TIMON
```
remove the conda environment
```bash
conda env remove -n timon_conda
```

remove the conda environment on jupyter notebook
```bash
jupyter kernelspec uninstall monet_conda
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
