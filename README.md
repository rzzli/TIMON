[![python-version](https://img.shields.io/badge/python-3.7+-blue.svg)](https://www.python.org/downloads/release/python-360/)

# MONet 
Motif Co-Occurrence Network analysis (MONet) is a tool to analyze transcription factor (TF) motif co-occurrences from epigenetic datasets (open chromatin, transcription factor binding sites, enhancers, etc.). The first step of MONet is to identify non-overlapping motifs to identify the best matching TFs in the sequences of interest.  Next, MONet builds co-occurrences matrices between TFs. Through comparing the TF co-occurrences in the cell type of interest with background cell types (50+ cell types from ENCODE), MONet identifies significantly enriched TF pair in the cell type of interest and thus construct the TF co-occurrence networks. 

For questions on installation or usage, please open an issue, submit a pull request, or contact Rick Z. Li (zhl022@eng.ucsd.edu)
<p align="center">
<img src="https://github.com/rzzli/MONet/blob/main/image/monet.jpg" width="900" height="512">
</p>

### Installation
```bash
git clone https://github.com/rzzli/MONet.git
cd MONet
pip install -e .
```
#### optional: install in a new conda environment
```bash
conda create -n monet_conda python=3.7.7
conda activate monet_conda
git clone https://github.com/rzzli/MONet.git
cd MONet
pip install -e .
```
add the new conda environment to jupyter notebook
```bash
conda activate monet_conda
python -m ipykernel install --user --name monet_conda
```

### Uninstall MONet/conda environment
uninstall MONet package
```bash
pip uninstall MONet
```
remove the conda environment
```bash
conda env remove -n monet_conda
```

remove the conda environment on jupyter notebook
```bash
jupyter kernelspec uninstall monet_conda
```

### import MONet 
```python
from MONet.motifFreq import *
```

### Tutorials

This section is actively under development.

- [Fetal and Postnatal microglia motif co-occurrence analysis.](tutorials/microglia_cooccurrence.ipynb)

### License
[This project is licensed under MIT](https://github.com/rzzli/MONet/blob/main/LICENSE)

### Contribution
MONet was developed primarily by Rick Z Li, with contributions and suggestions from Claudia Han and supervision from Christopher Glass.
