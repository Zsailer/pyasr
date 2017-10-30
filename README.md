# PyASR

**Ancestral Sequence Reconstruction in Python**

PyASR provides a simple Python interface for Ancestral Sequence Reconstruction (ASR).
Reconstruct ancestral sequences from a phylogenetic tree and multiple sequence alignment. 
Paired with [PhyloPandas](https://github.com/Zsailer/phylopandas), PyASR makes ASR
simple and human readable. 

## Basic Usage

```python
import phylopandas as pd
import dendropy as d
import pyasr

# Use phylopandas to read a set of ancestor.s
df_seqs = pd.read_fasta('test.fasta')

# Use dendropy to read in tree.
tree = d.Tree.get(path='tree.newick', schema='newick')

# Reconstruct nodes in tree.
tree, df_seqs, df_anc = pyasr.reconstruct(df_seqs, tree, working_dir='test', alpha=1.235)

# Write out ancestor dataframe to a CSV file.
df_anc.to_csv('ancestors.csv')
```

We can visualize the ancestors side-by-side with the tree using inside of JupyterLab
thanks to the ToyTree library.

<img src="docs/jlab-example.png" align="middle">

## Install

```
git clone 
cd 
pip install -e .
```

## Dependencies

The following Python dependencies are required for PyASR to work.

- Pandas
- Biopython
- PhyloPandas
- DendroPy
