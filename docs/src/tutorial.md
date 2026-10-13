# Tutorial

This tutorial walks through building, visualizing, and analyzing cyclic peptides step by step.

## 1. Install Package

Clone the repository and install with [uv](https://docs.astral.sh/uv/):

```bash
git clone https://github.com/partrita/cyclicpeptide.git
cd cyclicpeptide
uv sync
```

Then activate the environment and run Python:

```bash
uv run python
```

In Python:

## 2. Create a Cyclic Peptide by Sequence2Structure

Convert an essential amino acid sequence to a cyclic structure and visualize it using `IOManager`:

```python
from cyclicpeptide import Sequence2Structure, IOManager

# Generate cyclic peptide SMILES and RDKit mol object
smiles, peptide = Sequence2Structure.seq2stru_essentialAA(sequence='APG', cyclic=True)
print('SMILES:', smiles)

# Plot SMILES structure using IOManager
IOManager.plot_smiles(smiles, w=300, h=300, isdisplay=True)
```

## 3. Structure2Sequence for Essential Amino Acid Peptides

Derive peptide sequences from a molecular structure:

```python
from rdkit import Chem
from cyclicpeptide import Structure2Sequence

m = Chem.MolFromSmiles(smiles)
m_renum, seq = Structure2Sequence.detect_order_essentialAA(m)
print('Detected Sequence:', seq)
```

### Interactive HTML Report with Structure2Sequence.transform

`transform` creates an interactive HTML report containing structure highlights and residue breakdown:

```python
from IPython.display import HTML
from cyclicpeptide import Structure2Sequence

html = Structure2Sequence.transform(smiles, monomers_path='') # use default monomer DB
HTML(html)
```

## 4. Cyclic Peptide with Non-Natural Amino Acids or Custom Components

Generate and inspect structures with non-canonical monomers:

```python
from cyclicpeptide import IOManager

smiles = 'C[C@@H]1NC(=O)CNC(=O)[C@@H]2CC(CCCCC)CN2C1=O'
IOManager.plot_smiles(smiles, w=300, h=300, isdisplay=True)
```

## 5. Sequence Extraction for Custom Cyclic Peptides

Extract sequence components:

```python
from IPython.display import HTML
from cyclicpeptide import Structure2Sequence

html = Structure2Sequence.transform(smiles, monomers_path='')
HTML(html)
```

## 6. Extending the Monomer Database

You can extend the monomer database with new residue definitions:

```python
import pandas as pd
from cyclicpeptide import IOManager, Structure2Sequence

monomers_path = 'cyclicpeptide/states/monomer.tsv'
monomer_db = pd.read_csv(monomers_path, sep='\t', index_col=0)

# Inspect an existing monomer (e.g., Proline)
print(monomer_db.loc[monomer_db['Code'] == 'Pro', :])

# Create a customized monomer entry
new_monomer = monomer_db.loc[monomer_db['Code'] == 'Pro', :].values[0].copy()
new_monomer[0] = 'Pro**'
new_monomer[1] = 'Proline*'

# Preview monomer structure
IOManager.plot_smiles(new_monomer[6], w=300, h=300, isdisplay=True)

# Append to DB and save
monomer_db.loc[len(monomer_db), :] = new_monomer
monomer_db.to_csv('monomer_new.tsv', sep='\t')

# Run Structure2Sequence using the custom monomer database
smiles = 'C[C@@H]1NC(=O)CNC(=O)[C@@H]2CC(CCCCC)CN2C1=O'
html = Structure2Sequence.transform(smiles, monomers_path='monomer_new.tsv')
```
