# Sequence2Structure

`cyclicpeptide.Sequence2Structure` builds 2D and 3D molecular structures from linear or cyclic peptide sequences.

## Key Functions

### `seq2stru_essentialAA(sequence, cyclic=True)`
Generates a peptide structure composed of essential (canonical) amino acids.
- **Parameters**:
  - `sequence` (*str*): Amino acid sequence in one-letter format (e.g., `'APG'`) or three-letter dash-separated format (e.g., `'Ala-Pro-Gly'`).
  - `cyclic` (*bool*): If `True`, cyclizes the N-to-C terminus (default: `True`).
- **Returns**: `(smiles, peptide_mol)`

```python
from cyclicpeptide import Sequence2Structure, IOManager

smiles, peptide = Sequence2Structure.seq2stru_essentialAA('APG', cyclic=True)
print('SMILES:', smiles)
IOManager.plot_smiles(smiles, w=300, h=300, isdisplay=True)
```

### `seq2stru_no_essentialAA(sequence, references, cyclic=True)`
Generates peptide structures containing non-standard, modified, or D-amino acid residues using monomer definitions.
- **Parameters**:
  - `sequence` (*str*): Monomer codes separated by `--` (e.g., `'Aad--NMe-Ala--4OH-Thr--3Me-Pro'`).
  - `references` (*dict*): Monomer library (obtained from `reference_aa_monomer()`).
  - `cyclic` (*bool*): Whether to cyclize (default: `True`).
- **Returns**: `(smiles, peptide_mol)`

```python
from cyclicpeptide import Sequence2Structure, IOManager

references = Sequence2Structure.reference_aa_monomer(None)
smiles, peptide = Sequence2Structure.seq2stru_no_essentialAA(
    'Aad--NMe-Ala--4OH-Thr--3Me-Pro', references, cyclic=True
)
IOManager.plot_smiles(smiles, w=300, h=300, isdisplay=True)
```

### `reference_aa_monomer(monomers_path=None)`
Loads monomer library from a TSV file. If `monomers_path=None`, the default `states/monomer.tsv` is used.
