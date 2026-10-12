# Structure2Sequence

`cyclicpeptide.Structure2Sequence` deconstructs cyclic peptide molecular structures back into amino acid residue sequences.

## Key Functions

### `detect_order_essentialAA(mol)`
Detects the cyclic peptide backbone for standard amino acid peptides, splits residues, and derives the sequential order of amino acids.
- **Parameters**: `mol` (*rdkit.Chem.Mol*): RDKit Mol object.
- **Returns**: `(m_renum, seq)`

```python
from rdkit import Chem
from cyclicpeptide import Structure2Sequence

m = Chem.MolFromSmiles('C[C@H](N)C(=O)N1CCC[C@H]1C(=O)NCC(=O)O')
m_renum, seq = Structure2Sequence.detect_order_essentialAA(m)
print(seq)
```

### `detect_backbone(mol)`
Detects peptide backbone atoms and returns atom indices.

```python
from cyclicpeptide import Structure2Sequence

backbone_idx = Structure2Sequence.detect_backbone(m)
```

### `order_backbone(mol, backbone_idx)`
Orders and renumbers the backbone atoms sequentially along the peptide cycle.

```python
m_renum, ordered_idx = Structure2Sequence.order_backbone(m, backbone_idx)
```

### `highlight_atom(m, atoms, transparency=1, isdisplay=True)`
Generates an SVG depicting the molecule with specified atom indices highlighted.

```python
from cyclicpeptide import Structure2Sequence

svg = Structure2Sequence.highlight_atom(m, [backbone_idx])
```

### `transform(smiles, monomers_path=None)`
High-level function performing end-to-end residue mapping and generating a full HTML visual summary report.
- **Parameters**:
  - `smiles` (*str*): Input SMILES string.
  - `monomers_path` (*str* or *None*): Path to custom monomer TSV, or `None` / `""` for built-in database.
- **Returns**: HTML string.

```python
from IPython.display import HTML
from cyclicpeptide import Structure2Sequence

html = Structure2Sequence.transform(smiles, monomers_path='')
HTML(html)
```
