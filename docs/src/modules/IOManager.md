# IOManager

`cyclicpeptide.IOManager` handles visualization and file export for molecules, sequences, and molecular graphs.

## Key Functions

### `plot_smiles(smiles, output_file='output.svg', w=600, h=600, isdisplay=False)`
Render a 2D depiction of a molecule given its SMILES string.
- **Parameters**:
  - `smiles` (*str*): The SMILES representation of the molecule.
  - `output_file` (*str*): Output SVG file path (default: `'output.svg'`).
  - `w` (*int*): Image width in pixels (default: `600`).
  - `h` (*int*): Image height in pixels (default: `600`).
  - `isdisplay` (*bool*): If `True`, displays SVG within Jupyter notebooks.

```python
from cyclicpeptide import IOManager

smiles = 'C[C@H](N)C(=O)N1CCC[C@H]1C(=O)NCC(=O)O'
IOManager.plot_smiles(smiles, w=300, h=300, isdisplay=True)
```

### `save_mol(molblock, output_file='output.mol')`
Save an MDL Molfile string to a file.

```python
IOManager.save_mol(molblock, 'peptide.mol')
```

### `save_pdb(pdbblock, output_file='output.pdb')`
Save a PDB formatted block string to a file.

```python
IOManager.save_pdb(pdbblock, 'peptide.pdb')
```

### `plot_graph(G, output_file='output.pdf', dpi=60)`
Plot a NetworkX graph representation of a cyclic peptide and save to PDF.

```python
IOManager.plot_graph(G, output_file='graph.pdf')
```

### `graph2svg(G, output_file='output.svg')`
Convert a NetworkX graph to an SVG image.

```python
IOManager.graph2svg(G, output_file='graph.svg')
```
