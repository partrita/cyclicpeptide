# StructureTransformer

`cyclicpeptide.StructureTransformer` converts molecular representations and generates 3D conformations.

## Key Functions

### `output_molecule(mol)`
Converts an RDKit Mol object into common cheminformatics text formats.
- **Returns**: Dictionary with:
  - `'smiles'`: SMILES string
  - `'inchi'`: InChI identifier
  - `'inchikey'`: InChIKey
  - `'molblock'`: MDL Mol block
  - `'pdbblock'`: PDB block

```python
from rdkit import Chem
from cyclicpeptide import StructureTransformer as trans

mol = Chem.MolFromSmiles('N[C@@H](C)C(=O)N1CCC[C@H]1C(=O)')
blocks = trans.output_molecule(mol)
print(blocks['inchikey'])
```

### `predict_3d_conformation(mol)`
Adds hydrogens and generates 3D atomic coordinates using ETKDG embedding.

```python
mol_3d = trans.predict_3d_conformation(mol)
```

### `mol_optimize(mol)`
Optimizes molecular geometry using the Universal Force Field (UFF).

```python
mol_opt = trans.mol_optimize(mol_3d)
```
