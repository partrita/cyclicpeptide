# PropertyAnalysis

`cyclicpeptide.PropertyAnalysis` computes physicochemical properties, drug-likeness filters, and molecular fingerprints for cyclic peptides.

## Key Functions

### `chemial_physical_properties_from_smiles(smiles)`
Computes a comprehensive property dictionary from a SMILES string:
- **Physicochemical**: Exact mass, topological polar surface area (TPSA), Crippen LogP, H-bond donors/acceptors, rotatable bond count, formal charge, molar refractivity, ring count.
- **Rules & Filters**: Lipinski's Rule of 5, Veber's rule, Ghose filter.
- **Fingerprints**: RDKit topological, Daylight-like, Morgan (radius 2), MACCS keys.

```python
from cyclicpeptide import PropertyAnalysis as pa

smiles = 'CC(C1C(=O)NC(C(SSCC(C(=O)NC(C(=O)NC(C(=O)NC(C(=O)N1)CCCN)CC2=CNC3=CC=CC=C32)CC4=CC=C(C=C4)O)NC(=O)C(CC5=CC=CC=C5)N)(C)C)C(=O)NC(C(C)O)C(=O)N)O'
props = pa.chemial_physical_properties_from_smiles(smiles)

print('Exact Mass:', props['Exact_Mass'])
print('Rule of Five:', props['Rule_of_Five'])
print('TPSA:', props['TPSA'])
```

### `calculate_amino_acid_composition(sequence)`
Calculates the relative percentage and counts of constituent amino acids in a sequence.

```python
composition = pa.calculate_amino_acid_composition('ACDEFGHIKLMNPQRSTVWY')
```
