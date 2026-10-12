# SequenceGeneration

`cyclicpeptide.SequenceGeneration` generates variant libraries of cyclic peptide sequences based on residue mutation and substitution rules.

## Key Functions

### `generate_similar_sequences(sequence, Replacement_ratio=0.25, rules=None)`
Generates modified sequences by substituting residues according to property-based substitution rules.
- **Parameters**:
  - `sequence` (*str*): Input sequence in amino acid chain format (e.g., `'Glu--Lys--Leu--Phe'`).
  - `Replacement_ratio` (*float*): Fraction of residues eligible for replacement (default: `0.25`).
  - `rules` (*dict* or *None*): Custom substitution matrix/dictionary. If `None`, built-in rules are used.
- **Returns**: List of generated similar sequence strings.

```python
from cyclicpeptide import SequenceGeneration as sg

similar = sg.generate_similar_sequences('Glu--Lys--Leu--Phe', Replacement_ratio=0.25)
print(f"Generated {len(similar)} variant sequences.")
```

### `save_sequences_to_csv(sequences, output_file)`
Saves generated sequences to a CSV file.

```python
sg.save_sequences_to_csv(similar, 'similar_sequences.csv')
```
