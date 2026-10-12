# SequenceTransformer

`cyclicpeptide.SequenceTransformer` detects, parses, and converts peptide sequence representations between standard notations.

## Supported Formats

- **IUPAC condensed**: `cyclo(Ala-Tyr-Gly-...)`
- **Amino acid chain**: `Ala(1)--Tyr--Gly--Val(1)`
- **Graph presentation**: `Ala,Tyr,Gly @0,2`
- **One-letter peptide**: `AYG...`

## Key Functions

### `read_sequence(sequence)`
Auto-detects the sequence format and extracts graph components (nodes and connectivity edges).
- **Parameters**: `sequence` (*str*): Peptide sequence in any supported format.
- **Returns**: `(seq_format, nodes, edges)`

```python
from cyclicpeptide import SequenceTransformer as st

seq_format, nodes, edges = st.read_sequence(
    'aThr,Tyr,dhAbu,bOH-Gln,Gly,Gln,His,Dab,C13:2(t4.t6)-OH(2.3),Lyx,dhAbu @1,5 @6,10 @0,8'
)
print("Detected format:", seq_format)
```

### `create_sequence(nodes, edges)`
Converts graph nodes and edges into all four supported string notations.
- **Parameters**: `nodes` (*list*), `edges` (*list*)
- **Returns**: Dictionary with keys:
  - `'iupac_condensed'`
  - `'amino_acid_chain'`
  - `'graph_presentation'`
  - `'one_letter_peptide'`

```python
formats = st.create_sequence(nodes, edges)
print(formats['iupac_condensed'])
print(formats['amino_acid_chain'])
```
