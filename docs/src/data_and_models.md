# Data Files & Pretrained Models

## Data Files (`states/`)

The `states/` directory supplies reference data used across `cyclicpeptide`:

| File | Content |
| --- | --- |
| `AminoAcids.txt` | Standard amino acid reference definitions (name, symbol, SMILES). |
| `aas.txt` | Essential amino acid names used for substructure ordering. |
| `aa_smiles.txt` | Mapping of essential amino acid names to SMILES representations. |
| `monomer.tsv` | Comprehensive monomer library (Code, SMILES, MW, Symbol, Flags). |
| `Chemical_P.csv` | Reference table describing physicochemical property descriptors. |

## Graph Neural Network Model (`GCN/`)

The repository includes a Graph Convolutional Network (GCN) model designed for learning structural representations of cyclic peptides and computing graph-level similarities.

- Model architecture and checkpoint files are located in `cyclicpeptide/GCN/`.
- See [`cyclicpeptide/GCN/readme.md`](https://github.com/partrita/cyclicpeptide/blob/main/cyclicpeptide/GCN/readme.md) for details on training and inference.
