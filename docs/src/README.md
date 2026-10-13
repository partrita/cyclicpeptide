# cyclicpeptide: A Python Package for Cyclic Peptide Drug Design

Welcome to the documentation for **cyclicpeptide**.

## Overview

**cyclicpeptide** provides computational building blocks for cyclic peptide drug design:

- Convert cyclic peptide **sequences to structures** (SMILES/mol) and **structures back to sequences**, including peptides built from non-canonical / modified amino acids;
- Interconvert between common sequence notations (**IUPAC condensed**, amino acid chain, graph presentation, one-letter code);
- Compute **molecular properties**, drug-likeness rules, and **molecular fingerprints**;
- Measure cyclic peptide **similarity** via exact graph alignment or a pretrained GCN model.

## Installation

Clone the repository and install with [uv](https://docs.astral.sh/uv/) (recommended):

```bash
git clone https://github.com/partrita/cyclicpeptide.git
cd cyclicpeptide
uv sync        # creates .venv, installs pinned dependencies, generates uv.lock
uv run python  # run anything inside the managed environment
```

> **Requirements**: Python >=3.9, <3.12.

## Quick Tour of Modules

| Module | Description |
| --- | --- |
| [IOManager](modules/IOManager.md) | Save/plot molecules and graphs (SVG, PDF, MOL, PDB) |
| [Sequence2Structure](modules/Sequence2Structure.md) | Build peptide molecules (mol objects) from sequences |
| [Structure2Sequence](modules/Structure2Sequence.md) | Identify amino acid units in a structure and derive sequence(s) |
| [SequenceTransformer](modules/SequenceTransformer.md) | Sequence format detection and interconversion |
| [StructureTransformer](modules/StructureTransformer.md) | Structure format interconversion and 3D conformation handling |
| [PropertyAnalysis](modules/PropertyAnalysis.md) | Chemical/physical properties, drug-likeness rules, fingerprints |
| [GraphAlignment](modules/GraphAlignment.md) | Graph-based cyclic peptide alignment and similarity |
| [SequenceGeneration](modules/SequenceGeneration.md) | Generate similar sequences by residue substitution rules |
| [Data & Models](data_and_models.md) | Monomer libraries and pretrained GCN similarity model |
