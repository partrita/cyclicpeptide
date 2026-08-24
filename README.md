# cyclicpeptide: A Python Package for Cyclic Peptide Drug Design

Detailed documentation for cyclicpeptide is available at https://dfwlab.github.io/cyclicpeptide/.

## Overview

**cyclicpeptide** provides computational building blocks for cyclic peptide drug design:

- Convert cyclic peptide **sequences to structures** (SMILES/mol) and **structures back to sequences**, including peptides built from non-canonical / modified amino acids;
- Interconvert between common sequence notations (**IUPAC condensed**, amino acid chain, graph presentation, one-letter code);
- Compute **molecular properties**, drug-likeness rules, and **molecular fingerprints**;
- Measure cyclic peptide **similarity** via exact graph alignment or a pretrained GCN model.

## Installation

Install from PyPI with pip:

    pip install cyclicpeptide

Or set up a development environment from source with [uv](https://docs.astral.sh/uv/):

```bash
git clone https://github.com/Willow0316/cyclicpeptide.git
cd cyclicpeptide
uv sync        # creates .venv, installs pinned dependencies, generates uv.lock
uv run python  # run anything inside the managed environment
```

Requires Python >=3.9,<3.12.

## Modules

| Module | Description |
|--------|-------------|
| `PropertyAnalysis` | Chemical/physical properties, drug-likeness rules, molecular fingerprints |
| `SequenceTransformer` | Sequence format detection and interconversion |
| `Sequence2Structure` | Build peptide molecules (mol objects) from sequences |
| `Structure2Sequence` | Identify amino acid units in a structure and derive sequence(s) |
| `StructureTransformer` | Structure format interconversion and 3D conformation handling |
| `GraphAlignment` | Graph-based cyclic peptide alignment and similarity |
| `SequenceGeneration` | Generate similar sequences by residue substitution rules |
| `IOManager` | Save/plot molecules and graphs (SVG, PDF, MOL, PDB) |

### PropertyAnalysis

`chemial_physical_properties_from_smiles` computes properties such as exact mass, topological polar surface area, Crippen LogP, H-bond donor/acceptor counts, rotatable bonds, formal charge, refractivity and ring count; evaluates **Lipinski's Rule of Five**, **Veber's rule**, and the **Ghose filter**; and returns four fingerprints: RDKit, Daylight-like topological, Morgan (radius 2), and MACCS keys.

```python
from cyclicpeptide import PropertyAnalysis as pa

smiles = 'CC(C1C(=O)NC(C(SSCC(C(=O)NC(C(=O)NC(C(=O)NC(C(=O)N1)CCCN)CC2=CNC3=CC=CC=C32)CC4=CC=C(C=C4)O)NC(=O)C(CC5=CC=CC=C5)N)(C)C)C(=O)NC(C(C)O)C(=O)N)O'
props = pa.chemial_physical_properties_from_smiles(smiles)
print(props['Exact_Mass'], props['Rule_of_Five'])
```

Amino acid composition of a chain can be obtained with `calculate_amino_acid_composition(sequence)`.

### SequenceTransformer

`read_sequence` auto-detects the input notation and returns `(seq_format, nodes, edges)`; supported formats are *graph presentation*, *IUPAC condensed*, *amino acid chain*, and *one-letter code*. `create_sequence` converts nodes/edges back into all formats at once.

```python
from cyclicpeptide import SequenceTransformer as st

seq_format, nodes, edges = st.read_sequence(
    'aThr,Tyr,dhAbu,bOH-Gln,Gly,Gln,His,Dab,C13:2(t4.t6)-OH(2.3),Lyx,dhAbu @1,5 @6,10 @0,8')

seq_list = st.create_sequence(nodes, edges)
# {'iupac_condensed': ..., 'amino_acid_chain': ..., 'graph_presentation': ..., 'one_letter_peptide': ...}
```

### Sequence2Structure

Build peptide mol objects from sequences. Essential-amino-acid sequences accept one-letter codes or three-letter chains; non-essential residues are assembled from a monomer reference library (`states/monomer.tsv`).

```python
from cyclicpeptide import Sequence2Structure as s2s
from cyclicpeptide import IOManager

smiles, peptide = s2s.seq2stru_essentialAA('Ala-Ala-Cys-Asp', cyclic=True)
IOManager.plot_smiles(smiles, w=300, h=300, isdisplay=True)

references = s2s.reference_aa_monomer(None)  # default monomer library
smiles, peptide = s2s.seq2stru_no_essentialAA('Aad--NMe-Ala--4OH-Thr--3Me-Pro', references, cyclic=True)
```

### Structure2Sequence

Convert a cyclic peptide SMILES back into sequence information: detect the backbone, split amino acid units, match them against the monomer reference library, and enumerate possible chains. `transform` produces a complete HTML report (`output.html`) with structure highlights, per-residue structures, mappings, and locations.

```python
from IPython.display import HTML
from cyclicpeptide import Structure2Sequence as struc2seq

html = struc2seq.transform(smiles, monomers_path=None)  # default monomer library
HTML(html)
```

### StructureTransformer

Interconvert SMILES / InChI / InChIKey / Molblock / PDBblock and handle 3D conformations (hydrogen addition + ETKDG embedding, UFF optimization).

```python
from rdkit import Chem
from cyclicpeptide import StructureTransformer as trans

mol = Chem.MolFromSmiles('N[C@@H](C)C(=O)N1CCC[C@H]1C(=O)')
blocks = trans.output_molecule(mol)   # {'smiles', 'inchi', 'inchikey', 'molblock', 'pdbblock'}
mol_3d = trans.predict_3d_conformation(mol)
mol_3d = trans.mol_optimize(mol_3d)
```

### GraphAlignment

Represent a cyclic peptide as a graph (nodes = amino acids, edges = backbone and cross-links) and compute similarities: maximum common subgraph (MCS) similarity, graph-edit-distance similarity, and composition-frequency similarity.

```python
from cyclicpeptide import GraphAlignment as ga
from cyclicpeptide import SequenceTransformer as st

_, q_nodes, q_edges = st.read_sequence('Ala(1)--Ala--Gly--Phe--Pro--Val--Phe--Phe(1)')
_, r_nodes, r_edges = st.read_sequence('Cys(1)(2)--Cys--OH-DL-Val(2)--4OH-Leu--OH-Ile(1)')

query = ga.create_graph(q_nodes, q_edges)
reference = ga.create_graph(r_nodes, r_edges)

n_common, mcs_sim = ga.mcs_similarity(query, reference)
ged, ged_sim = ga.graph_similarity(query, reference)
```

### SequenceGeneration

Generate sequences similar to an input by substituting residues under property-based default rules (or a custom rule file), limited by a replacement ratio, and save results to CSV.

```python
from cyclicpeptide import SequenceGeneration as sg

similar = sg.generate_similar_sequences('Glu--Lys--Leu--Phe', Replacement_ratio=0.25)
sg.save_sequences_to_csv(similar, 'similar_sequences.csv')
```

### IOManager

Plot/save molecules and graphs:

```python
from cyclicpeptide import IOManager as io_m

io_m.plot_smiles(smiles, output_file='peptide.svg', w=600, h=600)   # molecule SVG
io_m.save_mol(molblock, 'peptide.mol')                              # MDL molfile
io_m.save_pdb(pdbblock, 'peptide.pdb')                              # PDB file
io_m.plot_graph(G, output_file='graph.pdf')                         # networkx graph -> PDF
io_m.graph2svg(G, output_file='graph.svg')                          # networkx graph -> SVG
```

## Data Files (`states/`)

| File | Content |
|------|---------|
| `AminoAcids.txt` | Amino acid reference tuples (name, one-letter code, ...) |
| `aas.txt` | Essential amino acid names (order matters for substructure matching) |
| `aa_smiles.txt` | `NAME:SMILES` mapping for essential amino acids |
| `monomer.tsv` | Default monomer reference library (code, SMILES, weight, symbol, flags) |
| `Chemical_P.csv` | Property name/description table used in the documentation |

## Graph Similarity Model (`GCN/`)

A GCN-based graph-alignment similarity model trained on cyclic peptide graphs. See [`GCN/readme.md`](GCN/readme.md).

## License

MIT — see [LICENSE.md](LICENSE.md).
