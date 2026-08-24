# GA-GCN: GCN-Based Graph Alignment Similarity Model

This folder contains a graph neural network model that predicts the **similarity between two cyclic peptide graphs**, serving as a learned alternative to exact graph alignment (`cyclicpeptide.GraphAlignment`).

## Contents

| File | Description |
|------|-------------|
| `GA_GCN_modeling.ipynb` | Model definition and training pipeline (PyTorch Lightning) |
| `GA_GCN_usage.ipynb` | Inference utilities: `load_ga_gcn`, `ga_prediction`, sequence-to-graph conversion |
| `GA_GCN_metric.ipynb` | Evaluation of the trained model |
| `GA_GCN.pth` | Pretrained model weights |

## Architecture

- **Input**: cyclic peptide graphs built from amino-acid-chain sequences via `sequence_to_node_edge` (nodes = amino acids, edges = backbone and cross-link connections).
- **Node features**: one-hot encoding over 28 residue labels (20 essential AAs plus common non-canonical ones such as `Orn`, `Dap`, `DL`, `D`, `4OH`, ...) with an extra "unknown" flag (29 dimensions in total).
- **Encoder** (`GCNEmbedding`): two `GCNConv` layers (29 → 64 → 32) with ReLU, followed by global mean pooling over nodes.
- **Similarity head** (`GraphSimilarityModel`): the embeddings of two graphs are combined as `[h1, h2, |h1 − h2|]` and passed through an MLP (96 → 64 → 1) with a sigmoid output, yielding a similarity score in [0, 1].
- **Training**: MSE loss against alignment-based similarity targets (see the modeling notebook).

## Usage

```python
import torch
# define/copy GraphSimilarityModel and helpers from GA_GCN_usage.ipynb first

model = load_ga_gcn(model_path='GA_GCN.pth')

seq1 = 'Ala(1)--Ala--Gly--Phe--Pro--Val--Phe--Phe(1)'
seq2 = 'Pro(1)--Val--Phe--Phe--Ala--Ala--Gly--Phe(1)'   # same ring, rotated
print(ga_prediction(model, seq1, seq2))
```

To score a query against a reference library in batch:

```python
similarities = [ga_prediction(model, query, ref) for ref in ref_seqs]
```

## Dependencies

The notebooks require packages that are **not** part of the core `cyclicpeptide` dependencies:

```bash
uv pip install torch torch_geometric pytorch_lightning torchmetrics pandas openpyxl
```
