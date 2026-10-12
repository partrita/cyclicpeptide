# GraphAlignment

`cyclicpeptide.GraphAlignment` treats cyclic peptides as annotated graphs (nodes = amino acid residues, edges = peptide bonds and cross-links) to compute structural similarities.

## Key Functions

### `create_graph(nodes, edges)`
Constructs a NetworkX graph from sequence node attributes and edge connections.

```python
from cyclicpeptide import GraphAlignment as ga
from cyclicpeptide import SequenceTransformer as st

_, q_nodes, q_edges = st.read_sequence('Ala(1)--Ala--Gly--Phe--Pro--Val--Phe--Phe(1)')
_, r_nodes, r_edges = st.read_sequence('Cys(1)(2)--Cys--OH-DL-Val(2)--4OH-Leu--OH-Ile(1)')

query = ga.create_graph(q_nodes, q_edges)
reference = ga.create_graph(r_nodes, r_edges)
```

### `mcs_similarity(graph1, graph2)`
Computes Maximum Common Subgraph (MCS) similarity between two cyclic peptide graphs.
- **Returns**: `(n_common_nodes, similarity_score)`

```python
n_common, mcs_sim = ga.mcs_similarity(query, reference)
print("MCS Similarity:", mcs_sim)
```

### `graph_similarity(graph1, graph2)`
Computes Graph Edit Distance (GED) similarity between two graphs.
- **Returns**: `(ged, similarity_score)`

```python
ged, ged_sim = ga.graph_similarity(query, reference)
print("GED Similarity:", ged_sim)
```
