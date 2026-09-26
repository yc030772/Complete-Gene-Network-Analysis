# Gene Network & Biological Enrichment

Gene interaction network analysis that ranks genes by how much their interaction pattern changes between CADM1-positive and CADM1-negative tumor samples (329 ANOVA-selected genes; 27 vs 28 patients). Full write-up, including STRING/KEGG enrichment: [report.pdf](report.pdf).

## Quick start

```bash
python3 -m venv .venv && .venv/bin/pip install -r requirements.txt
.venv/bin/python gene_network.py AIIM_HW2_ANOVA_grouped_CADM1.npz
```

The course dataset `AIIM_HW2_ANOVA_grouped_CADM1.npz` is **not included**. It must contain `x_pos` (samples × genes), `x_neg`, and `ANOVA_genes`. Place it in the repo root, or pass its path as the first argument.

Tests (synthetic data, no dataset needed): `.venv/bin/python test_gene_network.py`

## Method

**Interaction network.** Each gene's expression is modeled as a linear combination of all other genes:

$$x_{i}[n]=\sum_{j \neq i} a_{ij}x_{j}[n]+\epsilon_{i}[n]$$

$A$ is solved row by row with Ridge regression ($\alpha=1$, no intercept), so $a_{ii}=0$. This is done separately for each group, giving $A^{(+)}$ and $A^{(-)}$.

**Relevance Value.** How much gene $j$'s influence on the network changes between groups:

$$RV_j = \sum_{i} |a_{ij}^{(+)} - a_{ij}^{(-)}|$$

## Top 15 genes

| Rank | Gene | RV | | Rank | Gene | RV |
|---|---|---|---|---|---|---|
| 1 | PLAC8 | 4.9781 | | 9 | LOC389765 | 3.8162 |
| 2 | MMP9 | 4.8135 | | 10 | S100A16 | 3.7446 |
| 3 | FAM198B | 4.5533 | | 11 | NDUFAF2 | 3.6928 |
| 4 | PRICKLE1 | 4.0398 | | 12 | SCIN | 3.5561 |
| 5 | ZFP3 | 4.0019 | | 13 | ANKRD12 | 3.5531 |
| 6 | SNORD89 | 3.9309 | | 14 | HOXB2 | 3.4953 |
| 7 | LOC202181 | 3.9125 | | 15 | ADAMTS5 | 3.4848 |
| 8 | ZNF383 | 3.8391 | | | | |

## Biological findings

STRING network of the top-15 RV genes plus the provided biomarkers is enriched for:

- **Central carbon metabolism in cancer** (KEGG hsa05230): `PKM` and `SLC2A1` drive the Warburg effect.
- **HIF-1 signaling**: hypoxia response promoting angiogenesis and metastasis.
- **Cell adhesion / ECM remodeling**: `CADM1`, `ALCAM`, `EPCAM` junction loss plus `MMP9` and `ADAMTS5` matrix degradation.

A second network built from CADM1's physical interactors (BioGRID) adds renal cell carcinoma (hsa05211) and thyroid hormone signaling (hsa04919). See [report.pdf](report.pdf) for details.

## Files

| File | Purpose |
|---|---|
| `gene_network.py` | The analysis (network solve + RV ranking) |
| `network_analysis.ipynb` | Notebook entry point; calls `gene_network.main()` |
| `test_gene_network.py` | Sanity checks on synthetic data |
| `report.pdf` | Full homework report (AIIM 114-1, HW2) |
