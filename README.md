# Gene Network & Biological Enrichment

## 1. Project Overview
This project implements a gene interaction network analysis to identify key biomarkers distinguishing different cancer-related states (CADM1 positive vs. negative groups).

## 2. Methodology
### Interaction Network Modeling
The interaction between genes is modeled as:
$$x_{i}[n]=\sum_{j \in G, j \neq i} a_{ij}x_{j}[n]+\epsilon_{i}[n]$$
We solved the interaction matrix $A$ using **Ridge Regression** ($\alpha=1.0$) to handle potential multicollinearity in high-dimensional genetic data.

### Relevance Value (RV)
To identify influential genes, we calculated the **Relevance Value**:
$$RV_j = \sum_{i} |a_{ij}^{(+)} - a_{ij}^{(-)}|$$
This identifies genes whose interaction patterns change most significantly between groups.

## 3. Top 15 Genes Identified
| Rank | Gene | RV |
|---|---|---|
| 1 | PLAC8 | 4.9781 |
| 2 | MMP9 | 4.8135 |
| ... | ... | ... |

## 4. Biological Findings
- **Metabolism:** Identified key roles of `PKM` and `SLC2A1` in the Warburg effect.
- **Microenvironment:** `MMP9` and `ADAMTS5` highlight significant ECM remodeling activity.