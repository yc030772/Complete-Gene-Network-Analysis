"""Gene interaction network + Relevance Value (RV) ranking for CADM1+ vs CADM1- groups.

Model: x_i[n] = sum_{j != i} a_ij x_j[n] + eps_i[n], solved row-by-row with Ridge.
RV_j  = sum_i |a_ij(+) - a_ij(-)|
"""
import sys

import numpy as np
import pandas as pd
from sklearn.linear_model import Ridge


def solve_interaction_matrix(X: np.ndarray, alpha: float = 1.0) -> np.ndarray:
    """Return A (p x p) where A[i, j] is gene j's effect on gene i; diagonal is 0.

    X has shape (samples, genes). No intercept: the model has no constant term.
    """
    p = X.shape[1]
    A = np.zeros((p, p))
    model = Ridge(alpha=alpha, fit_intercept=False)
    for i in range(p):
        model.fit(np.delete(X, i, axis=1), X[:, i])
        A[i] = np.insert(model.coef_, i, 0)
    return A


def relevance_values(A_pos: np.ndarray, A_neg: np.ndarray, genes) -> pd.DataFrame:
    """Column sums of |A+ - A-|, sorted descending."""
    rv = np.abs(A_pos - A_neg).sum(axis=0)
    df = pd.DataFrame({"Gene": genes, "RV": rv})
    df = df.sort_values("RV", ascending=False, ignore_index=True)
    df.index += 1
    return df


def main(path: str = "AIIM_HW2_ANOVA_grouped_CADM1.npz", top: int = 15) -> pd.DataFrame:
    data = np.load(path, allow_pickle=True)
    X_pos, X_neg, genes = data["x_pos"], data["x_neg"], data["ANOVA_genes"]
    print(f"X+ {X_pos.shape}, X- {X_neg.shape}, {len(genes)} genes")

    rv = relevance_values(solve_interaction_matrix(X_pos), solve_interaction_matrix(X_neg), genes)
    print(rv.head(top).to_string(float_format="%.4f"))
    return rv


if __name__ == "__main__":
    main(*sys.argv[1:2])
