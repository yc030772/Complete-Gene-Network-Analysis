import numpy as np

from gene_network import relevance_values, solve_interaction_matrix


def test_recovers_known_interaction():
    rng = np.random.default_rng(0)
    X = rng.normal(size=(500, 3))
    X[:, 0] = 2 * X[:, 1] + 0.01 * rng.normal(size=500)
    A = solve_interaction_matrix(X, alpha=1e-6)
    assert np.all(np.diag(A) == 0)
    assert abs(A[0, 1] - 2) < 0.05


def test_rv_ranks_changed_gene_first():
    A_pos = np.zeros((3, 3))
    A_neg = np.zeros((3, 3))
    A_pos[0, 2] = 5  # gene C's influence changes most
    A_pos[1, 0] = 1
    rv = relevance_values(A_pos, A_neg, ["A", "B", "C"])
    assert list(rv["Gene"]) == ["C", "A", "B"]
    assert rv.loc[1, "RV"] == 5


if __name__ == "__main__":
    test_recovers_known_interaction()
    test_rv_ranks_changed_gene_first()
    print("ok")
