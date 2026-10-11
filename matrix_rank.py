"""Matrix Rank Determination in Database Queries (mini project)

Each row of A is a query, each column is a record.
rank(A) = number of independent queries.
"""
import numpy as np


def rref(M, tol=1e-9):
    """Gaussian elimination with partial pivoting. Returns (R, pivot_columns)."""
    R = np.array(M, dtype=float)
    m, n = R.shape
    pivots, row = [], 0
    for c in range(n):
        if row >= m:
            break
        best = row + np.argmax(np.abs(R[row:, c]))
        if abs(R[best, c]) <= tol:
            continue
        R[[row, best]] = R[[best, row]]
        R[row] /= R[row, c]
        for i in range(m):
            if i != row:
                R[i] -= R[i, c] * R[row]
        pivots.append(c)
        row += 1
    return R, pivots


def null_space(M, tol=1e-9):
    R, piv = rref(M, tol)
    n = R.shape[1]
    basis = []
    for f in (c for c in range(n) if c not in piv):
        v = np.zeros(n)
        v[f] = 1.0
        for k, pc in enumerate(piv):
            v[pc] = -R[k, f]
        basis.append(v)
    return basis


def redundant_queries(M, tol=1e-9):
    """Return indices (1-based) of queries that depend on earlier ones."""
    M = np.array(M, dtype=float)
    out, prev = [], 0
    for q in range(len(M)):
        r = len(rref(M[: q + 1], tol)[1])
        if r == prev:
            out.append(q + 1)
        prev = r
    return out


def report(M, tol=1e-9):
    M = np.array(M, dtype=float)
    R, piv = rref(M, tol)
    print("A =\n", M)
    print("\nRREF(A) =\n", np.round(R, 6))
    print("\nrank (elimination) =", len(piv))
    print("rank (numpy SVD)   =", np.linalg.matrix_rank(M))
    print("redundant queries  =", redundant_queries(M, tol))
    print("null space basis   =", [np.round(v, 6).tolist() for v in null_space(M, tol)])
    print("singular values    =", np.linalg.svd(M, compute_uv=False))


def stability_demo(A):
    print("\n--- Stability experiment: perturb A[2,3] ---")
    for eps in (0.0, 1e-16, 1e-10):
        B = np.array(A, dtype=float)
        B[2, 3] = eps
        s = np.linalg.svd(B, compute_uv=False)
        print(f"eps={eps:.0e}  rank={np.linalg.matrix_rank(B)}  smallest singular value={s[-1]:.2e}")


if __name__ == "__main__":
    A = [
        [1, 1, 0, 0],   # Q1
        [0, 1, 1, 0],   # Q2
        [1, 2, 1, 0],   # Q3 = Q1 + Q2
        [0, 0, 1, 1],   # Q4
        [1, 1, 1, 1],   # Q5 = Q1 + Q4
    ]
    report(A)
    stability_demo(A)
