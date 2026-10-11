# Matrix Rank Determination in Database Queries

![Python](https://img.shields.io/badge/Python-3.8%2B-blue)
![NumPy](https://img.shields.io/badge/NumPy-required-lightgrey)
![Topic](https://img.shields.io/badge/Topic-Linear%20Algebra-green)

A mini project that models a set of database queries as a matrix and uses its **rank** to find redundant queries and to understand what the queries can (and cannot) reveal about the data.

---
**Live demo:** https://loading-02.github.io/Matrix-rank-project/

## Table of Contents

1. [Overview](#overview)
2. [The Idea](#the-idea)
3. [Features](#features)
4. [Project Structure](#project-structure)
5. [Getting Started](#getting-started)
6. [Example](#example)
7. [Mathematical Background](#mathematical-background)
8. [Numerical Stability](#numerical-stability)
9. [Interactive Explorer](#interactive-explorer)
10. [Verification with MATLAB and Excel](#verification-with-matlab-and-excel)
11. [Author](#author)

---

## Overview

Suppose a user runs several SUM queries over a table. Each query can be written as a row of 0s and 1s (or weights), one column per record. Stacking the rows gives a **query matrix A**.

The **rank of A** is the number of linearly independent queries. This tells us:

- **Which queries are redundant**: they can be answered from earlier ones, so a cache or optimizer can skip them.
- **Whether the data can be recovered**: if rank is less than the number of records, the queries cannot pin down every value.
- **What stays hidden**: the null space lists the data changes that no query can detect.

## The Idea

| Concept | In the database | In the matrix |
|---|---|---|
| A query | `SUM(x1 + x2)` | A row, e.g. `[1 1 0 0]` |
| A record | `x1`, `x2`, ... | A column |
| Independent queries | New information | Pivot rows |
| Redundant query | Derivable from others | Row that reduces to zero |
| Hidden data changes | Invisible to all queries | Null space of A |

## Features

- Gaussian elimination with **partial pivoting** to reach reduced row echelon form (RREF)
- Rank by elimination, cross-checked against NumPy's SVD-based rank
- Detection of **redundant queries**
- **Null space** basis (changes the queries cannot detect)
- Singular values and a **stability experiment** showing how a tiny perturbation changes the computed rank
- A browser-based interactive tool, no installation needed

## Project Structure

```
Matrix-rank-project/
├── matrix_rank.py              # Python implementation
├── matrix-rank-explorer.html   # Interactive rank calculator
└── README.md
```

## Getting Started

### Requirements

- Python 3.8 or newer
- NumPy

### Install and run

```bash
git clone https://github.com/Loading-02/Matrix-rank-project.git
cd Matrix-rank-project
pip install numpy
python matrix_rank.py
```

On Mac/Linux, use `python3` instead of `python`.

## Example

Five queries over four records:

| Query | Meaning | Row of A |
|---|---|---|
| Q1 | x1 + x2 | `[1 1 0 0]` |
| Q2 | x2 + x3 | `[0 1 1 0]` |
| Q3 | x1 + 2·x2 + x3 | `[1 2 1 0]` |
| Q4 | x3 + x4 | `[0 0 1 1]` |
| Q5 | x1 + x2 + x3 + x4 | `[1 1 1 1]` |

Here **Q3 = Q1 + Q2** and **Q5 = Q1 + Q4**, so two queries are redundant.

### Output

```
rank (elimination) = 3
rank (numpy SVD)   = 3
redundant queries  = [3, 5]
null space basis   = [[-1.0, 1.0, -1.0, 1.0]]
```

The null space vector `(1, -1, 1, -1)` means that adding `+1, -1, +1, -1` to the four record values changes none of the five query answers.

## Mathematical Background

- **Rank** = number of nonzero rows in the row echelon form = number of pivots.
- **Elementary row operations** (swap, scale, add a multiple of one row to another) do not change the rank.
- **Rank-nullity theorem**: `rank(A) + nullity(A) = number of columns`. In the example, 3 + 1 = 4.
- **Row rank = column rank**, and `rank(A) <= min(rows, columns)`.
- A second check uses **minors**: the rank is the size of the largest square submatrix with a nonzero determinant.

## Numerical Stability

Rank is decided by a tolerance, so floating-point results need care. The script perturbs one entry of Q3 and recomputes:

| Perturbation | Smallest singular value | Computed rank |
|---|---|---|
| 0 (original) | about 2e-16 | 3 |
| 1e-16 | about 6e-18 | 3 |
| 1e-10 | about 3e-11 | 4 |

Takeaways:

- A tiny singular value means the matrix is **nearly** rank-deficient (large condition number), so the rank depends on the tolerance.
- Pivoting and SVD-based rank are the standard ways to make the decision robust.
- With exact 0/1 integer entries, exact arithmetic is possible, so instability mostly appears with weighted or measured entries.

## Interactive Explorer

Open `matrix-rank-explorer.html` in any browser (double-click the file, or run `start matrix-rank-explorer.html` on Windows).

You can edit the query matrix, add or remove queries and records, and see the rank, RREF, redundant queries and null space. Change the **zero tolerance** (for example from `1e-9` to `1e-13`) on the near-dependent example to see the stability effect.

## Verification with MATLAB and Excel

**MATLAB**

```matlab
A = [1 1 0 0; 0 1 1 0; 1 2 1 0; 0 0 1 1; 1 1 1 1];
rank(A)      % 3
rref(A)
null(A)
svd(A)
```

**Excel** has no built-in matrix rank function, so use `MDETERM` on minors:

- The 3×3 minor from rows Q1, Q2, Q4 (columns 1-3) has determinant 1, so rank is at least 3.
- Every 4×4 minor has determinant 0, so rank is at most 3.

## Author

**Loading-02**
GitHub: [@Loading-02](https://github.com/Loading-02)
