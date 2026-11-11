# SpringRank-R-sparse

R package `springrank`: a sparse, scalable R implementation of **SpringRank** (De Bacco, Larremore & Moore, 2018) for ranking the nodes of large directed, weighted networks — for example worker moves between employers in full-population register data.

It is a port of the Python package [LarremoreLab/SpringRank](https://github.com/LarremoreLab/SpringRank). An R script also exists in the authors' reference repository [cdebacco/SpringRank](https://github.com/cdebacco/SpringRank), but it fills the matrix diagonals in loops and builds dense N×N matrices when `alpha = 0`, which does not scale to hundreds of thousands of nodes. This package keeps everything sparse.

This is an **unofficial** port. It is not affiliated with or endorsed by the original authors; see [Credits](#credits).

## Installation

```r
# install.packages("remotes")
remotes::install_github("GITHUB-USER/SpringRank-R-sparse")
```

Only depends on `Matrix` (ships with R). `igraph` is optional and only used to warn about disconnected networks when `alpha = 0`.

## Usage

```r
library(springrank)

# Worker moves between employers (one row per move, or aggregated with n)
moves <- data.frame(origin      = c("P", "P", "M", "M", "C"),
                    destination = c("M", "C", "C", "P", "M"),
                    n           = c(6, 1, 5, 2, 1))

# Revealed preference: the destination ranks above the origin
fit <- springrank_moves(moves$origin, moves$destination, moves$n, alpha = 0.5)
fit
as.data.frame(fit)        # node, rank, out_weight (moves in), in_weight (moves out)

sr_beta(fit)              # steepness of the hierarchy (beta_L)
sr_violations(fit)        # share of moves going down the hierarchy
sr_predict(fit, "P", "C") # P(a move between P and C goes from C to P, i.e. P "wins")
sr_rescale(fit, 0.75)     # 1 unit = 75% chance of winning a move
sr_test(fit, n_null = 1000)          # is the network significantly hierarchical?
sr_bootstrap(fit, R = 200)           # standard errors per employer
```

For a training/placement logic (origin employers rank above the employers they send people to, as in academic hiring), use `higher = "origin"`. Reversing the convention simply flips the sign of the ranks, so which end is "top" must be justified substantively (e.g. validated against external measures).

General networks: `springrank(A, alpha)` takes an adjacency matrix where `A[i, j]` counts observations of `i` above `j`; `springrank_edgelist(from, to, weight)` takes an edge list.

## Functions

| Function | Purpose |
|---|---|
| `springrank()` | Ranks from an adjacency matrix (sparse CG solver, or `method = "direct"`) |
| `springrank_edgelist()` | Ranks from an edge list (`from` above `to`) |
| `springrank_moves()` | Ranks from origin–destination moves, with explicit direction convention |
| `sr_beta()` | Inverse temperature: `"global"` (β_L, max. conditional likelihood) or `"local"` (β_a, max. accuracy σ_a) |
| `sr_predict()`, `sr_rescale()` | Edge-direction probabilities; rescaled ranks |
| `sr_violations()` | Edge weight against the ranking |
| `sr_energy()`, `sr_test()` | SpringRank energy; permutation test for hierarchy |
| `sr_bootstrap()` | Poisson-bootstrap standard errors (memory constant in R) |
| `sr_generate()` | Simulate from the SpringRank generative model (small N) |

## Choosing `alpha`

- `alpha = 0`: unregularized; ranks are mean-centred and only comparable within one weakly connected component.
- `alpha > 0`: every node is tied to 0 by a weak spring. Gives a unique solution for any network and shrinks nodes with few edges towards 0 — recommended for sparse data. Report results for a few values.

## Validation

`validation/compare_python.R` compares against the Python package (`pip install springrank`). When the port was written, on simulated networks (400 and 3,000 nodes, `alpha` 0 and 0.5) the ranks agreed to < 1e-7, and β_L and edge predictions agreed to all printed digits. Re-run the script to reproduce this on your own installation.

Independently of Python, the unit tests check the solvers against a dense solve, and β_L, σ_a and the energy against the paper's formulas evaluated on dense matrices.

Intentional differences from the Python package:

- **β_a (local):** σ_a (Eq. 12 in the paper) sums over all ordered node pairs. The Python code visits only stored edges, so reciprocated pairs get twice the weight of one-directional pairs. This package follows the paper, so β_a can differ. β_a is often very large or at the cap on sparse data; use β_L as the steepness measure.
- **β_L at the boundary:** the estimating equation is monotone in β, so if its root is at 0 this package returns 0 (the Python version returns the cap of 20 whenever root-finding on [0.01, 20] fails).
- **Scalability:** the Python check for a perfect hierarchy converts the matrix to dense; here all post-estimation functions work on sparse triplets.

## Performance

On simulated mobility networks with skewed employer sizes, on an Apple M3 laptop with R 4.2 (`validation/benchmark.R`):

| Employers | Moves | Build + fit | β_L | Violations |
|---|---|---|---|---|
| 100,000 | 1 million | 1.6 s | 1.5 s | 0.1 s |
| 500,000 | 10 million | 13 s | 15 s | 0.7 s |

Most of the fit time is building the matrix from the edge list; the CG solve takes about 20 iterations. Each `sr_test()` / `sr_bootstrap()` replicate costs about one fit.

## Data

This repository contains **code only**. No data are included, and none should be committed. Examples use simulated networks.

## Citation

If you use this package, please cite the original method:

> De Bacco, C., Larremore, D. B., & Moore, C. (2018). A physical model for efficient ranking in networks. *Science Advances*, 4(7), eaar8260. https://doi.org/10.1126/sciadv.aar8260

## Credits

SpringRank was developed by **Caterina De Bacco, Daniel B. Larremore and Cristopher Moore**. All credit for the method belongs to them.

This R package follows two upstream implementations:

- [LarremoreLab/SpringRank](https://github.com/LarremoreLab/SpringRank), the sparse Python package. Code by Dan Larremore, K. Hunter Wapman, Apara Venkateswaran and Ben Aoki-Sherwood.
- [cdebacco/SpringRank](https://github.com/cdebacco/SpringRank), the authors' reference implementations (Python, MATLAB, R, SAS/IML).

The R port is by Jingying Wang. It was written with the assistance of an AI coding assistant (Anthropic's Claude), working from the upstream code and the paper, and is checked by the unit tests in `tests/` and the comparison script in `validation/`. Errors in this port are not the original authors' responsibility; please report them in this repository's issue tracker, not upstream.

## License

Copyright (C) 2026 Jingying Wang.

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. It is distributed without any warranty. See `LICENSE.md`.

The upstream repositories both carry MIT license files (Copyright (c) 2024 Daniel B Larremore; Copyright (c) 2017 cdebacco), reproduced verbatim in `inst/COPYRIGHTS`. The Python package's metadata additionally declares `GPL-3.0-or-later`; GPL-3-or-later was chosen here because it is compatible with both.
