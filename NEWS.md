# springrank 0.1.0

* First public release.
* Sparse SpringRank fit (`springrank()`, `springrank_edgelist()`,
  `springrank_moves()`) with a Jacobi-preconditioned conjugate-gradient solver
  and an optional sparse direct solver.
* Post-estimation: `sr_beta()`, `sr_predict()`, `sr_rescale()`,
  `sr_violations()`, `sr_energy()`, `sr_test()`, `sr_bootstrap()`,
  `sr_generate()`.
* Ported from the Python package
  [LarremoreLab/SpringRank](https://github.com/LarremoreLab/SpringRank)
  (v0.0.10); see the README for intentional differences.
