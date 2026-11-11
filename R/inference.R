#' Edges that go against the ranking
#'
#' Counts the edge weight pointing from a lower- to a higher-ranked node,
#' i.e. "violations" of the hierarchy. For [springrank_moves()] with
#' `higher = "destination"` these are moves from a higher- to a lower-ranked
#' employer (downward moves). Self-loops are ignored.
#'
#' @param fit A `"springrank"` object.
#' @return A list with `rate` (share of edge weight against the ranking, ties
#'   counted as half), `against`, `tied`, `along` and `total` edge weight.
#' @export
sr_violations <- function(fit) {
  .sr_check_fit(fit)
  t <- .sr_triplets(fit$A)
  keep <- t$i != t$j
  s <- unname(fit$ranks)
  d <- s[t$i[keep]] - s[t$j[keep]]
  x <- t$x[keep]
  against <- sum(x[d < 0]); tied <- sum(x[d == 0]); along <- sum(x[d > 0])
  total <- against + tied + along
  list(rate = if (total > 0) (against + tied / 2) / total else NA_real_,
       against = against, tied = tied, along = along, total = total)
}

#' SpringRank energy
#'
#' Energy of the interaction springs, \eqn{H(s) = \frac12 \sum_{i \ne j}
#' A_{ij}(s_i - s_j - 1)^2}. Lower energy per edge means the network is closer
#' to a perfect hierarchy.
#'
#' @param A Adjacency matrix (or a `"springrank"` object, whose `A` and
#'   `ranks` are then used).
#' @param s Rank vector (ignored if `A` is a `"springrank"` object).
#' @return A single number.
#' @export
sr_energy <- function(A, s = NULL) {
  if (inherits(A, "springrank")) { s <- A$ranks; A <- A$A }
  A <- .sr_as_sparse(A)
  s <- unname(s)
  if (length(s) != nrow(A)) stop("`s` must have length nrow(A).")
  t <- .sr_triplets(A)
  keep <- t$i != t$j
  0.5 * sum(t$x[keep] * (s[t$i[keep]] - s[t$j[keep]] - 1)^2)
}

#' Significance test for hierarchy
#'
#' Permutation test of De Bacco et al. (2018): the observed SpringRank energy
#' is compared with networks in which the direction of every edge is
#' randomized while the total weight between each pair of nodes,
#' \eqn{A_{ij} + A_{ji}}, is kept fixed. Ranks are re-estimated (same
#' `alpha`) for each null network. A small p-value means the network is more
#' hierarchical than expected given who interacts with whom.
#'
#' @param fit A `"springrank"` object. Edge weights must be integer counts.
#' @param n_null Number of randomized networks.
#' @param seed Optional random seed.
#' @param progress Print progress every 10% of draws.
#' @param ... Passed to [springrank()] for the null fits (e.g. `tol`).
#' @return A list with `p_value` (one-sided, with +1 correction),
#'   `energy` (observed), `energy_null` (vector), and the same quantities per
#'   unit of edge weight (`energy_per_edge`, `energy_null_per_edge`).
#' @export
sr_test <- function(fit, n_null = 1000L, seed = NULL, progress = FALSE, ...) {
  .sr_check_fit(fit)
  if (!is.null(seed)) set.seed(seed)
  A <- fit$A
  N <- nrow(A)
  dn <- dimnames(A)
  p <- .sr_pairs(A)
  abar <- p$a_ij + p$a_ji
  if (length(abar) == 0L) stop("Network has no edges between distinct nodes.")
  if (any(abs(A@x - round(A@x)) > 1e-8))
    stop("sr_test() randomizes the direction of individual edges and therefore ",
         "needs integer edge weights (counts).", call. = FALSE)
  abar <- round(abar)
  dg <- Matrix::diag(A)
  loops <- which(dg != 0)
  M <- sum(abar)
  H_obs <- sr_energy(A, fit$ranks)
  H_null <- numeric(n_null)
  step <- max(1L, n_null %/% 10L)
  for (r in seq_len(n_null)) {
    fw <- stats::rbinom(length(abar), size = abar, prob = 0.5)
    bw <- abar - fw
    Ar <- Matrix::sparseMatrix(i = c(p$i, p$j, loops), j = c(p$j, p$i, loops),
                               x = c(fw, bw, dg[loops]), dims = c(N, N))
    Ar <- Matrix::drop0(Ar)
    sr <- suppressWarnings(springrank(Ar, alpha = fit$alpha, method = fit$method, ...))
    H_null[r] <- sr_energy(Ar, sr$ranks)
    if (progress && r %% step == 0L) message(sprintf("sr_test: %d / %d", r, n_null))
  }
  list(p_value = (1 + sum(H_null <= H_obs)) / (n_null + 1),
       energy = H_obs, energy_null = H_null,
       energy_per_edge = H_obs / M, energy_null_per_edge = H_null / M)
}

#' Bootstrap standard errors for SpringRank scores
#'
#' Poisson bootstrap: each edge weight \eqn{A_{ij}} is replaced by a draw
#' from \eqn{Poisson(A_{ij})}, which approximates resampling individual
#' edges (moves) with replacement and scales to large networks. Ranks are
#' re-estimated for each replicate (warm-started from the original fit), and
#' running means and variances are accumulated so memory use does not grow
#' with `R` unless `keep = TRUE`.
#'
#' @param fit A `"springrank"` object.
#' @param R Number of bootstrap replicates.
#' @param seed Optional random seed.
#' @param keep If `TRUE`, also return the full `nodes x R` matrix of
#'   replicate ranks (memory: 8 * nodes * R bytes).
#' @param progress Print progress every 10% of replicates.
#' @param ... Passed to [springrank()] (e.g. `tol`).
#' @return A data frame with `node`, `rank`, `boot_mean`, `se`, sorted by
#'   rank; with `keep = TRUE` a list with that data frame (`summary`) and the
#'   matrix (`replicates`).
#' @note With `alpha = 0` ranks are only identified up to a shift, so each
#'   replicate is mean-centred; standard errors then refer to positions
#'   relative to the mean.
#' @export
sr_bootstrap <- function(fit, R = 200L, seed = NULL, keep = FALSE, progress = FALSE, ...) {
  .sr_check_fit(fit)
  if (!is.null(seed)) set.seed(seed)
  A <- fit$A
  N <- nrow(A)
  t <- .sr_triplets(A)
  mean_ <- numeric(N); m2 <- numeric(N)
  reps <- if (keep) matrix(NA_real_, N, R, dimnames = list(names(fit$ranks), NULL)) else NULL
  step <- max(1L, R %/% 10L)
  for (r in seq_len(R)) {
    x <- stats::rpois(length(t$x), t$x)
    Ab <- Matrix::drop0(Matrix::sparseMatrix(i = t$i, j = t$j, x = x, dims = c(N, N)))
    sb <- suppressWarnings(springrank(Ab, alpha = fit$alpha, method = fit$method,
                                      x0 = if (fit$method == "cg") unname(fit$ranks) else NULL,
                                      ...))$ranks
    delta <- sb - mean_
    mean_ <- mean_ + delta / r
    m2 <- m2 + delta * (sb - mean_)
    if (keep) reps[, r] <- sb
    if (progress && r %% step == 0L) message(sprintf("sr_bootstrap: %d / %d", r, R))
  }
  out <- data.frame(node = names(fit$ranks), rank = unname(fit$ranks),
                    boot_mean = mean_, se = sqrt(m2 / max(1L, R - 1L)),
                    stringsAsFactors = FALSE)
  out <- out[order(out$rank, decreasing = TRUE), , drop = FALSE]
  rownames(out) <- NULL
  if (keep) list(summary = out, replicates = reps) else out
}

#' Simulate a network from the SpringRank generative model
#'
#' Draws ranks \eqn{s_i \sim N(0, 1/(\alpha\beta))} and edge counts
#' \eqn{A_{ij} \sim Poisson(c \exp[-\frac{\beta}{2}(s_i - s_j - 1)^2])},
#' with \eqn{c} set so the mean out-degree is `K`. Uses dense \eqn{N \times N}
#' computations, so it is intended for simulation studies with up to a few
#' thousand nodes. Self-loops are not generated.
#'
#' @param N Number of nodes.
#' @param beta Inverse temperature (steepness).
#' @param alpha Prior precision of the ranks.
#' @param K Mean (weighted) out-degree.
#' @param seed Optional random seed.
#' @return A list with `A` (sparse adjacency matrix) and `ranks` (true ranks).
#' @export
sr_generate <- function(N, beta = 1, alpha = 1, K = 5, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  s <- stats::rnorm(N, 0, 1 / sqrt(alpha * beta))
  D <- outer(s, s, "-")
  E <- exp(-0.5 * beta * (D - 1)^2)
  diag(E) <- 0
  cc <- K * N / sum(E)
  A <- matrix(stats::rpois(N * N, cc * E), N, N)
  nm <- as.character(seq_len(N))
  A <- Matrix::Matrix(A, sparse = TRUE)
  A <- methods::as(methods::as(methods::as(A, "CsparseMatrix"), "generalMatrix"), "dMatrix")
  dimnames(A) <- list(nm, nm)
  list(A = A, ranks = stats::setNames(s, nm))
}
