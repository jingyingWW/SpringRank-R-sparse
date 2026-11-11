#' Fit SpringRank to a directed network
#'
#' Computes SpringRank scores (De Bacco, Larremore & Moore 2018) by solving the
#' sparse linear system
#' \deqn{[D^{out} + D^{in} - (A + A^T) + \alpha I] s = (D^{out} - D^{in}) 1.}
#' An edge \eqn{i \to j} with weight \eqn{A_{ij}} is evidence that node
#' \eqn{i} ranks *above* node \eqn{j}. Only observed edges enter the system, so
#' the many zero cells of a sparse network cost nothing.
#'
#' @param A A square adjacency matrix: a base `matrix` or any sparse matrix from
#'   the \pkg{Matrix} package. `A[i, j]` is the (non-negative) weight of edges
#'   from `i` to `j`, i.e. the number of times `i` was observed above `j`.
#' @param alpha Regularization strength (`>= 0`). With `alpha = 0` the ranks
#'   are only defined up to a constant and are mean-centred (as in the Python
#'   package); this requires a weakly connected network. With `alpha > 0` each
#'   node is tied to zero by a spring of stiffness `alpha`, giving a unique
#'   solution and shrinking nodes with few edges towards zero.
#' @param method `"cg"` (default) uses a Jacobi-preconditioned conjugate
#'   gradient solver, which scales to millions of nodes. `"direct"` uses a
#'   sparse factorization, which is exact but may run out of memory on very
#'   large networks.
#' @param tol Relative residual tolerance for `"cg"`.
#' @param maxit Maximum number of `"cg"` iterations.
#' @param x0 Optional starting vector for `"cg"` (e.g. ranks from an earlier
#'   fit).
#'
#' @return An object of class `"springrank"`: a list with `ranks` (named
#'   numeric vector), `alpha`, the sparse adjacency matrix `A`, `out_weight`
#'   and `in_weight` (weighted out-/in-degree), and solver diagnostics.
#'
#' @references De Bacco, C., Larremore, D. B., & Moore, C. (2018). A physical
#'   model for efficient ranking in networks. *Science Advances*, 4(7),
#'   eaar8260.
#' @examples
#' A <- matrix(c(0, 3, 1,
#'               0, 0, 2,
#'               0, 1, 0), 3, byrow = TRUE,
#'             dimnames = list(c("a", "b", "c"), c("a", "b", "c")))
#' fit <- springrank(A, alpha = 0)
#' fit$ranks
#' @export
springrank <- function(A, alpha = 0, method = c("cg", "direct"),
                       tol = 1e-10, maxit = 10000L, x0 = NULL) {
  method <- match.arg(method)
  A <- .sr_as_sparse(A)
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) || alpha < 0)
    stop("`alpha` must be a single non-negative number.")
  N <- nrow(A)

  kout <- Matrix::rowSums(A)
  kin <- Matrix::colSums(A)
  b <- kout - kin
  L <- .sr_operator(A, kout, kin, alpha)

  if (alpha == 0) .sr_check_connected(A)

  if (method == "cg") {
    d <- kout + kin + alpha
    minv <- ifelse(d > 0, 1 / d, 1)
    if (is.null(x0)) x0 <- numeric(N)
    if (length(x0) != N) stop("`x0` must have length nrow(A).")
    sol <- .sr_pcg(L, b, x0, minv, tol, maxit)
    if (!sol$converged)
      warning(sprintf("CG did not converge in %d iterations (relative residual %.2e).",
                      sol$iter, sol$relres))
    s <- sol$x
    diag_info <- sol[c("iter", "converged", "relres")]
  } else {
    s <- .sr_direct(L, b, alpha)
    diag_info <- list(iter = NA_integer_, converged = TRUE,
                      relres = sqrt(sum((b - as.vector(L %*% s))^2)) / max(sqrt(sum(b^2)), 1e-300))
  }

  if (alpha == 0) s <- s - mean(s)
  names(s) <- .sr_names(A)

  structure(list(ranks = s, alpha = alpha, A = A,
                 out_weight = stats::setNames(kout, names(s)),
                 in_weight = stats::setNames(kin, names(s)),
                 method = method, iterations = diag_info$iter,
                 converged = diag_info$converged, relres = diag_info$relres),
            class = "springrank")
}

#' Fit SpringRank from an edge list
#'
#' Convenience wrapper around [springrank()] for data stored as an edge list.
#' Repeated `(from, to)` rows are summed.
#'
#' @param from,to Vectors of node identifiers (any type coercible to
#'   character). An edge `from -> to` means `from` ranks above `to`.
#' @param weight Edge weights (counts); recycled to the length of `from`.
#' @param nodes Optional vector of all node identifiers, to include nodes that
#'   have no edges (they get rank 0 when `alpha > 0`).
#' @param ... Passed to [springrank()] (e.g. `alpha`, `method`).
#' @return A `"springrank"` object.
#' @examples
#' fit <- springrank_edgelist(from = c("a", "a", "b"), to = c("b", "c", "c"),
#'                            alpha = 0.1)
#' fit$ranks
#' @export
springrank_edgelist <- function(from, to, weight = 1, nodes = NULL, ...) {
  A <- .sr_edgelist_matrix(from, to, weight, nodes)
  springrank(A, ...)
}

#' Fit SpringRank to worker moves between employers
#'
#' Builds the SpringRank network from job-to-job moves. By default a move
#' from `origin` to `destination` is read as revealed preference: the
#' destination ranks above the origin (edge `destination -> origin`). Set
#' `higher = "origin"` for a placement/training logic in which origin
#' employers rank above the employers they send workers to (as in
#' academic hiring networks).
#'
#' @param origin,destination Employer identifiers of each move (or of each
#'   aggregated origin-destination cell).
#' @param n Number of moves in each row (default 1, one row per move).
#' @param higher Which side of a move is treated as higher-ranked:
#'   `"destination"` (default) or `"origin"`.
#' @param nodes Optional vector of all employer identifiers.
#' @param ... Passed to [springrank()] (e.g. `alpha`).
#' @return A `"springrank"` object. `out_weight` is then the number of moves
#'   *into* an employer (hires from other employers) when
#'   `higher = "destination"`, and `in_weight` the number of moves out.
#' @examples
#' moves <- data.frame(origin = c("P", "P", "M", "M", "C"),
#'                     destination = c("M", "C", "C", "P", "M"),
#'                     n = c(6, 1, 5, 2, 1))
#' fit <- springrank_moves(moves$origin, moves$destination, moves$n, alpha = 0.5)
#' sort(fit$ranks, decreasing = TRUE)
#' @export
springrank_moves <- function(origin, destination, n = 1,
                             higher = c("destination", "origin"),
                             nodes = NULL, ...) {
  higher <- match.arg(higher)
  if (higher == "destination") {
    fit <- springrank_edgelist(destination, origin, n, nodes, ...)
  } else {
    fit <- springrank_edgelist(origin, destination, n, nodes, ...)
  }
  fit$higher <- higher
  fit
}

#' @export
print.springrank <- function(x, n = 5L, ...) {
  A <- x$A
  M <- sum(A@x) - sum(Matrix::diag(A))
  cat("SpringRank fit\n")
  cat(sprintf("  nodes: %s   edges (weighted, excl. self-loops): %s   alpha: %g\n",
              format(nrow(A), big.mark = ","), format(M, big.mark = ","), x$alpha))
  if (!is.null(x$higher))
    cat(sprintf("  moves read as: %s ranks higher\n", x$higher))
  if (identical(x$method, "cg"))
    cat(sprintf("  solver: cg, %d iterations, converged = %s (rel. residual %.1e)\n",
                x$iterations, x$converged, x$relres))
  else cat("  solver: direct\n")
  o <- order(x$ranks, decreasing = TRUE)
  k <- min(n, length(o))
  cat(sprintf("  top %d:\n", k))
  print(round(x$ranks[o[seq_len(k)]], 3))
  invisible(x)
}

#' @export
as.data.frame.springrank <- function(x, row.names = NULL, optional = FALSE, ...) {
  out <- data.frame(node = names(x$ranks), rank = unname(x$ranks),
                    out_weight = unname(x$out_weight),
                    in_weight = unname(x$in_weight),
                    stringsAsFactors = FALSE)
  out[order(out$rank, decreasing = TRUE), , drop = FALSE]
}
