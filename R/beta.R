#' Inverse temperature (steepness) of a SpringRank hierarchy
#'
#' Estimates \eqn{\beta}, which controls how strictly edges follow the ranks:
#' the probability that an edge between \eqn{i} and \eqn{j} points
#' \eqn{i \to j} is \eqn{P_{ij} = 1 / (1 + e^{-2\beta(s_i - s_j)})}.
#' Larger \eqn{\beta} means a steeper hierarchy.
#'
#' * `type = "global"` returns \eqn{\hat\beta_L}, the maximizer of the
#'   conditional log-likelihood (root of Eq. S39 in De Bacco et al. 2018).
#'   The estimating equation is monotone in \eqn{\beta}, so the root is found
#'   on `[0, max_beta]`; if it lies at 0 the function returns 0, and if it
#'   exceeds `max_beta` it is capped with a warning.
#' * `type = "local"` returns \eqn{\hat\beta_a}, the maximizer of the edge
#'   direction accuracy \eqn{\sigma_a} (Eq. 12), found by bounded 1-d
#'   optimization over `[0, 30]` and capped at `max_beta`.
#'
#' If no edge points against the ranking, the optimum is infinite and
#' `max_beta` is returned with a warning. If all connected nodes are tied,
#' \eqn{\beta} is not identified and 0 is returned with a warning.
#'
#' @param fit A `"springrank"` object.
#' @param type `"global"` (default) or `"local"`.
#' @param max_beta Upper cap for \eqn{\beta} (default 20, as in the Python
#'   package).
#' @return A single number, with attributes `sigma_a` (edge-direction
#'   accuracy) and `loglik_per_edge` (mean conditional log-likelihood per
#'   edge) evaluated at the returned \eqn{\beta}.
#' @export
sr_beta <- function(fit, type = c("global", "local"), max_beta = 20) {
  .sr_check_fit(fit)
  type <- match.arg(type)
  p <- .sr_pairs(fit$A)
  s <- unname(fit$ranks)
  d <- s[p$i] - s[p$j]
  abar <- p$a_ij + p$a_ji

  if (length(d) == 0L) stop("Network has no edges between distinct nodes.")

  viol <- sr_violations(fit)
  if (viol$against == 0 && viol$along == 0) {
    warning("All connected nodes have the same rank: beta is not identified; returning 0.",
            call. = FALSE)
    beta <- 0
  } else if (viol$against == 0) {
    warning(sprintf("No edge points against the ranking: optimal beta is infinite; returning max_beta = %g.",
                    max_beta), call. = FALSE)
    beta <- max_beta
  } else if (type == "global") {
    f <- function(beta) sum(d * (p$a_ij - abar * stats::plogis(2 * beta * d)))
    if (f(0) <= 0) {
      beta <- 0
    } else if (f(max_beta) > 0) {
      warning(sprintf(paste0("Estimating equation for beta_L has no root below max_beta = %g ",
                             "(steep or over-fitted hierarchy). Returning max_beta; ",
                             "consider larger alpha."), max_beta), call. = FALSE)
      beta <- max_beta
    } else {
      beta <- stats::uniroot(f, c(0, max_beta), tol = 1e-12)$root
    }
  } else {
    g <- function(beta) sum(abs(p$a_ij - abar * stats::plogis(2 * beta * d)))
    grid <- seq(0, 30, length.out = 6L)
    best <- NULL
    for (k in seq_len(length(grid) - 1L)) {
      r <- stats::optimize(g, c(grid[k], grid[k + 1L]))
      if (is.null(best) || r$objective < best$objective) best <- r
    }
    beta <- best$minimum
    if (beta > max_beta) {
      warning(sprintf("beta_a exceeds max_beta (possible over-fitting); capping at %g.", max_beta),
              call. = FALSE)
      beta <- max_beta
    }
  }
  acc <- .sr_accuracy(p, d, abar, beta)
  structure(beta, sigma_a = acc$sigma_a, loglik_per_edge = acc$loglik_per_edge, type = type)
}

.sr_accuracy <- function(p, d, abar, beta) {
  P <- stats::plogis(2 * beta * d)
  M <- sum(abar)
  # sigma_a (Eq. 12): 1 - (1/2M) sum over ordered pairs |A_ij - Abar_ij P_ij|;
  # both orderings of a pair contribute the same amount, hence 2 * sum / 2M.
  sigma_a <- 1 - sum(abs(p$a_ij - abar * P)) / M
  ll <- sum(p$a_ij * stats::plogis(2 * beta * d, log.p = TRUE) +
              p$a_ji * stats::plogis(-2 * beta * d, log.p = TRUE))
  list(sigma_a = sigma_a, loglik_per_edge = ll / M)
}

#' Predict edge direction between two nodes
#'
#' Probability that an edge between `i` and `j` points from `i` to `j`
#' (i.e. that `i` is observed above `j`), \eqn{1/(1+e^{-2\beta(s_i - s_j)})}.
#'
#' @param fit A `"springrank"` object.
#' @param i,j Node names or indices (vectorized).
#' @param beta Inverse temperature; estimated with [sr_beta()] if `NULL`.
#' @param ... Passed to [sr_beta()] when `beta` is `NULL`.
#' @return Numeric vector of probabilities.
#' @export
sr_predict <- function(fit, i, j, beta = NULL, ...) {
  .sr_check_fit(fit)
  if (is.null(beta)) beta <- as.numeric(sr_beta(fit, ...))
  ii <- .sr_index(fit, i); jj <- .sr_index(fit, j)
  stats::plogis(2 * beta * (fit$ranks[ii] - fit$ranks[jj]))
}

#' Rescale SpringRank scores to an interpretable unit
#'
#' Rescales ranks so that a difference of one unit corresponds to a
#' probability `target` that the edge points from the higher- to the
#' lower-ranked node (e.g. `target = 0.75`: a firm one unit higher wins a
#' move between the two firms 75% of the time).
#'
#' @param fit A `"springrank"` object.
#' @param target Probability in (0.5, 1).
#' @param beta Inverse temperature; estimated with [sr_beta()] if `NULL`.
#' @param ... Passed to [sr_beta()] when `beta` is `NULL`.
#' @return Named numeric vector of rescaled ranks.
#' @export
sr_rescale <- function(fit, target = 0.75, beta = NULL, ...) {
  .sr_check_fit(fit)
  if (!(target > 0.5 && target < 1)) stop("`target` must be in (0.5, 1).")
  if (is.null(beta)) beta <- as.numeric(sr_beta(fit, ...))
  fit$ranks * (2 * beta / stats::qlogis(target))
}
