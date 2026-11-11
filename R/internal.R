# Internal helpers (not exported)

.sr_as_sparse <- function(A) {
  if (inherits(A, "springrank")) return(A$A)
  if (is.data.frame(A))
    stop("`A` is a data frame; use springrank_edgelist() or springrank_moves() for edge lists.")
  if (!inherits(A, "Matrix")) {
    A <- as.matrix(A)
    if (!is.numeric(A) && !is.logical(A)) stop("`A` must be numeric.")
    dn <- dimnames(A)
    A <- Matrix::Matrix(A * 1, sparse = TRUE)
    if (!is.null(dn)) dimnames(A) <- dn
  }
  A <- methods::as(A, "CsparseMatrix")
  A <- methods::as(A, "generalMatrix")
  A <- methods::as(A, "dMatrix")
  if (nrow(A) != ncol(A)) stop("Adjacency matrix must be square.")
  if (any(!is.finite(A@x))) stop("Adjacency matrix contains non-finite entries.")
  if (any(A@x < 0)) stop("Adjacency matrix cannot contain negative entries.")
  Matrix::drop0(A)
}

.sr_names <- function(A) {
  rn <- rownames(A)
  if (is.null(rn)) rn <- colnames(A)
  if (is.null(rn)) rn <- as.character(seq_len(nrow(A)))
  rn
}

.sr_operator <- function(A, kout, kin, alpha) {
  Matrix::Diagonal(x = kout + kin + alpha) - (A + Matrix::t(A))
}

# Triplets (i, j, x) of a dgCMatrix, 1-based
.sr_triplets <- function(A) {
  list(i = A@i + 1L,
       j = rep.int(seq_len(ncol(A)), diff(A@p)),
       x = A@x)
}

# Unordered node pairs with at least one edge (self-loops dropped).
# Returns i < j, a_ij = weight i -> j, a_ji = weight j -> i.
.sr_pairs <- function(A) {
  N <- nrow(A)
  t <- .sr_triplets(A)
  keep <- t$i != t$j
  i <- t$i[keep]; j <- t$j[keep]; x <- t$x[keep]
  lo <- pmin(i, j); hi <- pmax(i, j)
  fwd <- ifelse(i == lo, x, 0)
  bwd <- x - fwd
  key <- (as.numeric(lo) - 1) * N + hi
  o <- order(key)
  key <- key[o]; lo <- lo[o]; hi <- hi[o]; fwd <- fwd[o]; bwd <- bwd[o]
  n <- length(key)
  dup <- if (n > 1L) c(FALSE, key[-1L] == key[-n]) else rep(FALSE, n)
  grp <- cumsum(!dup)
  first <- !dup
  a_ij <- fwd[first]; a_ji <- bwd[first]
  if (any(dup)) {
    g <- grp[dup]
    a_ij[g] <- a_ij[g] + fwd[dup]
    a_ji[g] <- a_ji[g] + bwd[dup]
  }
  list(i = lo[first], j = hi[first], a_ij = a_ij, a_ji = a_ji)
}

# Jacobi-preconditioned conjugate gradient for symmetric PSD systems
.sr_pcg <- function(L, b, x0, minv, tol, maxit) {
  bn <- sqrt(sum(b^2))
  if (bn == 0) return(list(x = numeric(length(b)), iter = 0L, converged = TRUE, relres = 0))
  x <- x0
  r <- b - as.vector(L %*% x)
  rn <- sqrt(sum(r^2))
  if (rn <= tol * bn) return(list(x = x, iter = 0L, converged = TRUE, relres = rn / bn))
  z <- minv * r
  p <- z
  rz <- sum(r * z)
  k <- 0L
  conv <- FALSE
  while (k < maxit) {
    k <- k + 1L
    Lp <- as.vector(L %*% p)
    a <- rz / sum(p * Lp)
    x <- x + a * p
    r <- r - a * Lp
    rn <- sqrt(sum(r^2))
    if (rn <= tol * bn) { conv <- TRUE; break }
    z <- minv * r
    rz_new <- sum(r * z)
    p <- z + (rz_new / rz) * p
    rz <- rz_new
  }
  list(x = x, iter = k, converged = conv, relres = rn / bn)
}

.sr_direct <- function(L, b, alpha) {
  N <- length(b)
  if (alpha > 0) {
    return(as.vector(Matrix::solve(Matrix::forceSymmetric(L), b)))
  }
  # alpha = 0: Lagrange multiplier for the constraint sum(s) = 0
  one <- Matrix::sparseMatrix(i = seq_len(N), j = rep(1L, N), x = 1, dims = c(N, 1L))
  K <- rbind(cbind(L, one),
             cbind(Matrix::t(one), Matrix::sparseMatrix(i = integer(0), j = integer(0),
                                                        x = numeric(0), dims = c(1L, 1L))))
  K <- methods::as(methods::as(K, "generalMatrix"), "CsparseMatrix")
  sol <- tryCatch(Matrix::solve(K, c(b, 0)),
                  error = function(e) stop("Direct solve failed (is the network weakly connected?). ",
                                           "Use alpha > 0 or method = 'cg'. Original error: ",
                                           conditionMessage(e), call. = FALSE))
  as.vector(sol)[seq_len(N)]
}

# Warn if alpha = 0 and the network has several weak components.
# Uses igraph if installed; otherwise skipped.
.sr_check_connected <- function(A) {
  if (!requireNamespace("igraph", quietly = TRUE)) return(invisible(NULL))
  t <- .sr_triplets(A)
  g <- igraph::graph_from_edgelist(cbind(t$i, t$j), directed = FALSE)
  g <- igraph::add_vertices(g, max(0L, nrow(A) - igraph::vcount(g)))
  nc <- igraph::components(g)$no
  if (nc > 1L)
    warning(sprintf(paste0("Network has %d weakly connected components (incl. isolated nodes). ",
                           "With alpha = 0 ranks are not comparable across components; ",
                           "consider alpha > 0."), nc), call. = FALSE)
  invisible(nc)
}

.sr_id <- function(x) {
  # Node identifiers as character. as.character() on doubles gives "1e+05" for
  # 100000 and keeps only 15 significant digits, which breaks matching against
  # integer ids and can merge distinct long ids; whole numbers are therefore
  # printed in full.
  if (is.double(x) && !inherits(x, c("Date", "POSIXt", "difftime"))) {
    ok <- is.finite(x)
    if (all(x[ok] == floor(x[ok])) && all(abs(x[ok]) < 2^53)) {
      out <- rep(NA_character_, length(x))
      out[ok] <- sprintf("%.0f", x[ok])
      return(out)
    }
  }
  as.character(x)
}

.sr_edgelist_matrix <- function(from, to, weight = 1, nodes = NULL) {
  if (length(from) != length(to)) stop("`from` and `to` must have the same length.")
  weight <- rep_len(as.numeric(weight), length(from))
  ok <- !(is.na(from) | is.na(to) | is.na(weight))
  if (any(!ok)) warning(sprintf("Dropping %d rows with missing values.", sum(!ok)), call. = FALSE)
  from <- .sr_id(from[ok]); to <- .sr_id(to[ok]); weight <- weight[ok]
  if (any(!is.finite(weight))) stop("Edge weights must be finite.")
  if (any(weight < 0)) stop("Edge weights cannot be negative.")
  lev <- if (is.null(nodes)) sort(unique(c(from, to))) else unique(.sr_id(nodes))
  fi <- match(from, lev); ti <- match(to, lev)
  if (anyNA(fi) || anyNA(ti)) stop("Some edge endpoints are not in `nodes`.")
  Matrix::sparseMatrix(i = fi, j = ti, x = weight, dims = c(length(lev), length(lev)),
                       dimnames = list(lev, lev))
}

.sr_index <- function(fit, v) {
  if (is.character(v)) {
    idx <- match(v, names(fit$ranks))
    if (anyNA(idx)) stop("Unknown node name(s): ", paste(v[is.na(idx)], collapse = ", "))
    idx
  } else {
    v <- as.integer(v)
    if (any(v < 1L | v > length(fit$ranks))) stop("Node index out of bounds.")
    v
  }
}

.sr_check_fit <- function(fit) {
  if (!inherits(fit, "springrank")) stop("`fit` must be a 'springrank' object from springrank().")
}
