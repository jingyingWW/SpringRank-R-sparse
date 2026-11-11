# Cross-check the R port against the Python package (pip install springrank).
# Usage: Rscript validation/compare_python.R   (writes and reads files in tempdir())
library(springrank)
dir <- tempdir()
set.seed(11)
cases <- list(gen = sr_generate(400, beta = 1.5, alpha = 1, K = 5)$A)
# a sparse, skewed "mobility-like" network
N <- 3000; M <- 20000
true <- rnorm(N); w <- rexp(N)^2
o <- sample(N, M, TRUE, w); d <- sample(N, M, TRUE, w)
up <- runif(M) < plogis(2 * (true[d] - true[o]))
from <- ifelse(up, d, o); to <- ifelse(up, o, d)
keep <- from != to
A <- Matrix::sparseMatrix(from[keep], to[keep], x = 1, dims = c(N, N))
cases$mobility <- A

py <- '
import sys, numpy as np, scipy.sparse as sp
from springrank import SpringRank
f, alpha, out = sys.argv[1], float(sys.argv[2]), sys.argv[3]
e = np.loadtxt(f, ndmin=2); N = int(open(f + ".n").read())
A = sp.csr_matrix((e[:,2], (e[:,0].astype(int), e[:,1].astype(int))), shape=(N, N))
m = SpringRank(alpha=alpha, rtol=1e-12, inverse_temp_fit_warning=False); m.fit(A)
bg = m.get_beta()
m2 = SpringRank(alpha=alpha, rtol=1e-12, inverse_temp_type="local", inverse_temp_fit_warning=False); m2.fit(A)
bl = m2.get_beta()
np.savetxt(out, m.ranks); open(out + ".beta", "w").write(f"{bg} {bl} {m.predict([0,1])}")
'
pyfile <- file.path(dir, "sr.py"); writeLines(py, pyfile)

for (nm in names(cases)) for (alpha in c(0, 0.5)) {
  A <- cases[[nm]]
  if (alpha == 0) {  # alpha = 0 needs a connected graph: restrict to the largest component
    g <- igraph::graph_from_adjacency_matrix(A, mode = "max", weighted = TRUE)
    cmp <- igraph::components(g); keep <- which(cmp$membership == which.max(cmp$csize))
    A <- A[keep, keep]
  }
  t <- Matrix::summary(as(A, "generalMatrix"))
  f <- file.path(dir, "edges.txt"); out <- file.path(dir, "ranks.txt")
  write.table(cbind(t$i - 1, t$j - 1, t$x), f, row.names = FALSE, col.names = FALSE)
  writeLines(as.character(nrow(A)), paste0(f, ".n"))
  system2("python3", c(pyfile, f, alpha, out))
  rp <- scan(out, quiet = TRUE); bp <- scan(paste0(out, ".beta"), quiet = TRUE)
  fit <- springrank(A, alpha = alpha)
  bg <- suppressWarnings(as.numeric(sr_beta(fit)))
  bl <- suppressWarnings(as.numeric(sr_beta(fit, "local")))
  cat(sprintf("%-9s alpha=%-4s N=%5d  max|rank diff|=%.1e  beta_L R=%.5f Py=%.5f  beta_a R=%.3f Py=%.3f  P(0->1) R=%.5f Py=%.5f\n",
              nm, alpha, nrow(A), max(abs(fit$ranks - rp)), bg, bp[1], bl, bp[2],
              sr_predict(fit, 1, 2, beta = bg), bp[3]))
}
