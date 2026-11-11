# Timing on simulated mobility-like networks (skewed employer sizes).
library(springrank)
bench <- function(N, M, alpha = 1, seed = 1) {
  set.seed(seed)
  true <- rnorm(N); w <- (rexp(N))^2
  o <- sample.int(N, M, TRUE, w); d <- sample.int(N, M, TRUE, w)
  up <- runif(M) < plogis(2 * (true[d] - true[o]))
  t0 <- proc.time()[["elapsed"]]
  fit <- springrank_moves(ifelse(up, o, d), ifelse(up, d, o), alpha = alpha,
                          nodes = seq_len(N))
  t1 <- proc.time()[["elapsed"]]
  b <- sr_beta(fit); t2 <- proc.time()[["elapsed"]]
  v <- sr_violations(fit); t3 <- proc.time()[["elapsed"]]
  act <- (fit$out_weight + fit$in_weight) >= 10
  cat(sprintf("N=%9s moves=%11s pairs=%11s | build+fit %5.1fs (%d CG it.) | beta %4.1fs | violations %4.1fs | cor(true) firms>=10 moves: %.2f\n",
              format(N, big.mark = ",", scientific = FALSE), format(M, big.mark = ",", scientific = FALSE),
              format(length(fit$A@x), big.mark = ","), t1 - t0, fit$iterations,
              t2 - t1, t3 - t2, cor(fit$ranks[act], true[act])))
}
bench(1e5, 1e6)
bench(5e5, 1e7)
