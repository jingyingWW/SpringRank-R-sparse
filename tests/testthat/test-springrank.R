toy <- function() {
  # moves origin -> destination: P->M 6, M->P 2, M->C 5, C->M 1, P->C 1
  data.frame(origin = c("P", "M", "M", "C", "P"),
             destination = c("M", "P", "C", "M", "C"),
             n = c(6, 2, 5, 1, 1))
}

dense_solution <- function(A, alpha) {
  A <- as.matrix(A)
  kout <- rowSums(A); kin <- colSums(A)
  L <- diag(kout + kin + alpha) - (A + t(A))
  b <- kout - kin
  if (alpha > 0) return(solve(L, b))
  N <- nrow(A)
  K <- rbind(cbind(L, 1), c(rep(1, N), 0))
  solve(K, c(b, 0))[seq_len(N)]
}

test_that("cg and direct solvers match a dense solve", {
  set.seed(1)
  g <- sr_generate(60, beta = 2, alpha = 1, K = 6)
  for (a in c(0, 0.3, 2)) {
    ref <- unname(dense_solution(g$A, a))
    expect_equal(unname(springrank(g$A, alpha = a)$ranks), ref, tolerance = 1e-6)
    expect_equal(unname(springrank(g$A, alpha = a, method = "direct")$ranks), ref, tolerance = 1e-8)
  }
})

test_that("alpha = 0 ranks are mean-centred", {
  g <- sr_generate(40, beta = 1, K = 8, seed = 2)
  expect_equal(mean(springrank(g$A)$ranks), 0, tolerance = 1e-10)
})

test_that("toy mobility example gives the expected ordering and values", {
  m <- toy()
  fit <- springrank_moves(m$origin, m$destination, m$n, alpha = 0)
  expect_equal(unname(round(fit$ranks[c("P", "M", "C")], 2)), c(-0.54, -0.05, 0.59))
  # reversing the convention flips the ranking
  fit2 <- springrank_moves(m$origin, m$destination, m$n, alpha = 0, higher = "origin")
  expect_equal(fit2$ranks[c("P", "M", "C")], -fit$ranks[c("P", "M", "C")], tolerance = 1e-8)
  # destination convention: out_weight = moves in
  expect_equal(unname(fit$out_weight["C"]), 6)
  expect_equal(unname(fit$in_weight["C"]), 1)
})

test_that("edge list duplicates are summed and isolated nodes kept", {
  fit <- springrank_edgelist(c("a", "a", "b"), c("b", "b", "c"), weight = 1,
                             nodes = c("a", "b", "c", "z"), alpha = 1)
  expect_equal(fit$A["a", "b"], 2)
  expect_equal(unname(fit$ranks["z"]), 0)
  expect_true(fit$ranks["a"] > fit$ranks["b"] && fit$ranks["b"] > fit$ranks["c"])
})

test_that("predictions are complementary and rescaling works", {
  g <- sr_generate(50, beta = 2, K = 6, seed = 3)
  fit <- springrank(g$A, alpha = 0.1)
  b <- as.numeric(sr_beta(fit))
  p1 <- sr_predict(fit, 1:10, 11:20, beta = b)
  p2 <- sr_predict(fit, 11:20, 1:10, beta = b)
  expect_equal(p1 + p2, rep(1, 10), ignore_attr = TRUE)
  r <- sr_rescale(fit, 0.75, beta = b)
  # one rescaled unit = 75% probability
  expect_equal(stats::plogis(2 * b * (1 / (2 * b / stats::qlogis(0.75)))), 0.75)
  expect_length(r, 50)
})

test_that("beta estimates recover the generating steepness roughly", {
  g <- sr_generate(300, beta = 2, alpha = 1, K = 10, seed = 4)
  fit <- springrank(g$A, alpha = 0)
  b <- as.numeric(sr_beta(fit))
  expect_gt(b, 1)
  expect_lt(b, 4)
  bl <- sr_beta(fit, type = "local")
  expect_true(attr(bl, "sigma_a") > 0.5)
})

test_that("perfect hierarchy returns max_beta with a warning", {
  A <- matrix(0, 3, 3); A[1, 2] <- 2; A[2, 3] <- 1; A[1, 3] <- 1
  fit <- springrank(A, alpha = 0)
  expect_warning(b <- sr_beta(fit), "infinite")
  expect_equal(as.numeric(b), 20)
})

test_that("violations count edges against the ranking", {
  m <- toy()
  fit <- springrank_moves(m$origin, m$destination, m$n, alpha = 0)
  v <- sr_violations(fit)
  # downward moves: M->P (2), C->M (1)
  expect_equal(v$against, 3)
  expect_equal(v$total, 15)
})

test_that("significance test separates hierarchy from random directions", {
  g <- sr_generate(80, beta = 3, K = 8, seed = 5)
  fit <- springrank(g$A, alpha = 0.5)
  expect_lt(sr_test(fit, n_null = 99, seed = 1)$p_value, 0.05)
  # symmetric network: no hierarchy
  set.seed(6)
  B <- matrix(rpois(80 * 80, 0.1), 80); B <- B + t(B); diag(B) <- 0
  fitB <- springrank(B, alpha = 0.5)
  expect_gt(sr_test(fitB, n_null = 99, seed = 1)$p_value, 0.2)
})

test_that("bootstrap SEs are smaller for nodes with more edges", {
  g <- sr_generate(100, beta = 2, K = 6, seed = 7)
  fit <- springrank(g$A, alpha = 1)
  bs <- sr_bootstrap(fit, R = 60, seed = 1)
  deg <- (fit$out_weight + fit$in_weight)[bs$node]
  expect_true(all(bs$se >= 0))
  expect_lt(stats::cor(deg, bs$se, method = "spearman"), 0)
})

test_that("input validation", {
  expect_error(springrank(matrix(-1, 2, 2)), "negative")
  expect_error(springrank(matrix(1, 2, 3)), "square")
  expect_error(springrank(diag(2), alpha = -1), "alpha")
})

test_that("numeric node ids are not mangled", {
  fit <- springrank_edgelist(c(1e5, 2e5), c(2e5, 3e5),
                             nodes = c(100000L, 200000L, 300000L), alpha = 1)
  expect_equal(names(fit$ranks), c("100000", "200000", "300000"))
  big <- c(1234567890123456, 1234567890123457)
  fit2 <- springrank_edgelist(big, c(1, 1), alpha = 1)
  expect_length(fit2$ranks, 3)
  expect_true(all(c("1234567890123456", "1234567890123457") %in% names(fit2$ranks)))
})

test_that("sr_test rejects non-integer weights", {
  fit <- springrank(matrix(c(0, 1.5, 0.5, 0), 2), alpha = 1)
  expect_error(sr_test(fit, n_null = 5), "integer")
})

test_that("beta is 0 when all ranks are tied", {
  fit <- springrank(matrix(c(0, 2, 2, 0), 2), alpha = 1)
  expect_warning(b <- sr_beta(fit), "not identified")
  expect_equal(as.numeric(b), 0)
  expect_warning(b2 <- sr_beta(fit, type = "local"), "not identified")
  expect_equal(as.numeric(b2), 0)
})

test_that("global beta solves the paper's estimating equation (Eq. S39)", {
  g <- sr_generate(120, beta = 1.5, alpha = 1, K = 8, seed = 8)
  fit <- springrank(g$A, alpha = 0.2)
  b <- as.numeric(sr_beta(fit))
  A <- as.matrix(fit$A); s <- unname(fit$ranks)
  D <- outer(s, s, "-")
  eq <- sum(D * (A - (A + t(A)) * stats::plogis(2 * b * D)))
  expect_lt(abs(eq), 1e-6)
  # sigma_a of Eq. 12, summed over all ordered pairs
  M <- sum(A)
  sig <- 1 - sum(abs(A - (A + t(A)) * stats::plogis(2 * b * D))) / (2 * M)
  expect_equal(attr(sr_beta(fit), "sigma_a"), sig, tolerance = 1e-10)
  # energy against the dense formula
  expect_equal(sr_energy(fit), 0.5 * sum(A * (D - 1)^2), tolerance = 1e-10)
})
