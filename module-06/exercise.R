# ============================================================================
# Module 6 Exercise: Causal Forest - Honest Splitting and Pointwise CIs
# ============================================================================
# Each question states a task; the solution follows. Try before reading on.
# Goal: (1) build a partition-based CATE estimator from scratch and use it to
# expose the honest-vs-adaptive winner's curse, then (2) verify that grf's
# pointwise confidence intervals cover the truth near their nominal 95% rate.
# Runtime target: under ~5 minutes. Progress is printed as it goes.

suppressMessages({
  library(tidyverse)
  library(grf)
})
set.seed(1)

# ===== Q1. A partition-based CATE estimator, from scratch =====
# Task: write a recursive binary-split tree (depth 2-3) that maximizes the
# squared difference in leaf treatment effects, plus helpers to (a) push a
# fresh sample down a grown tree and (b) estimate leaf effects on any sample.
# Keeping "grow the structure" separate from "estimate the leaf effects" is
# exactly what will let us run the SAME tree adaptively vs honestly.

leaf_effect <- function(y, w) {                 # difference in means within a leaf
  if (sum(w == 1) < 2 || sum(w == 0) < 2) return(NA_real_)
  mean(y[w == 1]) - mean(y[w == 0])
}

grow_tree <- function(X, y, w, depth, grid, min_leaf = 50) {
  n <- length(y)
  node <- list(leaf = TRUE, effect = leaf_effect(y, w), n = n)
  if (depth == 0 || n < 2 * min_leaf) return(node)
  best <- list(crit = -Inf)
  for (v in 1:ncol(X)) {                        # search variable x threshold
    for (c in grid) {
      L <- X[, v] <= c
      if (sum(L) < min_leaf || sum(!L) < min_leaf) next
      dl <- leaf_effect(y[L], w[L]); dr <- leaf_effect(y[!L], w[!L])
      if (is.na(dl) || is.na(dr)) next
      crit <- sum(L) * sum(!L) * (dl - dr)^2    # heterogeneity-maximizing split
      if (crit > best$crit) best <- list(crit = crit, v = v, c = c)
    }
  }
  if (!is.finite(best$crit)) return(node)
  L <- X[, best$v] <= best$c
  list(leaf = FALSE, v = best$v, c = best$c,
       left  = grow_tree(X[L, , drop = FALSE],  y[L],  w[L],  depth - 1, grid, min_leaf),
       right = grow_tree(X[!L, , drop = FALSE], y[!L], w[!L], depth - 1, grid, min_leaf))
}

assign_leaf <- function(tree, X, path = "r") {  # portable: apply structure to any X
  if (tree$leaf) return(rep(path, nrow(X)))
  L <- X[, tree$v] <= tree$c
  out <- character(nrow(X))
  out[L]  <- assign_leaf(tree$left,  X[L, , drop = FALSE],  paste0(path, "L"))
  out[!L] <- assign_leaf(tree$right, X[!L, , drop = FALSE], paste0(path, "R"))
  out
}

leaf_estimates <- function(tree, X, y, w) {     # estimate effects on a GIVEN sample
  tibble(leaf = assign_leaf(tree, X), y = y, w = w) |>
    group_by(leaf) |>
    summarise(effect = leaf_effect(y, w), n = n(), .groups = "drop")
}

cat("Q1. Partition estimator built (grow_tree / assign_leaf / leaf_estimates).\n")

# ===== Q2. The winner's curse: adaptive vs honest on a CONSTANT effect =====
# Task: simulate a 2D DGP with tau = 2 EVERYWHERE (no real heterogeneity).
# Grow a depth-2 tree two ways over ~200 Monte Carlo reps at n = 2000:
#   - ADAPTIVE: grow the tree AND estimate leaf effects on the same data.
#   - HONEST: grow the tree on half the data, estimate leaf effects on the
#     OTHER half.
# Report the average leaf-effect spread (max leaf minus min leaf). Adaptive
# manufactures heterogeneity out of noise; honest reports ~0.

toy_dgp <- function(n, seed, hetero = FALSE) {
  set.seed(seed)
  x1 <- runif(n); x2 <- runif(n)
  tau <- if (hetero) 1 + 2 * (x1 > 0.5) else rep(2, n)   # constant unless hetero
  w <- rbinom(n, 1, 0.5)
  y <- x1 + x2 + tau * w + rnorm(n, 0, 1)
  list(X = cbind(x1, x2), y = y, w = w, tau = tau)
}

grid  <- seq(0.2, 0.8, by = 0.05)
reps  <- 200
n_toy <- 2000

wc <- map_dfr(seq_len(reps), function(r) {
  d <- toy_dgp(n_toy, seed = r, hetero = FALSE)
  # ADAPTIVE: structure and estimates from the same n rows
  tr_a <- grow_tree(d$X, d$y, d$w, depth = 2, grid = grid)
  est_a <- leaf_estimates(tr_a, d$X, d$y, d$w)
  # HONEST: structure from half A, estimates from held-out half B
  idx <- sample(n_toy, n_toy / 2)
  A <- rep(FALSE, n_toy); A[idx] <- TRUE; B <- !A
  tr_h <- grow_tree(d$X[A, ], d$y[A], d$w[A], depth = 2, grid = grid)
  est_h <- leaf_estimates(tr_h, d$X[B, ], d$y[B], d$w[B])
  if (r %% 50 == 0) cat(sprintf("   ...winner's-curse rep %d/%d\n", r, reps))
  tibble(
    adaptive_spread = max(est_a$effect) - min(est_a$effect),
    honest_spread   = max(est_h$effect) - min(est_h$effect),
    adaptive_rmse   = sqrt(mean((est_a$effect - 2)^2)),   # truth is 2 in every leaf
    honest_rmse     = sqrt(mean((est_h$effect - 2)^2))
  )
})

cat(sprintf("\nQ2. Constant effect tau = 2 everywhere (NO real heterogeneity):\n"))
cat(sprintf("    ADAPTIVE mean leaf spread = %.3f  | leaf RMSE vs truth = %.3f\n",
            mean(wc$adaptive_spread), mean(wc$adaptive_rmse)))
cat(sprintf("    HONEST   mean leaf spread = %.3f  | leaf RMSE vs truth = %.3f\n",
            mean(wc$honest_spread), mean(wc$honest_rmse)))
cat("    -> adaptive invents heterogeneity from noise; honest is centered.\n")
stopifnot(mean(wc$adaptive_spread) > 1.4 * mean(wc$honest_spread))
stopifnot(mean(wc$adaptive_rmse)   > mean(wc$honest_rmse))

# ===== Q3. Honesty does not destroy REAL heterogeneity =====
# Task: rerun on a DGP with a genuine step (tau = 1 for x1 <= 0.5, tau = 3
# above). Confirm both methods recover the true gap of 2.0: when the signal is
# strong and the split is clean, there is little noise to capitalize on, so the
# winner's curse is negligible and honesty costs nothing. Honesty is not
# conservatism; it removes bias where bias exists (Q2), not real signal.

het <- map_dfr(seq_len(100), function(r) {
  d <- toy_dgp(n_toy, seed = 500 + r, hetero = TRUE)
  tr_a <- grow_tree(d$X, d$y, d$w, depth = 1, grid = grid)   # one split suffices
  est_a <- leaf_estimates(tr_a, d$X, d$y, d$w)
  idx <- sample(n_toy, n_toy / 2)
  A <- rep(FALSE, n_toy); A[idx] <- TRUE; B <- !A
  tr_h <- grow_tree(d$X[A, ], d$y[A], d$w[A], depth = 1, grid = grid)
  est_h <- leaf_estimates(tr_h, d$X[B, ], d$y[B], d$w[B])
  tibble(adaptive_gap = max(est_a$effect) - min(est_a$effect),
         honest_gap   = max(est_h$effect) - min(est_h$effect))
})
cat(sprintf("\nQ3. Real step effect (true gap = 3 - 1 = 2.0):\n"))
cat(sprintf("    ADAPTIVE mean estimated gap = %.3f\n", mean(het$adaptive_gap)))
cat(sprintf("    HONEST   mean estimated gap = %.3f\n", mean(het$honest_gap)))
cat("    -> both recover the true gap ~2.0: honesty preserves real signal.\n")
stopifnot(mean(het$honest_gap) > 1.5, mean(het$adaptive_gap) > 1.5)  # both recover it

# ===== Q4. grf pointwise CI coverage over many refits =====
# Task: fit modest causal forests (num.trees = 500, n = 1500) many times on
# fresh data, and at a few FIXED query points check how often the 95%
# interval tau_hat +/- 1.96*se covers the true tau(x). Coverage should sit
# near 0.95. This is the empirical face of the Wager-Athey pointwise CLT.

mk_grf_data <- function(n, seed) {
  set.seed(seed)
  X <- cbind(x1 = runif(n), x2 = runif(n))
  tau <- 1 + 2 * (X[, 1] > 0.5)                 # true CATE: 1 below, 3 above
  w <- rbinom(n, 1, 0.5)
  y <- X[, 1] + X[, 2] + tau * w + rnorm(n, 0, 1)
  list(X = X, y = y, w = w)
}

query    <- rbind(c(0.25, 0.5), c(0.60, 0.5), c(0.85, 0.5))
truth_q  <- c(1, 3, 3)                           # tau at each query point
K        <- 60
hits     <- matrix(NA, K, nrow(query))
for (k in seq_len(K)) {
  d <- mk_grf_data(1500, seed = 2000 + k)
  cf <- causal_forest(d$X, d$y, d$w, W.hat = 0.5, num.trees = 500, seed = k)
  p  <- predict(cf, query, estimate.variance = TRUE)
  se <- sqrt(p$variance.estimates)
  hits[k, ] <- abs(p$predictions - truth_q) <= 1.96 * se
  if (k %% 20 == 0) cat(sprintf("   ...coverage fit %d/%d\n", k, K))
}
cover <- colMeans(hits)
cat(sprintf("\nQ4. grf pointwise 95%% CI coverage over %d refits:\n", K))
for (j in seq_len(nrow(query))) {
  cat(sprintf("    x1 = %.2f (true tau = %d):  coverage = %.2f\n",
              query[j, 1], truth_q[j], cover[j]))
}
cat(sprintf("    mean coverage across points = %.2f (nominal 0.95)\n", mean(cover)))
stopifnot(mean(cover) > 0.80)   # near nominal; honest trees make this hold

cat("\nAll checks passed: honest splitting removes the winner's curse, and\n")
cat("grf's pointwise intervals cover near their nominal rate.\n")
