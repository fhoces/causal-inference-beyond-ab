# ============================================================================
# Module 7 Exercise: Policy Learning - Trees, DR Values, and the Optimism Gap
# ============================================================================
# Each question states a task; the solution follows. Try before reading on.
# Goal: learn a cost-aware depth-2 policy tree from causal-forest DR scores,
# then hand-build a 5-fold cross-fitted DR value estimator and use it to show
# the optimism of same-data evaluation and to rank policies honestly.
#
# Runtime note: the 5-fold loop fits 10 forests; expect ~20-30s. Progress is
# printed per fold.

suppressMessages({
  library(tidyverse)
  library(grf)
  library(policytree)
})

# Shared HTE DGP (identical to Module 6 and the slides)
make_hte_data <- function(n = 6000, seed = 42) {
  set.seed(seed)
  tibble(
    density  = runif(n),                        # city density percentile
    tenure   = pmin(rexp(n, 1 / 18), 90),       # months on platform
    peak_shr = rbeta(n, 2, 3),                  # share of hours in peak
    rating   = pmin(pmax(rnorm(n, 4.7, 0.2), 3.5), 5),
    w        = rbinom(n, 1, 0.5),
    tau      = 1.5 * density - 0.02 * tenure + 2.0 * density * peak_shr,
    y        = 10 + 5 * density + 0.05 * tenure + 3 * peak_shr +
               tau * w + rnorm(n, 0, 2)
  )
}

dat  <- make_hte_data(6000)
cost <- 0.6
e    <- 0.5
covs <- c("density", "tenure", "peak_shr", "rating")
X <- as.matrix(dat[, covs]); Y <- dat$y; W <- dat$w

# ---- helpers ----------------------------------------------------------------
# manual AIPW (doubly-robust) scores for held-out rows, from training nuisances
aipw_scores <- function(cf, yf, newX, newY, newW, e = 0.5) {
  tau <- predict(cf, newX)$predictions        # effect from the causal forest
  m   <- predict(yf, newX)$predictions        # marginal outcome E[Y|X]
  mu1 <- m + (1 - e) * tau                     # implied E[Y|X, W=1]
  mu0 <- m - e * tau                           # implied E[Y|X, W=0]
  cbind(control = mu0 + (1 - newW) / (1 - e) * (newY - mu0),
        treated = mu1 +      newW  /      e  * (newY - mu1))
}
# DR value of a policy (actions in {1,2}), net of cost, with influence-fn SE
dr_value <- function(action, G, cost) {
  chosen <- ifelse(action == 2, G[, "treated"] - cost, G[, "control"])
  c(value = mean(chosen), se = sd(chosen) / sqrt(length(chosen)))
}

# ===== Q1. Learn a cost-aware depth-2 policy tree =====
# Task: fit a causal forest, build doubly-robust scores, subtract the per-push
# cost from the treated column, and learn a depth-2 policy_tree. Read the rule
# and compare its treated share to the true optimal share (tau > cost).

cf      <- causal_forest(X, Y, W, W.hat = e, num.trees = 2000, seed = 42)
tau_hat <- predict(cf)$predictions
Gamma   <- double_robust_scores(cf)            # n x 2: control, treated

Gamma_net <- Gamma
Gamma_net[, "treated"] <- Gamma_net[, "treated"] - cost
tree <- policy_tree(X, Gamma_net, depth = 2)
act  <- predict(tree, X)                        # 1 = control, 2 = treat

cat("Q1. Cost-aware depth-2 policy tree:\n")
print(tree)
cat(sprintf("\n    true optimal share (tau > %.1f)   = %.3f\n", cost,
            mean(dat$tau  > cost)))
cat(sprintf("    forest plug-in share (tauhat > c) = %.3f\n",
            mean(tau_hat  > cost)))
cat(sprintf("    policy-tree treated share         = %.3f\n", mean(act == 2)))
stopifnot(mean(act == 2) > 0.3, mean(act == 2) < 0.8)   # sensible, not degenerate

# ===== Q2. The 5-fold cross-fitted DR value estimator =====
# Task: implement cross-fitting yourself. For each fold, fit the nuisances
# (causal forest + outcome forest) on the TRAINING folds, learn a policy tree
# on the training folds, then score the HELD-OUT fold with DR scores built
# from the training nuisances and evaluate each policy there. Pool the
# held-out contributions. This keeps every policy off the data that trained it.

set.seed(42)
K     <- 5
folds <- sample(rep(1:K, length.out = nrow(X)))
acc   <- tibble()                               # pooled held-out influence terms

for (k in 1:K) {
  tr <- folds != k; te <- !tr
  cf_k <- causal_forest(X[tr, ], Y[tr], W[tr], W.hat = e, num.trees = 2000, seed = 1)
  yf_k <- regression_forest(X[tr, ], Y[tr], num.trees = 2000, seed = 1)

  # learn the tree on training-fold DR scores (net of cost)
  G_tr <- double_robust_scores(cf_k); G_tr[, "treated"] <- G_tr[, "treated"] - cost
  tree_k <- policy_tree(X[tr, ], G_tr, depth = 2)

  # held-out DR scores from the TRAINING nuisances
  G_te   <- aipw_scores(cf_k, yf_k, X[te, ], Y[te], W[te], e)
  tau_te <- predict(cf_k, X[te, ])$predictions

  # candidate policies, evaluated on the held-out fold
  a_none   <- rep(1L, sum(te))
  a_all    <- rep(2L, sum(te))
  a_oracle <- ifelse(dat$tau[te] > cost, 2L, 1L)
  a_plugin <- ifelse(tau_te      > cost, 2L, 1L)
  a_tree   <- predict(tree_k, X[te, ])

  score <- function(a) ifelse(a == 2, G_te[, "treated"] - cost, G_te[, "control"])
  acc <- bind_rows(acc, tibble(
    none   = score(a_none),   all    = score(a_all),
    oracle = score(a_oracle), plugin = score(a_plugin), tree = score(a_tree)))
  cat(sprintf("    fold %d done | tree treated share = %.3f\n",
              k, mean(a_tree == 2)))
}

# naive same-data value of the learned tree (fit and scored on all the data)
Gamma_full <- Gamma                              # scores from the Q1 full-data forest
a_naive    <- act
naive_gain <- mean(ifelse(a_naive == 2, Gamma_full[, "treated"] - cost,
                          Gamma_full[, "control"])) - mean(Gamma_full[, "control"])
honest_gain <- mean(acc$tree - acc$none)

cat(sprintf("\nQ2. optimism gap:\n    naive same-data tree gain  = %.3f\n",
            naive_gain))
cat(sprintf("    honest cross-fit tree gain = %.3f\n", honest_gain))
cat(sprintf("    optimism (naive - honest)  = %.3f\n", naive_gain - honest_gain))
stopifnot(naive_gain > honest_gain)              # same-data evaluation is optimistic

# ===== Q3. The honest league table, with SEs and paired tests =====
# Task: report the cross-fitted DR value of each policy as a gain over
# treat-none, with influence-function SEs. Then assert two facts for THIS DGP,
# verified numerically above before asserting: (a) the tree beats treat-all net
# of cost, and (b) the tree's value is within noise of the (infeasible) oracle.

summ <- function(v) c(value = mean(v), se = sd(v) / sqrt(length(v)))
league <- map_dfr(c("none", "all", "oracle", "plugin", "tree"), function(nm) {
  s <- summ(acc[[nm]] - acc$none)                # gain over treat-none
  tibble(policy = nm, gain = s["value"], se = s["se"])
})
cat("\nQ3. Cross-fitted league table (gain over treat-none, net of cost):\n")
print(league |> mutate(across(where(is.numeric), ~ round(.x, 4))))

# paired comparisons (nuisance noise cancels in the differences)
paired <- function(a, b) {
  d <- acc[[a]] - acc[[b]]
  c(diff = mean(d), se = sd(d) / sqrt(length(d)),
    t = mean(d) / (sd(d) / sqrt(length(d))))
}
p_ta <- paired("tree", "all")
p_to <- paired("tree", "oracle")
cat(sprintf("\n    tree - all    : diff = %+.3f  se = %.3f  t = %+.2f\n",
            p_ta["diff"], p_ta["se"], p_ta["t"]))
cat(sprintf("    tree - oracle : diff = %+.3f  se = %.3f  t = %+.2f\n",
            p_to["diff"], p_to["se"], p_to["t"]))

# (a) tree beats treat-all net of cost, decisively
stopifnot(p_ta["diff"] > 0, p_ta["t"] > 3)
# (b) tree is within noise of the oracle: not statistically distinguishable
stopifnot(abs(p_to["t"]) < 3, honest_gain > 0.8 * league$gain[league$policy == "oracle"])

# ===== Q4. Cost sensitivity of the rule =====
# Task: reuse the full-data DR scores and refit the tree at several costs.
# Confirm the treated share falls monotonically as the per-push cost rises,
# tracking the true optimal share (tau > c).

cost_grid <- c(0.0, 0.3, 0.6, 0.9, 1.2)
cost_tbl <- map_dfr(cost_grid, function(cc) {
  G <- Gamma; G[, "treated"] <- G[, "treated"] - cc
  tt <- policy_tree(X, G, depth = 2)
  tibble(cost = cc,
         tree_share = mean(predict(tt, X) == 2),
         opt_share  = mean(dat$tau > cc))
})
cat("\nQ4. Cost sensitivity (treated share vs cost):\n")
print(cost_tbl |> mutate(across(where(is.numeric), ~ round(.x, 3))))
stopifnot(all(diff(cost_tbl$tree_share) <= 0.02))   # weakly decreasing in cost

cat("\nAll checks passed: cost-aware tree learned, cross-fit removes the\n")
cat("optimism, and the interpretable rule matches the oracle within noise.\n")
