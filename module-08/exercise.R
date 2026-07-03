# ============================================================================
# Module 8 Exercise: Matrix Completion - Hand-Coding Soft-Impute
# ============================================================================
# Each question states a task; the solution follows. Try before reading on.
# Goal: (a) hand-code the soft-impute algorithm and match gsynth's own MC
# estimator; (b) show the extrapolation-bias structural fact on a noiseless
# corner; (c) run the placebo-imputation tournament from scratch; (d) bridge
# to the SDID block-panel setting from Module 5.

suppressMessages({
  library(tidyverse)
  library(fixest)
  library(did)
  library(gsynth)
  library(synthdid)
})

two_way_fe <- function(Y, obs) {
  a <- rep(0, nrow(Y)); b <- rep(0, ncol(Y)); m <- 0
  for (k in 1:300) {
    m_old <- m
    m <- mean(Y[obs] - outer(a, rep(1, ncol(Y)))[obs] - outer(rep(1, nrow(Y)), b)[obs])
    R <- Y - m - outer(rep(1, nrow(Y)), b)
    a <- sapply(1:nrow(Y), function(i) mean(R[i, obs[i, ]]))
    R <- Y - m - outer(a, rep(1, ncol(Y)))
    b <- sapply(1:ncol(Y), function(t) mean(R[obs[, t], t]))
    if (abs(m - m_old) < 1e-10 && k > 3) break
  }
  m + outer(a, rep(1, ncol(Y))) + outer(rep(1, nrow(Y)), b)
}

soft_impute <- function(Yr, obs, lambda, tol = 1e-8, maxit = 2000) {
  L <- matrix(0, nrow(Yr), ncol(Yr)); d <- 0
  for (it in 1:maxit) {
    Z <- ifelse(obs, Yr, L)
    s <- svd(Z); d <- pmax(s$d - lambda, 0)
    L_new <- s$u %*% diag(d) %*% t(s$v)
    if (sum((L_new - L)^2) < tol * max(sum(L^2), 1)) { L <- L_new; break }
    L <- L_new
  }
  list(L = L, rank = sum(d > 1e-10))
}

mc_impute <- function(Y, obs, lambdas = 2^seq(-2, 4), cv_frac = 0.1, seed = 1) {
  FEmat <- two_way_fe(Y, obs); Yr <- Y - FEmat
  set.seed(seed)
  oi <- which(obs); hold <- sample(oi, round(cv_frac * length(oi)))
  obs_cv <- obs; obs_cv[hold] <- FALSE
  FE_cv <- two_way_fe(Y, obs_cv); Yr_cv <- Y - FE_cv
  errs <- sapply(lambdas, function(l) {
    fit <- soft_impute(Yr_cv, obs_cv, l)
    sqrt(mean((Yr_cv[hold] - fit$L[hold])^2))
  })
  lam <- lambdas[which.min(errs)]
  fit <- soft_impute(Yr, obs, lam)
  list(Y0 = FEmat + fit$L, lambda = lam, rank = fit$rank)
}

gs_impute_rmse <- function(Y, mask, estimator = "ife") {
  d <- expand_grid(city = 1:nrow(Y), t = 1:ncol(Y)) |>
    mutate(y = Y[cbind(city, t)], D = as.integer(mask[cbind(city, t)]))
  tryCatch({
    co <- capture.output(
      g <- if (estimator == "ife")
        gsynth(y ~ D, data = as.data.frame(d), index = c("city", "t"),
               force = "two-way", CV = TRUE, r = c(0, 4), se = FALSE)
      else
        gsynth(y ~ D, data = as.data.frame(d), index = c("city", "t"),
               force = "two-way", estimator = "mc", se = FALSE))
    err <- t(g$eff)
    sqrt(mean(err[mask]^2, na.rm = TRUE))
  }, error = function(e) NA_real_)
}

# ===== Q1. The Demo A panel: staggered launch, selection on growth =====
# Task: build the panel and confirm its treatment structure matches the
# slides (25 treated cities, at least 10 distinct launch cohorts, 35
# never-treated).

make_att_panel <- function(n_cities = 60, n_t = 40, n_treated = 25,
                           tau = 2.0, seed = 13, noise = 0.5) {
  set.seed(seed)
  f1 <- 0.12 * (1:n_t) + as.numeric(arima.sim(list(ar = 0.85), n_t, sd = 0.5))
  f2 <- 1.2 * sin(2 * pi * (1:n_t) / 12)
  L1 <- runif(n_cities, 0.5, 2.0); L2 <- runif(n_cities, 0.2, 1.5)
  mu <- runif(n_cities, 20, 40)
  score <- L1 + rnorm(n_cities, 0, 0.4)
  treated <- rank(-score) <= n_treated
  g <- rep(Inf, n_cities)
  g[treated] <- round(34 - (L1[treated] - 0.5) / 1.5 * 19) +
    sample(-2:2, n_treated, replace = TRUE)
  g[treated] <- pmax(pmin(g[treated], 36), 13)
  expand_grid(city = 1:n_cities, t = 1:n_t) |>
    mutate(g = g[city], D = as.integer(t >= g),
           y0 = mu[city] + L1[city] * f1[t] + L2[city] * f2[t] + rnorm(n(), 0, noise),
           y = y0 + tau * D)
}

att_dat <- make_att_panel()
n_treated_cities <- n_distinct(att_dat$city[att_dat$D == 1])
n_cohorts <- n_distinct(att_dat$g[is.finite(att_dat$g)])
n_never <- n_distinct(att_dat$city) - n_treated_cities
cat(sprintf("Q1. %d treated cities, %d launch cohorts, %d never-treated\n",
            n_treated_cities, n_cohorts, n_never))
stopifnot(n_treated_cities == 25, n_cohorts >= 10, n_never == 35)

# ===== Q2. Hand-code soft-impute; show the extrapolation-bias fact =====
# Task: on a NOISELESS rank-2 matrix, a random 20% mask (no selection) should
# be trivial to impute (RMSE near zero): the missing cells' rows are still
# observed elsewhere, so their loadings are easy to pin down. The SAME total
# amount of missingness, but concentrated where the TOP-loading rows miss
# their last third of columns (selection on the loading, exactly how Demo A's
# treated corner is built), should be much harder: that gap IS the
# extrapolation lesson of the module, not "corners are hard" but "missingness
# correlated with an extreme, under-identified part of the loading space is
# hard".

set.seed(21)
N <- 30; Tt <- 30
loadings1 <- runif(N, 0.5, 2.0); loadings2 <- runif(N, 0.2, 1.5)
v1 <- 0.15 * (1:Tt); v2 <- sin(2 * pi * (1:Tt) / 8)
Ytest <- outer(loadings1, v1) + outer(loadings2, v2)   # exactly rank 2, no noise

sel <- rank(-loadings1) <= 10                          # top-loading rows
mask_selection <- matrix(FALSE, N, Tt); mask_selection[sel, (Tt - 9):Tt] <- TRUE
set.seed(22)
mask_random <- matrix(FALSE, N, Tt)
mask_random[sample(N * Tt, sum(mask_selection))] <- TRUE  # same missing count

fit_random <- soft_impute(Ytest, !mask_random, lambda = 0.05)
rmse_random <- sqrt(mean((Ytest[mask_random] - fit_random$L[mask_random])^2))

fit_selection <- soft_impute(Ytest, !mask_selection, lambda = 0.05)
rmse_selection <- sqrt(mean((Ytest[mask_selection] - fit_selection$L[mask_selection])^2))

cat(sprintf("\nQ2. random-mask RMSE = %.4f | selection-mask RMSE = %.4f\n",
            rmse_random, rmse_selection))
stopifnot(rmse_random < 0.05, rmse_selection > 1)

# ===== Q3. ATT on Demo A: soft-impute vs gsynth's own MC estimator =====
# Task: they implement the same estimator with different CV/optimization
# details, so expect close agreement, not machine precision. Both should
# land near TWFE, not near the truth (the conservatism fact from the
# slides).

twfe_att <- coef(feols(y ~ D | city + t, data = att_dat))[["D"]]

Y  <- att_dat |> select(city, t, y) |> pivot_wider(names_from = t, values_from = y) |>
  column_to_rownames("city") |> as.matrix()
Dm <- att_dat |> select(city, t, D) |> pivot_wider(names_from = t, values_from = D) |>
  column_to_rownames("city") |> as.matrix()

mc <- mc_impute(Y, Dm == 0)
si_att <- mean(Y[Dm == 1] - mc$Y0[Dm == 1])

co <- capture.output(
  g_mc <- gsynth(y ~ D, data = as.data.frame(att_dat), index = c("city", "t"),
                force = "two-way", estimator = "mc", se = FALSE))
gsynth_mc_att <- g_mc$att.avg

cat(sprintf("\nQ3. TWFE = %.3f | soft-impute ATT = %.3f | gsynth-MC ATT = %.3f\n",
            twfe_att, si_att, gsynth_mc_att))
stopifnot(
  abs(si_att - gsynth_mc_att) < 0.2,           # same estimator, different optimizer
  abs(si_att - twfe_att) < 0.6,                # both land near TWFE
  abs(gsynth_mc_att - twfe_att) < 0.6
)

# ===== Q4. gsynth IFE recovers the truth; TWFE and CS do not =====
# Task: fit gsynth with CV over the rank, and Callaway-Sant'Anna. IFE
# should land close to the true effect (2.0); TWFE and CS should not.

co <- capture.output(
  g_ife <- gsynth(y ~ D, data = as.data.frame(att_dat), index = c("city", "t"),
                  force = "two-way", CV = TRUE, r = c(0, 4), se = FALSE))

d_cs <- att_dat |> mutate(g0 = ifelse(is.finite(g), g, 0)) |>
  select(city, t, y, g0) |> as.data.frame()
cs <- suppressWarnings(att_gt("y", tname = "t", idname = "city", gname = "g0",
      data = d_cs, control_group = "nevertreated", base_period = "universal",
      est_method = "dr", bstrap = FALSE, cband = FALSE))
cs_att <- suppressWarnings(aggte(cs, type = "simple")$overall.att)

truth <- 2.0
cat(sprintf("\nQ4. gsynth-IFE = %.3f (r=%d) | CS = %.3f | TWFE = %.3f | truth = %.1f\n",
            g_ife$att.avg, g_ife$r.cv, cs_att, twfe_att, truth))
stopifnot(
  abs(g_ife$att.avg - truth) < 0.4,
  twfe_att - truth > 1,
  cs_att - truth > 0.8
)

# ===== Q5. Cross-validate lambda for the soft-impute fit =====
# Task: mask 10% of observed cells, refit at each candidate lambda, pick
# the one minimizing held-out RMSE. Confirm the chosen value is interior
# to the search grid, not pinned at either endpoint.

lambdas <- 2^seq(-2, 4)
cat(sprintf("\nQ5. lambda grid: %s\n", paste(round(lambdas, 2), collapse = ", ")))
cat(sprintf("    chosen lambda = %.2f (rank %d)\n", mc$lambda, mc$rank))
stopifnot(mc$lambda > min(lambdas), mc$lambda < max(lambdas))

# ===== Q6. The placebo-imputation tournament =====
# Task: on a no-treatment factor panel, hide observed cells under four
# patterns and score each method's imputation RMSE against the known
# truth. Soft-impute should beat FE-only everywhere, including the one
# pattern (random holes) where gsynth cannot even run.

make_factor_panel <- function(n_cities = 40, n_t = 40, seed = 11, noise = 0.5) {
  set.seed(seed)
  f1 <- 0.12 * (1:n_t) + as.numeric(arima.sim(list(ar = 0.85), n_t, sd = 0.5))
  f2 <- 1.2 * sin(2 * pi * (1:n_t) / 12)
  L1 <- runif(n_cities, 0.5, 2.0); L2 <- runif(n_cities, 0.2, 1.5)
  mu <- runif(n_cities, 20, 40)
  outer(mu, rep(1, n_t)) + outer(L1, f1) + outer(L2, f2) +
    matrix(rnorm(n_cities * n_t, 0, noise), n_cities)
}

set.seed(11)
f1p <- 0.12 * (1:40) + as.numeric(arima.sim(list(ar = 0.85), 40, sd = 0.5))
L1p <- runif(40, 0.5, 2.0)
Y0 <- make_factor_panel()
Np <- nrow(Y0); Tp <- ncol(Y0)

tourn_masks <- list()
m <- matrix(FALSE, Np, Tp); m[which.max(L1p), 31:40] <- TRUE
tourn_masks[["one unit, block"]] <- m
set.seed(2); tr <- order(-L1p)[1:10]
m <- matrix(FALSE, Np, Tp); m[tr, 31:40] <- TRUE
tourn_masks[["10 units, block"]] <- m
set.seed(3); tr <- order(-L1p + rnorm(Np, 0, 0.4))[1:20]
gg <- sample(12:36, 20, replace = TRUE)
m <- matrix(FALSE, Np, Tp); for (i in seq_along(tr)) m[tr[i], gg[i]:Tp] <- TRUE
tourn_masks[["20 units, staggered"]] <- m
set.seed(4)
m <- matrix(FALSE, Np, Tp); m[sample(Np * Tp, round(0.15 * Np * Tp))] <- TRUE
tourn_masks[["15% random holes"]] <- m

tourn_res <- imap_dfr(tourn_masks, function(mask, nm) {
  obs <- !mask
  fe_hat <- two_way_fe(Y0, obs)
  fe_rmse <- sqrt(mean((Y0[mask] - fe_hat[mask])^2))
  mc_t <- mc_impute(Y0, obs)
  mc_rmse <- sqrt(mean((Y0[mask] - mc_t$Y0[mask])^2))
  absorbing <- all(apply(mask, 1, function(r) all(diff(r) >= 0)))
  ife_rmse <- if (absorbing) gs_impute_rmse(Y0, mask, "ife") else NA_real_
  tibble(pattern = nm, fe_only = fe_rmse, gsynth_ife = ife_rmse, soft_impute = mc_rmse)
})
cat("\nQ6. Placebo-imputation tournament (RMSE):\n")
print(tourn_res |> mutate(across(where(is.numeric), ~ round(.x, 3))))

structured <- tourn_res |> filter(pattern != "15% random holes")
holes <- tourn_res |> filter(pattern == "15% random holes")
stopifnot(
  all(structured$soft_impute < structured$fe_only),
  is.na(holes$gsynth_ife),
  holes$soft_impute < holes$fe_only
)

# ===== Q7. The SDID bridge: block panel, one simultaneous cohort =====
# Task: same factor structure as Demo A/B, but 10 cities treated together
# in a simple block (last 10 weeks), for direct comparison with Module 5's
# synthdid machinery. Compare DiD, SDID, gsynth-IFE, and soft-impute.

make_block_panel <- function(n_cities = 60, n_t = 40, n_treated = 10,
                             tau = 2.0, seed = 8, noise = 0.5) {
  set.seed(seed)
  f1 <- 0.12 * (1:n_t) + as.numeric(arima.sim(list(ar = 0.85), n_t, sd = 0.5))
  f2 <- 1.2 * sin(2 * pi * (1:n_t) / 12)
  L1 <- runif(n_cities, 0.5, 2.0); L2 <- runif(n_cities, 0.2, 1.5)
  mu <- runif(n_cities, 20, 40)
  score <- L1 + rnorm(n_cities, 0, 0.4)
  treated <- rank(-score) <= n_treated
  g_block <- n_t - 9
  expand_grid(city = 1:n_cities, t = 1:n_t) |>
    mutate(is_tr = treated[city], D = as.integer(is_tr & t >= g_block),
           y0 = mu[city] + L1[city] * f1[t] + L2[city] * f2[t] + rnorm(n(), 0, noise),
           y = y0 + tau * D)
}

bp <- make_block_panel()
twfe_b <- coef(feols(y ~ D | city + t, data = bp))[["D"]]
bp_setup <- panel.matrices(as.data.frame(bp), unit = "city", time = "t",
                           outcome = "y", treatment = "D")
sdid_b <- as.numeric(synthdid_estimate(bp_setup$Y, bp_setup$N0, bp_setup$T0))

co <- capture.output(
  g_ife_b <- gsynth(y ~ D, data = as.data.frame(bp), index = c("city", "t"),
                    force = "two-way", CV = TRUE, r = c(0, 4), se = FALSE))

Yb  <- bp |> select(city, t, y) |> pivot_wider(names_from = t, values_from = y) |>
  column_to_rownames("city") |> as.matrix()
Dmb <- bp |> select(city, t, D) |> pivot_wider(names_from = t, values_from = D) |>
  column_to_rownames("city") |> as.matrix()
mc_b <- mc_impute(Yb, Dmb == 0)
si_att_b <- mean(Yb[Dmb == 1] - mc_b$Y0[Dmb == 1])

cat(sprintf("\nQ7. Block panel (10 treated, truth 2.0): DiD %.2f | SDID %.2f | gsynth-IFE %.2f | soft-impute %.2f\n",
            twfe_b, sdid_b, g_ife_b$att.avg, si_att_b))
stopifnot(
  abs(sdid_b - truth) < 0.5,
  abs(g_ife_b$att.avg - truth) < 0.5,
  abs(twfe_b - truth) > abs(sdid_b - truth),
  abs(twfe_b - truth) > abs(g_ife_b$att.avg - truth)
)

cat("\nAll checks passed: soft-impute matches gsynth-MC, the extrapolation\n")
cat("gap is confirmed on a noiseless corner, the tournament favors matrix\n")
cat("completion only on the irregular pattern, and the SDID bridge holds.\n")
