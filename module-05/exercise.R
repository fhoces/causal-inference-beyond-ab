# ============================================================================
# Module 5 Exercise: Synthetic DiD - Reconstructing the Dual-Weight Estimator
# ============================================================================
# Each question states a task; the solution follows. Try before reading on.
# Goal: (a) rebuild the SDID point estimate EXACTLY from the package's own
# weights via the weighted double difference; (b) re-solve BOTH weight
# programs yourself with quadprog and compare to the package's Frank-Wolfe.

suppressMessages({
  library(tidyverse)
  library(synthdid)
  library(quadprog)
})

# ===== Q1. The empirical comparison (California Prop 99) =====
# Task: reshape the canonical panel, then compute DiD, SC, and SDID with the
# synthdid package plus placebo SEs. Reproduce the AER paper's headline
# ordering: DiD most negative and least precise, SDID between and sharpest.

data("california_prop99", package = "synthdid")
setup <- panel.matrices(california_prop99)      # long -> Y (block), N0, T0, W
Y <- setup$Y; N0 <- setup$N0; T0 <- setup$T0
N <- nrow(Y); T <- ncol(Y); N1 <- N - N0; T1 <- T - T0

tau_sdid <- synthdid_estimate(Y, N0, T0)
tau_sc   <- sc_estimate(Y, N0, T0)
tau_did  <- did_estimate(Y, N0, T0)
ests <- list(SDID = tau_sdid, SC = tau_sc, DID = tau_did)

set.seed(1)   # the placebo variance resamples control units; seed for reproducibility
comp <- imap_dfr(ests, ~ tibble(
  Estimator = .y, Estimate = as.numeric(.x),
  SE = sqrt(vcov(.x, method = "placebo"))))
cat("Q1. Prop 99: DiD vs SC vs SDID (packs per capita)\n")
print(comp |> mutate(across(where(is.numeric), ~ round(.x, 2))))
sdid_se <- comp$SE[comp$Estimator == "SDID"]
did_se  <- comp$SE[comp$Estimator == "DID"]
stopifnot(
  as.numeric(tau_did) < as.numeric(tau_sc),          # DiD most negative (deterministic)
  as.numeric(tau_sc)  < as.numeric(tau_sdid),        # SDID least negative (deterministic)
  sdid_se < did_se                                   # SDID much sharper than DiD
)

# ===== Q2. Reconstruct SDID exactly from the package's own weights =====
# Task: pull omega and lambda off the estimate object and apply the weighted
# double-difference formula by hand. Assert equality to ~1e-10.

w <- attr(tau_sdid, "weights")
omega <- w$omega                      # length N0, on the simplex
lambda <- w$lambda                    # length T0, on the simplex

Y_co <- Y[1:N0, , drop = FALSE]; Y_tr <- Y[(N0 + 1):N, , drop = FALSE]
pre <- 1:T0; post <- (T0 + 1):T

double_diff <- function(omega, lambda) {
  tr_post <- mean(colMeans(Y_tr[, post, drop = FALSE]))               # treated, post (uniform)
  tr_pre  <- sum(lambda * colMeans(Y_tr[, pre, drop = FALSE]))        # treated, lambda-pre
  co_post <- mean(as.numeric(omega %*% Y_co[, post, drop = FALSE]))   # omega-ctrl, post
  co_pre  <- sum(lambda * as.numeric(omega %*% Y_co[, pre, drop = FALSE]))  # omega,lambda-ctrl,pre
  (tr_post - tr_pre) - (co_post - co_pre)
}

tau_rebuilt <- double_diff(omega, lambda)
cat(sprintf("\nQ2. package tau = %.10f | rebuilt = %.10f | diff = %.2e\n",
            as.numeric(tau_sdid), tau_rebuilt, abs(tau_rebuilt - as.numeric(tau_sdid))))
stopifnot(abs(tau_rebuilt - as.numeric(tau_sdid)) < 1e-10)

# ===== Q3. Compute the regularization parameter zeta from scratch =====
# Task: zeta = (N1 * T1)^(1/4) * sigma_hat, where sigma_hat is the SD of
# first-differenced control outcomes over the pre-period (paper's formula,
# denominator N0*(T0-1)).

delta   <- Y_co[, 2:T0, drop = FALSE] - Y_co[, 1:(T0 - 1), drop = FALSE]  # first diffs, pre
delta_v <- as.numeric(delta)
sigma_hat <- sqrt(sum((delta_v - mean(delta_v))^2) / (N0 * (T0 - 1)))
zeta <- (N1 * T1)^(1 / 4) * sigma_hat
cat(sprintf("\nQ3. sigma_hat = %.4f | zeta = %.4f\n", sigma_hat, zeta))
stopifnot(sigma_hat > 0, zeta > sigma_hat)   # (N1*T1)^(1/4) > 1 here

# ===== Q4. Solve the UNIT-weight program yourself with quadprog =====
# Task: minimize ||intercept + Y_ctrl' w - treated_pre_avg||^2 + zeta^2*T0*||w||^2
# over w >= 0, sum(w) = 1, with a free intercept. Compare to the package.

solve_weights <- function(Xctrl, target, eta) {
  # Xctrl: (n obs) x (k candidate weights). Free intercept, ridge eta on w only.
  k <- ncol(Xctrl)
  X <- cbind(1, Xctrl)
  P <- diag(c(0, rep(1, k)))                          # do not penalize intercept
  Dmat <- 2 * (t(X) %*% X + eta * P) + diag(1e-8, k + 1)  # jitter -> positive definite
  dvec <- 2 * as.numeric(t(X) %*% target)
  Amat <- cbind(c(0, rep(1, k)), rbind(0, diag(k)))   # sum(w)=1 (eq), then w>=0
  bvec <- c(1, rep(0, k))
  solve.QP(Dmat, dvec, Amat, bvec, meq = 1)$solution[-1]  # drop intercept
}

target_omega <- colMeans(Y_tr[, pre, drop = FALSE])   # treated pre-period path (length T0)
X_omega      <- t(Y_co[, pre, drop = FALSE])          # T0 x N0
omega_qp     <- solve_weights(X_omega, target_omega, eta = zeta^2 * T0)

cat(sprintf("\nQ4. unit weights: cor(QP, package) = %.5f | sum = %.6f | min = %.2e\n",
            cor(omega_qp, omega), sum(omega_qp), min(omega_qp)))
stopifnot(cor(omega_qp, omega) > 0.99, abs(sum(omega_qp) - 1) < 1e-6,
          min(omega_qp) > -1e-8)

# ===== Q5. Solve the TIME-weight program yourself with quadprog =====
# Task: minimize ||intercept + Y_ctrl_pre lambda - ctrl_post_avg||^2 over the
# simplex (time program uses only a negligible ridge). Compare to the package.

target_lambda <- rowMeans(Y_co[, post, drop = FALSE])  # each control's post avg (length N0)
X_lambda      <- Y_co[, pre, drop = FALSE]             # N0 x T0
lambda_qp     <- solve_weights(X_lambda, target_lambda, eta = (1e-6 * sigma_hat)^2 * N0)

cat(sprintf("Q5. time weights: cor(QP, package) = %.5f | sum = %.6f\n",
            cor(lambda_qp, lambda), sum(lambda_qp)))
stopifnot(cor(lambda_qp, lambda) > 0.99, abs(sum(lambda_qp) - 1) < 1e-6)

# ===== Q6. Your weights -> tau, vs the package =====
# Task: feed YOUR quadprog weights through the same double-difference formula.
# The package uses Frank-Wolfe, not a QP solver, so agreement is CLOSE BUT NOT
# EXACT. Report the gap honestly (expect ~1e-3, not machine precision).

tau_qp <- double_diff(omega_qp, lambda_qp)
cat(sprintf("\nQ6. package tau = %.5f | quadprog tau = %.5f | abs diff = %.2e\n",
            as.numeric(tau_sdid), tau_qp, abs(tau_qp - as.numeric(tau_sdid))))
cat("    (Frank-Wolfe vs QP: same program, different optimizer -> close, not identical.)\n")
stopifnot(abs(tau_qp - as.numeric(tau_sdid)) < 0.05)   # close, but NOT 1e-10

# ===== Q7. Factor structure: SDID recovers truth, DiD does not =====
# Task: simulate a 21-city x 48-month panel with a 2-factor interactive-FE
# structure (parallel trends fails by construction) and a policy at month 36
# with true effect 4.0. Show DiD is badly biased while SC/SDID recover truth.

make_rideshare <- function(n_donors = 20, n_t = 48, g = 36, seed = 7) {
  set.seed(seed)
  N  <- n_donors + 1                                              # city 1 = treated
  f1 <- seq(0, 6, length.out = n_t) + cumsum(rnorm(n_t, 0, 0.3))  # trending factor
  f2 <- 2 * sin(2 * pi * (1:n_t) / 12) + cumsum(rnorm(n_t, 0, 0.2))  # seasonal factor
  L1 <- c(1.5, runif(n_donors, 0.2, 1.7))                        # treated loads HIGH on f1
  L2 <- c(1.0, runif(n_donors, 0.2, 1.6))
  mu <- c(30, runif(n_donors, 20, 40))                           # city levels
  expand_grid(city = 1:N, t = 1:n_t) |>
    mutate(post = t >= g, D = (city == 1) & post,
           y = mu[city] + L1[city] * f1[t] + L2[city] * f2[t] +
               if_else(D, 4.0, 0) + rnorm(n(), 0, 1))            # true effect = 4.0
}

rs <- make_rideshare()
rs_setup <- panel.matrices(as.data.frame(mutate(rs, w = as.integer(D))),
                           unit = "city", time = "t", outcome = "y", treatment = "w")
rs_ests <- list(
  SDID = synthdid_estimate(rs_setup$Y, rs_setup$N0, rs_setup$T0),
  SC   = sc_estimate(rs_setup$Y, rs_setup$N0, rs_setup$T0),
  DID  = did_estimate(rs_setup$Y, rs_setup$N0, rs_setup$T0))

set.seed(1)   # placebo variance resampling
rs_tab <- imap_dfr(rs_ests, ~ tibble(
  Estimator = .y, Estimate = as.numeric(.x),
  SE = sqrt(vcov(.x, method = "placebo")),
  abs_error = abs(as.numeric(.x) - 4.0)))
cat("\nQ7. Factor-structure panel (truth = 4.0):\n")
print(rs_tab |> mutate(across(where(is.numeric), ~ round(.x, 3))))

err <- setNames(rs_tab$abs_error, rs_tab$Estimator)
stopifnot(
  err["DID"]  > 1.5,          # DiD badly biased by the factor structure
  err["SDID"] < 0.5,          # SDID recovers truth
  err["SDID"] < err["DID"],   # SDID closer than DiD
  err["SC"]   < err["DID"]    # SC also beats DiD here
)

cat("\nAll checks passed: SDID rebuilt from its weights, and re-solved from scratch.\n")
