# ============================================================================
# Module 3 Exercise: Honest DiD - Roth (2022) Pretrend Power and Breakdown Values
# ============================================================================
# Each question states a task; the solution follows. Try before reading on.
# Goal: reproduce Roth's two results (pre-tests have low power; conditioning on
# passing makes bias worse) analytically from one event-study fit, then compute
# HonestDiD breakdown values on the same design.

suppressMessages({
  library(tidyverse)
  library(fixest)
  library(HonestDiD)
  library(MASS)          # mvrnorm for the analytic draws
})
select <- dplyr::select  # MASS masks dplyr::select

# Shared DGP for modules 1-3 (identical in the slides)
make_panel <- function(scenario = c("constant", "dynamic", "heterogeneous"),
                       n_cities = 30, n_t = 60, seed = 42) {
  scenario <- match.arg(scenario)
  set.seed(seed)
  cohorts <- tibble(
    city = 1:n_cities,
    g = c(rep(15, 8), rep(25, 8), rep(35, 8), rep(Inf, 6))
  )
  expand_grid(city = 1:n_cities, t = 1:n_t) |>
    left_join(cohorts, by = "city") |>
    mutate(
      treated = t >= g,
      eff = case_when(
        !treated ~ 0,
        scenario == "constant"      ~ 1.0,
        scenario == "dynamic"       ~ 0.4 * (1 + 0.10 * (t - g)),
        scenario == "heterogeneous" ~ 0.4 * (1 + 0.10 * (t - g)) *
                                      (1 + 0.8 * (g == 15) - 0.8 * (g == 35))
      ),
      y = 5 + 0.05 * t + city * 0.1 + eff + rnorm(n(), 0, 0.3)
    )
}

# ===== Q1. A clean single-cohort event study and its variance matrix =====
# Task: from the constant-effects panel, build the event study for cohort
# g = 25 versus the six never-treated cities over the window -5..+5, clustered
# by city. Extract the coefficient vector and vcov in HonestDiD's expected
# order (earliest pre .. -2, then 0 .. K), and CONFIRM constructOriginalCS
# reproduces the feols CI for the first post period before trusting anything.

event_fit <- function(scenario = "constant", confound = 0) {
  d <- make_panel(scenario) |>
    filter(g == 25 | !is.finite(g)) |>
    mutate(rel_time = if_else(is.finite(g), as.integer(t - 25), -1L)) |>
    filter(t >= 20, t <= 30)                       # rel_time in [-5, 5]
  if (confound != 0)                               # differential linear trend
    d <- d |> mutate(y = y + if_else(g == 25, confound * (t - 25), 0))
  feols(y ~ i(rel_time, ref = -1) | city + t, data = d, cluster = ~city)
}

fit <- event_fit("constant")
ev  <- as.integer(sub("rel_time::", "", names(coef(fit))))
ord <- order(ev)                                   # ascending event time
betahat  <- unname(coef(fit)[ord])
Sigma    <- vcov(fit)[ord, ord]
evt      <- ev[ord]
numPre   <- sum(evt < 0)                           # 4
numPost  <- sum(evt >= 0)                          # 6
post1    <- which(evt == 0)                        # first post-period index
pre_idx  <- which(evt < 0)

# validate: HonestDiD's "original" CS uses NORMAL critical values, so it must
# match beta +/- 1.96*SE (feols confint uses a small-sample t, hence is a touch
# wider - the only difference is the degrees-of-freedom correction).
se_post1 <- sqrt(Sigma[post1, post1])
ci_norm  <- betahat[post1] + c(-1, 1) * qnorm(0.975) * se_post1
ocs <- constructOriginalCS(betahat, Sigma, numPre, numPost,
                           l_vec = basisVector(1, numPost))
cat(sprintf("Q1. first-post coef = %.3f (SE %.3f)\n", betahat[post1], se_post1))
cat(sprintf("    normal CI (1.96*SE) = [%.3f, %.3f]\n", ci_norm[1], ci_norm[2]))
cat(sprintf("    HonestDiD CI        = [%.3f, %.3f]  (must match)\n", ocs$lb, ocs$ub))
cat(sprintf("    feols t-based CI    = [%.3f, %.3f]  (wider: small-sample df)\n",
            confint(fit)["rel_time::0", ][[1]], confint(fit)["rel_time::0", ][[2]]))
stopifnot(abs(ocs$lb - ci_norm[1]) < 1e-6, abs(ocs$ub - ci_norm[2]) < 1e-6)

# ===== Q2. Roth's pretrend-power simulation (the analytic short-cut) =====
# Task: instead of re-simulating panels, draw betahat ~ N(beta_true, Sigma)
# using the Sigma you just estimated (Roth's approach). Impose a purely LINEAR
# parallel-trends violation of slope delta and NO true effect (tau = 0), so
# beta_true[s] = delta * (s + 1) relative to the base period -1. For a grid of
# delta, compute the power of the joint pre-trend Wald test. Show that a
# violation big enough to bias the estimate is detected well under 80% of the
# time.

crit  <- qchisq(0.95, df = numPre)                 # 5% joint-test threshold
Spre_inv <- solve(Sigma[pre_idx, pre_idx])

roth_sim <- function(delta, R = 4000, tau = 0) {
  set.seed(100 + round(1000 * delta))
  beta_true <- tau * (evt >= 0) + delta * (evt + 1)      # effect + linear trend
  draws <- MASS::mvrnorm(R, mu = beta_true, Sigma = Sigma)
  wald  <- apply(draws[, pre_idx, drop = FALSE], 1,
                 function(b) as.numeric(t(b) %*% Spre_inv %*% b))
  reject <- wald > crit
  bhat_post1 <- draws[, post1]
  tibble(
    delta       = delta,
    power       = mean(reject),                    # P(pre-test rejects)
    bias_uncond = mean(bhat_post1 - tau),          # naive estimate bias
    bias_cond   = mean(bhat_post1[!reject] - tau)  # bias | passed the pre-test
  )
}

grid <- c(0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12)
power_tbl <- map_dfr(grid, roth_sim)
cat("\nQ2. Pretrend-test power and naive bias by violation slope:\n")
print(power_tbl |> mutate(across(where(is.numeric), ~ round(.x, 3))))

# size is ~5% at delta = 0; power stays below 80% even at delta = 0.12
stopifnot(abs(power_tbl$power[power_tbl$delta == 0] - 0.05) < 0.03)
stopifnot(power_tbl$power[power_tbl$delta == 0.10] < 0.80)
stopifnot(all(power_tbl$power[power_tbl$delta > 0 & power_tbl$delta <= 0.12] < 0.80))
cat(sprintf("    A slope-0.10 trend biases the estimate yet is caught only %.0f%% of the time.\n",
            100 * power_tbl$power[power_tbl$delta == 0.10]))

# ===== Q3. Pre-testing does not clean up the survivors (it is worse) =====
# Task: compare the naive bias to the bias CONDITIONAL on passing the pre-test.
# Because leads and lags share the base period, their errors are positively
# correlated; passing selects draws whose noise pushed the leads (and thus the
# lags) the same way, so conditioning AMPLIFIES the bias.

corr_pre_post <- mean(cov2cor(Sigma)[pre_idx, post1])
cat(sprintf("\nQ3. mean lead-vs-first-lag correlation = %.2f (positive)\n",
            corr_pre_post))
worse <- power_tbl |> filter(delta > 0) |>
  mutate(ratio = bias_cond / bias_uncond)
cat("    unconditional vs conditional-on-passing bias:\n")
print(worse |> select(delta, bias_uncond, bias_cond, ratio) |>
        mutate(across(where(is.numeric), ~ round(.x, 3))))
stopifnot(all(worse$bias_cond >= worse$bias_uncond - 1e-6))   # conditioning never helps
cat(sprintf("    At delta = 0.10 pre-testing turns a %.2f bias into a %.2f bias.\n",
            power_tbl$bias_uncond[power_tbl$delta == 0.10],
            power_tbl$bias_cond[power_tbl$delta == 0.10]))

# ===== Q4. The breakdown value on the clean design =====
# Task: on the real (unconfounded) fit, report the relative-magnitudes
# breakdown value for the first post-period: the largest Mbar at which the
# robust CI still excludes zero.

bd_value <- function(sens) {
  col   <- if ("M" %in% names(sens)) "M" else "Mbar"
  incl0 <- sens$lb <= 0 & sens$ub >= 0
  if (all(!incl0)) return(Inf)                     # never breaks in the grid
  min(sens[[col]][incl0])                          # first value that covers 0
}

rm_clean <- createSensitivityResults_relativeMagnitudes(
  betahat, Sigma, numPre, numPost,
  Mbarvec = seq(0.25, 2, 0.25), l_vec = basisVector(1, numPost), gridPoints = 80)
Mbar_star <- bd_value(rm_clean)   # first Mbar at which the robust CI covers zero
cat(sprintf("\nQ4. Delta^RM breakdown Mbar* = %.2f: the robust CI covers zero here\n",
            Mbar_star))
cat(sprintf("    and excludes it for every smaller Mbar. Reading: the post-period\n"))
cat(sprintf("    violation would have to be %.1fx the worst pre-period one before\n",
            Mbar_star))
cat("    the effect becomes insignificant.\n")
stopifnot(Mbar_star >= 1.25)

# ===== Q5. The restriction must match the threat =====
# Task: inject a differential linear trend of slope 0.10 into the treated
# cohort. Recompute both breakdown values. Delta^RM reacts (the pre-period
# violation is now large); Delta^SD does NOT, because a linear trend has zero
# second difference and Delta^SD extrapolates it away. A smoothness restriction
# is blind to a linear confound by construction.

fit_c <- event_fit("constant", confound = -0.10)
evc   <- as.integer(sub("rel_time::", "", names(coef(fit_c))))
ordc  <- order(evc)
betac <- unname(coef(fit_c)[ordc]); Sigmac <- vcov(fit_c)[ordc, ordc]

sd_of <- function(b, S) createSensitivityResults(b, S, numPre, numPost,
  Mvec = seq(0.05, 0.6, 0.05), l_vec = basisVector(1, numPost), method = "C-LF")
rm_of <- function(b, S) createSensitivityResults_relativeMagnitudes(b, S,
  numPre, numPost, Mbarvec = seq(0.25, 2, 0.25),
  l_vec = basisVector(1, numPost), gridPoints = 80)

M_clean  <- bd_value(sd_of(betahat, Sigma))
M_conf   <- bd_value(sd_of(betac, Sigmac))
RM_clean <- bd_value(rm_of(betahat, Sigma))
RM_conf  <- bd_value(rm_of(betac, Sigmac))

cat("\nQ5. Breakdown values, clean vs a slope-0.10 linear confound:\n")
print(tibble(
  restriction = c("Delta^SD (smoothness)", "Delta^RM (relative magnitudes)"),
  clean       = c(M_clean, RM_clean),
  confounded  = c(M_conf, RM_conf)
))
stopifnot(abs(M_clean - M_conf) < 1e-9)     # SD blind to the linear confound
stopifnot(RM_conf < RM_clean)               # RM flags it
cat("    Delta^SD is unchanged (it de-trends linear violations); Delta^RM drops.\n")

cat("\nAll checks passed: pre-tests are underpowered, conditioning is worse,\n")
cat("and the breakdown value depends on which violation you are honest about.\n")
