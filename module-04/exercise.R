# ============================================================================
# Module 4 Exercise: Synthetic Control from Scratch vs Synth, Placebo Inference
# ============================================================================
# Each question states a task; the solution follows. Try before reading on.
# Goal: hand-code the SC quadratic program on California Prop 99, match it to
# Synth::synth(), then run the full placebo-in-space and placebo-in-time tests.

suppressMessages({
  library(tidyverse)
  library(quadprog)
  library(Synth)
})

# ---- Data: California Prop 99 (39 states, 1970-2000, CA treated 1989) --------
data("california_prop99", package = "synthdid")
cp <- california_prop99 |> mutate(State = as.character(State))

# Wide outcome matrix: rows = years, columns = states
wide <- cp |>
  select(State, Year, PacksPerCapita) |>
  pivot_wider(names_from = State, values_from = PacksPerCapita) |>
  arrange(Year)
years    <- wide$Year
Ymat     <- as.matrix(wide[, -1])          # 31 x 39
pre_idx  <- which(years <= 1988)           # 19 pre-treatment years
post_idx <- which(years >  1988)           # 12 post years
donors   <- setdiff(colnames(Ymat), "California")

# ===== Q1. Build the synthetic control from scratch (V = I) =====
# Task: solve the SC quadratic program on the pre-period outcome lags. Minimize
#   || Y1_pre - Y0_pre W ||^2   s.t.  w_j >= 0,  sum_j w_j = 1
# using solve.QP. Report the donor weights, pre-RMSPE, and post-period ATT.
# Note: with 38 donors and 19 pre-years the design is rank-deficient, so add a
# tiny ridge to keep the QP's Dmat strictly positive definite (solve.QP needs it).

sc_weights <- function(y1, Y0, idx, lambda = 1e-4) {
  A <- Y0[idx, , drop = FALSE]; b <- y1[idx]; J <- ncol(Y0)
  Dmat <- 2 * (t(A) %*% A + lambda * diag(J))   # ridge -> PD, deterministic
  dvec <- 2 * as.numeric(t(A) %*% b)
  Amat <- cbind(rep(1, J), diag(J))             # 1st col: equality (sum = 1)
  bvec <- c(1, rep(0, J))                        # then w_j >= 0
  solve.QP(Dmat, dvec, Amat, bvec, meq = 1)$solution
}

y1 <- Ymat[, "California"]
Y0 <- Ymat[, donors, drop = FALSE]
w_qp <- sc_weights(y1, Y0, pre_idx); names(w_qp) <- donors
gap_qp   <- y1 - as.numeric(Y0 %*% w_qp)
pre_rmspe_qp  <- sqrt(mean(gap_qp[pre_idx]^2))
att_qp   <- mean(gap_qp[post_idx])

cat("Q1. From-scratch SC (V = I):\n")
print(round(sort(w_qp[w_qp > 0.01], decreasing = TRUE), 3))
cat(sprintf("    pre-RMSPE = %.3f | post-period ATT = %.2f packs\n",
            pre_rmspe_qp, att_qp))
stopifnot(abs(sum(w_qp) - 1) < 1e-6, all(w_qp > -1e-8), att_qp < -10)

# ===== Q2. Replicate with Synth::synth() and compare =====
# Task: fit Synth::synth() with one outcome-lag special.predictor per pre-year,
# so it optimizes V over the same information. Compare weights side by side,
# both pre-RMSPEs, and the correlation of the two gap series.

ids <- cp |> distinct(State) |> arrange(State) |> mutate(id = row_number())
d   <- cp |> left_join(ids, by = "State") |> as.data.frame()  # Synth needs data.frame
ca_id     <- ids$id[ids$State == "California"]
donor_ids <- setdiff(ids$id, ca_id)
pre_years <- 1970:1988

dp <- dataprep(
  foo = d, dependent = "PacksPerCapita",
  unit.variable = "id", time.variable = "Year", unit.names.variable = "State",
  special.predictors = lapply(pre_years, function(y) list("PacksPerCapita", y, "mean")),
  treatment.identifier = ca_id, controls.identifier = donor_ids,
  time.predictors.prior = pre_years, time.optimize.ssr = pre_years,
  time.plot = 1970:2000)
invisible(capture.output(so <- synth(dp, verbose = FALSE)))   # hide verbose optimizer

w_synth <- as.numeric(so$solution.w)
names(w_synth) <- ids$State[match(rownames(so$solution.w), ids$id)]
synth_path <- as.numeric(dp$Y0plot %*% so$solution.w)
gap_synth  <- as.numeric(dp$Y1plot) - synth_path
pre_rmspe_synth <- sqrt(mean(gap_synth[pre_idx]^2))

cmp <- tibble(donor = donors,
              w_scratch = round(w_qp, 3),
              w_synth   = round(w_synth[donors], 3)) |>
  filter(w_scratch > 0.01 | w_synth > 0.01) |>
  arrange(desc(w_scratch))
cat("\nQ2. Weight comparison (scratch vs Synth):\n"); print(cmp)
gap_cor <- cor(gap_qp, gap_synth)
cat(sprintf("    pre-RMSPE: scratch %.3f | Synth %.3f\n", pre_rmspe_qp, pre_rmspe_synth))
cat(sprintf("    gap-series correlation = %.4f\n", gap_cor))
stopifnot(gap_cor > 0.95)   # the two estimators agree on the counterfactual path

# ===== Q3. Placebo-in-space and California's exact p-value =====
# Task: assign treatment to each state in turn, refit SC against the rest,
# and compute the post/pre RMSPE ratio. Rank California; the exact p-value is
# rank / N. Confirm California is among the top few ratios.
# Watch the tibble trap: DON'T name a column `pre` and then index gap[-pre] in
# the same tibble() call -- the new column masks the index vector. Compute the
# scalars first, then assemble.

placebo_ratio <- function(nm) {
  yy <- Ymat[, nm]; DD <- Ymat[, setdiff(colnames(Ymat), nm), drop = FALSE]
  w  <- sc_weights(yy, DD, pre_idx)
  g  <- yy - as.numeric(DD %*% w)
  pre_r  <- sqrt(mean(g[pre_idx]^2))
  post_r <- sqrt(mean(g[post_idx]^2))
  tibble(state = nm, pre_rmspe = pre_r, post_rmspe = post_r, ratio = post_r / pre_r)
}
placebos <- map_dfr(colnames(Ymat), placebo_ratio)
ca_ratio <- placebos$ratio[placebos$state == "California"]
rank_ca  <- sum(placebos$ratio >= ca_ratio)
p_value  <- rank_ca / nrow(placebos)

cat("\nQ3. Placebo-in-space (RMSPE-ratio test):\n")
print(placebos |> arrange(desc(ratio)) |> head(5) |>
      mutate(across(where(is.numeric), ~ round(.x, 2))))
cat(sprintf("    California ratio = %.2f | rank %d of %d | exact p = %.3f\n",
            ca_ratio, rank_ca, nrow(placebos), p_value))
stopifnot(rank_ca <= 3, ca_ratio > 8)   # California is in the top few

# ===== Q4. Trimming bad pre-fitters =====
# Task: apply Abadie's screen -- drop placebo states whose pre-RMSPE exceeds
# 2x California's. Show it removes only poorly-fit states (which carry low
# ratios anyway), so California's rank is unchanged here: the two states above
# it (Missouri, Virginia) are themselves well fit.

ca_pre <- placebos$pre_rmspe[placebos$state == "California"]
kept   <- placebos |> filter(pre_rmspe <= 2 * ca_pre)
rank_trim <- sum(kept$ratio >= ca_ratio)
cat(sprintf("\nQ4. Trim pre-RMSPE > 2x CA: kept %d of %d states | CA rank %d | p = %.3f\n",
            nrow(kept), nrow(placebos), rank_trim, rank_trim / nrow(kept)))
cat("    States above California (well-fit placebos that drift post-period):\n")
print(placebos |> filter(ratio > ca_ratio) |> arrange(desc(ratio)) |>
      mutate(across(where(is.numeric), ~ round(.x, 2))))

# ===== Q5. Placebo-in-time (backdating) =====
# Task: pretend treatment happened in 1980. Refit SC on 1970-1979 only, then
# look at the 1980-1988 gap. A credible design shows a small pseudo-effect
# before the real 1989 date, far smaller than the real post-1989 effect.

fake_pre <- which(years <= 1979)
w_bd  <- sc_weights(y1, Y0, fake_pre)
gap_bd <- y1 - as.numeric(Y0 %*% w_bd)
pseudo_effect <- mean(gap_bd[years >= 1980 & years <= 1988])
real_effect   <- mean(gap_bd[years > 1988])
cat(sprintf("\nQ5. Backdate to 1980 (fit 1970-1979):\n"))
cat(sprintf("    pseudo-effect 1980-1988 = %+.2f | real effect 1989-2000 = %+.2f\n",
            pseudo_effect, real_effect))
stopifnot(abs(pseudo_effect) < 0.5 * abs(real_effect))  # real effect dominates

cat("\nAll checks passed: SC from scratch matches Synth, and the placebo",
    "\ntests localize California's effect to the post-1989 period.\n")
