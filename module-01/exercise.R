# ============================================================================
# Module 1 Exercise: TWFE Diagnosed - Hand-Coding the Bacon Decomposition
# ============================================================================
# Each question states a task; the solution follows. Try before reading on.
# Goal: replicate bacondecomp::bacon() from scratch, to machine precision.

suppressMessages({
  library(tidyverse)
  library(fixest)
  library(bacondecomp)
})

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
      cohort = case_when(g == 15 ~ "Early (g=15)", g == 25 ~ "Mid (g=25)",
                         g == 35 ~ "Late (g=35)", TRUE ~ "Never"),
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

# ===== Q1. TWFE vs the truth =====
# Task: simulate the dynamic scenario; compute the true ATT (mean effect over
# treated cells) and the TWFE estimate. How large is the attenuation?

panel <- make_panel("dynamic")

true_att <- mean(panel$eff[panel$treated])
twfe <- coef(feols(y ~ treated | city + t, data = panel))[["treatedTRUE"]]

cat(sprintf("Q1. true ATT = %.4f | TWFE = %.4f | bias = %.1f%%\n",
            true_att, twfe, 100 * (twfe / true_att - 1)))
stopifnot(twfe < 0.75 * true_att)  # the attenuation is not subtle

# ===== Q2. Enumerate all nine 2x2 estimates by hand =====
# Task: for the three cohorts (g = 15, 25, 35) and the never-treated group,
# compute every pairwise 2x2 DiD estimate from subsample means alone:
#   - each cohort vs never-treated (full panel, post = t >= g)
#   - earlier vs later (window t < g_later, control not yet treated)
#   - later vs earlier (window t >= g_earlier, control ALREADY treated)

T_ <- max(panel$t)

mean_y <- function(cities, t_lo, t_hi) {
  panel |> filter(city %in% cities, t >= t_lo, t <= t_hi) |> pull(y) |> mean()
}
cities_of <- function(gv) {
  panel |> distinct(city, g) |> filter(g == gv) |> pull(city)
}
cities_never <- panel |> distinct(city, g) |> filter(!is.finite(g)) |> pull(city)

groups <- panel |> distinct(city, g) |> count(g, name = "n_units") |>
  mutate(Dbar = ifelse(is.finite(g), (T_ - g + 1) / T_, 0))
gs <- groups |> filter(is.finite(g))
U  <- groups |> filter(!is.finite(g))

res <- list()

# -- cohort k vs never-treated --
for (i in seq_len(nrow(gs))) {
  gk <- gs$g[i]; nk <- gs$n_units[i]; Dk <- gs$Dbar[i]; nU <- U$n_units
  ck <- cities_of(gk)
  b <- (mean_y(ck, gk, T_) - mean_y(ck, 1, gk - 1)) -
       (mean_y(cities_never, gk, T_) - mean_y(cities_never, 1, gk - 1))
  nkU <- nk / (nk + nU)
  s <- (nk + nU)^2 * nkU * (1 - nkU) * Dk * (1 - Dk)   # Bacon weight (raw)
  res[[length(res) + 1]] <- tibble(
    g_treat = gk, g_ctrl = Inf, type = "Treated vs Untreated",
    estimate = b, s_raw = s)
}

# -- timing pairs: earlier k vs later l --
pairs <- expand_grid(i = seq_len(nrow(gs)), j = seq_len(nrow(gs))) |>
  filter(gs$g[i] < gs$g[j])
for (r in seq_len(nrow(pairs))) {
  gk <- gs$g[pairs$i[r]]; gl <- gs$g[pairs$j[r]]
  nk <- gs$n_units[pairs$i[r]]; nl <- gs$n_units[pairs$j[r]]
  Dk <- gs$Dbar[pairs$i[r]]; Dl <- gs$Dbar[pairs$j[r]]
  ck <- cities_of(gk); cl <- cities_of(gl)
  nkl <- nk / (nk + nl)

  # earlier vs later: window t < g_l, l is a clean not-yet-treated control
  b_k <- (mean_y(ck, gk, gl - 1) - mean_y(ck, 1, gk - 1)) -
         (mean_y(cl, gk, gl - 1) - mean_y(cl, 1, gk - 1))
  s_k <- ((nk + nl) * (1 - Dl))^2 * nkl * (1 - nkl) *
         ((Dk - Dl) / (1 - Dl)) * ((1 - Dk) / (1 - Dl))
  res[[length(res) + 1]] <- tibble(
    g_treat = gk, g_ctrl = gl, type = "Earlier vs Later Treated",
    estimate = b_k, s_raw = s_k)

  # later vs earlier: window t >= g_k, k is treated the WHOLE window (forbidden)
  b_l <- (mean_y(cl, gl, T_) - mean_y(cl, gk, gl - 1)) -
         (mean_y(ck, gl, T_) - mean_y(ck, gk, gl - 1))
  s_l <- ((nk + nl) * Dk)^2 * nkl * (1 - nkl) *
         (Dl / Dk) * ((Dk - Dl) / Dk)
  res[[length(res) + 1]] <- tibble(
    g_treat = gl, g_ctrl = gk, type = "Later vs Earlier Treated",
    estimate = b_l, s_raw = s_l)
}

hand <- bind_rows(res) |> mutate(weight = s_raw / sum(s_raw))
cat("\nQ2. Nine hand-coded 2x2s:\n")
print(hand |> select(-s_raw) |> mutate(across(where(is.numeric), ~ round(.x, 4))),
      n = 9)
stopifnot(abs(sum(hand$weight) - 1) < 1e-12)

# ===== Q3. The decomposition identity =====
# Task: verify that the weighted sum of your nine 2x2s equals the TWFE
# coefficient exactly.

hand_sum <- sum(hand$estimate * hand$weight)
cat(sprintf("\nQ3. weighted sum = %.10f | TWFE = %.10f | diff = %.2e\n",
            hand_sum, twfe, abs(hand_sum - twfe)))
stopifnot(abs(hand_sum - twfe) < 1e-8)

# ===== Q4. Replicate bacondecomp::bacon() =====
# Task: run the package on the same panel and match your estimates and
# weights, cell by cell.

df <- panel |>
  transmute(city, t, y, treat = as.numeric(treated)) |>  # avoid `treated` name clash
  as.data.frame()
bd <- bacon(y ~ treat, data = df, id_var = "city", time_var = "t")

cmp <- bd |>
  select(g_treat = treated, g_ctrl = untreated, w_pkg = weight, est_pkg = estimate) |>
  left_join(hand |> mutate(g_ctrl = ifelse(is.infinite(g_ctrl), 99999, g_ctrl)) |>
              select(g_treat, g_ctrl, w_hand = weight, est_hand = estimate),
            by = c("g_treat", "g_ctrl"))

cat(sprintf("\nQ4. max |weight diff| = %.2e | max |estimate diff| = %.2e\n",
            max(abs(cmp$w_pkg - cmp$w_hand)),
            max(abs(cmp$est_pkg - cmp$est_hand))))
stopifnot(max(abs(cmp$w_pkg - cmp$w_hand)) < 1e-12,
          max(abs(cmp$est_pkg - cmp$est_hand)) < 1e-10)

# ===== Q5. Where does the bias live? =====
# Task: aggregate the decomposition by comparison type. Which species drags
# the TWFE coefficient down, and how much weight does it carry?

by_type <- hand |>
  group_by(type) |>
  summarise(avg_estimate = weighted.mean(estimate, weight),
            weight = sum(weight), .groups = "drop")
cat("\nQ5. Decomposition by comparison type:\n")
print(by_type |> mutate(across(where(is.numeric), ~ round(.x, 4))))

forbidden_w <- by_type |> filter(type == "Later vs Earlier Treated") |> pull(weight)
cat(sprintf("\n    Forbidden comparisons carry %.0f%% of the weight.\n",
            100 * forbidden_w))

# ===== Q6. Scenario sweep: when is TWFE fine? =====
# Task: run TWFE on all three scenarios. Confirm that the constant scenario
# is unbiased and that heterogeneity makes things worse than dynamics alone.

sweep <- map_dfr(c("constant", "dynamic", "heterogeneous"), function(s) {
  p <- make_panel(s)
  tibble(scenario = s,
         true_att = mean(p$eff[p$treated]),
         twfe = coef(feols(y ~ treated | city + t, data = p))[[1]])
}) |> mutate(bias_pct = 100 * (twfe / true_att - 1))
cat("\nQ6. Scenario sweep:\n")
print(sweep |> mutate(across(where(is.numeric), ~ round(.x, 3))))
stopifnot(abs(sweep$bias_pct[sweep$scenario == "constant"]) < 5)

# ===== Q7. Sign reversal =====
# Task: first try steep dynamics (0.2 + 0.35 per week) on the standard panel
# and observe that TWFE stays positive: the 6 never-treated cities anchor
# half the weight on clean comparisons. Then drop the never-treated group
# (everyone eventually adopts) and show that even MILD dynamics
# (0.2 + 0.15 per week) flip the TWFE sign while every effect is positive.

# (a) steep dynamics, never-treated group present: attenuated but positive
extreme <- make_panel("dynamic") |>
  mutate(eff = if_else(treated, 0.2 + 0.35 * (t - g), 0),
         y = 5 + 0.05 * t + city * 0.1 + eff + rnorm(n(), 0, 0.3))
twfe_ext <- coef(feols(y ~ treated | city + t, data = extreme))[[1]]
cat(sprintf("\nQ7a. with never-treated: true ATT = %.2f | TWFE = %.3f (still > 0)\n",
            mean(extreme$eff[extreme$treated]), twfe_ext))
stopifnot(twfe_ext > 0)

# (b) no never-treated units, mild dynamics: the flip
set.seed(42)
flip <- expand_grid(city = 1:30, t = 1:60) |>
  left_join(tibble(city = 1:30, g = rep(c(15, 25, 35), each = 10)),
            by = "city") |>
  mutate(treated = t >= g,
         eff = if_else(treated, 0.2 + 0.15 * (t - g), 0),
         y = 5 + 0.05 * t + city * 0.1 + eff + rnorm(n(), 0, 0.3))
twfe_flip <- coef(feols(y ~ treated | city + t, data = flip))[[1]]
cat(sprintf("Q7b. no never-treated:  all effects > 0: %s | true ATT = %.2f | TWFE = %.3f\n",
            all(flip$eff[flip$treated] > 0),
            mean(flip$eff[flip$treated]), twfe_flip))
stopifnot(all(flip$eff[flip$treated] > 0), twfe_flip < 0)

cat("\nAll checks passed: you have replicated bacon() from scratch.\n")
