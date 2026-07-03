# ============================================================================
# Module 2 Exercise: Heterogeneity-Robust DiD
# ============================================================================
# Each question states a task; the solution follows. Try before reading on.
# Centerpiece: implement Callaway-Sant'Anna from scratch and match did::att_gt
# to machine precision, then line up all four estimators against their true
# targets.

suppressMessages({
  library(tidyverse)
  library(fixest)
  library(did)
  library(didimputation)
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

panel <- make_panel("heterogeneous")

# ===== Q1. The true targets, computed from `eff` =====
# Task: three aggregations answer three questions. Compute each TRUE target
# directly from the effect column, so you have something to validate against.
#   simple  = mean effect over every treated cell
#   group   = average over cohorts of each cohort's mean post-effect
#   dynamic = average over event times of the mean effect at that event time

true_simple <- mean(panel$eff[panel$treated])
true_group <- panel |> filter(treated, is.finite(g)) |>
  group_by(g) |> summarise(a = mean(eff), .groups = "drop") |> pull(a) |> mean()
true_dynamic <- panel |> filter(treated, is.finite(g)) |>
  mutate(e = t - g) |> group_by(e) |> summarise(a = mean(eff), .groups = "drop") |>
  pull(a) |> mean()

cat(sprintf("Q1. true targets: simple = %.4f | group = %.4f | dynamic = %.4f\n",
            true_simple, true_group, true_dynamic))
stopifnot(true_dynamic > true_simple, true_simple > true_group)  # they differ

# ===== Q2. Callaway-Sant'Anna FROM SCRATCH =====
# Task: compute ATT(g,t) for every (g, t >= g) as a clean 2x2 DiD, using
# not-yet-treated controls and a UNIVERSAL base period (g - 1):
#   ATT(g,t) = [Ybar_{g,t} - Ybar_{g,g-1}] - [Ybar_{C,t} - Ybar_{C,g-1}]
# where C = cohorts not yet treated by t (those with G > t, incl. never).
# Pool the control cohorts weighting by their unit counts.

gm <- panel |> group_by(g, t) |> summarise(ybar = mean(y), .groups = "drop")
ybar <- function(gv, tt) gm$ybar[gm$g == gv & gm$t == tt]

cohorts_g   <- sort(unique(panel$g[is.finite(panel$g)]))   # 15, 25, 35
all_g       <- sort(unique(panel$g))                       # + Inf (never)
periods     <- sort(unique(panel$t))
unit_counts <- panel |> distinct(city, g) |> count(g)

hand <- map_dfr(cohorts_g, function(gv) {
  map_dfr(periods[periods >= gv], function(tt) {
    ctrl_g <- all_g[all_g > tt]                            # not-yet-treated by t
    w <- unit_counts$n[match(ctrl_g, unit_counts$g)]
    cbar_t    <- weighted.mean(vapply(ctrl_g, ybar, 0, tt = tt),     w)
    cbar_base <- weighted.mean(vapply(ctrl_g, ybar, 0, tt = gv - 1), w)
    tibble(g = gv, t = tt, e = tt - gv,
           att_hand = (ybar(gv, tt) - ybar(gv, gv - 1)) - (cbar_t - cbar_base))
  })
})

# Validate cell by cell against the package
dat <- panel |> mutate(g0 = ifelse(is.finite(g), g, 0)) |>  # did wants 0 = never
  select(city, t, y, g0) |> as.data.frame()
cs_nyt <- att_gt(yname = "y", tname = "t", idname = "city", gname = "g0",
                 data = dat, control_group = "notyettreated",
                 base_period = "universal", est_method = "dr")
pkg <- tibble(g = cs_nyt$group, t = cs_nyt$t, att_pkg = cs_nyt$att) |> filter(t >= g)

cmp <- hand |> inner_join(pkg, by = c("g", "t"))
cat(sprintf("\nQ2. hand-coded CS vs did::att_gt: %d cells, max |diff| = %.2e\n",
            nrow(cmp), max(abs(cmp$att_hand - cmp$att_pkg))))
stopifnot(max(abs(cmp$att_hand - cmp$att_pkg)) < 1e-10)  # DR = 2x2 means exactly

# ===== Q3. Aggregate your ATT(g,t), match aggte() =====
# Task: aggregate the hand-coded surface to group and dynamic effects and
# confirm they match did::aggte(). group = equal weight per (equal-size)
# cohort; dynamic = equal weight per event time e >= 0.

hand_group <- hand |> group_by(g) |> summarise(a = mean(att_hand), .groups = "drop") |>
  pull(a) |> mean()
hand_dynamic <- hand |> group_by(e) |> summarise(a = mean(att_hand), .groups = "drop") |>
  filter(e >= 0) |> pull(a) |> mean()

agg_group <- aggte(cs_nyt, type = "group")$overall.att
agg_dyn   <- aggte(cs_nyt, type = "dynamic")$overall.att
cat(sprintf("Q3. group:   hand = %.4f | aggte = %.4f | diff = %.2e\n",
            hand_group, agg_group, abs(hand_group - agg_group)))
cat(sprintf("    dynamic: hand = %.4f | aggte = %.4f | diff = %.2e\n",
            hand_dynamic, agg_dyn, abs(hand_dynamic - agg_dyn)))
stopifnot(abs(hand_group - agg_group) < 1e-8, abs(hand_dynamic - agg_dyn) < 1e-8)

# ===== Q4. Sun-Abraham equals CS (never-treated, same base) =====
# Task: fit the interaction-weighted estimator with fixest::sunab (never-treated
# cohort recoded to a large finite value), aggregate to an overall ATT, and show
# it equals the CS SIMPLE aggregation with never-treated controls.

panel_sa <- panel |> mutate(g_sa = ifelse(is.finite(g), g, 10000))
fit_sa <- feols(y ~ sunab(g_sa, t) | city + t, data = panel_sa)
sa_att <- summary(fit_sa, agg = "att")$coeftable["ATT", "Estimate"]

cs_nev <- att_gt(yname = "y", tname = "t", idname = "city", gname = "g0",
                 data = dat, control_group = "nevertreated",
                 base_period = "universal", est_method = "dr")
cs_simple <- aggte(cs_nev, type = "simple")$overall.att
cat(sprintf("\nQ4. Sun-Abraham att = %.5f | CS-never simple = %.5f | diff = %.2e\n",
            sa_att, cs_simple, abs(sa_att - cs_simple)))
stopifnot(abs(sa_att - cs_simple) < 1e-6)  # algebraically identical

# ===== Q5. BJS imputation: efficiency under parallel trends =====
# Task: run the imputation estimator (untreated cells pin the FE model, impute
# Y(0), average the residuals). Show it nails the true SIMPLE ATT and does so
# with a far smaller SE than CS, because it uses every untreated cell.

bjs <- did_imputation(dat, yname = "y", gname = "g0", tname = "t", idname = "city")
cs_simple_se <- aggte(cs_nev, type = "simple")$overall.se
cat(sprintf("Q5. BJS = %.4f (se %.3f) vs true simple %.4f | CS se = %.3f\n",
            bjs$estimate, bjs$std.error, true_simple, cs_simple_se))
cat(sprintf("    BJS is %.1fx more precise than CS on the same panel.\n",
            cs_simple_se / bjs$std.error))
stopifnot(abs(bjs$estimate - true_simple) < 0.05, bjs$std.error < cs_simple_se)

# ===== Q6. dCDH switcher estimator (DID_M), hand-coded =====
# Task: DIDmultiplegt does not load here (X11). Hand-code the instantaneous
# DID_M: at each period where units switch on, take switchers' first difference
# minus stayers'-at-zero first difference, weight by switcher count. Validate on
# the DYNAMIC scenario, whose switch-time (e = 0) effect is 0.4 by construction.

did_m_instant <- function(df) {
  df <- df |> arrange(city, t) |> group_by(city) |>
    mutate(D = as.integer(treated), D_lag = lag(D), dy = y - lag(y)) |> ungroup()
  switchers <- df |> filter(!is.na(D_lag), D == 1, D_lag == 0)  # 0 -> 1 at t
  stayers0  <- df |> filter(!is.na(D_lag), D == 0, D_lag == 0)  # 0 -> 0 at t
  map_dfr(sort(unique(switchers$t)), function(tt) {
    tibble(t = tt, n_sw = sum(switchers$t == tt),
           didm = mean(switchers$dy[switchers$t == tt]) -
                  mean(stayers0$dy[stayers0$t == tt]))
  })
}

dm <- did_m_instant(make_panel("dynamic"))
dm_overall <- weighted.mean(dm$didm, dm$n_sw)
cat("\nQ6. DID_M by switch time (dynamic scenario):\n")
print(dm |> mutate(didm = round(didm, 3)))
cat(sprintf("    DID_M overall = %.4f | true instantaneous effect = 0.4000\n",
            dm_overall))
stopifnot(abs(dm_overall - 0.4) < 0.1)  # instantaneous switch effect recovered

# ===== Q7. Four estimators, one table =====
# Task: line up TWFE and all four robust estimators on the heterogeneous panel,
# each next to the target it is actually estimating. TWFE should be far below
# its target; the robust estimators cluster near theirs. simple/group/dynamic
# differ because they are different estimands.

twfe <- coef(feols(y ~ treated | city + t, data = panel))[["treatedTRUE"]]
final <- tribble(
  ~estimator,                 ~estimate,                              ~true_target,
  "TWFE",                      twfe,                                   true_simple,
  "CS (simple)",               cs_simple,                              true_simple,
  "Sun-Abraham",               sa_att,                                 true_simple,
  "BJS imputation",            bjs$estimate,                           true_simple,
  "CS (group)",                aggte(cs_nev, "group")$overall.att,     true_group,
  "CS (dynamic)",              aggte(cs_nev, "dynamic")$overall.att,   true_dynamic
) |> mutate(gap = estimate - true_target)
cat("\nQ7. Four estimators vs their true targets (heterogeneous):\n")
print(final |> mutate(across(where(is.numeric), ~ round(.x, 3))))
stopifnot(twfe < 0.5 * true_simple)                 # TWFE badly attenuated
stopifnot(abs(bjs$estimate - true_simple) < 0.05)   # BJS on target

cat("\nAll checks passed: CS reproduced from scratch, four estimators validated.\n")
