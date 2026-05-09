# Causal Inference Beyond A/B Tests — Learning Plan

A deep-dive companion to `experimentation-refresher/module-08`. Status:
**planned, not built.** This document is the spec; modules will follow.

## Audience

You've completed the experimentation-refresher Module 8 tour and want the
formal treatment of each method — the kind of detail an applied
econometrician or research scientist needs to make tooling decisions, not
just recognize names.

## Module 1 — TWFE Diagnosed
- The Goodman-Bacon (2021) decomposition: full derivation, weights formula.
- Bacon-decomposition diagnostic in practice (`bacondecomp` package).
- The "forbidden comparison" (already-treated as control) — visualizing it.
- When TWFE is *unbiased*: constant effects, no staggering.
- **Application:** ride-sharing zone-notification staggered rollout (extends
  the M8 example to 30 cities and 60 periods, with three different effect
  heterogeneity scenarios).
- **Exercise:** hand-code Bacon weights, replicate the package output.

## Module 2 — Heterogeneity-Robust DiD
- Callaway & Sant'Anna (2021) — group-time ATT, full estimator, not-yet-treated
  vs never-treated controls, doubly-robust adjustment, influence-function SEs.
- Sun & Abraham (2021) — saturated event-study, when it equals CS.
- Borusyak, Jaravel & Spiess (2024) — imputation, efficiency under PT.
- de Chaisemartin & D'Haultfœuille (2020) — multi-period switchers.
- Aggregation choices: simple average vs cohort-weighted vs event-time
  aggregation. What each answers.
- **Application:** same staggered rollout, all four estimators side-by-side.
- **Exercise:** implement CS from scratch using `Mostly_Harmless`-style data,
  validate against `did::att_gt()`.

## Module 3 — Honest DiD
- Rambachan & Roth (2023) — the formal sensitivity framework.
- Restrictions: smoothness ($M$), relative magnitudes, sign restrictions.
- The breakdown value as a sufficient statistic for robustness.
- Practical workflow with the `HonestDiD` package.
- **Application:** event-study estimates from M2, robustified.
- **Exercise:** replicate the Roth (2022) pretrend-power simulation.

## Module 4 — Synthetic Control: Full Treatment
- Abadie/Diamond/Hainmueller (2010, 2015) — original estimator, V-matrix.
- Donor pool selection: how to choose, what to avoid.
- Inference — placebo-in-space, placebo-in-time, conformal SCM.
- Augmented SC (Ben-Michael et al. 2021) — bias correction via outcome model.
- Generalized SC / `gsynth` (Xu 2017) — interactive fixed effects model.
- **Application:** California Prop 99 (canonical), then ride-sharing
  one-city policy change.
- **Exercise:** build SC from scratch with `quadprog`; compare to `Synth::synth()`.

## Module 5 — Synthetic DiD
- Arkhangelsky et al. (2021) — formal estimator, dual weights, regularization.
- The bridge: when SDID = DiD, when SDID = SC.
- Empirical comparison from the AER paper.
- `synthdid` package: estimation, jackknife SE, the comparison plot.
- **Application:** California Prop 99 again, now with SDID overlay.
- **Exercise:** implement the dual-weight optimization; compare to package.

## Module 6 — Causal Forest: Honest Splitting and Asymptotics
- Wager & Athey (2018) — honest splitting, asymptotic normality.
- Athey, Tibshirani & Wager (2019) — generalized random forests.
- The local moment-condition view: forests as adaptive nearest neighbors.
- Tuning `grf` — `min.node.size`, `honesty.fraction`, `tune.parameters`.
- HTE diagnostics — calibration, BLP test, omnibus test.
- **Application:** zone-notification HTE by city density, tenure, time-of-day.
- **Exercise:** implement honest splitting on a 2D toy DGP; verify pointwise CIs.

## Module 7 — Policy Learning
- Athey & Wager (2021) — regret bounds for learned policies.
- `grf::policy_tree()` — depth-$L$ decision trees over $X$ for treatment rules.
- Welfare maximization vs welfare estimation — different objectives.
- Cost-aware policies: $\pi(x) = \mathbb{1}\{\hat\tau(x) > c(x)\}$.
- Off-policy evaluation: doubly-robust value estimates.
- **Application:** zone-notification — when to push, to whom, given per-push cost.
- **Exercise:** policy-tree on causal-forest output; OPE with cross-fitting.

## Module 8 — Matrix Completion and the Modern Panel
- Athey et al. (2021) — matrix completion estimator for panel causal effects.
- The connection to SC, SDID, and interactive FE.
- `gsynth`, `MCPanel`, and the convergence of the panel literature.
- When to reach for matrix completion vs DiD vs SC.
- **Application:** sparse panel — many cities, some treated at irregular times.
- **Exercise:** simulate a low-rank panel; compare MC, SDID, and CS RMSEs.

## Stack

- **R** — primary.
- **fixest** — TWFE and Sun-Abraham.
- **did** — Callaway-Sant'Anna.
- **didimputation** — Borusyak et al.
- **DIDmultiplegt** — de Chaisemartin & D'Haultfœuille.
- **HonestDiD** — Rambachan & Roth.
- **Synth** — original SC.
- **gsynth** — generalized SC, interactive FE.
- **synthdid** — Synthetic DiD.
- **augsynth** — augmented SC.
- **grf** — generalized random forests, causal forest, policy tree.
- **MCPanel** — matrix completion.

## Key References

- Goodman-Bacon (2021), J. Econometrics.
- Callaway & Sant'Anna (2021), J. Econometrics.
- Sun & Abraham (2021), J. Econometrics.
- Borusyak, Jaravel & Spiess (2024).
- de Chaisemartin & D'Haultfœuille (2020), AER.
- Rambachan & Roth (2023), ReStud.
- Roth (2022), AER:Insights.
- Abadie, Diamond & Hainmueller (2010), JASA; (2015), AJPS.
- Ben-Michael, Feller & Rothstein (2021), JASA — augmented SC.
- Xu (2017), Political Analysis — generalized SC / interactive FE.
- Arkhangelsky et al. (2021), AER — Synthetic DiD.
- Wager & Athey (2018), JASA.
- Athey, Tibshirani & Wager (2019), Annals of Statistics.
- Athey & Wager (2021), Econometrica — policy learning.
- Athey et al. (2021), JASA — matrix completion.
