# Module 2: Heterogeneity-Robust DiD

Module 1 diagnosed the disease: under staggered adoption TWFE is a
variance-weighted average of every 2x2, and the forbidden comparisons
(already-treated units as controls) subtract other cohorts' effect
dynamics. This module is the cure. Four estimators that use only **clean**
comparisons (treated vs never-treated or not-yet-treated), plus the
aggregation choices that decide *which* average you end up reporting.

The running application is the same staggered zone-notification rollout: 30
cities, 60 weeks, cohorts at weeks 15/25/35 (8 cities each) plus 6
never-treated cities, parallel trends true by construction. Every number
below is on the `heterogeneous` scenario unless stated otherwise.

## Callaway and Sant'Anna (2021): the workhorse

### The estimand: group-time ATT

TWFE collapses a two-dimensional object into one number. CS refuses to
collapse. The building block is the **group-time average treatment effect**,

$$\text{ATT}(g, t) = \mathbb{E}\!\left[Y_t(g) - Y_t(0) \mid G = g\right],$$

the effect at calendar time $t$ on the cohort first treated at $g$ (where
$G=g$ indexes the cohort and $Y_t(0)$ is the never-treated potential
outcome). Two indices, deliberately kept apart: a cohort dimension $g$ and a
time dimension $t$. Everything else in CS is a choice about how to identify
each $\text{ATT}(g,t)$ and how to average them back down.

### Identification: which controls, which base period

For $t \geq g$, under parallel trends and no anticipation, $\text{ATT}(g,t)$
is a clean 2x2 DiD between cohort $g$ and a comparison group $C$, anchored at
a pre-period base:

$$\text{ATT}(g, t) = \big[\bar Y_{g,t} - \bar Y_{g,g-1}\big]
- \big[\bar Y_{C,t} - \bar Y_{C,g-1}\big].$$

Two identification levers:

- **Comparison group $C$.** *Never-treated* uses only the 6 cities that never
  adopt: cleanest, but throws away information and is high-variance when the
  never-treated group is small. *Not-yet-treated* adds, for each $(g,t)$, any
  cohort not treated by $t$ (cohorts with $G > t$, plus the never-treated):
  more controls, lower variance, valid under the same parallel-trends
  assumption extended to the not-yet-treated.
- **Base period.** The *universal* base fixes the anchor at $g-1$ (cohort
  $g$'s last pre-treatment period) for every $t$, so each $\text{ATT}(g,t)$ is
  a direct test of "did the gap move from $g-1$ to $t$?" The *varying* base
  uses $t-1$, which re-anchors as $t$ moves and is the natural choice for
  building an event study from long differences. Universal base is what makes
  CS coincide with Sun-Abraham.

### The doubly-robust estimator

The formula above is the unconditional version. With covariates $X$ you want
robustness to functional-form error. CS uses the **doubly-robust** moment
(Sant'Anna and Zhao 2020): combine an outcome-regression model
$m(X) = \mathbb{E}[\Delta Y \mid X, C]$ for the control trend with a
propensity model $p(X) = \Pr(G = g \mid X)$ for group membership,

$$\widehat{\text{ATT}}^{dr}(g,t) = \mathbb{E}\!\left[ \left(\frac{\mathbb{1}\{G=g\}}{\mathbb{E}[\mathbb{1}\{G=g\}]} - \frac{p(X)(1-\mathbb{1}\{G=g\})\,/\,(1-p(X))}{\mathbb{E}[p(X)(1-\mathbb{1}\{G=g\})\,/\,(1-p(X))]}\right)\big(\Delta Y - m(X)\big)\right].$$

Consistent if *either* $m$ or $p$ is correctly specified, not both. With no
covariates both models are trivial and the estimator collapses back to the
2x2 of subsample means: `est_method = "dr"` with no `xformla` returns exactly
the hand-coded difference. That identity is the exercise's validation target
and it holds to machine precision (3e-15 here).

### Inference: influence functions and uniform bands

CS derives the **influence function** of each $\widehat{\text{ATT}}(g,t)$,
the observation-level score whose sample variance is the estimator's
variance. That buys two things. First, analytic standard errors that account
for the estimated nuisance models. Second, because the influence functions
for all $(g,t)$ are jointly available, a **multiplier bootstrap** delivers
*uniform* confidence bands: intervals that cover the entire ATT path
simultaneously, not one point at a time. The simultaneous critical value
exceeds the pointwise 1.96 (here 2.90 for the dynamic aggregation), so honest
event-study plots use the wider band. `did` clusters at the unit level by
default and turns on `bstrap` and `cband` automatically.

## Aggregation: simple, group, dynamic

$\text{ATT}(g,t)$ is a surface. A headline number is a weighted average over
it, and the weights are a modeling choice that changes the answer.

| Aggregation | Averages over | Answers | True target here |
|---|---|---|---|
| **simple** | every treated $(g,t)$ cell equally | "average effect across all treated place-weeks" | 1.407 |
| **group** | cohorts (weighted by cohort size), each cohort's post-effects averaged | "average effect per treated cohort" | 1.207 |
| **dynamic** | event time $e = t - g$, equal weight per $e$ | "how does the effect evolve with exposure?" | 1.792 |

These are genuinely different estimands and they differ numerically. On this
panel the true targets are 1.41, 1.21, and 1.79. The dynamic average is
largest because it equal-weights event times, and at long exposure only the
early high-effect cohort survives, so long-$e$ periods (dominated by the big
cohort) get full weight. The group average is smallest because it gives the
low-effect late cohort the same cohort weight as the high-effect early one.
Report the aggregation that matches the decision: simple for a rollout's
average lift, group to see which cohorts drove it, dynamic to judge whether
the effect is still growing.

A practical caution specific to this panel: CS with only 6 never-treated
cities is high-variance. The point estimates sit above the true targets by a
fixed amount (about +0.25 for simple, +0.20 for group, +0.30 for dynamic),
and that offset is *identical across all three effect scenarios*, which is the
signature of a shared finite-sample draw in the small never-treated anchor,
not bias. Standard errors are correspondingly large (0.10 to 0.25). This is
exactly the setting where imputation wins.

## Sun and Abraham (2021): fixing the event study

### The contamination problem

The reflex diagnostic for DiD is the TWFE event study: regress the outcome on
leads and lags of treatment, $y_{it} = \alpha_i + \lambda_t + \sum_{e \neq -1}
\mu_e \, \mathbb{1}\{t - g_i = e\} + \varepsilon_{it}$, and read $\mu_e$ as the
effect $e$ periods after adoption. Sun and Abraham show each $\mu_e$ is *not*
the average effect at relative time $e$. It is a weighted sum of cohort-specific
effects at $e$ **and at other relative times**, with weights that can be
negative and need not sum to one within $e$. The contamination is worst under
cohort heterogeneity: a lead coefficient $\mu_{-3}$ can pick up post-treatment
effects of a different cohort, producing a **nonzero pre-trend even when
parallel trends holds exactly**.

You can see it with the noise switched off. On the deterministic panel (no
sampling error, parallel trends exact by construction) the naive TWFE event
study still shows pre-period coefficients growing to 0.145 and post
coefficients attenuated below the truth. Every wiggle is contamination, not
noise. The Sun-Abraham and CS event studies on the same deterministic panel
return exactly zero pre and the exact post path.

### The interaction-weighted estimator

The fix is to saturate: interact every cohort with every relative time,

$$y_{it} = \alpha_i + \lambda_t + \sum_{g}\sum_{e \neq -1}
\delta_{g,e}\, \mathbb{1}\{G_i = g\}\,\mathbb{1}\{t - g = e\} + \varepsilon_{it},$$

so $\delta_{g,e}$ is a clean cohort-specific effect. Then average the
$\delta_{g,e}$ across cohorts at each $e$ using **each cohort's share of
treated units** as weights: the interaction-weighted (IW) estimator. Because
the weights are cohort shares, they are non-negative and sum to one, and no
cohort's effect leaks into another relative time. `fixest::sunab()` builds the
saturated design and the aggregation in one call; the never-treated cohort is
coded as a large finite value so it acts as the excluded control.

### When SA equals CS

The IW estimator equals CS **when both use never-treated controls and the same
base period** ($g-1$). Under those conditions the two are algebraically the
same clean comparison, differing only in the software's variance estimator.
On this panel the equality is exact: SA's overall ATT and CS's simple
aggregation both return 1.6606 (difference 1e-11), and the two event-study
paths are identical at every $e$. The standard errors differ (SA 0.10, CS
0.25) because they estimate the control-group sampling variance differently.

## Borusyak, Jaravel and Spiess (2024): imputation

The imputation estimator is the most efficient of the four under parallel
trends, and the most transparent. Three steps:

1. Fit the two-way fixed-effects model $Y_{it} = \alpha_i + \lambda_t +
   \varepsilon_{it}$ on **untreated cells only** (all not-yet-treated and
   never-treated observations, including pre-periods of treated units).
2. Impute the untreated potential outcome $\hat Y_{it}(0) = \hat\alpha_i +
   \hat\lambda_t$ for every treated cell.
3. The estimated effect is $\hat Y_{it} - \hat Y_{it}(0)$; average with
   whatever weights define your target (`horizon = TRUE` gives the event
   study).

Because step 1 uses *every* untreated cell to pin down the fixed effects, the
imputed counterfactual is far lower variance than a CS 2x2 that leans on a
6-city never-treated group. Under parallel trends and homoskedasticity BJS
show this estimator is efficient (it attains the semiparametric variance
bound). On this panel the payoff is stark: BJS returns 1.406 against a true
simple ATT of 1.407, with a standard error of 0.033, roughly one-seventh of
CS's 0.25. The framing is the mirror image of Module 1's forbidden
comparison: never extrapolate a treated cell's counterfactual from other
treated cells, only from untreated ones.

The cost is the flip side of the efficiency: imputation leans hard on the
functional form of the untreated model. If parallel trends is shaky, the
imputed $\hat Y(0)$ inherits the misspecification with no propensity-model
insurance. That is why the doubly-robust CS estimator and imputation are
complements, not substitutes.

## de Chaisemartin and D'Haultfoeuille (2020): switchers

dCDH attack the problem one *switch* at a time. Their $\text{DID}_M$ estimator
averages, over every period $t$ where some unit changes treatment status, the
difference between the outcome change of the **switchers** and the outcome
change of units whose treatment is **stable** across $t-1, t$. In our
no-exit staggered design the only switchers are units turning on at their
adoption date, and the stable comparison group is the not-yet-treated (plus
never-treated). The instantaneous estimator is

$$\widehat{\text{DID}}_M = \sum_{t} \frac{N_t^{sw}}{\sum_s N_s^{sw}}
\left[\overline{\Delta Y}^{\,sw}_t - \overline{\Delta Y}^{\,stable\,0}_t\right],$$

the switcher-count-weighted average of "switchers' first difference minus
stayers' first difference" at each switch time. It targets the average
*instantaneous* (event-time-zero) effect. On the `dynamic` scenario, whose
switch-time effect is 0.4 by construction, the hand-coded estimator returns
0.43, within sampling noise. The `DIDmultiplegt` package computes the full
dynamic version and, critically, a **negative-weights diagnostic**: it reports
how much weight a TWFE specification places on cohort-time cells with the
wrong sign, the quantity that lets TWFE flip a uniformly positive effect
negative (Module 1's sign-reversal demonstration).

## Choosing among the four

| Estimator | Package | Reach for it when | Watch out for |
|---|---|---|---|
| Callaway-Sant'Anna | `did` | you want group-time ATTs and flexible aggregation with DR robustness | high variance with a small never-treated group |
| Sun-Abraham | `fixest` | you live in a regression workflow and want a clean event study | equals CS, so no robustness gain over it |
| Borusyak et al. | `didimputation` | parallel trends is credible and you want efficiency | leans on the untreated model; no DR insurance |
| dCDH | `DIDmultiplegt` | treatment turns on and off, or you want the negative-weights audit | instantaneous estimand differs from a full dynamic ATT |

In practice: pick one workhorse (CS or BJS), report its event study with a
uniform band, and check robustness with a second estimator. If they agree,
the staggering was handled; if they diverge, parallel trends or the control
group is doing the work and you investigate before reporting.

## References

- Callaway, B. and Sant'Anna, P. (2021). Difference-in-differences with
  multiple time periods. *Journal of Econometrics*, 225(2), 200-230.
- Sun, L. and Abraham, S. (2021). Estimating dynamic treatment effects in
  event studies with heterogeneous treatment effects. *Journal of
  Econometrics*, 225(2), 175-199.
- Borusyak, K., Jaravel, X. and Spiess, J. (2024). Revisiting event-study
  designs: robust and efficient estimation. *Review of Economic Studies*,
  91(6), 3253-3285.
- de Chaisemartin, C. and D'Haultfoeuille, X. (2020). Two-way fixed effects
  estimators with heterogeneous treatment effects. *American Economic
  Review*, 110(9), 2964-2996.
- Sant'Anna, P. and Zhao, J. (2020). Doubly robust difference-in-differences
  estimators. *Journal of Econometrics*, 219(1), 101-122.
- Roth, J., Sant'Anna, P., Bilinski, A. and Poe, J. (2023). What's trending in
  difference-in-differences? *Journal of Econometrics*, 235(2), 2218-2244.
