# Module 8: Matrix Completion and the Modern Panel Toolbox

The course closes where panel causal inference is heading: treat the whole
problem as missing data. Athey, Bayati, Doudchenko, Imbens and Khosravi
(2021) show that TWFE (Module 1), the heterogeneity-robust DiD estimators
(Module 2), synthetic control (Module 4) and synthetic DiD (Module 5) are
all, underneath, rules for filling in a matrix of untreated potential
outcomes. Matrix completion (MC-NNM) is the member of that family that
imposes the least structure: it works on any missingness pattern, including
the irregular ones none of the earlier estimators can handle.

## The imputation view

Every panel causal problem has a matrix of potential outcomes under
control, $Y_{it}(0)$, for units $i = 1, \dots, N$ and periods
$t = 1, \dots, T$. You observe $Y_{it}(0)$ for control cells; for treated
cells it is missing, not caused by treatment, simply unobserved. Once you
have an estimate $\hat Y_{it}(0)$ for every missing cell, the ATT falls out
directly:

$$\widehat{ATT} = \frac{1}{|\mathcal{M}|}\sum_{(i,t) \in \mathcal{M}}
\big(Y_{it} - \hat Y_{it}(0)\big)$$

where $\mathcal{M}$ is the set of treated (missing) cells. Every course
estimator differs only in *how* it fills those cells:

- **TWFE / DiD** fills a missing cell with $\hat\alpha_i + \hat\beta_t$: a
  pure two-way fixed-effects imputation, no interactive structure at all.
- **Synthetic control** fills a treated unit's missing cells with a convex
  combination of donor units' outcomes at the same period.
- **Synthetic DiD** fills them with a weighted two-way FE fit, unit weights
  matching pre-period shape, time weights matching post-period relevance.
- **`gsynth` / interactive fixed effects (IFE)** fills them with a
  factor model of fixed rank $r$, chosen by cross-validation.
- **Matrix completion (MC-NNM)** fills them with a low-rank matrix whose
  effective rank is chosen continuously, by nuclear-norm shrinkage.

The missingness pattern itself varies by application: one unit going
missing in a block at the end of the panel (single-market policy change),
many units in one simultaneous block (a coordinated rollout), staggered
absorbing missingness (a rollout ops schedule city by city), or genuinely
irregular missingness (data gaps plus treatment, or treatment on-again
off-again). Only matrix completion handles the last case; the others
require some form of block or absorbing structure.

## The MC-NNM estimator

Let $\mathcal{O}$ denote the observed (untreated) cells. The estimator
solves for a low-rank matrix $L$ and unpenalized two-way fixed effects
jointly:

$$\min_{L,\, \alpha,\, \beta}\; \frac{1}{|\mathcal{O}|}
\sum_{(i,t) \in \mathcal{O}} \big(Y_{it} - L_{it} - \alpha_i - \beta_t\big)^2
+ \lambda \lVert L \rVert_*$$

Two design choices matter as much as the objective itself:

1. **The fixed effects are not penalized.** $\alpha_i$ and $\beta_t$ absorb
   level and trend structure every unit or every period shares. Leaving
   them inside $L$ would force the nuclear-norm penalty to pay for that
   shared structure, wasting rank budget that should go to the genuinely
   interactive (unit-by-time) part of the residual.
2. **$\lambda$ is chosen by cross-validation on the observed cells only.**
   Hide a random subset of *observed* control cells, refit, and pick the
   $\lambda$ that best predicts the held-out subset. The treated cells are
   never part of this CV loop, which is exactly why the estimator can
   under-extrapolate into them (see below).

Given $\hat L, \hat\alpha, \hat\beta$, the counterfactual for every missing
cell is $\hat Y_{it}(0) = \hat L_{it} + \hat\alpha_i + \hat\beta_t$, and the
ATT is the average gap over the treated cells.

## Why the nuclear norm

Minimizing the rank of $L$ directly is combinatorial and NP-hard in
general. The nuclear norm

$$\lVert L \rVert_* = \sum_k \sigma_k(L)$$

(the sum of $L$'s singular values) is its convex relaxation, playing the
same role for a matrix that the lasso's $\ell_1$ penalty plays for a
vector: it shrinks every singular value toward zero and sets the smallest
ones exactly to zero, choosing the *effective* rank of $L$ continuously
rather than fixing it in advance.

### The soft-impute algorithm

Mazumder, Hastie and Tibshirani (2010) give a simple fixed-point algorithm.
Starting from $L = 0$: fill the missing cells of the (FE-residualized)
outcome matrix with the current $L$, take an SVD of the filled matrix,
soft-threshold every singular value by $\lambda$ (subtract $\lambda$, floor
at zero), reassemble $L$ from the thresholded SVD, and repeat to a fixed
point. Each step is an exact solution to a proximal problem,
$\min_L \tfrac12 \lVert Z - L \rVert_F^2 + \lambda \lVert L \rVert_*$, so the
iteration is a majorize-minimize scheme that provably converges to the
minimizer of the nuclear-norm objective over the observed cells. The
exercise hand-codes this loop and checks it against `gsynth(estimator =
"mc")`.

## One family: SC, SDID, `gsynth`, MC-NNM

| Method | Imputation model | Regularization | Pattern it needs |
|---|---|---|---|
| SC (Module 4) | convex combination of donors matching pre-period levels | simplex constraint, no intercept | block, one or a few treated units |
| SDID (Module 5) | weighted two-way FE regression | simplex-plus-ridge weights, local parallel trends | block |
| `gsynth` / IFE | factor model, rank fixed by CV, loadings by least squares | none beyond the hard rank cutoff | staggered adoption is fine; needs enough pre-periods per treated unit |
| MC-NNM (this module) | low-rank matrix plus two-way FE | continuous nuclear-norm shrinkage | any pattern, including irregular missingness |

Athey et al.'s own framing is useful for interviews: SC-style estimators
run "horizontal" regressions, exploiting similarity across *units* at a
given time; DiD-style estimators run "vertical" regressions, exploiting
similarity across *time* for a given unit. Matrix completion is the member
of the family that uses both directions simultaneously, through the
two-way fixed effects plus a shared low-rank residual structure.

## What the demos show

**Demo A (a staggered launch with selection on growth).** 60 cities, 40
weeks, two latent factors (a mildly trending, persistent market factor and
a seasonal factor) with city-specific loadings. Ops launch the
fastest-growing (highest-loading) cities first, so launch timing correlates
with the trending factor by construction; parallel trends is false. 25
treated cities across 12 distinct launch weeks, 35 never-treated. True
effect 2.0.

| Estimator | Estimate | Error vs truth |
|---|---|---|
| TWFE | 3.49 | +1.49 |
| Callaway-Sant'Anna | 3.15 | +1.15 |
| `gsynth` (IFE, rank 2 by CV) | 2.07 | +0.07 |
| `gsynth` (MC-NNM) | 3.66 | +1.66 |
| soft-impute (hand-coded) | 3.67 | +1.67 |

TWFE and CS both inherit the launch-timing selection bias. Fixed-rank IFE,
because the true factor structure is exactly rank 2 and CV finds it,
recovers the truth almost exactly. Matrix completion, whether via
`gsynth`'s own `mc` estimator or the hand-coded soft-impute (the two agree
to within 0.01), lands close to the TWFE answer, not the truth.

**Why.** Cross-validation for $\lambda$ only ever sees observed control
cells; it never sees the treated corner it is asked to predict. Shrinking
singular values flattens exactly the growth-factor divergence that needs
to be extrapolated into that corner. A noiseless version of the same DGP
(zero noise, no treatment effect, so any imputation error is pure
extrapolation bias) makes the point starkly: FE-only imputation has corner
RMSE about 1.96 against a corner scale (SD) of about 6.3; nuclear-norm
shrinkage at several values of $\lambda$ gets RMSE 1.85-1.93, barely
better than ignoring the factor structure entirely; a rank-2 `gsynth`/IFE
fit recovers the corner to machine precision (RMSE 0). This is structural,
not a tuning failure: shrinkage is a conservative extrapolation policy, and
the corner is precisely where it is most conservative.

**Demo B (the placebo-imputation tournament).** Athey et al.'s own
validation recipe, applied here to a 40-city, no-treatment version of the
same factor panel (noise floor 0.5). Hide cells you actually observe under
four patterns, impute with each method, score RMSE against ground truth
you happen to know.

| Pattern | FE only | `gsynth` IFE | `gsynth` MC | soft-impute |
|---|---|---|---|---|
| one unit, block | 2.22 | 0.41 | 0.88 | 1.18 |
| 10 units, block | 2.19 | 0.55 | 1.82 | 1.81 |
| 20 units, staggered | 2.25 | 0.97 | 2.12 | 2.16 |
| 15% random holes | 1.03 | cannot run | cannot run | 0.59 |

Two robust conclusions. First, FE-only imputation (exactly the
Borusyak-Jaravel-Spiess counterfactual of Module 2) is dominated by every
factor-aware method on every pattern: parallel-trends imputation pays for
ignoring the factor structure everywhere, not just under selection. Second,
fixed-rank IFE wins whenever it can run (it exploits the exact rank), but
`gsynth` requires an *absorbing* treatment pattern (once missing, always
missing within a unit), so it cannot even attempt the random-holes case.
Matrix completion is the only method in the table that runs on every
pattern, and on the one pattern that rules out everything else, it clearly
beats FE-only imputation.

## Shrinkage is an extrapolation policy you choose

The two demos together are the module's central lesson. Regularization
choice is not a technical footnote, it is a bet about how far you are
willing to extrapolate:

- **Fixed-rank factor estimation (`gsynth`/IFE)** extrapolates aggressively
  once it commits to a rank: no shrinkage inside that rank, so a
  well-specified low-rank structure is recovered almost exactly, even deep
  into a treated corner far from any observed cell.
- **Nuclear-norm shrinkage (MC-NNM)** extrapolates conservatively: it
  never fully trusts a factor structure it cannot validate on observed
  data, so it degrades gracefully when the rank is unknown, the pattern is
  irregular, or the panel is thin, at the cost of under-correcting a
  confound that a correctly-specified fixed-rank model would catch.

Neither dominates. The practical rule: a strongly trending confound plus a
long extrapolation horizon, with enough pre-periods per treated unit,
favors `gsynth`/IFE. An irregular missingness pattern, a panel with real
data gaps, or genuine uncertainty about the rank favors matrix completion.
Either way, validate with Demo B's recipe on your own panel before trusting
the answer: hide cells you actually observe, impute, and see which method's
errors you can live with.

## Staggered adoption across the toolbox

How the block-treatment tools of Modules 4-5 extend, or fail to, under
staggered timing:

| Pattern | What runs |
|---|---|
| One block, all units treated together | SC, SDID, CS, BJS, `gsynth`, MC-NNM |
| Staggered, absorbing | CS, BJS, `gsynth`/IFE, MC-NNM (not plain SC or SDID, which need a single block) |
| Irregular or genuinely missing cells | MC-NNM only |

The imputation view unifies the entire course. Borusyak-Jaravel-Spiess
imputation (Module 2) is *exactly* the FE-only imputer in the tournament
above, MC-NNM with $L$ constrained to zero: no low-rank term, pure additive
fixed effects. `gsynth` is MC-NNM with a hard, CV-selected rank in place of
continuous nuclear-norm shrinkage. Every difference-in-differences
estimator in this course sits on the same regularization dial, from "no
low-rank term at all" (BJS) through "hard-cutoff rank" (`gsynth`) to
"continuous shrinkage" (MC-NNM).

## When to reach for matrix completion

- **Irregular adoption or genuinely missing panel cells.** `gsynth`, SC,
  and SDID all require a structured, absorbing pattern. MC-NNM is the only
  estimator here that runs on arbitrary missingness, and Demo B shows it
  wins there.
- **You are not confident in the rank.** Continuous shrinkage degrades
  gracefully if the true factor count is misjudged; a hard-rank model can
  be badly mis-specified with the wrong $r$.
- **Many treated units, a thin control pool.** This is the paper's own
  motivating case (a hard-rank factor model is harder to estimate with few
  controls), but treat it as a reason to *test*, not a guarantee: in a
  single check on the tournament DGP with 35 of 40 units treated,
  fixed-rank IFE still beat soft-impute (RMSE 0.72 vs 1.15). Run Demo B's
  recipe on your own panel rather than assuming.

**Caveat.** A strongly trending confound plus a long extrapolation horizon,
with enough pre-periods per treated unit, favors `gsynth`/IFE (Demo A).
Nuclear-norm shrinkage is the safer default when the pattern rules out
fixed-rank methods entirely, or when you are unsure the rank is right.
Always validate with the placebo-imputation exercise before trusting either
estimator on a new panel.

## Long-run outcomes from short panels: the surrogate index

A different problem entirely, but one that shares the "the data you have
is not the data you need" theme: the experimental panel runs 8 weeks, but
the decision depends on 12-month retention. Waiting a year to evaluate
every launch is not viable.

Athey, Chetty, Imbens and Kang (2019) combine two data sources:

1. **The experiment itself**, which records treatment and short-run
   surrogates $S$ (trips, active days, cancellation rate) but not the
   long-run outcome $Y_{\text{long}}$.
2. **Historical data**, from units observed for the full horizon, linking
   those same surrogates to $Y_{\text{long}}$.

Fit the surrogate index $\hat h(S, X) = \hat E[Y_{\text{long}} \mid S, X]$
on the historical sample, then use $\hat h(S_i, X_i)$ as the outcome inside
the experiment. The resulting treatment-effect estimate targets the
long-run effect without ever observing the long-run panel for the
experimental units.

**Assumptions.** Surrogacy: treatment affects the long-run outcome only
through the observed surrogates, $Y_i(\text{long}) \perp W_i \mid S_i,
X_i$. Comparability: the surrogate-to-long-run relationship estimated
historically also holds in the experimental sample. Overlap: the
historical sample covers the surrogate values the experiment actually
produces. A useful bonus when the assumptions hold: the index also reduces
noise, since it averages away idiosyncratic long-run variation the
surrogates do not explain, so the experiment's effective sample size for
the long-run question is larger than running the full panel would give
directly.

**Failure modes.** A treatment channel that bypasses the surrogates
entirely, for example a fee change that leaves short-run usage flat but
slowly erodes trust and long-run retention, breaks surrogacy silently: the
index would report no effect while the true long-run effect is nonzero.
Distribution shift between the historical and experimental samples (a
changed product, a changed market) breaks comparability. Neither failure
is directly testable from the experiment alone; both need domain judgment
about the mechanism connecting treatment to the long-run outcome.

## Choosing the method: the course's capstone

- **Can you randomize? Then randomize.** Everything below answers "what do
  you do when the experiment you wanted is not available": a
  platform-wide launch, one treated market, a legal or PR constraint, or a
  retroactive question.
- **Staggered rollout with clean never- or late-treated controls** &rarr;
  Callaway-Sant'Anna or BJS imputation (Module 2), audited with a
  Goodman-Bacon decomposition (Module 1) and stress-tested with HonestDiD
  sensitivity (Module 3).
- **One treated market with a good donor pool** &rarr; synthetic control,
  augmented SC, or synthetic DiD (Modules 4-5), with placebo or conformal
  inference.
- **Irregular adoption, holes in the panel, or many treated units** &rarr;
  matrix completion or `gsynth` (this module), validated by hiding observed
  cells before trusting either.
- **Need to know who benefits** &rarr; causal forest (Module 6), with
  calibration, best-linear-projection, and RATE diagnostics before any
  subgroup claim.
- **Need to decide who gets the treatment** &rarr; policy tree with
  cross-fitted off-policy evaluation (Module 7).
- **Long-run outcome, short panel** &rarr; the surrogate index, layered on
  top of whichever design above supplies the experiment or
  quasi-experiment.

## Practitioner guidance

- Start from the pattern, not the method. Draw the missingness pattern
  first (one block, many blocks, staggered, irregular); it rules out most
  of the toolbox before you touch an estimator.
- Cross-validate $\lambda$ (or $r$, for `gsynth`) only against cells you
  actually observe, never against treated cells, and report sensitivity to
  the tuning choice rather than a single point.
- Validate any panel estimator, before trusting its ATT, by hiding cells
  you already know and scoring the imputation. This is cheap, requires no
  new data, and directly answers whether the estimator's extrapolation
  policy fits your setting.
- Treat "irregular pattern" and "thin control pool" as two different
  reasons to consider matrix completion, and check the second empirically:
  it does not always favor MC over a fixed-rank alternative.
- For long-horizon outcomes, reach for the surrogate index only when you
  can defend surrogacy and comparability qualitatively; neither is testable
  from the experiment alone.

## References

- Athey, S., Bayati, M., Doudchenko, N., Imbens, G. and Khosravi, K.
  (2021). Matrix completion methods for causal panel data models. *Journal
  of the American Statistical Association*, 116(536), 1716-1730.
- Athey, S., Chetty, R., Imbens, G. W. and Kang, H. (2019). The surrogate
  index: combining short-term proxies to estimate long-term treatment
  effects more rapidly and precisely. NBER Working Paper 26463.
- Xu, Y. (2017). Generalized synthetic control method: causal inference
  with interactive fixed effects models. *Political Analysis*, 25(1),
  57-76.
- Arkhangelsky, D., Athey, S., Hirshberg, D. A., Imbens, G. W. and Wager,
  S. (2021). Synthetic difference-in-differences. *American Economic
  Review*, 111(12), 4088-4118.
- Borusyak, K., Jaravel, X. and Spiess, J. (2024). Revisiting event-study
  designs: robust and efficient estimation. *Review of Economic Studies*.
- Mazumder, R., Hastie, T. and Tibshirani, R. (2010). Spectral
  regularization algorithms for learning large incomplete matrices.
  *Journal of Machine Learning Research*, 11, 2287-2322.
- Rambachan, A. and Roth, J. (2023). A more credible approach to
  parallel trends. *Review of Economic Studies*, 90(5), 2555-2591.
- Liu, L., Wang, Y. and Xu, Y. (2024). A practical guide to
  counterfactual estimators for causal inference with time-series
  cross-sectional data. *American Journal of Political Science*, 68(1),
  160-176.
