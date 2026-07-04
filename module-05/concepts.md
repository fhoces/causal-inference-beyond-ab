# Module 5: Synthetic DiD

The bridge from synthetic control to difference-in-differences. Module 4
built synthetic control (SC) as a weighted match on donor levels; Module 1
diagnosed TWFE/DiD as a global parallel-trends bet. Arkhangelsky, Athey,
Hirshberg, Imbens and Wager (2021) show these are two corners of one
estimator. Synthetic DiD (SDID) sits between them: it reweights *units*
like SC and reweights *time periods* like nobody else, then runs a
weighted two-way fixed effects regression. The payoff is double robustness
in the panel sense: SDID is consistent if *either* the unit weights or the
time weights succeed, and it is more precise than both parents on the
canonical data.

## The estimator, formally

Block treatment structure: units $i = 1, \dots, N$, periods
$t = 1, \dots, T$. The last $N_1 = N - N_0$ units are treated in the last
$T_1 = T - T_0$ periods; $W_{it} = 1$ marks a treated cell. The SDID point
estimate solves a *weighted* two-way fixed effects problem:

$$(\hat\tau, \hat\mu, \hat\alpha, \hat\beta) =
\arg\min_{\tau, \mu, \alpha, \beta}
\sum_{i=1}^{N}\sum_{t=1}^{T}
\left(Y_{it} - \mu - \alpha_i - \beta_t - \tau W_{it}\right)^2
\hat\omega_i\, \hat\lambda_t$$

The unit weights $\hat\omega_i$ and time weights $\hat\lambda_t$ are
computed first, in two separate convex programs, then plugged in. Given the
weights, $\hat\tau$ has a closed form: a weighted double difference (see
below). The whole estimator is therefore "solve for weights, then run
weighted DiD".

### The unit-weight program

The unit weights are chosen so that the weighted average of control units
tracks the treated units' *trajectory* across the pre-period, up to a
constant level shift:

$$(\hat\omega_0, \hat\omega) =
\arg\min_{\omega_0 \in \mathbb{R},\, \omega \in \Omega}
\sum_{t=1}^{T_0}
\left(\omega_0 + \sum_{i=1}^{N_0}\omega_i Y_{it}
- \frac{1}{N_1}\sum_{i=N_0+1}^{N} Y_{it}\right)^2
+ \zeta^2\, T_0\, \lVert\omega\rVert_2^2,$$

over the simplex $\Omega = \{\omega \in \mathbb{R}^{N_0}_{\geq 0} :
\sum_i \omega_i = 1\}$. Two features distinguish this from the SC program of
Module 4:

1. **An intercept $\omega_0$.** It is free (not on the simplex). Adding it
   means the donors only have to match the *shape* of the treated
   pre-trend, not its level. Level differences are absorbed by $\omega_0$
   (and, in the final regression, by the unit fixed effect $\alpha_i$).
2. **A ridge penalty $\zeta^2 T_0 \lVert\omega\rVert_2^2$.** It spreads
   weight across donors and guarantees a unique solution when donors are
   collinear. SC uses no such penalty.

### The regularization parameter

The ridge magnitude is not a free tuning knob; the paper pins it to the
noise scale of the data:

$$\zeta = (N_1 T_1)^{1/4}\, \hat\sigma,
\qquad
\hat\sigma^2 = \frac{1}{N_0 (T_0 - 1)}
\sum_{i=1}^{N_0}\sum_{t=1}^{T_0 - 1}
\left(\Delta_{it} - \bar\Delta\right)^2,
\quad \Delta_{it} = Y_{i, t+1} - Y_{it}.$$

Here $\hat\sigma$ is the standard deviation of *first-differenced* control
outcomes over the pre-period (a scale-free estimate of period-to-period
noise), and $\bar\Delta$ is their mean. On California Prop 99,
$\hat\sigma \approx 5.5$ and $\zeta \approx 10.2$. The $(N_1 T_1)^{1/4}$
factor makes the penalty vanish at the right rate as the panel grows.

### The time-weight program

Symmetrically, the time weights make the pre-periods, weighted, resemble
the post-period for the control units:

$$(\hat\lambda_0, \hat\lambda) =
\arg\min_{\lambda_0 \in \mathbb{R},\, \lambda \in \Lambda}
\sum_{i=1}^{N_0}
\left(\lambda_0 + \sum_{t=1}^{T_0}\lambda_t Y_{it}
- \frac{1}{T_1}\sum_{t=T_0+1}^{T} Y_{it}\right)^2,$$

over $\Lambda = \{\lambda \in \mathbb{R}^{T_0}_{\geq 0} :
\sum_t \lambda_t = 1\}$. The time program uses only a negligible
regularizer (numerical, not statistical): with $T_0$ typically small there
is no collinearity problem to fix. The time weights concentrate on the
handful of pre-periods most predictive of the post-period, so the "before"
baseline is not a flat average over the whole history but a tailored match.
On Prop 99, only three of the nineteen pre-periods carry meaningful weight.

### Local parallel trends

DiD assumes parallel trends *globally*: the raw average of all controls
would have moved parallel to the treated units absent treatment. SDID
weakens this to a *local* condition. The unit weights construct one
synthetic control whose pre-trend already matches the treated trajectory;
parallel trends need only hold for *that* weighted comparison, over the
$\hat\lambda$-weighted periods. You are no longer betting that Utah is a
good counterfactual for California, only that a specific convex combination
of donor states is, and only over the periods the time weights select. This
is why SDID tolerates the interactive-fixed-effects (factor) structures that
sink plain DiD: it matches on the factor loadings implicitly.

## The bridge: DiD, SC, SDID

Both parents are corners of the SDID objective.

| | uniform time weights $\lambda_t = 1/T$ | optimized time weights $\hat\lambda$ |
|---|---|---|
| **uniform unit weights** $\omega_i = 1/N$ | **DiD** | time-weighted DiD |
| **optimized unit weights** $\hat\omega$ | **SC** (drop unit FE) | **SDID** |

- **SDID = DiD** when both weight vectors are uniform. The weighted TWFE
  objective collapses to ordinary TWFE, and $\hat\tau$ is the plain
  two-way-fixed-effects DiD coefficient.
- **SDID = SC** when the time weights put all mass where SC's implicit
  comparison does *and* there is no unit intercept. SC matches the treated
  unit's pre-period levels with a simplex of donors and has no unit fixed
  effect; SDID with $\omega_0 = 0$ and $\alpha_i \equiv 0$ recovers exactly
  that.

### Why the intercept changes everything

The single most consequential difference between SC and SDID is the unit
intercept $\omega_0$ in the weight program (equivalently the unit fixed
effects $\alpha_i$ in the final regression). Pure SC has to match *levels*:
the synthetic control must sit on top of the treated unit throughout the
pre-period, which forces the weights to chase level as well as shape. If no
convex combination of donors reaches the treated unit's level, SC is
biased. SDID absorbs the level gap into $\omega_0$, so the donors only have
to reproduce the treated unit's *trend*. This is exactly the
difference-in-differences logic (levels cancel in the double difference)
imported into the synthetic-control weighting step. It is why SDID needs no
perfect pre-treatment fit, only a parallel one.

## Inference

SDID ships with two variance estimators. Both treat the weights as fixed
(the paper shows the estimation of $\hat\omega, \hat\lambda$ is
first-order negligible).

- **Placebo variance.** Reassign the treatment to control units, one at a
  time (or in blocks matching the treated count), re-estimate on the pure
  control panel where the true effect is zero, and take the variance of the
  placebo estimates. This works with a single treated unit, which is the
  common SC/SDID case, and is the default on Prop 99. It requires a
  homoskedasticity-style assumption across units.
- **Jackknife variance.** Delete one unit at a time, re-estimate, and form
  the usual jackknife variance. It is cheaper and assumption-lighter but
  **requires at least two treated units**: with one treated unit, deleting
  it leaves no treated variation and the estimator is undefined (the
  package returns `NA`). Use it when several units are treated; fall back to
  the placebo estimator otherwise.

A third option, bootstrap, is available for many treated units but is the
most expensive.

## Empirical comparison: California Prop 99

Reproducing the AER paper's headline table live (all three computed with
`synthdid`, none hardcoded):

| Estimator | Effect (packs per capita) | Placebo SE |
|---|---|---|
| DiD | -27.3 | 16.3 |
| SC | -19.6 | 11.0 |
| SDID | -15.6 | 10.0 |

The point estimates are exact; the placebo standard errors are
simulation-based, so they shift by a point or so from run to run (the table
above fixes `set.seed(1)`). DiD is the most negative and least precise: it
forces all 38 donor states
to weight equally and bets on global parallel trends, which the raw donor
average violates (California's cigarette consumption was already on a
different path). SC tightens this by reweighting donors but pays for a
perfect level match. SDID lands between the two in point estimate and has
the smallest standard error: the time weights concentrate the pre-period
comparison and the ridge-regularized unit weights spread risk across donors.
The signature `synthdid_plot()` overlays the treated trajectory, the
$\hat\omega$-weighted synthetic trajectory, the $\hat\lambda$-weighted
pre-period band (the shaded region marking which periods anchor the
baseline), and the parallelogram whose vertical side is $\hat\tau$.

## Application: ride-share factor-structure panel

A cleaner demonstration of *why* SDID beats DiD uses a simulated panel with
an interactive fixed-effects structure, where parallel trends fails by
construction. One treated city and 20 donor cities over 48 months; outcomes
are driven by two latent factors with city-specific loadings, plus a policy
(airport-pickup rule) at month 36 with a true effect of 4.0. The treated
city loads heavily on the trending factor, so the raw donor average drifts
away from its counterfactual and DiD is badly biased upward:

| Estimator | Effect | Placebo SE | Error vs truth (4.0) |
|---|---|---|---|
| SDID | 3.94 | 0.62 | 0.06 |
| SC | 4.31 | 1.63 | 0.31 |
| DiD | 6.29 | 1.76 | 2.29 |

DiD misreads more than half the trending-factor divergence as treatment
effect. SC and SDID both reweight donors toward the high-loading cities and
recover the truth; SDID's time weights and regularization make it the
sharpest of the three. This is the interactive-fixed-effects regime where
SDID and generalized SC (`gsynth`, Module 4) shine and DiD should not be
trusted.

## Practitioner guidance

- Reach for SDID when you have a **panel with a clear pre-period and a
  handful of treated units**, and you suspect the raw control average is not
  a parallel counterfactual (visible pre-trend divergence, factor
  structure, level mismatch).
- SDID **degrades gracefully to DiD**: if the unit weights come out near
  uniform and the time weights near flat, you are in the DiD-is-fine regime
  and SDID confirms it.
- With a **single treated unit**, use the placebo variance estimator; with
  **multiple**, prefer the jackknife.
- SDID does not fix everything: it still assumes a low-rank-plus-noise
  structure and a stable-effect block. For staggered adoption you combine it
  with the timing machinery of Modules 1 and 2, and for irregular missingness
  you reach for matrix completion (Module 8).

## The same problem at an online retailer

When a handful of metros receive next-day delivery upgrades in the same quarter, the single-unit synthetic control no longer applies directly and a full DiD requires parallel trends that the endogenous rollout order makes suspect. Synthetic DiD reweights the donor metros (those still on two-day delivery) to match the pre-launch order trend of the treated group, relaxing the strict parallel-trends assumption that DiD requires. Simultaneously, it downweights remote pre-periods where the common-trend assumption is least credible, concentrating identification on the recent pre-launch window. The unit weights that SDID produces reveal which donor metros anchor the counterfactual: metros with order trajectories and demand characteristics similar to the treated cohort receive the most weight, while dissimilar metros receive near-zero weight. In contrast to pure synthetic control, SDID accommodates multiple treated units naturally, and in contrast to pure DiD it does not assume that untreated metros follow the same trend as treated ones unconditionally. When the treated group is small but larger than one, SDID occupies the right position in the method-choice space.

## References

- Arkhangelsky, D., Athey, S., Hirshberg, D. A., Imbens, G. W. and Wager, S.
  (2021). Synthetic difference-in-differences. *American Economic Review*,
  111(12), 4088-4118.
- Abadie, A., Diamond, A. and Hainmueller, J. (2010). Synthetic control
  methods for comparative case studies. *Journal of the American
  Statistical Association*, 105(490), 493-505.
- Doudchenko, N. and Imbens, G. W. (2016). Balancing, regression,
  difference-in-differences and synthetic control methods: a synthesis.
  NBER Working Paper 22791.
- Xu, Y. (2017). Generalized synthetic control method: causal inference with
  interactive fixed effects models. *Political Analysis*, 25(1), 57-76.
- Clarke, D., Pailañir, D., Athey, S. and Imbens, G. (2023). Synthetic
  difference-in-differences estimation. *IZA Discussion Paper* / Stata and R
  `synthdid` implementation notes.
