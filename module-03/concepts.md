# Module 3: Honest DiD

Sensitivity bounds for parallel trends. The experimentation-refresher's
Module 8 showed the *idea* that a pre-trend test is weak evidence and that
Rambachan and Roth (2023) bound the damage. This module gives the formal
treatment: the restriction sets $\Delta$, the partial-identification
machinery, the robust confidence sets, the breakdown value, and the
`HonestDiD` workflow on the event study from Modules 1 and 2.

## Why a passing pre-trend test is weak evidence

The event-study workhorse regresses the outcome on leads and lags of
treatment relative to a base period (here event time $-1$),

$$y_{it} = \alpha_i + \lambda_t + \sum_{s \neq -1} \beta_s\, \mathbb{1}\{t - g_i = s\} + \varepsilon_{it},$$

and the "pre-trend test" is a joint Wald test that the lead coefficients
$\beta_s$, $s < 0$, are zero. Failing to reject is read as "parallel trends
looks fine." Two problems, both from Roth (2022).

### Low power against the violations that matter

A differential trend that is roughly linear is exactly the shape a pre-trend
test struggles to see, because with a handful of pre-periods and realistic
standard errors the lead coefficients are individually small. In the
single-cohort event study on these slides (four pre-periods, eight treated
cities), a linear differential trend of slope $\delta = 0.10$ log-trips per
week biases the naive first-post estimate and yet is detected only about
58% of the time. Even $\delta = 0.12$ is caught under 80% of the time. The
violations large enough to overturn a conclusion are missed the majority of
the time.

### Conditioning on passing distorts inference

Worse, dropping the analyses that fail the pre-test (a "pre-test estimator")
does not clean up the survivors. Lead and lag coefficients share the base
period, so their sampling errors are positively correlated (about $+0.49$ in
the design here). A sample passes the pre-test when noise happens to pull the
lead coefficients toward zero; the same noise pulls the lag coefficients in
the same direction, so conditional on passing, the post-period estimate is
*more* biased, not less. In the simulation the naive bias at $\delta = 0.10$
is $0.10$ unconditionally but $0.21$ conditional on passing: pre-testing more
than doubles it. This is a selective-inference / pre-test-bias problem, not a
finite-sample curiosity.

The lesson is not "never look at pre-trends." It is: a binary pass/fail hides
both how underpowered the test was and how much bias a violation you could not
detect would create. Report a continuous robustness summary instead.

## The Rambachan-Roth framework

Stack the event-study coefficients as $\hat\beta = (\hat\beta_{pre},
\hat\beta_{post})$ with $\hat\beta \sim \mathcal{N}(\beta, \Sigma)$, and
decompose the estimand into the causal effect and a differential-trend
nuisance,

$$\beta_{post} = \underbrace{\tau_{post}}_{\text{effect}} + \underbrace{\delta_{post}}_{\text{PT violation}}, \qquad \beta_{pre} = \delta_{pre}.$$

Exact parallel trends is $\delta = 0$. Rambachan and Roth drop it and instead
assume only that the full violation vector $\delta$ lies in a researcher-chosen
set $\Delta$,

$$\delta \in \Delta.$$

The pre-period coefficients $\beta_{pre}$ estimate $\delta_{pre}$ directly, so
$\Delta$ *links* the observed pre-period violation to the unobserved
post-period violation $\delta_{post}$. The target parameter is a scalar
$\theta = \ell' \tau_{post}$ (for a weight vector $\ell$; e.g. the first
post-period effect, or an average). Because only the sum $\tau_{post} +
\delta_{post}$ is observed, $\theta$ is **partially identified**: the identified
set is

$$\mathcal{S}(\beta, \Delta) = \left\{ \ell'(\beta_{post} - \delta_{post}) : \delta \in \Delta,\ \delta_{pre} = \beta_{pre} \right\}.$$

Inference targets this set. A robust confidence set $\mathcal{C}_{1-\alpha}$ is
constructed to cover the true $\theta$ with probability at least $1 - \alpha$
uniformly over $\delta \in \Delta$. As $\Delta$ grows, $\mathcal{S}$ widens and
so does $\mathcal{C}$: honesty about the parallel-trends assumption is paid for
in interval width, not in a point estimate that silently moves.

## The restriction menu

The whole method reduces to choosing $\Delta$. Each choice encodes an economic
story about *what kind* of confounding trend you are worried about.

### Smoothness: $\Delta^{SD}(M)$

Bound the discrete second difference of the violation by $M$:

$$\Delta^{SD}(M) = \left\{ \delta : \left| (\delta_{s+1} - \delta_s) - (\delta_s - \delta_{s-1}) \right| \leq M \ \ \forall s \right\}.$$

The story: the differential trend does not *accelerate* by more than $M$ per
period. $M = 0$ forces the violation to be exactly linear, extrapolated from
the pre-period slope; larger $M$ allows curvature. A key and often-missed
consequence: **any linear trend has zero second difference, so it lives in
$\Delta^{SD}(M)$ for every $M \geq 0$.** Under $\Delta^{SD}$ the estimator
extrapolates the pre-period linear trend into the post-period and nets it out.
A purely linear confound is therefore something $\Delta^{SD}$ *corrects for*,
not something it flags. This is the right restriction when your worry is that a
smooth secular trend (adoption momentum, a slow demand shift) contaminates the
comparison.

### Relative magnitudes: $\Delta^{RM}(\bar M)$

Bound the largest post-period change in the violation by $\bar M$ times the
largest pre-period change:

$$\Delta^{RM}(\bar M) = \left\{ \delta : \left| \delta_{s+1} - \delta_s \right| \leq \bar M \cdot \max_{s' < 0} \left| \delta_{s'+1} - \delta_{s'} \right| \ \ \forall s \geq 0 \right\}.$$

The story: whatever process moved the groups apart before treatment cannot
suddenly move them much faster after. $\bar M = 1$ says the post-period
violation is no larger, period for period, than the worst pre-period wiggle you
already see in the data. This restriction keys directly on the *magnitude* of
the estimated pre-period violation, so unlike $\Delta^{SD}$ it does react to a
linear confound: a steeper pre-trend inflates the benchmark and widens the
bounds. It is the natural default when you do not want to assume the trend is
smooth, only that it does not change character at the treatment date.

### Sign and monotonicity restrictions

Additional shape information tightens $\Delta$ further:

- $\Delta^{SDPB}$, $\Delta^{RMB}$ (bias sign): impose $\delta_{post} \geq 0$
  (or $\leq 0$) when you know the direction of the confound (e.g. treated
  cities were on a known upswing).
- $\Delta^{SDI}$, $\Delta^{SDM}$ (monotonicity): impose that the violation is
  monotone in event time.

These combine with $\Delta^{SD}$ or $\Delta^{RM}$ (the "combined" sets), and
each added restriction shrinks the identified set. Impose only what you can
defend from institutional knowledge; a sign restriction you cannot justify buys
a narrower interval you should not trust.

## The breakdown value

Rather than pick one $M$, trace the robust CI as $M$ (or $\bar M$) grows and
report the threshold at which the conclusion flips:

$$M^{\ast} = \sup \left\{ M : 0 \notin \mathcal{C}_{1-\alpha}\big(\Delta^{SD}(M)\big) \right\},$$

the largest smoothness allowance under which the robust CI still excludes zero
(the analogous $\bar M^{\ast}$ for $\Delta^{RM}$). The breakdown value is a
one-number summary of robustness: it converts "is the effect significant?" into
"how large a parallel-trends violation would it take to make it
insignificant?" You then judge whether a violation that large is plausible,
using the observed pre-period wiggle as the yardstick (this is why $\bar M^\ast$
is especially interpretable: $\bar M^\ast = 1.5$ means the post-period violation
would have to be 1.5 times the worst pre-period one).

## Robust confidence sets: how they are built

Two families, from Rambachan and Roth:

- **Fixed-length CI (FLCI).** A CI of fixed length centered on an affine
  estimator, optimal when $\Delta$ is convex and centrosymmetric (like
  $\Delta^{SD}(M)$). It solves a min-max problem for the shortest interval with
  correct worst-case coverage.
- **Conditional and hybrid (ARP).** The moment-inequality approach of Andrews,
  Roth and Pakes: invert a conditional test over the identified set. The
  conditional least-favorable hybrid (**C-LF**) is the recommended default for
  non-centrosymmetric sets such as $\Delta^{RM}(\bar M)$, where FLCI does not
  apply.

On the machine used to build this module the FLCI solver is unavailable, so all
runs here use **C-LF**, which is valid for both $\Delta^{SD}$ and $\Delta^{RM}$.
In practice, for $\Delta^{SD}$ the FLCI is usually a touch shorter; the
breakdown-value logic is identical.

## The running application

We take the shared 30-city rollout DGP from Modules 1 and 2 and carve out a
clean single-cohort event study: the eight cities that adopt at week $g = 25$
versus the six never-treated cities, event window $-5$ to $+5$, estimated with
`fixest::feols(y ~ i(rel_time, ref = -1) | city + t)`, clustered by city. Never-
treated cities are pinned to the base period so they serve as controls at every
calendar week through the time fixed effects.

A practical constraint fixes the window: a clustered variance matrix has rank at
most (number of clusters $- 1$). With 14 clusters (8 treated $+$ 6 never), an
event study with 20 coefficients ($-10$ to $+10$) yields a rank-deficient
$\Sigma$ and the robust bounds explode. The $-5$ to $+5$ window has 10
coefficients, comfortably below the rank ceiling. Watch the coefficient count
against the cluster count whenever you feed an event study into `HonestDiD`.

### The workflow and a mandatory sanity check

`HonestDiD` takes two inputs: `betahat`, the event-study coefficient vector
ordered (earliest pre $\dots$ $-2$, then $0 \dots K$) with the base period
omitted, and `sigma`, its clustered vcov. The `fixest` `i()` coefficients come
out in exactly that order. Before trusting any sensitivity output, confirm that
`constructOriginalCS` reproduces the conventional `feols` confidence interval
for the target period. Here the first-post `feols` estimate is $0.888$ (SE
$0.220$, CI $[0.457, 1.319]$) and `constructOriginalCS` returns the same
$[0.457, 1.319]$: the coefficient ordering is right.

### Clean case: parallel trends holds by construction

Parallel trends holds by design (city fixed effects, a common weekly trend, no
differential slope), and the pre-test passes ($p = 0.16$). Targeting the first
post-period effect:

| Restriction | Breakdown value | Reading |
|---|---|---|
| $\Delta^{SD}$ | $M^{\ast} = 0.45$ | violation could accelerate by up to 0.45/period before significance is lost |
| $\Delta^{RM}$ | $\bar M^{\ast} = 1.50$ | post-period violation could be 1.5x the worst pre-period one |

Both are large relative to anything the pre-period suggests, so the effect is
robust and you can say so with a number instead of "the pre-trends looked flat."

### The target vector $\ell$ matters

Switching the target from the first post-period to the average of all six
post-periods collapses robustness: the average's breakdown values are
$M^{\ast} = 0.05$ and $\bar M^{\ast} = 0.50$. Under $\Delta^{SD}$ the worst-case
bias at post-period $s$ grows like $M \cdot s^2$, so an average that loads on
late post-periods is exposed to compounding extrapolation. The first post-period
is the most robust target; a long-run average is the least. Choose $\ell$ to
match the estimand you actually care about, and know that later horizons are
inherently more fragile to trend extrapolation.

### Confounded case: the restriction must match the threat

Inject a differential linear trend of slope $0.10$ into the treated cohort. Now
the pre-test rejects sharply ($p < 10^{-3}$). The breakdown values split:

- $\Delta^{RM}$: $\bar M^{\ast}$ falls from $1.50$ to $1.00$ - the steeper
  pre-period violation inflates the benchmark and the bounds widen.
- $\Delta^{SD}$: $M^{\ast}$ stays at $0.45$ - a linear trend has zero second
  difference, so $\Delta^{SD}$ extrapolates and nets it out. The smoothness
  restriction is, by construction, blind to a linear confound.

This is the central practitioner point. If the confounding trend you fear is
smooth and linear, $\Delta^{SD}$ *assumes it continues and corrects for it* (a
strong assumption). If you are unwilling to assume the post-period trend behaves
like the pre-period one, $\Delta^{RM}$ is the honest choice. The restriction is
not a knob to tune until you like the answer; it is a statement about the
economics of the confounder.

### Do not over-read a single rejection

In the confounded sample the pre-test rejected, so here it "worked." But that is
one draw. Whether the test *reliably* catches a violation is a repeated-sampling
property, and the simulation above shows a slope-$0.10$ trend is missed roughly
40% of the time. The gap between "this sample rejected" and "the test reliably
rejects" is the whole Roth (2022) argument, and it is why the deliverable is a
breakdown value, not a pre-test $p$.

## Practitioner workflow

1. Estimate the event study; confirm the coefficient count is below the cluster
   rank ceiling.
2. Validate `constructOriginalCS` against the `feols` CI for the target period.
3. Choose $\Delta$ to match the confounder you fear: $\Delta^{SD}$ for a smooth
   trend, $\Delta^{RM}$ for a "no worse than pre-period" bound; add sign or
   monotonicity only if defensible.
4. Report the breakdown value $M^{\ast}$ / $\bar M^{\ast}$, and calibrate it
   against the observed pre-period violation.
5. Report the sensitivity plot, not just the pass/fail of a pre-test.

## References

- Rambachan, A. and Roth, J. (2023). A more credible approach to parallel
  trends. *Review of Economic Studies*, 90(5), 2555-2591.
- Roth, J. (2022). Pretest with caution: event-study estimates after testing
  for parallel trends. *American Economic Review: Insights*, 4(3), 305-322.
- Andrews, I., Roth, J. and Pakes, A. (2023). Inference for linear conditional
  moment inequalities. *Review of Economic Studies*, 90(6), 2763-2791.
- Roth, J., Sant'Anna, P., Bilinski, A. and Poe, J. (2023). What's trending in
  difference-in-differences? A synthesis of the recent econometrics literature.
  *Journal of Econometrics*, 235(2), 2218-2244.
- Manski, C. (2003). *Partial Identification of Probability Distributions*.
  Springer. (The partial-identification tradition the framework sits in.)
