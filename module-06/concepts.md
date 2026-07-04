# Module 6: Causal Forest

Honest splitting, asymptotics, and HTE diagnostics. The
experimentation-refresher tour showed *that* a causal forest estimates
$\tau(x)$ with honest trees. This module gives the formal treatment: why
honesty plus subsampling delivers a pointwise central limit theorem, the
generalized-random-forest view of a forest as an adaptive kernel, the
R-learner orthogonalization that `grf` runs under the hood, and the
diagnostics that tell you whether the heterogeneity you found is real and
targetable.

## From the ATE to the CATE function

The average treatment effect answers "should we ship?". The conditional
average treatment effect (CATE)

$$\tau(x) = \mathbb{E}[Y(1) - Y(0) \mid X = x]$$

answers "to whom?". Estimating $\tau(x)$ as a smooth function of covariates
is a nonparametric regression problem, but with a twist: $\tau(x)$ is never
observed, not even noisily, because no unit reveals both potential outcomes.
A causal forest sidesteps this by estimating $\tau(x)$ through *local*
treatment-effect comparisons rather than by regressing an observed label.

The prize is not the point estimate alone. It is a **pointwise confidence
interval**: an interval for $\tau(x)$ at a fixed query point $x$ with valid
asymptotic coverage. That is what separates a causal forest from a black-box
CATE model such as a two-model T-learner or a gradient boosting machine,
which give predictions but no honest uncertainty.

## Honest trees (Wager and Athey 2018)

### The honesty split

A regression tree that uses the same data to (a) choose where to split and
(b) estimate the response inside each leaf is *adaptive*, and adaptive leaf
estimates are biased. The split is chosen to make the leaves look as
different as possible, so the estimate inside a leaf capitalizes on the same
noise that placed the split there. This is a within-sample winner's curse:
apparent heterogeneity appears even when the true effect is constant.

**Honesty** breaks the dependence. Each tree draws a subsample, then randomly
partitions it into two disjoint halves:

- the **splitting half** ($\mathcal{J}$) decides the tree structure: which
  covariate, which threshold, at every node;
- the **estimation half** ($\mathcal{I}$) is dropped down the finished tree,
  and the leaf treatment effect is computed only from these held-out points.

Because $\mathcal{I}$ played no role in placing the splits, the leaf estimate
is (conditionally) unbiased for the average effect over that leaf's region.
The forest averages honest trees over many subsamples.

### Subsampling, not bootstrap

Wager and Athey build the forest from **subsamples** of size $s$ drawn
without replacement, not bootstrap resamples. The reason is inferential: a
forest built from subsamples is a *U-statistic* (an average of a symmetric
kernel over subsets), and the Hajek projection of a U-statistic is
asymptotically normal. Bootstrap resampling breaks the clean U-statistic
structure needed for the variance estimator.

The subsample size must grow, but not too fast:
$s = n^{\beta}$ with $\beta \in (0, 1)$ and, for the CLT, $\beta$ bounded
below (roughly $\beta > 1 - (\log \text{something})^{-1}$ in their Theorem 3;
the practical reading is "let $s/n \to 0$ slowly"). Trees must also be
**$\alpha$-regular**: every split leaves at least a fraction $\alpha$ of the
parent's observations on each side, and leaves keep shrinking, so no leaf
locks onto a vanishing neighborhood.

### Asymptotic normality: the result and the conditions

**Theorem (Wager and Athey 2018, informal).** Under honesty,
$\alpha$-regularity, a subsample rate $s = n^\beta$ in the admissible range,
overlap ($0 < \eta \le e(x) \le 1 - \eta$), and Lipschitz-type smoothness of
$\tau(\cdot)$ and the nuisance functions, the causal-forest estimate is
pointwise consistent and asymptotically Gaussian:

$$\frac{\hat\tau(x) - \tau(x)}{\sqrt{\operatorname{Var}[\hat\tau(x)]}}
\;\xrightarrow{d}\; \mathcal{N}(0, 1),$$

and the variance is consistently estimated by an **infinitesimal jackknife**
(the Bayesian-bootstrap-of-little-jackknives estimator), which `grf` returns
as `variance.estimates`. Honesty is what makes the bias asymptotically
negligible relative to the standard error, so the ratio centers at zero and
the interval $\hat\tau(x) \pm 1.96\,\widehat{\text{se}}(x)$ covers.

### What honesty buys and what it does not

You get **valid pointwise inference**: a CI at each fixed $x$. You do **not**
get **uniform inference**: the theorem says nothing about
$\sup_x |\hat\tau(x) - \tau(x)|$, so you cannot read a confidence *band* over
the whole covariate space off the pointwise intervals, and you cannot use
them for multiplicity-correct statements like "the effect is positive
everywhere". For questions about the shape of $\tau(\cdot)$ as a whole, use
the aggregate diagnostics below (calibration test, best linear projection,
RATE), not a naive union of pointwise CIs. Honesty also costs statistical
efficiency: each tree uses only half its subsample to estimate, so honest
forests have higher variance than adaptive ones at fixed $n$. The trade is
bias (which invalidates inference) for variance (which only widens it).

## Generalized random forests (Athey, Tibshirani and Wager 2019)

### Forests as adaptive nearest neighbors

GRF reframes the forest as a way to produce **data-adaptive kernel weights**.
For a target point $x$, tree $b$ drops $x$ to a leaf $L_b(x)$; a training
point $i$ gets weight $1/|L_b(x)|$ if it shares that leaf and 0 otherwise.
Averaging over $B$ trees gives

$$\alpha_i(x) = \frac{1}{B}\sum_{b=1}^{B}
\frac{\mathbb{1}\{X_i \in L_b(x)\}}{|L_b(x)|},
\qquad \sum_i \alpha_i(x) = 1.$$

The forest is then a **weighted nearest-neighbor** method whose metric is
learned from the data: points that repeatedly land in $x$'s leaf are treated
as its neighbors, and the splitting rule decides "near" in the directions
that matter for the effect. `grf::get_forest_weights` returns the full matrix
of $\alpha_i(x)$.

### The local moment condition

GRF estimates a parameter $\theta(x)$ defined by a **local estimating
equation** (moment condition):

$$\mathbb{E}\!\left[\psi_{\theta(x)}(O_i) \mid X_i = x\right] = 0,$$

solved locally with the forest weights,

$$\hat\theta(x) = \arg\min_\theta
\left\| \sum_i \alpha_i(x)\, \psi_\theta(O_i) \right\|.$$

For the causal forest the score is the residual-on-residual moment
$\psi_\tau = (Y_i - \hat m(X_i) - \tau (W_i - \hat e(X_i)))(W_i - \hat e(X_i))$,
so the local solution is a **weighted residualized regression** of outcome on
treatment, exactly a local Robinson (1988) partially linear estimator. The
same machinery with a different $\psi$ gives quantile forests, instrumental
forests (score = the IV moment), and survival forests: this is why GRF is
"generalized".

### Gradient-based splitting

Solving the moment exactly at every candidate split would be too slow, so GRF
splits on **influence-function pseudo-outcomes**. It computes the gradient of
the score at the parent node's estimate, forms a one-step (Newton)
pseudo-outcome $\rho_i$ for each observation, and then runs a fast standard
regression-tree split that maximizes the between-child variance of $\rho_i$.
This approximates the split that would most increase heterogeneity in
$\theta(x)$ while costing about the same as a CART split.

## Local centering and the R-learner

Before growing the forest, `grf` **orthogonalizes** (local centering,
Robinson-style): it fits $\hat m(x) = \mathbb{E}[Y \mid X = x]$ and
$\hat e(x) = \mathbb{E}[W \mid X = x]$ with separate regression forests, then
grows the causal forest on the residuals $Y - \hat m(X)$ and $W - \hat e(X)$.
This is the **R-learner** objective (Nie and Wager 2021): $\tau(\cdot)$
minimizes

$$\sum_i \big[(Y_i - \hat m(X_i)) - \tau(X_i)(W_i - \hat e(X_i))\big]^2.$$

Orthogonalization makes the estimate **Neyman-orthogonal**: first-order
insensitive to small errors in $\hat m$ and $\hat e$. Two payoffs. First,
confounding by a smooth prognostic signal (a covariate that shifts $Y$ and is
correlated with $W$) is projected out through $\hat m$, so the forest spends
its splits on effect *heterogeneity* rather than re-discovering the main
effect. Second, even in a randomized experiment where $e(x) \equiv 0.5$ by
design, chance covariate imbalance means the finite-sample association between
$W$ and $X$ is not exactly zero; centering on $\hat m(x)$ removes the
prognostic-covariate variance from the outcome and sharpens $\hat\tau(x)$. In
an RCT you should set the propensity to its known value,
`W.hat = 0.5`, so the forest does not waste data estimating $e(x)$ and cannot
introduce propensity noise.

## Tuning grf

The parameters that actually move results:

- `honesty.fraction` (default 0.5): share of the subsample used to *place*
  splits; the rest estimates leaves. Lower it toward 0 and you recover an
  adaptive forest (better fit, invalid CIs).
- `min.node.size`: floor on leaf size; larger means more smoothing and more
  stable, wider-support estimates.
- `sample.fraction` (default 0.5): subsample size $s/n$. Governs the
  bias-variance and the CLT rate.
- `mtry`: covariates tried per split. Matters when there are many noise
  covariates.
- `honesty.prune.leaves` (default TRUE): prunes leaves that end up empty of
  treated or control units in the estimation half, avoiding undefined leaf
  effects.

`tune.parameters = "all"` cross-validates these against the R-learner
"debiased error". In practice `num.trees` (use enough for stable variance
estimates, 2000 or more), `min.node.size`, and `sample.fraction` move results
the most; `mtry` matters mainly with many junk covariates; and turning
honesty off changes the *inference*, not just the fit.

## HTE diagnostics

Run these in order; each answers a sharper question than the last.

### Calibration test

`test_calibration` fits the "best linear predictor" of the true effect on two
constructed regressors: the **mean** forest prediction and the
**differential** (demeaned) forest prediction. A `mean.forest.prediction`
coefficient near 1 says the forest is calibrated on average (its ATE is
right); a `differential.forest.prediction` coefficient near 1 and significant
says the forest's *ranking* of who has large vs small effects carries real
signal. It is an omnibus test for the presence of heterogeneity.

### Best linear projection

`best_linear_projection` regresses the (AIPW-scored) CATE on chosen
covariates and reports the coefficients of the best linear approximation to
$\tau(x)$, with heteroskedasticity-robust standard errors. It is the
interpretable summary you show stakeholders: "the effect rises with density
and falls with tenure". A covariate with a coefficient indistinguishable from
zero is not driving heterogeneity, which is how a true nuisance covariate
reveals itself.

### RATE and the TOC curve

The **Targeting Operator Characteristic** (TOC) curve plots, against the
fraction $q$ of the population you treat when you treat the highest
$\hat\tau$ units first, the average benefit among that top fraction $q$ minus
the overall ATE. The **Rank-Weighted Average Treatment Effect**
(`rank_average_treatment_effect`) is the area under the TOC (AUTOC, or the
Qini weighting), estimated on held-out data with a standard error. A RATE
significantly above zero is the direct answer to "is there heterogeneity I
can *target*?". Crucially, RATE evaluated with a **bad** priority (a nuisance
covariate) returns a null, so it also validates *which* signal is targetable.

### Validation against truth (simulation only)

In a simulation you know $\tau(x)$, so you can report the CATE RMSE and the
empirical coverage of the pointwise intervals directly. On real data these
are unavailable, which is exactly why the three model-free diagnostics above
exist. Always show that the honest forest's intervals cover near their
nominal rate and that turning honesty off both worsens RMSE and blows up the
variance estimates.

## The running application

A cross-sectional driver experiment (the shared HTE DGP): 6000 drivers, a
randomized push notification $w$, weekly completed trips $y$. The true effect
$\tau = 1.5\,\text{density} - 0.02\,\text{tenure} + 2.0\,\text{density}
\cdot\text{peak share}$ rises with city density, falls slowly with tenure, and
is amplified for peak-hour drivers in dense cities; `rating` is a pure
nuisance covariate. The average true effect is about 0.79 trips per week.

Fitting `causal_forest` with `W.hat = 0.5` and 2000 trees: the AIPW ATE is
0.81 (se 0.05), recovering the truth; the CATE RMSE is about 0.28; the
calibration test returns a mean coefficient near 1.0 and a differential
coefficient near 1.05 (both highly significant, so the heterogeneity is real
and correctly signed); the best linear projection recovers density (about
2.2), tenure (about -0.016), peak share (about 0.97), and a `rating`
coefficient indistinguishable from zero; the AUTOC RATE with the forest's own
$\hat\tau$ priority is 0.59 (se 0.05) but is a null -0.05 when the useless
`rating` is the priority. Refitting with `honesty = FALSE` doubles the CATE
RMSE (to about 0.60) and inflates the pointwise standard errors roughly
fourfold, so the intervals become uninformative even where they nominally
cover.

## Practitioner checklist

1. Set `W.hat` to the known propensity in an experiment; estimate it
   otherwise. Keep local centering on.
2. Fit with enough trees for stable variance (2000+); tune
   `min.node.size` / `sample.fraction` if the calibration test flags
   miscalibration.
3. Run `test_calibration` first: if the differential coefficient is not
   significant, there is no heterogeneity to model, and you should report the
   ATE and stop.
4. If heterogeneity exists, summarize it with `best_linear_projection` and
   quantify targetability with `rank_average_treatment_effect` (RATE).
5. Report pointwise CIs from `estimate.variance`, but never as a uniform
   band; for whole-function claims use the aggregate diagnostics.
6. In simulation, validate RMSE and coverage, and show the honesty
   ablation.

## The same problem at an online retailer

The randomized signup discount provides a clean experimental setting for estimating heterogeneous treatment effects. Customers are randomized to receive a discount offer or to see the standard price; the outcome is 12-month spend. The causal forest takes account tenure, pre-period spend, primary purchase category, and device type as covariates and estimates the conditional average treatment effect of the signup offer across the covariate space. Honest splitting prevents the winner's-curse inflation that arises when the same data selects the subgroup boundaries and evaluates the treatment effect within them: separate subsamples for tree construction and leaf-level estimation give asymptotically valid pointwise intervals. The forest reveals that the signup offer's effect on spend varies substantially: customers with high pre-period spend and long tenure respond little, because heavy members are likely to join regardless of the discount, while newer or lighter customers show larger effects at the margin. These conditional estimates feed directly into the policy-learning problem in Module 7, where the goal is to deploy a targeting rule that allocates the discount efficiently across the customer base.

## References

- Wager, S. and Athey, S. (2018). Estimation and inference of heterogeneous
  treatment effects using random forests. *Journal of the American
  Statistical Association*, 113(523), 1228-1242.
- Athey, S., Tibshirani, J. and Wager, S. (2019). Generalized random forests.
  *Annals of Statistics*, 47(2), 1148-1178.
- Nie, X. and Wager, S. (2021). Quasi-oracle estimation of heterogeneous
  treatment effects. *Biometrika*, 108(2), 299-319.
- Yadlowsky, S., Fleming, S., Shah, N., Brunskill, E. and Wager, S. (2021).
  Evaluating treatment prioritization rules via rank-weighted average
  treatment effects. *arXiv:2111.07966*.
- Robinson, P. (1988). Root-N-consistent semiparametric regression.
  *Econometrica*, 56(4), 931-954.
