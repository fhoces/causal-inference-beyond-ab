# Module 1: TWFE Diagnosed

Goodman-Bacon and the zoo of 2x2s. The experimentation-refresher's Module 8
showed *that* two-way fixed effects breaks under staggered adoption; this
module shows *exactly how*, with the full decomposition, the weight
formulas, and a hand-coded replication of `bacondecomp::bacon()`.

## The estimator and the question it dodges

The two-way fixed effects (TWFE) regression is

$$y_{it} = \alpha_i + \lambda_t + \beta^{DD} D_{it} + \varepsilon_{it}$$

with unit effects $\alpha_i$, period effects $\lambda_t$, and a single
treatment dummy $D_{it}$. In the canonical 2x2 (two groups, two periods,
one adoption date) this *is* the difference-in-differences estimator and,
under parallel trends, it is unbiased for the ATT.

With **staggered adoption** (units switch on at different dates) the
regression still returns one number. The question is which number. By
Frisch-Waugh-Lovell,

$$\hat\beta^{DD} = \frac{\widehat{Cov}(y_{it}, \tilde D_{it})}{\widehat{Var}(\tilde D_{it})}$$

where $\tilde D_{it} = D_{it} - \bar D_i - \bar D_t + \bar{\bar D}$ is the
double-demeaned treatment. Three facts follow directly:

1. TWFE is a fixed linear combination of the outcome cells, so an exact
   algebraic decomposition exists.
2. Always-treated (and never-treated) units have $\tilde D \approx 0$: they
   identify the fixed effects, not $\beta$.
3. The identifying variance from a cohort treated a share $\bar D_k$ of the
   panel scales with $\bar D_k(1 - \bar D_k)$: units treated near the middle
   of the panel dominate.

## The Goodman-Bacon decomposition

**Theorem (Goodman-Bacon 2021).** With timing groups $k = 1, \dots, K$
(adoption dates $g_k$) and possibly a never-treated group $U$, the TWFE
coefficient is exactly

$$\hat\beta^{DD} = \sum_{k \neq U} s_{kU}\, \hat\beta_{kU}
+ \sum_{k} \sum_{l > k} \left[ s_{kl}^{k}\, \hat\beta_{kl}^{k}
+ s_{kl}^{l}\, \hat\beta_{kl}^{l} \right]$$

a weighted average of every pairwise 2x2 DiD in the data, with non-negative
weights summing to one.

### The three species of 2x2

| 2x2 | Treatment | Control | Window |
|---|---|---|---|
| $\hat\beta_{kU}$ | cohort $k$ | never-treated | full panel |
| $\hat\beta_{kl}^{k}$ | earlier cohort $k$ | later cohort $l$, not yet treated | $t < g_l$ |
| $\hat\beta_{kl}^{l}$ | later cohort $l$ | earlier cohort $k$, already treated | $t \geq g_k$ |

The first two are legitimate DiDs: their controls are untreated throughout
the comparison window. The third, the **forbidden comparison**, uses
already-treated units as controls.

### The weights

Let $n_k$ be group sizes, $n_{kl} = n_k/(n_k + n_l)$, and $\bar D_k$ the
share of the panel that cohort $k$ spends treated. Then

$$s_{kU} \propto (n_k + n_U)^2\, n_{kU}(1-n_{kU})\, \bar D_k(1-\bar D_k)$$

$$s_{kl}^{k} \propto \left[(n_k+n_l)(1-\bar D_l)\right]^2 n_{kl}(1-n_{kl})\,
\frac{\bar D_k - \bar D_l}{1-\bar D_l}\cdot\frac{1-\bar D_k}{1-\bar D_l}$$

$$s_{kl}^{l} \propto \left[(n_k+n_l)\bar D_k\right]^2 n_{kl}(1-n_{kl})\,
\frac{\bar D_l}{\bar D_k}\cdot\frac{\bar D_k-\bar D_l}{\bar D_k}$$

normalized to sum to one. Each weight is a product of a subsample-size term,
a group-balance term, and a treatment-share-variance term. Two design
lessons: mid-panel adopters get the most weight, and the weights depend only
on sizes and calendar timing, never on effect magnitudes or precision. Two
rollouts with identical treatment effects but different calendars estimate
different TWFE targets.

### Why the forbidden comparison bites

For the pair earlier $k$, later $l$, window $t \geq g_k$, split the window at
$g_l$ into $W_1 = [g_k, g_l)$ and $W_2 = [g_l, T]$. Under parallel trends,

$$\hat\beta_{kl}^{l} \xrightarrow{p} \text{ATT}_l(W_2)
- \left[\text{ATT}_k(W_2) - \text{ATT}_k(W_1)\right].$$

The bracketed term is the control cohort's *effect drift* across the window
split. Constant effects kill it; dynamic effects (effects growing with time
since adoption) make it positive, so it is subtracted: attenuation first,
sign reversal when dynamics are strong enough.

### The estimand decomposition

Taking plims under parallel trends,

$$\text{plim}\;\hat\beta^{DD} = \text{VWATT} + \text{VWCT} - \Delta\text{ATT}$$

- **VWATT**: variance-weighted ATT, the "good" part, though note the
  aggregation weights are the Bacon $s$, not treated-cell shares.
- **VWCT**: variance-weighted common trends, zero when parallel trends holds.
- **$\Delta$ATT**: the accumulated within-cohort effect changes subtracted by
  forbidden comparisons; zero when effects are constant in event time.

Dynamics enter only through $\Delta$ATT and always with a minus sign: TWFE
bias from staggering is toward zero and beyond, never away from it. PT
violations (VWCT) and heterogeneity bias ($\Delta$ATT) are separate terms:
fixing one does not fix the other. Module 2's estimators fix $\Delta$ATT;
Module 3 stress-tests VWCT.

## The running application

Staggered rollout of a driver zone-notification feature: 30 cities, 60
weeks, adoption cohorts at weeks 15/25/35 (8 cities each) plus 6
never-treated cities. City fixed effects and a common weekly trend are built
in, so parallel trends holds by construction and every bias is the
estimator's fault. Three treatment-effect scenarios on one DGP:

| Scenario | Effect | True ATT | TWFE | Bias |
|---|---|---|---|---|
| constant | 1.0 flat | 1.000 | 0.997 | none |
| dynamic | $0.4[1 + 0.10(t-g)]$ | 1.137 | 0.686 | -40% |
| heterogeneous | dynamic $\times$ cohort scaling | 1.407 | 0.423 | -70% |

The decomposition on the dynamic scenario puts roughly half the weight on
vs-never comparisons (average estimate near 1.1), a fifth on legitimate
earlier-vs-later comparisons, and a third on forbidden comparisons whose
average estimate is 0.08. Same weights on the constant scenario, but every
2x2 estimates 1.0, so the weighting is harmless.

Sign reversal needs one more ingredient: **no never-treated anchor**. With
the 6 never-treated cities in place, even steep dynamics leave TWFE positive
because half the weight sits on clean vs-never comparisons. Drop them
(every city eventually adopts, the default at a company that ships
everywhere) and even mild dynamics of $0.2 + 0.15(t-g)$ produce a *negative*
TWFE coefficient against a true ATT near 3, with every cohort-time effect
strictly positive.

## When TWFE is fine

1. **Non-staggered designs.** One adoption date plus a clean control group
   is the 2x2 world; nothing here applies.
2. **Constant, homogeneous effects.** $\Delta\text{ATT} = 0$ and every 2x2
   estimates the same quantity; the weighting is irrelevant.
3. **Negligible forbidden weight.** With a large never-treated group and
   compressed adoption timing, the forbidden 2x2s may carry trivial weight;
   TWFE is then approximately VWATT, and you should say so rather than
   "the ATT".

## Practitioner checklist

1. Map the cohorts: adoption dates, sizes, never-treated share.
2. Run `bacondecomp::bacon()`; report weight by comparison type.
3. Plot weight vs estimate per 2x2; look for high-weight outliers.
4. Re-estimate with a heterogeneity-robust estimator (Module 2); report the
   gap.
5. If the design is not staggered, say so and skip the ceremony.

## The same problem at an online retailer

A large online retailer upgrades metros from two-day to next-day delivery as new fulfillment centers open, staggering the rollout over several quarters. The resulting dataset is a metro-week orders panel with multiple adoption cohorts and a slice of metros that never receive the upgrade. This is exactly the staggered-adoption setting the module analyzes: TWFE on the panel uses already-upgraded metros as controls for late-adopting metros during the late cohorts' post periods, the forbidden comparison. The mechanism is habit formation: customers update their ordering behavior as they learn the reliability of the faster delivery promise, so effects grow with time since adoption rather than jumping to a fixed level at the upgrade date. Dynamic effects are precisely the condition under which the forbidden comparison biases TWFE. Early cohorts' growing effect gets subtracted from late cohorts' estimates, attenuating TWFE toward zero and potentially reversing its sign if the never-treated slice is small. Preserving a set of metros that never receive next-day delivery gives every estimator in this course a clean comparison group.

## References

- Goodman-Bacon, A. (2021). Difference-in-differences with variation in
  treatment timing. *Journal of Econometrics*, 225(2), 254-277.
- de Chaisemartin, C. and D'Haultfoeuille, X. (2020). Two-way fixed effects
  estimators with heterogeneous treatment effects. *American Economic
  Review*, 110(9), 2964-2996.
- Borusyak, K., Jaravel, X. and Spiess, J. (2024). Revisiting event-study
  designs: robust and efficient estimation. *Review of Economic Studies*.
- Baker, A., Larcker, D. and Wang, C. (2022). How much should we trust
  staggered difference-in-differences estimates? *Journal of Financial
  Economics*, 144(2), 370-395.
- Roth, J., Sant'Anna, P., Bilinski, A. and Poe, J. (2023). What's trending
  in difference-in-differences? A synthesis of the recent econometrics
  literature. *Journal of Econometrics*, 235(2), 2218-2244.
