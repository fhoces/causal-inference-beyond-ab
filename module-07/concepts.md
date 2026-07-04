# Module 7: Policy Learning

From HTE estimates to deployment rules. Module 6 gave us $\hat\tau(x)$, a
calibrated map from driver features to conditional treatment effects. That
answers "to whom is the push notification most effective?" It does *not* tell
us what to actually do. This module closes the loop: turn heterogeneous
effects into a deployable rule $\pi: X \to \{0, 1\}$, learn it with the
regret guarantees of Athey and Wager (2021), and evaluate a candidate policy
honestly with off-policy value estimates.

## The policy learning problem

A policy is a map $\pi: X \to \{0, 1\}$ assigning a treatment decision to
every covariate profile. Its **value** is the expected outcome if we deploy
it to the population:

$$V(\pi) = E\big[ Y(\pi(X)) \big] = E\big[ Y(0) + \pi(X)\,\tau(X) \big].$$

The second equality uses $Y(\pi(X)) = Y(0) + \pi(X)(Y(1) - Y(0))$ and the
tower property. We want $\pi^\star = \arg\max_{\pi \in \Pi} V(\pi)$ over a
restricted class $\Pi$ (shallow trees, linear-index rules, monotone rules).
The unconstrained optimum is the pointwise rule
$\pi^\star_{\text{full}}(x) = \mathbb{1}\{\tau(x) > 0\}$, or with a per-unit
cost $c$, $\pi^\star_{\text{full}}(x) = \mathbb{1}\{\tau(x) > c\}$. Restricting
to $\Pi$ trades a little value for interpretability, fairness, and
implementability.

### Welfare maximization is not welfare estimation

This is the single most important idea in the module. Estimating $\tau(x)$
well (the Module 6 objective) and choosing a good policy are *different*
objectives.

- A CATE estimator is judged by its error in **levels**: mean-squared error
  $E[(\hat\tau(x) - \tau(x))^2]$, pointwise coverage.
- A policy is judged only by whether it gets the **sign of $\tau(x) - c$
  right** where it matters, weighted by $|\tau(x) - c|$. A rule needs the
  *ranking* around the cost threshold, not unbiased levels.

Consequences:

1. A biased $\hat\tau$ that preserves the ordering of units relative to $c$
   yields the optimal policy. An unbiased-but-noisy $\hat\tau$ can yield a
   worse policy.
2. Regions where $\tau(x)$ is far from $c$ are easy: any reasonable estimate
   gets the decision right. All the action is near the decision boundary
   $\tau(x) = c$, exactly where CATE estimation is hardest and least
   consequential for MSE.
3. You should optimize the welfare objective directly, not a plug-in of a
   CATE fit. That is empirical welfare maximization.

## Doubly-robust scores

To optimize welfare from a finite RCT (or observational data) we need an
unbiased, low-variance estimate of $V(\pi)$ for any $\pi$. The
augmented-inverse-propensity (AIPW), or doubly-robust, score does this. For
each unit build a per-arm score that is an unbiased estimate of the potential
outcome under that arm:

$$\Gamma_i(w) = \hat\mu_w(X_i) + \frac{\mathbb{1}\{W_i = w\}}{\hat e_w(X_i)}
\big( Y_i - \hat\mu_w(X_i) \big), \qquad w \in \{0, 1\},$$

where $\hat\mu_w(x) = \hat E[Y \mid X = x, W = w]$ is an outcome model and
$\hat e_w(x) = \hat P(W = w \mid X = x)$ is the propensity. In our RCT the
propensity is known: $\hat e_1 = 0.5$. The doubly-robust name: $\Gamma_i(w)$
is unbiased for $E[Y(w) \mid X_i]$ if *either* $\hat\mu_w$ *or* $\hat e_w$ is
correct. `grf` builds these from a causal forest: with the forest's marginal
outcome model $\hat m(x)$ and effect $\hat\tau(x)$,
$\hat\mu_1 = \hat m + (1 - \hat e)\hat\tau$ and
$\hat\mu_0 = \hat m - \hat e\hat\tau$;
`policytree::double_robust_scores(forest)` returns the $n \times 2$ matrix
with columns `control` and `treated`.

### Why plug-in rules are inferior

The tempting shortcut is the plug-in rule
$\hat\pi(x) = \mathbb{1}\{\hat\tau(x) > c\}$. It is dominated by DR-score
learning for two reasons:

1. **Propensity / regularization bias.** $\hat\tau$ from a forest is a
   regularized, biased estimate. The plug-in inherits that bias directly at
   the decision boundary. The DR score corrects the outcome model with the
   residual term, so first-order errors in $\hat\mu$ do not propagate into the
   value estimate.
2. **It optimizes the wrong loss.** The forest split criterion targets CATE
   heterogeneity (MSE-like), not welfare. Learning the policy on the DR
   scores optimizes value directly, and inherits Neyman-orthogonality: the
   welfare objective is insensitive to first-order nuisance error.

## Empirical welfare maximization and regret

Athey and Wager (2021) define the estimated policy as the empirical welfare
maximizer over the DR scores:

$$\hat\pi = \arg\max_{\pi \in \Pi} \frac{1}{n} \sum_{i=1}^n
\Big[ \Gamma_i(1)\,\pi(X_i) + \Gamma_i(0)\,(1 - \pi(X_i)) \Big].$$

Performance is measured by **regret** against the best rule in the class:

$$R(\hat\pi) = V(\pi^\star_\Pi) - V(\hat\pi), \qquad
\pi^\star_\Pi = \arg\max_{\pi \in \Pi} V(\pi).$$

The main theorem: with doubly-robust scores and a class $\Pi$ of bounded
Vapnik-Chervonenkis (VC) dimension,

$$R(\hat\pi) = O_p\!\left( \sqrt{\frac{\text{VC}(\Pi)}{n}} \right).$$

Two features of this bound deserve emphasis.

- **Restricting $\Pi$ is a feature, not a bug.** The regret is against the
  best rule *in the class*, and the rate improves as $\text{VC}(\Pi)$ shrinks.
  A depth-2 tree has small VC dimension, so it converges fast; it is also
  interpretable, auditable, and implementable in a rules engine. You choose
  $\Pi$ for deployment constraints and get a tighter bound for free.
- **Utilitarian, not pointwise.** The guarantee is on average welfare, not on
  getting every unit's decision right. Pointwise-optimal rules are impossible
  from finite data (near $\tau(x) = c$ you cannot resolve the sign), but the
  units you misclassify there are exactly the ones where $|\tau(x) - c|$ is
  small, so they cost little welfare. Regret is the right yardstick precisely
  because it downweights the unavoidable errors.

## policytree in practice

`policytree::policy_tree(X, Gamma, depth = 2)` performs an **exact** search
over all depth-$L$ axis-aligned trees, maximizing the summed DR reward. This
is not greedy CART: it globally optimizes the welfare objective. The cost:
exact search is combinatorial in the number of covariates and split points,
so runtime grows fast with depth. Depth 2 is the practical default (four
leaves, three splits); depth 3 is often too slow on wide data and rarely buys
enough welfare to justify the loss of interpretability. `predict(tree, X)`
returns actions in $\{1, 2\}$ where 1 = control and 2 = treat.

Reading a tree: each internal node is a covariate threshold, each leaf a
treat / control decision. Because the tree is fit on the *welfare* objective,
its splits identify the covariate regions where treating clears the cost, not
where the CATE is merely large.

### Cost-aware policies

A per-push cost $c$ enters cleanly. The value of treating unit $i$ net of cost
is $\Gamma_i(1) - c$, so subtract $c$ from the treated column of $\Gamma$
before fitting:

$$\Gamma_i^{\text{net}} = \big( \Gamma_i(0),\; \Gamma_i(1) - c \big),
\qquad \pi(x) = \mathbb{1}\{\tau(x) > c\}.$$

Raising $c$ makes the treat leaves shrink: the tree only treats regions whose
effect clears the higher bar. In the application, $c = 0.6$ (about 0.6 weekly
trips of margin per push).

## Off-policy evaluation

Given a candidate policy $\pi$ (learned however), estimate its value from RCT
data with the doubly-robust value estimator:

$$\hat V(\pi) = \frac{1}{n} \sum_{i=1}^n
\Big[ \Gamma_i(0) + \pi(X_i)\big(\Gamma_i(1) - \Gamma_i(0)\big) \Big],$$

net of cost by using $\Gamma_i(1) - c$ in place of $\Gamma_i(1)$. This is an
average of i.i.d. terms, so the standard error follows from the
influence-function form: $\widehat{\text{SE}} = \hat\sigma_\psi / \sqrt{n}$
where $\psi_i$ is the summand and $\hat\sigma_\psi$ its sample standard
deviation.

### The naive trap and cross-fitting

If you evaluate a *learned* policy on the same data that trained it, $\hat
V(\hat\pi)$ is optimistically biased: the policy has adapted to the noise in
those scores. The fix is **cross-fitting** (or an honest split): the policy
must never be evaluated on the data that trained it.

$K$-fold cross-fit:

1. Split into $K$ folds.
2. For each fold $k$: fit the nuisances (forest, outcome model) on the other
   $K - 1$ folds, learn $\hat\pi_{-k}$ on those folds, then compute DR scores
   on fold $k$ from the *training-fold* nuisances and evaluate
   $\hat\pi_{-k}$ on fold $k$.
3. Pool the held-out contributions; the SE comes from the pooled influence
   terms.

Because each held-out unit is scored by a policy and nuisances that never saw
it, the estimate is honest.

## The application

Zone-notification push to 6000 drivers, randomized $W$ with propensity 0.5,
outcome = weekly completed trips, cost $c = 0.6$ per push. Effect
$\tau(x) = 1.5\,\text{density} - 0.02\,\text{tenure} +
2.0\,\text{density}\times\text{peak\_shr}$; `rating` is a pure nuisance
covariate. The true optimal share treated ($\tau > 0.6$) is about 57%.

Learned depth-2 tree (net of cost) splits on tenure then density: newer
drivers (tenure below about 28 months) are treated once density exceeds about
0.27; veteran drivers are treated only in the densest markets (density above
about 0.73). It treats roughly 53% of drivers. The interpretation matches the
DGP: tenure lowers $\tau$, so veterans need a denser market to clear the cost.

Honest off-policy values, as gains over treat-none, net of cost:

| Policy | Gain over treat-none | Note |
|---|---|---|
| treat all | 0.26 | ignores heterogeneity, pays cost everywhere |
| oracle ($\tau > c$) | 0.40 | infeasible upper bound |
| forest plug-in | 0.40 | $\hat\tau > c$, unrestricted |
| depth-2 policy tree | 0.41 | interpretable, three splits |

The tree is within noise of the oracle and beats treat-all by about 0.16 net
of cost. The naive same-data tree value overstates the gain (about 0.45
versus the honest 0.39 to 0.41): that gap is the optimism a cross-fit
removes.

## When this breaks

1. **Bad nuisances.** DR corrects first-order error, not gross
   misspecification. In observational data with poor overlap the propensity
   blows up the score variance and the value estimate becomes useless. Trim or
   check overlap first.
2. **Evaluating on training data.** The naive trap: always cross-fit or split.
3. **Over-rich $\Pi$.** A deep tree has high VC dimension, slow regret, and no
   interpretability. Depth 2 is usually the right default.
4. **Value near the boundary.** If most of the population sits near
   $\tau(x) = c$, no policy separates much value from treat-all, and the honest
   gain is small with wide SEs. Report the SE, not just the point estimate.

## The same problem at an online retailer

The membership experiment's heterogeneous treatment effects become actionable through policy learning. Given a fixed promotional budget, the retailer must decide which customers to offer the signup discount. The goal is to maximize incremental membership revenue: the policy should direct offers toward customers whose spending would increase most because of membership, not toward customers who would join anyway at full price or who would not benefit even as members. The `policy_tree` function approximates this welfare-maximizing rule as a shallow decision tree, splitting on covariates such as tenure and pre-period spend to produce a rule that can be deployed without per-customer score lookup. Naive targeting on in-sample predicted uplift exploits noise: the customers with the highest $\hat\tau(x)$ in the training data often reflect overfitting, and their realized uplift at deployment is smaller. The doubly-robust AIPW score corrects for this in-sample optimism and provides a reliable estimate of the targeting policy's true lift, evaluated on a held-out test set. A fixed promo budget introduces a cost constraint that the policy must respect: the policy value calculation determines how many offers can be sent profitably and which customer segments should receive priority.

## References

- Athey, S. and Wager, S. (2021). Policy learning with observational data.
  *Econometrica*, 89(1), 133-161.
- Kitagawa, T. and Tetenov, A. (2018). Who should be treated? Empirical
  welfare maximization methods for treatment choice. *Econometrica*, 86(2),
  591-616.
- Zhou, Z., Athey, S. and Wager, S. (2023). Offline multi-action policy
  learning: generalization and optimization. *Operations Research*, 71(1),
  148-183.
- Wager, S. and Athey, S. (2018). Estimation and inference of heterogeneous
  treatment effects using random forests. *JASA*, 113(523), 1228-1242.
- Athey, S., Tibshirani, J. and Wager, S. (2019). Generalized random forests.
  *Annals of Statistics*, 47(2), 1148-1178.
