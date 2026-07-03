# Causal Inference Beyond A/B Tests

A deep-dive companion to the [experimentation refresher](https://github.com/fhoces/experimentation-refresher),
specifically expanding the three methods covered as a tour in that course's
Module 8: modern difference-in-differences, synthetic control, and causal
forests for heterogeneous treatment effects.

**Live site:** [fhoces.github.io/causal-inference-beyond-ab](https://fhoces.github.io/causal-inference-beyond-ab/) (rendered slides and notes).

## Status

**Modules 1-8 built** (concepts + slides + exercises, rendered decks
committed). Course complete.

## Why this exists

The experimentation-refresher covers each of these methods in a single slide
deck: enough to spot the failure modes of TWFE, recognize when synthetic
control applies, and read a `grf` output in an interview. That's a tour, not
a treatment.

This course is the treatment. Each method gets a multi-module sequence with
the formal estimator definitions, the simulation studies that motivate the
modern alternatives, and the practitioner choices that don't fit on a slide.

## Modules

See [learning-plan.md](learning-plan.md) for the full breakdown.

| # | Module | Status |
|---|--------|--------|
| 1 | [TWFE Diagnosed: Goodman-Bacon and the Zoo of 2×2s](module-01/concepts.md) | done |
| 2 | [Heterogeneity-Robust DiD: CS, SA, BJS, dCDH in Detail](module-02/concepts.md) | done |
| 3 | [Honest DiD: Sensitivity Bounds for Parallel Trends](module-03/concepts.md) | done |
| 4 | [Synthetic Control: Estimator, Inference, Variants](module-04/concepts.md) | done |
| 5 | [Synthetic DiD and the Bridge from SC to DiD](module-05/concepts.md) | done |
| 6 | [Causal Forest: Honest Splitting and Asymptotics](module-06/concepts.md) | done |
| 7 | [Policy Learning: From τ̂(x) to Deployment Rules](module-07/concepts.md) | done |
| 8 | [Matrix Completion and the Modern Panel Toolbox](module-08/concepts.md) | done |

Module 8 closes with the course capstone: a method-selection decision tree mapping problem shape to estimator, the single most interview-useful artifact in the course.

## Domain

Following the convention of the sibling course collection, applications use
ride-sharing (Uber / Lyft style) data. Where formal methods need a specific
empirical setting (e.g., the original Card-Krueger DiD or the Abadie-Diamond
California Prop 99 SC), the canonical paper's data is used. The methods
transfer directly to the standard platform settings: staggered feature
rollouts, geo-level policy changes, fee or pricing changes in a single
market, and targeted incentives.

## Structure

Each module follows the standard concept → show → drill pattern from the
sibling courses: `concepts.md` for the written reference, `slides.Rmd` for
the xaringan deck, and `exercise.R` for runnable drills.

## Running the exercises

Each module's `exercise.R` is self-contained: run it directly with
`Rscript module-0N/exercise.R`. The R packages needed across the course are
`fixest`, `did`, `didimputation`, `HonestDiD`, `Synth`, `gsynth`, `synthdid`,
`augsynth`, `grf`, `quadprog`, `bacondecomp`, `policytree`, plus tidyverse.
Every script ends with a block of assertions, exiting 0 when all checks
pass.

## Stack

R, with `did`, `fixest`, `synthdid`, `Synth`, `gsynth`, `grf`, `HonestDiD`,
and the standard tidyverse.
