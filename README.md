# Causal Inference Beyond A/B Tests

A deep-dive companion to the [experimentation refresher](../experimentation-refresher/),
specifically expanding the three methods covered as a tour in that course's
Module 8: modern difference-in-differences, synthetic control, and causal
forests for heterogeneous treatment effects.

## Status

**In development.** Module structure and learning plan are sketched; modules
are not yet built.

## Why this exists

The experimentation-refresher covers each of these methods in a single slide
deck — enough to spot the failure modes of TWFE, recognize when synthetic
control applies, and read a `grf` output in an interview. That's a tour, not
a treatment.

This course is the treatment. Each method gets a multi-module sequence with
the formal estimator definitions, the simulation studies that motivate the
modern alternatives, and the practitioner choices that don't fit on a slide.

## Modules (planned)

See [learning-plan.md](learning-plan.md) for the full breakdown.

| # | Module | Status |
|---|--------|--------|
| 1 | TWFE Diagnosed: Goodman-Bacon and the Zoo of 2×2s | upcoming |
| 2 | Heterogeneity-Robust DiD: CS, SA, BJS, dCDH in Detail | upcoming |
| 3 | Honest DiD: Sensitivity Bounds for Parallel Trends | upcoming |
| 4 | Synthetic Control: Estimator, Inference, Variants | upcoming |
| 5 | Synthetic DiD and the Bridge from SC to DiD | upcoming |
| 6 | Causal Forest: Honest Splitting and Asymptotics | upcoming |
| 7 | Policy Learning: From τ̂(x) to Deployment Rules | upcoming |
| 8 | Matrix Completion and the Modern Panel Toolbox | upcoming |

## Domain

Following the convention of the sibling course collection, applications use
ride-sharing (Uber / Lyft style) data. Where formal methods need a specific
empirical setting (e.g., the original Card-Krueger DiD or the Abadie-Diamond
California Prop 99 SC), the canonical paper's data is used.

## Structure

Each module follows the standard concept → show → drill pattern from the
sibling courses: `concepts.md` for the written reference, `slides.Rmd` for
the xaringan deck, and `exercise.R` for runnable drills.

## Stack

R, with `did`, `fixest`, `synthdid`, `Synth`, `gsynth`, `grf`, `HonestDiD`,
and the standard tidyverse.
