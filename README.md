# mvRobustIV

<!-- badges: start -->
<!-- badges: end -->

**Uniform inference for multiple causal effects with invalid instruments.**

`mvRobustIV` builds confidence regions for a vector of causal effects in
multivariable instrumental-variable and Mendelian randomization (MR) analyses
when some instruments may be invalid and the valid set is unknown. The
**searching confidence region** is obtained by inverting scalar
residual-validity tests under a majority rule, without selecting a valid set,
and has coverage at least `1 - alpha` uniformly over a class of data-generating
processes that includes locally invalid configurations.

The key idea: with a scalar outcome, each instrument's direct effect on the
outcome is a single number regardless of the number of exposures `q`. For any
candidate causal vector `beta`, every instrument induces a scalar residual
contrast `Gamma_j - gamma_j' beta` that is asymptotically normal regardless of
`q`. Tests on these contrasts can be inverted into a region for `beta`.

## Installation

```r
# install.packages("remotes")
remotes::install_github("siyichen/mvRobustIV")
```

The package depends only on base R and the recommended packages `stats`,
`graphics` and `grDevices`.

## Quick start

```r
library(mvRobustIV)

# simulate two-sample summary data (q = 2 exposures, 30 instruments, 21 valid)
sim <- simulate_mvmr(q = 2, p = 30, n_valid = 21, seed = 1)

# with your own data:
#   obj <- mr_input(beta_Y, se_Y, beta_X, se_X, n_X = ..., n_Y = ...)

# first-stage diagnostics: strength, boundedness, surplus
iv_diagnostics(sim)

# the searching confidence region
reg <- searching_region(sim)
reg
plot(reg)            # q = 2 only

# does a candidate vector lie in the region?
searching_membership(sim, c(0.3, -0.2))$in_region
```

## What the region tells you

A coordinate is reported as **bounded** only when a majority of the selected
instruments are *individually* strongly associated with that exposure
(`boundedness_rule()`). When fewer than a majority are, the region is unbounded
along that direction: under unknown invalidity the data do not identify that
coordinate. This is an inferential finding, not a numerical failure, and it is
distinct from the conditional `F`, collinearity and surplus diagnostics, which
can all look healthy while a coordinate is still unbounded.

## Main functions

| Function | Purpose |
|---|---|
| `mr_input()` | Build and validate a summary-data object |
| `relevance_screen()` | Select the instrument set by the first-stage rule |
| `searching_region()` | Searching confidence region and its coordinate sections |
| `searching_membership()` | Test whether a `beta` lies in the region |
| `searching_section()` / `searching_projection()` | One coordinate's section / projection |
| `boundedness_rule()` | Per-coordinate boundedness from marginal relevance |
| `iv_diagnostics()` | Relevance, boundedness, strength and surplus diagnostics |
| `mv_ivw()` | Conventional inverse-variance-weighted estimate (comparator) |
| `tsht_valid_set()` | Iterated hard-thresholding valid-set estimate (comparator) |
| `sampling_region()` | Resampling refinement (secondary sensitivity check) |
| `simulate_mvmr()` | Generate example data |

## Scope and caveats

* The searching region is the procedure that carries the uniform guarantee.
  `sampling_region()` is a calibration-oriented sensitivity check and is
  typically wider.
* `searching_projection()` approximates the exact projection on a grid and is
  intended for small `q` (2 or 3).
* The conditional `F` in `iv_diagnostics()` is an approximate strength summary;
  for the formal Sanderson-Windmeijer conditional `F` use a dedicated MVMR
  package.
* The method assumes a scalar, linear outcome model, strong instruments, and a
  majority of valid instruments.

## Citation

```r
citation("mvRobustIV")
```

Chen, S. (2026). *Uniform inference for multiple causal effects with invalid
instruments.*

## License

MIT (c) 2026 Siyi Chen.
