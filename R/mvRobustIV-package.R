#' mvRobustIV: Uniform Inference for Multiple Causal Effects with Invalid Instruments
#'
#' Confidence regions for a vector of causal effects in multivariable
#' instrumental-variable and Mendelian randomization analyses when some
#' instruments may be invalid and the valid set is unknown. The package
#' implements the searching confidence region obtained by inverting scalar
#' residual-validity tests under a majority rule, without selecting a valid
#' set, together with first-stage relevance screening, boundedness and strength
#' diagnostics, an iterated hard-thresholding valid-set estimate, and an
#' optional resampling refinement.
#'
#' The central observation is that, with a scalar outcome, each instrument's
#' direct effect on the outcome is a single number irrespective of the number
#' of exposures \eqn{q}. For any candidate causal vector \eqn{\beta}, every
#' instrument induces a scalar residual contrast
#' \eqn{\hat\Gamma_j - \hat\gamma_j^\top \beta} that is asymptotically normal
#' regardless of \eqn{q}. Inverting validity tests on these contrasts under a
#' majority rule yields a region with coverage at least \eqn{1-\alpha}, uniform
#' over a class of data-generating processes that includes locally invalid
#' configurations.
#'
#' Main entry points: [mr_input()] builds a summary-data object,
#' [searching_region()] computes the region and its coordinate sections,
#' [iv_diagnostics()] reports relevance, boundedness and strength diagnostics,
#' and [simulate_mvmr()] generates example data.
#'
#' @references Chen, S. (2026). Uniform inference for multiple causal effects
#'   with invalid instruments.
#'
#' @keywords internal
"_PACKAGE"

## ---- internal helpers -------------------------------------------------------

`%||%` <- function(a, b) if (is.null(a)) b else a

## Per-comparison Bonferroni normal cutoff used throughout.
.zalpha <- function(alpha, m) {
  stats::qnorm(1 - alpha / (2 * m))
}

## Standardized residual contrasts and their thresholds at a candidate beta.
##  pi_j(beta)        = Gamma_j - gamma_j^T beta
##  R_jj(beta)/n      = se_Y_j^2 + sum_k beta_k^2 se_X_jk^2   (eq. Rscale)
.residuals <- function(object, beta) {
  beta <- as.numeric(beta)
  pi <- object$beta_Y - as.numeric(object$beta_X %*% beta)
  Rn <- object$se_Y^2 + as.numeric((object$se_X^2) %*% (beta^2))
  list(pi = pi, Rn = Rn)
}
