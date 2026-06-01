#' Resampling refinement of the searching region
#'
#' A coverage-oriented sensitivity refinement that parametrically resamples the
#' summary statistics and records, for each coordinate, the envelope of the
#' resampled searching-region centers. It is secondary to [searching_region()],
#' which carries the uniform guarantee; the resampling region is reported only
#' as a calibration check and is typically wider.
#'
#' @inheritParams searching_membership
#' @param M Number of resamples.
#' @param level Envelope level; the interval spans the central `level` fraction
#'   of resampled centers per coordinate.
#' @param seed Optional random seed.
#'
#' @return An object of class `"sampling_region"` with `intervals`
#'   (a \eqn{q \times 2} matrix), `centers` (the resampled point estimates),
#'   `M` and `level`.
#'
#' @examples
#' sim <- simulate_mvmr(seed = 1)
#' sampling_region(sim, M = 100, seed = 2)$intervals
#'
#' @export
sampling_region <- function(object, alpha = 0.05, S = NULL, M = 500,
                            level = 0.95, seed = NULL) {
  stopifnot(inherits(object, "mr_input"))
  if (!is.null(seed)) set.seed(seed)
  S <- S %||% relevance_screen(object)
  q <- object$q
  centers <- matrix(NA_real_, M, q)
  for (m in seq_len(M)) {
    rs <- object
    rs$beta_Y <- object$beta_Y + stats::rnorm(object$p) * object$se_Y
    rs$beta_X <- object$beta_X + matrix(stats::rnorm(object$p * q), object$p, q) * object$se_X
    centers[m, ] <- tsht_valid_set(rs, alpha, S)$estimate
  }
  probs <- c((1 - level) / 2, 1 - (1 - level) / 2)
  intervals <- t(apply(centers, 2, stats::quantile, probs = probs, names = FALSE))
  rownames(intervals) <- object$exposure
  colnames(intervals) <- c("lower", "upper")
  structure(list(intervals = intervals, centers = centers, M = M,
                 level = level, exposure = object$exposure),
            class = "sampling_region")
}

#' @export
print.sampling_region <- function(x, ...) {
  cat(sprintf("<sampling_region> resampling refinement (M = %d, level = %.2f)\n",
              x$M, x$level))
  cat("  secondary to searching_region(); reported as a sensitivity check.\n")
  print(round(x$intervals, 4))
  invisible(x)
}
