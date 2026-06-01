#' Membership test for the searching confidence region
#'
#' Tests whether a candidate causal vector \eqn{\beta} lies in the searching
#' region. At \eqn{\beta} each instrument in \eqn{\hat{\mathcal S}} is flagged
#' when its standardized residual exceeds the Bonferroni cutoff,
#' \eqn{|\hat\pi_j(\beta)| \ge z_\alpha \{R_{jj}(\beta)/n\}^{1/2}}, and
#' \eqn{\beta} belongs to the region when the number flagged is at most
#' \eqn{\lfloor |\hat{\mathcal S}|/2 \rfloor - 1}.
#'
#' @param object An [mr_input()] object.
#' @param beta Numeric vector of length \eqn{q}.
#' @param alpha Nominal level; the region has coverage at least \eqn{1-\alpha}.
#' @param S Optional selected set from [relevance_screen()]; computed if `NULL`.
#'
#' @return A list with `in_region` (logical), `n_flagged`, `max_flagged`
#'   (the threshold \eqn{\lfloor |\hat{\mathcal S}|/2 \rfloor - 1}) and the
#'   logical vector `flagged` over the selected set.
#'
#' @examples
#' sim <- simulate_mvmr(seed = 1)
#' searching_membership(sim, attr(sim, "true_beta"))$in_region
#'
#' @export
searching_membership <- function(object, beta, alpha = 0.05, S = NULL) {
  stopifnot(inherits(object, "mr_input"))
  S <- S %||% relevance_screen(object)
  if (length(S) == 0) stop("No instruments pass the relevance screen.", call. = FALSE)
  z <- .zalpha(alpha, length(S))
  r <- .residuals(object, beta)
  flagged <- abs(r$pi[S]) >= z * sqrt(r$Rn[S])
  n_flagged <- sum(flagged)
  max_flagged <- floor(length(S) / 2) - 1L
  list(in_region = n_flagged <= max_flagged,
       n_flagged = n_flagged,
       max_flagged = max_flagged,
       flagged = flagged)
}

#' Per-coordinate boundedness rule
#'
#' A coordinate's bisection section is bounded if and only if a majority of the
#' selected instruments are individually strongly associated with that
#' exposure, that is at least \eqn{\lfloor |\hat{\mathcal S}|/2 \rfloor} of them
#' satisfy \eqn{|\hat\gamma_{jk}|/\mathrm{se}_{\gamma,jk} > z_\alpha}. This rule
#' is computable before any region is formed and predicts which coordinates a
#' uniform region can bound.
#'
#' @inheritParams searching_membership
#'
#' @return A data frame with one row per exposure: `exposure`, `strong` (count
#'   of individually strong instruments), `threshold`
#'   (\eqn{\lfloor |\hat{\mathcal S}|/2 \rfloor}) and `bounded` (logical).
#'
#' @examples
#' boundedness_rule(simulate_mvmr(seed = 1))
#'
#' @export
boundedness_rule <- function(object, alpha = 0.05, S = NULL) {
  stopifnot(inherits(object, "mr_input"))
  S <- S %||% relevance_screen(object)
  z <- .zalpha(alpha, length(S))
  zmat <- abs(object$beta_X[S, , drop = FALSE]) / object$se_X[S, , drop = FALSE]
  strong <- colSums(zmat > z)
  threshold <- floor(length(S) / 2)
  data.frame(exposure = object$exposure,
             strong = as.integer(strong),
             threshold = as.integer(threshold),
             bounded = strong >= threshold,
             row.names = NULL, stringsAsFactors = FALSE)
}

## Internal: one-sided endpoint of a coordinate section by expand-then-bisect.
.section_endpoint <- function(f_in, x0, dir, max_abs, tol) {
  step <- max(abs(x0) * 0.5, 0.5)
  t_in <- x0
  repeat {
    t_try <- t_in + dir * step
    if (abs(t_try) > max_abs) return(dir * Inf)
    if (!f_in(t_try)) { t_out <- t_try; break }
    t_in <- t_try
    step <- step * 2
  }
  while (abs(t_out - t_in) > tol) {
    mid <- (t_in + t_out) / 2
    if (f_in(mid)) t_in <- mid else t_out <- mid
  }
  t_in
}

#' Coordinate section of the searching region
#'
#' Computes the one-dimensional slice of the searching region along exposure
#' `k`, holding the other coordinates fixed at `center`, by bisection. When the
#' coordinate is unbounded according to [boundedness_rule()] the corresponding
#' endpoints are returned as infinite.
#'
#' @inheritParams searching_membership
#' @param k Index of the exposure whose section is computed.
#' @param center Feasible reference point (length \eqn{q}) at which the other
#'   coordinates are held. Defaults to the [mv_ivw()] estimate.
#' @param max_abs Outer bracket for the bisection on each side.
#' @param tol Bisection tolerance.
#'
#' @return Numeric vector `c(lower, upper)`.
#'
#' @examples
#' sim <- simulate_mvmr(seed = 1)
#' searching_section(sim, k = 1)
#'
#' @export
searching_section <- function(object, k, alpha = 0.05, center = NULL, S = NULL,
                              max_abs = 50, tol = 1e-3) {
  stopifnot(inherits(object, "mr_input"))
  S <- S %||% relevance_screen(object)
  center <- center %||% .feasible_center(object, alpha, S)
  bnd <- boundedness_rule(object, alpha, S)$bounded[k]
  if (!isTRUE(bnd)) return(c(-Inf, Inf))

  f_in <- function(t) {
    b <- center
    b[k] <- t
    searching_membership(object, b, alpha, S)$in_region
  }
  if (!f_in(center[k])) {
    warning("Reference center is not in the region along coordinate ", k,
            "; returning NA section.", call. = FALSE)
    return(c(NA_real_, NA_real_))
  }
  lo <- .section_endpoint(f_in, center[k], -1, max_abs, tol)
  hi <- .section_endpoint(f_in, center[k], +1, max_abs, tol)
  c(lo, hi)
}

## Internal: find a feasible center, defaulting to mv_ivw with a small fallback.
.feasible_center <- function(object, alpha, S) {
  b0 <- mv_ivw(object, S)$estimate
  if (searching_membership(object, b0, alpha, S)$in_region) return(b0)
  ## coarse random search around the IVW point
  sds <- pmax(mv_ivw(object, S)$se, 1e-2)
  for (i in seq_len(500)) {
    cand <- b0 + stats::rnorm(object$q) * sds * (1 + i / 50)
    if (searching_membership(object, cand, alpha, S)$in_region) return(cand)
  }
  warning("No feasible center found; sections may be NA.", call. = FALSE)
  b0
}

#' Searching confidence region
#'
#' Computes the searching confidence region for the causal vector: the
#' per-coordinate bisection sections, the boundedness diagnostic, and the
#' reference center. The region carries coverage at least \eqn{1-\alpha},
#' uniform over the class of data-generating processes described in the
#' reference.
#'
#' @inheritParams searching_section
#'
#' @return An object of class `"searching_region"` with elements `sections`
#'   (a \eqn{q \times 2} matrix of lower and upper endpoints), `bounded`
#'   (logical per coordinate), `center`, `alpha`, `S` and `exposure`.
#'
#' @examples
#' sim <- simulate_mvmr(seed = 1)
#' searching_region(sim)
#'
#' @export
searching_region <- function(object, alpha = 0.05, center = NULL, S = NULL,
                             max_abs = 50, tol = 1e-3) {
  stopifnot(inherits(object, "mr_input"))
  S <- S %||% relevance_screen(object)
  center <- center %||% .feasible_center(object, alpha, S)
  bnd <- boundedness_rule(object, alpha, S)
  sections <- t(vapply(seq_len(object$q), function(k) {
    searching_section(object, k, alpha, center, S, max_abs, tol)
  }, numeric(2)))
  rownames(sections) <- object$exposure
  colnames(sections) <- c("lower", "upper")
  structure(
    list(sections = sections, bounded = bnd$bounded, center = center,
         alpha = alpha, S = S, exposure = object$exposure, q = object$q),
    class = "searching_region"
  )
}

#' Approximate coordinate projection of the searching region
#'
#' Returns the exact coordinate projection (the range of \eqn{\beta_k} over the
#' whole region) approximated on a grid. The projection is wider than the
#' section and can be unbounded under correlated exposures; it is feasible only
#' for small \eqn{q}.
#'
#' @inheritParams searching_section
#' @param grid Number of grid points per non-target coordinate.
#' @param span Half-width of the grid (in the units of `beta`) around the
#'   center for the non-target coordinates.
#'
#' @return Numeric vector `c(lower, upper)` for exposure `k`; infinite endpoints
#'   indicate an unbounded projection.
#'
#' @examples
#' sim <- simulate_mvmr(seed = 1)
#' searching_projection(sim, k = 1)
#'
#' @export
searching_projection <- function(object, k, alpha = 0.05, center = NULL, S = NULL,
                                 grid = 41, span = 5, max_abs = 50, tol = 1e-3) {
  stopifnot(inherits(object, "mr_input"))
  S <- S %||% relevance_screen(object)
  if (object$q > 3) {
    warning("Grid projection is intended for q <= 3; results may be coarse.",
            call. = FALSE)
  }
  center <- center %||% .feasible_center(object, alpha, S)
  if (!isTRUE(boundedness_rule(object, alpha, S)$bounded[k])) return(c(-Inf, Inf))

  others <- setdiff(seq_len(object$q), k)
  axes <- lapply(others, function(j) seq(center[j] - span, center[j] + span,
                                         length.out = grid))
  combos <- if (length(others)) as.matrix(expand.grid(axes)) else matrix(0, 1, 0)

  lo <- Inf; hi <- -Inf
  for (i in seq_len(nrow(combos))) {
    b <- center
    if (length(others)) b[others] <- combos[i, ]
    sec <- suppressWarnings(searching_section(object, k, alpha, b, S, max_abs, tol))
    if (any(is.infinite(sec))) return(c(-Inf, Inf))
    if (!any(is.na(sec))) { lo <- min(lo, sec[1]); hi <- max(hi, sec[2]) }
  }
  if (!is.finite(lo) || !is.finite(hi)) return(c(NA_real_, NA_real_))
  c(lo, hi)
}
