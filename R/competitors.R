#' Multivariable inverse-variance-weighted estimate
#'
#' Inverse-variance-weighted multivariable estimate over the selected set,
#' \eqn{\hat\beta = (B_X^\top W B_X)^{-1} B_X^\top W \hat\Gamma} with
#' \eqn{W = \mathrm{diag}(\mathrm{se}_{Y}^{-2})}. Provided as a reference center
#' for [searching_region()] and as a conventional comparator. It assumes all
#' selected instruments are valid and carries no robustness guarantee.
#'
#' @inheritParams searching_membership
#'
#' @return An object of class `"mv_ivw"` with `estimate`, `se`, `ci` (a
#'   \eqn{q \times 2} matrix of Wald intervals) and `exposure`.
#'
#' @examples
#' mv_ivw(simulate_mvmr(seed = 1))
#'
#' @export
mv_ivw <- function(object, S = NULL, alpha = 0.05) {
  stopifnot(inherits(object, "mr_input"))
  S <- S %||% relevance_screen(object)
  BX <- object$beta_X[S, , drop = FALSE]
  BY <- object$beta_Y[S]
  w <- 1 / object$se_Y[S]^2
  XtWX <- crossprod(BX, w * BX)
  XtWY <- crossprod(BX, w * BY)
  V <- solve(XtWX)
  est <- as.numeric(V %*% XtWY)
  se <- sqrt(diag(V))
  zc <- stats::qnorm(1 - alpha / 2)
  ci <- cbind(est - zc * se, est + zc * se)
  rownames(ci) <- object$exposure
  colnames(ci) <- c("lower", "upper")
  structure(list(estimate = stats::setNames(est, object$exposure), se = se,
                 ci = ci, alpha = alpha, exposure = object$exposure),
            class = "mv_ivw")
}

#' @export
print.mv_ivw <- function(x, ...) {
  cat("<mv_ivw> multivariable inverse-variance-weighted estimate\n")
  tab <- data.frame(estimate = round(x$estimate, 4),
                    se = round(x$se, 4),
                    lower = round(x$ci[, 1], 4),
                    upper = round(x$ci[, 2], 4))
  rownames(tab) <- x$exposure
  print(tab)
  invisible(x)
}

#' Iterated hard-thresholding valid-set estimate
#'
#' A practical two-stage hard-thresholding (TSHT) valid-set estimate used for
#' the surplus diagnostic and as a comparator. Starting from the [mv_ivw()] fit,
#' instruments whose standardized residual exceeds the cutoff are removed and
#' the estimate is refitted, iterating to convergence. The uniform searching
#' region does *not* use this selection; it is provided only for diagnostics
#' and comparison.
#'
#' @inheritParams searching_membership
#' @param max_iter Maximum number of refitting iterations.
#'
#' @return A list with `valid` (estimated valid-set indices), `estimate`
#'   (refitted on the valid set), `surplus` (\eqn{|\hat{\mathcal V}| -
#'   (\lceil |\hat{\mathcal S}|/2 \rceil + q)}) and `S`.
#'
#' @examples
#' tsht_valid_set(simulate_mvmr(seed = 1))$surplus
#'
#' @export
tsht_valid_set <- function(object, alpha = 0.05, S = NULL, max_iter = 10) {
  stopifnot(inherits(object, "mr_input"))
  S <- S %||% relevance_screen(object)
  z <- .zalpha(alpha, length(S))
  valid <- S
  beta <- mv_ivw(object, S)$estimate
  for (i in seq_len(max_iter)) {
    r <- .residuals(object, beta)
    keep <- S[abs(r$pi[S]) < z * sqrt(r$Rn[S])]
    if (length(keep) <= object$q) break
    new_beta <- mv_ivw(object, keep)$estimate
    if (setequal(keep, valid)) { valid <- keep; beta <- new_beta; break }
    valid <- keep
    beta <- new_beta
  }
  surplus <- length(valid) - (ceiling(length(S) / 2) + object$q)
  list(valid = valid, estimate = beta, surplus = surplus, S = S)
}
