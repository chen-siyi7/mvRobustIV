#' First-stage diagnostics
#'
#' Collects the diagnostics that govern whether the searching region is bounded
#' and well conditioned: the size of the selected set, per-exposure
#' marginal-relevance counts with the [boundedness_rule()], an approximate
#' conditional \eqn{F}, variance inflation factors, partial condition numbers,
#' and the hard-thresholding surplus.
#'
#' The conditional \eqn{F} reported here is the mean squared partial
#' \eqn{z}-value after regressing each exposure's associations on the others; it
#' is a strength summary read in one direction only, a small value flags weak
#' conditioning while a large value does not certify the per-coordinate strength
#' that boundedness needs. For the formal Sanderson-Windmeijer conditional
#' \eqn{F} use a dedicated multivariable MR package.
#'
#' @inheritParams searching_membership
#'
#' @return An object of class `"iv_diagnostics"`: a list with `n_instruments`,
#'   `per_exposure` (a data frame of `cond_F`, `vif`, `kappa`, `strong`,
#'   `bounded`), `majority_threshold`, `surplus` and `joint_bounded`.
#'
#' @examples
#' iv_diagnostics(simulate_mvmr(q = 3, p = 50, seed = 2))
#'
#' @export
iv_diagnostics <- function(object, alpha = 0.05, S = NULL) {
  stopifnot(inherits(object, "mr_input"))
  S <- S %||% relevance_screen(object)
  BX <- object$beta_X[S, , drop = FALSE]
  SE <- object$se_X[S, , drop = FALSE]
  q <- object$q
  m <- length(S)

  bnd <- boundedness_rule(object, alpha, S)

  ## Approximate conditional F and VIF from partialling each exposure on others.
  cond_F <- numeric(q)
  vif <- rep(NA_real_, q)
  R <- if (q > 1) suppressWarnings(stats::cor(BX)) else matrix(1, 1, 1)
  for (k in seq_len(q)) {
    if (q == 1) {
      resid <- BX[, 1]
    } else {
      Xk <- BX[, -k, drop = FALSE]
      fit <- stats::lm.fit(cbind(1, Xk), BX[, k])
      resid <- fit$residuals
      r2 <- 1 - sum(resid^2) / sum((BX[, k] - mean(BX[, k]))^2)
      vif[k] <- if (r2 < 1) 1 / (1 - r2) else Inf
    }
    cond_F[k] <- mean((resid / SE[, k])^2)
  }

  ## Partial condition numbers of the selected first-stage matrix.
  kappa <- rep(NA_real_, q)
  if (m >= q) {
    G <- BX
    GtG_inv <- tryCatch(solve(crossprod(G)), error = function(e) NULL)
    if (!is.null(GtG_inv)) {
      for (k in seq_len(q)) {
        ek <- numeric(q); ek[k] <- 1
        kappa[k] <- sqrt(sum((G %*% (GtG_inv %*% ek))^2))
      }
    }
  }

  per_exposure <- data.frame(
    exposure = object$exposure,
    cond_F = round(cond_F, 2),
    vif = round(vif, 3),
    kappa = round(kappa, 3),
    strong = bnd$strong,
    bounded = bnd$bounded,
    row.names = NULL, stringsAsFactors = FALSE
  )

  ts <- tsht_valid_set(object, alpha, S)

  structure(list(
    n_instruments = m,
    per_exposure = per_exposure,
    majority_threshold = floor(m / 2),
    surplus = ts$surplus,
    joint_bounded = all(bnd$bounded),
    alpha = alpha
  ), class = "iv_diagnostics")
}

#' @export
print.iv_diagnostics <- function(x, ...) {
  cat("<iv_diagnostics>\n")
  cat(sprintf("  selected instruments: %d   majority threshold: %d\n",
              x$n_instruments, x$majority_threshold))
  cat(sprintf("  hard-thresholding surplus: %+d\n", x$surplus))
  cat(sprintf("  joint region bounded (predicted): %s\n",
              if (x$joint_bounded) "yes" else "no"))
  cat("  per exposure:\n")
  print(x$per_exposure, row.names = FALSE)
  invisible(x)
}
