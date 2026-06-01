#' First-stage relevance screen
#'
#' Selects the instrument set \eqn{\hat{\mathcal S}} by the joint first-stage
#' rule \eqn{n\,\hat\gamma_j^\top \hat V_{G,jj}^{-1} \hat\gamma_j \ge q \log n}.
#' With diagonal SNP-exposure covariance the statistic reduces to
#' \eqn{\sum_k (\hat\gamma_{jk}/\mathrm{se}_{\gamma,jk})^2}, a sum of squared
#' per-exposure \eqn{z}-values.
#'
#' @param object An [mr_input()] object.
#' @param n Sample size for the \eqn{q \log n} cutoff. Defaults to `object$n_X`.
#'   If neither is available, a \eqn{\chi^2_q} cutoff at level `chisq_level` is
#'   used and a warning is issued.
#' @param chisq_level Tail level for the fallback \eqn{\chi^2_q} cutoff.
#'
#' @return Integer vector of selected instrument indices (into the rows of
#'   `object`).
#'
#' @examples
#' sim <- simulate_mvmr(seed = 1)
#' S <- relevance_screen(sim)
#' length(S)
#'
#' @export
relevance_screen <- function(object, n = NULL, chisq_level = 0.999) {
  stopifnot(inherits(object, "mr_input"))
  stat <- rowSums((object$beta_X / object$se_X)^2)
  n <- n %||% object$n_X
  if (!is.null(n)) {
    cutoff <- object$q * log(n)
  } else {
    cutoff <- stats::qchisq(chisq_level, df = object$q)
    warning("n_X not supplied; using a chi-squared(q) relevance cutoff. ",
            "Provide n_X for the q*log(n) screen.", call. = FALSE)
  }
  which(stat >= cutoff)
}
