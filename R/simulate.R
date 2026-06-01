#' Simulate two-sample multivariable summary data with invalid instruments
#'
#' Generates SNP-exposure and SNP-outcome summary statistics under a linear
#' model with a known valid set, for examples, tests and reproduction. The valid
#' instruments have zero direct effect; the remaining instruments carry direct
#' effects following the chosen invalidity pattern.
#'
#' @param q Number of exposures.
#' @param p Number of instruments.
#' @param n_valid Number of valid instruments (must satisfy the majority-plus-q
#'   rule for a bounded region).
#' @param beta True causal vector of length \eqn{q}; a sensible default is used
#'   when `NULL`.
#' @param strength Range `c(low, high)` for the magnitude of first-stage
#'   associations.
#' @param n_X,n_Y Exposure- and outcome-sample sizes; set the standard errors.
#' @param rho Correlation injected between exposures' first-stage vectors.
#' @param invalidity Direct-effect pattern: `"symmetric"`, `"plurality"` or
#'   `"directional"`.
#' @param tau Direct-effect scale.
#' @param seed Optional random seed.
#'
#' @return An [mr_input()] object with attributes `true_beta` and `valid_set`.
#'
#' @examples
#' sim <- simulate_mvmr(q = 2, p = 30, n_valid = 21, seed = 1)
#' attr(sim, "true_beta")
#'
#' @export
simulate_mvmr <- function(q = 2, p = 30, n_valid = 21, beta = NULL,
                          strength = c(0.02, 0.06), n_X = 5e4, n_Y = 5e4,
                          rho = 0,
                          invalidity = c("symmetric", "plurality", "directional"),
                          tau = 0.2, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  invalidity <- match.arg(invalidity)
  if (is.null(beta)) {
    beta <- if (q == 2) c(0.3, -0.2) else if (q == 3) c(0.3, -0.15, 0.1) else rep(0.3, q)
  }
  if (length(beta) != q) stop("beta must have length q.", call. = FALSE)
  if (n_valid > p) stop("n_valid cannot exceed p.", call. = FALSE)

  ## first-stage effects: signed magnitudes in the strength band
  G <- matrix(stats::runif(p * q, strength[1], strength[2]), p, q) *
    matrix(sample(c(-1, 1), p * q, TRUE), p, q)
  if (rho != 0 && q >= 2) {
    shared <- G[, 1]
    for (k in 2:q) G[, k] <- sqrt(1 - rho^2) * G[, k] + rho * shared
  }

  ## direct (invalidity) effects, zero on the valid set
  pi <- numeric(p)
  invalid <- setdiff(seq_len(p), seq_len(n_valid))
  if (length(invalid)) {
    pi[invalid] <- switch(
      invalidity,
      symmetric  = sample(c(-1, 1), length(invalid), TRUE) * tau,
      plurality  = rep_len(c(5 * tau, -5 * tau, 6 * tau), length(invalid)),
      directional = stats::runif(length(invalid), 0.5 * tau, 1.5 * tau)
    )
  }

  Gamma <- as.numeric(G %*% beta) + pi
  se_X <- matrix(1 / sqrt(n_X), p, q)
  se_Y <- rep(1 / sqrt(n_Y), p)

  beta_X <- G + matrix(stats::rnorm(p * q), p, q) * se_X
  beta_Y <- Gamma + stats::rnorm(p) * se_Y

  obj <- mr_input(beta_Y, se_Y, beta_X, se_X, n_X = n_X, n_Y = n_Y)
  attr(obj, "true_beta") <- beta
  attr(obj, "valid_set") <- seq_len(n_valid)
  obj
}
