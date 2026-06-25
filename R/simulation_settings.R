## ============================================================
##  simulation_settings.R
##
##  The six data-generating processes S1-S6. Parameters reproduce the
##  Supplementary Material specification (tau = 0.2, sigma_x = sigma_y = 1,
##  500 replications, seed 2024). The few magnitude/correlation choices that
##  the supplement leaves to the implementation are collected in DGP_PARS so
##  you can set them to match your validated repository exactly.
##
##  Each generator returns the population truth; draw_sample() then adds the
##  two-sample summary-statistic noise.
## ============================================================

DGP_PARS <- list(
  tau     = 0.2,        # invalidity magnitude
  sigma_x = 1,          # per-SNP exposure noise scale
  sigma_y = 1,          # per-SNP outcome noise scale
  ## first-stage effect magnitudes for the STRONG settings (S1-S5).
  ## Align these with your repo if you need the exact published numbers.
  gamma_lo = 0.10,
  gamma_hi = 0.40,
  s2_load  = 0.6,       # S2: gamma_{j,2} = s2_load * gamma_{j,1} + s2_noise * noise
  s2_noise = 0.8,
  s6_lo    = 0.005,     # S6: very weak first stage
  s6_hi    = 0.015
)

## helper: build a population object
.make_pop <- function(Gstar, pistar, beta_star, valid, n_exp, n_out) {
  Gammastar <- as.numeric(Gstar %*% beta_star) + pistar
  list(Gstar = Gstar, pistar = pistar, Gammastar = Gammastar,
       beta_star = beta_star, valid = valid, n_exp = n_exp, n_out = n_out,
       n = round((n_exp + n_out) / 2))
}

## ---- S1: majority, q = 2 ----
pop_S1 <- function(P = DGP_PARS) {
  set_dim <- list(pz = 30, nv = 21, q = 2, beta = c(0.3, -0.2),
                  n_exp = 50000, n_out = 50000)
  with(set_dim, {
    Gstar <- matrix(stats::runif(pz * q, P$gamma_lo, P$gamma_hi), pz, q)
    valid <- sort(sample.int(pz, nv))
    pistar <- numeric(pz)
    invalid <- setdiff(seq_len(pz), valid)
    pistar[invalid] <- stats::runif(length(invalid), -P$tau, P$tau) *
                       sample(c(-1, 1), length(invalid), TRUE)
    .make_pop(Gstar, pistar, beta, valid, n_exp, n_out)
  })
}

## ---- S2: correlated exposures, q = 2 ----
pop_S2 <- function(P = DGP_PARS) {
  pz <- 30; nv <- 21; q <- 2; beta <- c(0.3, -0.2)
  g1 <- stats::runif(pz, P$gamma_lo, P$gamma_hi)
  g2 <- P$s2_load * g1 + P$s2_noise * stats::runif(pz, P$gamma_lo, P$gamma_hi)
  Gstar <- cbind(g1, g2)
  valid <- sort(sample.int(pz, nv)); invalid <- setdiff(seq_len(pz), valid)
  pistar <- numeric(pz)
  pistar[invalid] <- stats::runif(length(invalid), -P$tau, P$tau) *
                     sample(c(-1, 1), length(invalid), TRUE)
  .make_pop(Gstar, pistar, beta, valid, 50000, 50000)
}

## ---- S3: three traits, q = 3 ----
pop_S3 <- function(P = DGP_PARS) {
  pz <- 50; nv <- 35; q <- 3; beta <- c(0.3, -0.15, 0.1)
  Gstar <- matrix(stats::runif(pz * q, P$gamma_lo, P$gamma_hi), pz, q)
  valid <- sort(sample.int(pz, nv)); invalid <- setdiff(seq_len(pz), valid)
  pistar <- numeric(pz)
  pistar[invalid] <- stats::runif(length(invalid), -P$tau, P$tau) *
                     sample(c(-1, 1), length(invalid), TRUE)
  .make_pop(Gstar, pistar, beta, valid, 100000, 80000)
}

## ---- S4: multimodal plurality, q = 2 (clusters: 5 at +tau, 5 at -tau, 3 at +2tau) ----
pop_S4 <- function(P = DGP_PARS) {
  pz <- 30; nv <- 17; q <- 2; beta <- c(0.3, -0.2)
  Gstar <- matrix(stats::runif(pz * q, P$gamma_lo, P$gamma_hi), pz, q)
  valid <- sort(sample.int(pz, nv)); invalid <- setdiff(seq_len(pz), valid)
  pistar <- numeric(pz)
  ## assign clusters among the invalid set
  cl <- invalid
  pistar[cl[1:5]]   <-  P$tau
  pistar[cl[6:10]]  <- -P$tau
  pistar[cl[11:13]] <-  2 * P$tau
  ## any remaining invalid get small uniform pleiotropy
  if (length(cl) > 13) pistar[cl[14:length(cl)]] <-
      stats::runif(length(cl) - 13, -P$tau, P$tau)
  .make_pop(Gstar, pistar, beta, valid, 50000, 50000)
}

## ---- S5: directional pleiotropy, q = 2 (all invalid pi > 0, violates InSIDE) ----
pop_S5 <- function(P = DGP_PARS) {
  pz <- 30; nv <- 21; q <- 2; beta <- c(0.3, -0.2)
  Gstar <- matrix(stats::runif(pz * q, P$gamma_lo, P$gamma_hi), pz, q)
  valid <- sort(sample.int(pz, nv)); invalid <- setdiff(seq_len(pz), valid)
  pistar <- numeric(pz)
  pistar[invalid] <- stats::runif(length(invalid), 0.5 * P$tau, 1.5 * P$tau)
  .make_pop(Gstar, pistar, beta, valid, 50000, 50000)
}

## ---- S6: many weak instruments, q = 2 (diagnostic counter-example) ----
pop_S6 <- function(P = DGP_PARS) {
  pz <- 100; nv <- 70; q <- 2; beta <- c(0.3, -0.2)
  Gstar <- matrix(stats::runif(pz * q, P$s6_lo, P$s6_hi), pz, q)
  valid <- sort(sample.int(pz, nv)); invalid <- setdiff(seq_len(pz), valid)
  pistar <- numeric(pz)
  pistar[invalid] <- stats::runif(length(invalid), -P$tau, P$tau) *
                     sample(c(-1, 1), length(invalid), TRUE)
  .make_pop(Gstar, pistar, beta, valid, 20000, 20000)
}

SETTINGS <- list(
  S1 = list(label = "S1: majority, $q{=}2$",            gen = pop_S1),
  S2 = list(label = "S2: correlated, $q{=}2$",          gen = pop_S2),
  S3 = list(label = "S3: three traits, $q{=}3$",        gen = pop_S3),
  S4 = list(label = "S4: plurality, $q{=}2$",           gen = pop_S4),
  S5 = list(label = "S5: directional, $q{=}2$",         gen = pop_S5),
  S6 = list(label = "S6: strong-IV violation, $q{=}2$", gen = pop_S6)
)

## ---- one summary-statistic draw from a population ----
draw_sample <- function(pop, P = DGP_PARS) {
  pz <- nrow(pop$Gstar); q <- ncol(pop$Gstar)
  seG     <- matrix(P$sigma_x / sqrt(pop$n_exp), pz, q)
  seGamma <- rep(P$sigma_y / sqrt(pop$n_out), pz)
  G     <- pop$Gstar + matrix(stats::rnorm(pz * q, 0, seG), pz, q)
  Gamma <- pop$Gammastar + stats::rnorm(pz, 0, seGamma)
  list(Gamma = Gamma, seGamma = seGamma, G = G, seG = seG, n = pop$n,
       valid = pop$valid, beta_star = pop$beta_star)
}
