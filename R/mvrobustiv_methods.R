## ============================================================
##  mvrobustiv_methods.R
##
##  Reference implementation of the searching/sampling region for
##  multivariable invalid-instrument Mendelian randomization, the
##  first-stage screen and diagnostics, and the competitor methods used
##  in the simulation and application.
##
##  Data are summary statistics in a list `dat` with components
##    Gamma   : length-p vector of SNP-outcome estimates  (hatGamma_j)
##    seGamma : length-p vector of their standard errors
##    G       : p x q matrix of SNP-exposure estimates     (hatgamma_jk)
##    seG     : p x q matrix of their standard errors
##    n       : common growth index used in the log n screen (scalar)
##  Effects are on the scale Gamma_j = gamma_j^T beta + pi_j.
##
##  Conventions match the manuscript: residual contrast
##    pihat_j(beta) = hatGamma_j - hatgamma_j^T beta,
##  variance R_jj(beta)/n estimated by seGamma_j^2 + sum_k seG_jk^2 beta_k^2
##  (independence across exposures and between the two samples), and the
##  Bonferroni residual test at level alpha over the screened set.
## ============================================================

## ---- first-stage screen: Shat = { j : n gamma_j^T V_{G,jj}^{-1} gamma_j >= q log n } ----
## With diagonal V_{G,jj} = diag(seG_j^2), the statistic is sum_k (G_jk/seG_jk)^2.
screen_instruments <- function(dat, thresh = NULL) {
  q <- ncol(dat$G)
  chi <- rowSums((dat$G / dat$seG)^2)              # joint first-stage chi-square
  if (is.null(thresh)) thresh <- q * log(dat$n)    # q log n
  which(chi >= thresh)
}

## ---- Bonferroni quantile for the residual test over |Shat| instruments ----
z_alpha <- function(nS, alpha = 0.05) stats::qnorm(1 - alpha / (2 * nS))

## ---- per-SNP qxq covariance of the SNP-exposure estimate vector ----
## Returns a list over `idx` of qxq covariance matrices Sigma_{gamma,j}, or
## NULL when only marginal standard errors are available (independence across
## exposures, the two-sample default). The off-diagonals are needed when the q
## exposure associations are estimated in the same or overlapping samples;
## supply them through dat$Sigma_g (a length-p list or p x q x q array of
## per-SNP covariances) or, as the standard same-sample approximation, through
## dat$gencov, a q x q phenotypic/observational correlation matrix of the
## exposures, in which case
##   Sigma_{gamma,j} = diag(seG_j) %*% gencov %*% diag(seG_j),
## that is cov(hatgamma_jk, hatgamma_jl) = gencov_kl * seG_jk * seG_jl.
## A diagonal gencov (the identity) reproduces the marginal-SE case exactly.
.snp_cov_list <- function(dat, idx) {
  q <- ncol(dat$G)
  if (!is.null(dat$Sigma_g)) {
    Sg <- dat$Sigma_g
    if (is.array(Sg) && length(dim(Sg)) == 3L)
      return(lapply(idx, function(j) Sg[j, , ]))
    return(Sg[idx])                                   # list of length p
  }
  if (!is.null(dat$gencov)) {
    R <- dat$gencov
    if (isTRUE(all.equal(R, diag(q)))) return(NULL)   # identity -> diagonal path
    seG <- dat$seG[idx, , drop = FALSE]
    return(lapply(seq_along(idx), function(i) {
      D <- diag(seG[i, ], q); D %*% R %*% D }))
  }
  NULL
}

## ---- variance of the residual contrast pi_j(beta) = Gamma_j - gamma_j^T beta ----
## Var{pi_j(beta)} = seGamma_j^2 + beta^T Sigma_{gamma,j} beta, using the
## two-sample independence of the exposure and outcome samples. With only
## marginal SEs this is the diagonal seGamma_j^2 + sum_k seG_jk^2 beta_k^2; when
## a per-SNP covariance is available the off-diagonal cross terms are added.
## Returns an s x m matrix for a candidate matrix B (m x q), s = length(idx).
.pi_var_mat <- function(B, dat, idx, covlist = NULL) {
  seGam2 <- dat$seGamma[idx]^2
  seGm2  <- dat$seG[idx, , drop = FALSE]^2
  V <- seGm2 %*% t(B^2)                               # s x m diagonal part
  V <- V + seGam2                                     # recycles down rows
  if (!is.null(covlist)) {
    for (i in seq_along(idx)) {
      O <- covlist[[i]]; diag(O) <- 0                 # off-diagonal cross terms
      if (any(O != 0)) V[i, ] <- V[i, ] + rowSums((B %*% O) * B)
    }
  }
  ## optional two-sample overlap term: cov(hatGamma_j, hatgamma_jk) carried in
  ## dat$cross (p x q). Under independence (default) dat$cross is NULL. The
  ## contribution to Var{pi_j(beta)} is -2 beta^T cov(hatGamma_j, hatgamma_j).
  if (!is.null(dat$cross)) {
    cr <- dat$cross[idx, , drop = FALSE]
    V <- V - 2 * (cr %*% t(B))
    if (any(V <= 0)) { pos <- V[V > 0]; V[V <= 0] <- if (length(pos)) min(pos) * 1e-3 else 1e-12 }
  }
  V
}

## ---- standardized residuals at a candidate beta ----
residual_stats <- function(beta, dat, idx = seq_along(dat$Gamma)) {
  G  <- dat$G[idx, , drop = FALSE]
  r  <- dat$Gamma[idx] - as.numeric(G %*% beta)
  v  <- as.numeric(.pi_var_mat(matrix(beta, nrow = 1L), dat, idx,
                               .snp_cov_list(dat, idx)))
  se <- sqrt(v)
  list(r = r, se = se, z = abs(r) / se)
}

## ---- membership in the searching region: fewer than a majority flagged ----
## flagged_j(beta) = { |z_j(beta)| > za }; accept beta if #flagged <= T - 1,
## T = floor(|Shat|/2).
n_flagged <- function(beta, dat, idx, za) {
  st <- residual_stats(beta, dat, idx)
  sum(st$z > za)
}
in_region <- function(beta, dat, idx, za, Tmaj) {
  n_flagged(beta, dat, idx, za) <= Tmaj - 1L
}

## ---- coordinate-wise SECTION via bisection at a reference beta0 ----
## Returns the interval of beta_k (with beta_{-k} fixed at beta0) for which
## beta is accepted; marks unbounded if the accepted set reaches +/- span.
searching_section <- function(dat, idx, za, Tmaj, beta0, k,
                              span = 50, tol = 1e-3) {
  q <- length(beta0)
  acc <- function(bk) { b <- beta0; b[k] <- bk; in_region(b, dat, idx, za, Tmaj) }
  ## scan a coarse grid to locate the accepted band around beta0[k]
  grid <- seq(beta0[k] - span, beta0[k] + span, length.out = 2001)
  ok <- vapply(grid, acc, logical(1))
  if (!any(ok)) return(c(lower = NA, upper = NA, unbounded = NA))   # empty section (rare)
  lo_idx <- which(ok)[1]; hi_idx <- rev(which(ok))[1]
  lower <- grid[lo_idx]; upper <- grid[hi_idx]
  unbounded <- (lo_idx == 1L) || (hi_idx == length(grid))
  ## refine the two boundaries by bisection unless unbounded on that side
  refine <- function(inside, outside) {
    for (i in 1:40) { mid <- (inside + outside)/2; if (acc(mid)) inside <- mid else outside <- mid
      if (abs(inside - outside) < tol) break }
    inside
  }
  if (lo_idx > 1L)            lower <- refine(grid[lo_idx], grid[lo_idx - 1L])
  if (hi_idx < length(grid))  upper <- refine(grid[hi_idx], grid[hi_idx + 1L])
  c(lower = lower, upper = upper, unbounded = as.numeric(unbounded))
}

## ---- vectorized acceptance over a matrix B (m x q) of candidate betas ----
## Returns a logical vector of length m; the workhorse for grid/bootstrap.
## `covlist` (per-SNP qxq covariances over idx) is computed once when NULL.
accept_many <- function(B, dat, idx, za, Tmaj, covlist = NULL) {
  if (is.null(covlist)) covlist <- .snp_cov_list(dat, idx)
  Gm <- dat$G[idx, , drop = FALSE]; Gam <- dat$Gamma[idx]
  r  <- Gam - Gm %*% t(B)                               # s x m
  z  <- abs(r) / sqrt(.pi_var_mat(B, dat, idx, covlist))
  colSums(z > za) <= Tmaj - 1L
}

## ---- vectorized flag COUNT over a matrix B (m x q) ----
flag_count_many <- function(B, dat, idx, za, covlist = NULL) {
  if (is.null(covlist)) covlist <- .snp_cov_list(dat, idx)
  Gm <- dat$G[idx, , drop = FALSE]; Gam <- dat$Gamma[idx]
  r <- Gam - Gm %*% t(B)
  z <- abs(r) / sqrt(.pi_var_mat(B, dat, idx, covlist))
  colSums(z > za)
}

## ---- least-flagged interior point ----
## The membership rule is reference-free, but the box placement and the
## section need a point that actually lies in the region. We take the
## fewest-flagged candidate among the robust estimators and, if needed,
## descend the flag count by a short coordinate search.
interior_point <- function(dat, idx, za, Tmaj, extra = NULL) {
  cands <- list(mv_tsht(dat, idx = idx)$beta, mv_ivw(dat, idx)$beta,
                mv_egger(dat, idx)$beta)
  if (!is.null(extra)) cands <- c(list(extra), cands)
  fl <- vapply(cands, function(b) flag_count_many(matrix(b, nrow = 1), dat, idx, za), numeric(1))
  b <- cands[[which.min(fl)]]; best <- min(fl)
  if (best <= Tmaj - 1L) return(b)                 # already interior
  ## coordinate descent on the integer flag count
  step <- 0.2
  for (sweep in 1:6) {
    improved <- FALSE
    for (k in seq_along(b)) for (s in c(step, -step)) {
      bb <- b; bb[k] <- bb[k] + s
      f <- flag_count_many(matrix(bb, nrow = 1), dat, idx, za)
      if (f < best) { b <- bb; best <- f; improved <- TRUE }
    }
    if (!improved) step <- step / 2
    if (best <= Tmaj - 1L) break
  }
  b
}

## ---- enclosing grid: bracket the whole bounded region ----
## Centre a product grid at `center` and expand each axis until no accepted
## point touches a boundary (the region is enclosed) or the cap is reached
## (that axis is unbounded). Resolution is set by a target cell size so a
## narrow region is resolved; `ng_max` caps cost. Returns the accepted grid.
enclosing_grid <- function(dat, idx, za, Tmaj, center, h0 = 0.25, cap = 8,
                           target = 0.006, ng_max = 161, iters = 10,
                           covlist = NULL) {
  if (is.null(covlist)) covlist <- .snp_cov_list(dat, idx)
  q <- length(center); h <- rep(h0, q); unb <- rep(FALSE, q)
  axes <- gridpts <- acc <- NULL
  build <- function() {
    ng <- pmin(pmax(ceiling(2 * h / target), 41L), ng_max)
    axes <<- lapply(seq_len(q), function(k)
      seq(center[k] - h[k], center[k] + h[k], length.out = ng[k]))
    gridpts <<- as.matrix(expand.grid(axes))
    acc <<- accept_many(gridpts, dat, idx, za, Tmaj, covlist)
  }
  build()
  tries <- 0
  while (!any(acc) && tries < iters) { h <- h * 2; tries <- tries + 1
    if (all(h >= cap)) break; build() }
  for (it in seq_len(iters)) {
    if (!any(acc)) break
    touch <- logical(q)
    for (k in seq_len(q)) { bk <- gridpts[acc, k]
      touch[k] <- (min(bk) <= axes[[k]][1] + 1e-9) ||
                  (max(bk) >= axes[[k]][length(axes[[k]])] - 1e-9) }
    if (!any(touch)) break
    grow <- FALSE
    for (k in which(touch)) {
      if (h[k] < cap - 1e-9) { h[k] <- min(h[k] * 1.7, cap); grow <- TRUE }
      else unb[k] <- TRUE
    }
    if (!grow) break
    build()
  }
  if (any(acc)) for (k in seq_len(q)) { bk <- gridpts[acc, k]
    if (((min(bk) <= axes[[k]][1] + 1e-9) ||
         (max(bk) >= axes[[k]][length(axes[[k]])] - 1e-9)) && h[k] >= cap - 1e-9)
      unb[k] <- TRUE }
  cell <- prod(vapply(axes, function(a) a[2] - a[1], numeric(1)))
  list(axes = axes, gridpts = gridpts, acc = acc, cell = cell, unbounded = unb)
}

## ---- full searching-region summary used by the simulation/application ----
## projection : CI^proj_k of Corollary 1, the shadow of the region onto axis k,
##              computed from a grid certified to enclose the region (so it can
##              never be narrower than the joint coverage allows).
## section    : the cheaper diagnostic slice through an interior point.
## volume     : Vol_q of the region on the enclosing grid.
## When `certify = TRUE`, every coordinate flagged unbounded at radius `cap` is
## re-tested on a grid of radius 2*cap at half the cell size, and bounded
## endpoints are re-read at double resolution; the returned `certified` and
## `proj_doubled` let the caller confirm that "unbounded" is not a grid artifact
## and that endpoints are stable. `Rmax` records the search radius (cap) used.
searching_region <- function(dat, alpha = 0.05, idx = NULL, beta0 = NULL,
                             span_section = 50, cap = 8,
                             target = NULL, ng_max = NULL, certify = FALSE) {
  if (is.null(idx)) idx <- screen_instruments(dat)
  nS <- length(idx); Tmaj <- floor(nS / 2); za <- z_alpha(nS, alpha)
  q <- ncol(dat$G)
  if (is.null(target)) target <- if (q >= 3) 0.02  else 0.006
  if (is.null(ng_max)) ng_max <- if (q >= 3) 61L   else 161L
  covlist <- .snp_cov_list(dat, idx)
  cstart <- interior_point(dat, idx, za, Tmaj, extra = beta0)
  gr <- enclosing_grid(dat, idx, za, Tmaj, center = cstart, cap = cap,
                       target = target, ng_max = ng_max, covlist = covlist)
  vol <- sum(gr$acc) * gr$cell
  proj <- matrix(NA_real_, q, 3, dimnames = list(NULL, c("lower","upper","unbounded")))
  csec <- cstart
  if (any(gr$acc)) {
    accpts <- gr$gridpts[gr$acc, , drop = FALSE]
    for (k in seq_len(q))
      proj[k, ] <- c(min(accpts[, k]), max(accpts[, k]), as.numeric(gr$unbounded[k]))
    ## section reference = region centre (projection midpoint where bounded),
    ## falling back to the fewest-flagged point for unbounded coordinates.
    nf <- flag_count_many(accpts, dat, idx, za, covlist)
    cmin <- accpts[which.min(nf), ]
    csec <- vapply(seq_len(q), function(k)
      if (gr$unbounded[k]) cmin[k] else (proj[k, "lower"] + proj[k, "upper"]) / 2,
      numeric(1))
    if (flag_count_many(matrix(csec, nrow = 1), dat, idx, za, covlist) > Tmaj - 1L) csec <- cmin
  }
  sec <- t(vapply(seq_len(q), function(k)
    searching_section(dat, idx, za, Tmaj, csec, k, span = span_section), numeric(3)))
  out <- list(idx = idx, nS = nS, Tmaj = Tmaj, za = za, beta0 = cstart, center = csec,
              section = sec, projection = proj, volume = vol, Rmax = cap,
              in_joint = function(beta) in_region(beta, dat, idx, za, Tmaj))
  if (certify && any(gr$acc)) {
    gr2 <- enclosing_grid(dat, idx, za, Tmaj, center = cstart, cap = 2 * cap,
                          target = target / 2, ng_max = ng_max, covlist = covlist)
    proj2 <- matrix(NA_real_, q, 3, dimnames = list(NULL, c("lower","upper","unbounded")))
    if (any(gr2$acc)) { ap2 <- gr2$gridpts[gr2$acc, , drop = FALSE]
      for (k in seq_len(q))
        proj2[k, ] <- c(min(ap2[, k]), max(ap2[, k]), as.numeric(gr2$unbounded[k])) }
    ## a coordinate is certified if its unbounded flag is stable across radius
    ## cap and 2*cap; a bounded endpoint is certified stable if it moves by less
    ## than a few cells of the finer grid or 2% of its magnitude.
    lin <- (gr2$cell)^(1 / q)
    certified <- vapply(seq_len(q), function(k) {
      if (isTRUE(gr$unbounded[k])) return(isTRUE(as.logical(proj2[k, "unbounded"])))
      tolk <- pmax(6 * lin, 0.02 * pmax(1, abs(proj2[k, 1:2])))
      all(abs(proj[k, 1:2] - proj2[k, 1:2]) < tolk)
    }, logical(1))
    out$certified <- certified; out$proj_doubled <- proj2; out$Rmax <- 2 * cap
  }
  out
}

## ---- sampling (resampling-calibrated) region ----
## CR^samp is the union over M parametric-bootstrap resamples of the searching
## region applied to each resample at threshold lambda*z_alpha. With lambda = 1
## (the default) it is the resampling-inflated searching region: wider than the
## single-sample region and nominal but uninformative, as Section 5 describes.
## Resolution is matched to the region scale so the projection is not collapsed.
## Returns union volume, per-coordinate projection, the centre, and, if
## `beta_eval` is given, whether that point is covered (accepted in some resample).
sampling_region <- function(dat, alpha = 0.05, M = 100, lambda = 1.0,
                            beta0 = NULL, span_grid = NULL, ng = NULL,
                            idx = NULL, beta_eval = NULL) {
  if (is.null(idx)) idx <- screen_instruments(dat)
  nS <- length(idx); Tmaj <- floor(nS / 2); za <- lambda * z_alpha(nS, alpha)
  q <- ncol(dat$G)
  if (is.null(beta0)) beta0 <- interior_point(dat, idx, z_alpha(nS, alpha), Tmaj)
  if (is.null(ng)) ng <- if (q >= 3) 21L else 81L
  if (is.null(span_grid)) span_grid <- if (q >= 3) 1.0 else 1.2
  axes <- lapply(seq_len(q), function(k)
    seq(beta0[k] - span_grid, beta0[k] + span_grid, length.out = ng))
  gridpts <- as.matrix(expand.grid(axes))
  cell <- prod(vapply(axes, function(a) a[2] - a[1], numeric(1)))
  covlist <- .snp_cov_list(dat, idx)
  accepted <- rep(FALSE, nrow(gridpts)); cover <- FALSE
  for (m in seq_len(M)) {
    d <- dat
    d$Gamma <- stats::rnorm(length(dat$Gamma), dat$Gamma, dat$seGamma)
    d$G     <- matrix(stats::rnorm(length(dat$G), dat$G, dat$seG), nrow = nrow(dat$G))
    accepted <- accepted | accept_many(gridpts, d, idx, za, Tmaj, covlist)
    if (!is.null(beta_eval) && !cover)
      cover <- n_flagged(beta_eval, d, idx, za) <= Tmaj - 1L
  }
  vol <- sum(accepted) * cell
  proj <- matrix(NA_real_, q, 3, dimnames = list(NULL, c("lower","upper","unbounded")))
  for (k in seq_len(q)) { bk <- gridpts[accepted, k]
    if (length(bk)) proj[k, ] <- c(min(bk), max(bk),
      as.numeric(min(bk) <= axes[[k]][1] + 1e-9 || max(bk) >= axes[[k]][ng] - 1e-9)) }
  list(volume = vol, projection = proj, center = beta0, cover = cover)
}

## ============================================================
##  Competitor methods
## ============================================================

## ---- multivariable IVW (iterated GLS), with per-coordinate Wald CIs ----
mv_ivw <- function(dat, idx = seq_along(dat$Gamma), alpha = 0.05, iter = 3) {
  G <- dat$G[idx, , drop = FALSE]; Gam <- dat$Gamma[idx]
  seGam <- dat$seGamma[idx]; seG <- dat$seG[idx, , drop = FALSE]
  beta <- rep(0, ncol(G))
  for (i in seq_len(iter)) {
    w <- 1 / (seGam^2 + as.numeric((seG^2) %*% (beta^2)))
    A <- crossprod(G, G * w); b <- crossprod(G, Gam * w)
    beta <- as.numeric(solve(A, b))
  }
  w <- 1 / (seGam^2 + as.numeric((seG^2) %*% (beta^2)))
  cov <- solve(crossprod(G, G * w))
  se <- sqrt(diag(cov)); zq <- stats::qnorm(1 - alpha/2)
  list(beta = beta, se = se, cov = cov,
       ci = cbind(beta - zq*se, beta + zq*se))
}

## ---- MV-Egger (intercept for directional pleiotropy) ----
mv_egger <- function(dat, idx = seq_along(dat$Gamma), alpha = 0.05) {
  G <- cbind(1, dat$G[idx, , drop = FALSE]); Gam <- dat$Gamma[idx]
  w <- 1 / dat$seGamma[idx]^2
  A <- crossprod(G, G * w); b <- crossprod(G, Gam * w)
  est <- as.numeric(solve(A, b)); cov <- solve(A)
  q <- ncol(dat$G); sl <- 2:(q+1)
  se <- sqrt(diag(cov))[sl]; zq <- stats::qnorm(1 - alpha/2)
  list(beta = est[sl], se = se, ci = cbind(est[sl]-zq*se, est[sl]+zq*se),
       intercept = est[1])
}

## ---- Oracle: MV-IVW restricted to the true valid set ----
mv_oracle <- function(dat, valid, alpha = 0.05) mv_ivw(dat, idx = valid, alpha = alpha)

## ---- TSHT-style: hard-threshold residual votes on Shat, refit IVW on
##      the retained majority. A practical multivariable analogue. ----
mv_tsht <- function(dat, alpha = 0.05, idx = NULL) {
  if (is.null(idx)) idx <- screen_instruments(dat)
  fit0 <- mv_ivw(dat, idx, alpha)
  za <- z_alpha(length(idx), alpha)
  st <- residual_stats(fit0$beta, dat, idx)
  keep <- idx[st$z <= za]
  if (length(keep) < ncol(dat$G) + 1) keep <- idx          # fall back if over-pruned
  fit <- mv_ivw(dat, keep, alpha)
  c(fit, list(valid = keep))
}

## ---- MV-Median: weighted least-absolute-deviations (IRLS), bootstrap SE ----
## The multivariable analogue of the weighted median: minimise
## sum_j w_j |Gamma_j - gamma_j^T beta| with IVW weights, solved by iteratively
## reweighted least squares. Robust to a minority of invalid instruments.
mv_median <- function(dat, idx = NULL, alpha = 0.05, nboot = 100) {
  if (is.null(idx)) idx <- screen_instruments(dat)
  G <- dat$G[idx, , drop = FALSE]; Gam <- dat$Gamma[idx]
  w0 <- 1 / dat$seGamma[idx]^2
  lad <- function(G, Gam, w0) {
    beta <- as.numeric(solve(crossprod(G, G * w0), crossprod(G, Gam * w0)))
    for (i in 1:50) {
      res <- abs(Gam - as.numeric(G %*% beta)); res[res < 1e-6] <- 1e-6
      w <- w0 / res
      bn <- as.numeric(solve(crossprod(G, G * w), crossprod(G, Gam * w)))
      if (max(abs(bn - beta)) < 1e-7) { beta <- bn; break }
      beta <- bn
    }
    beta
  }
  beta <- lad(G, Gam, w0)
  bb <- matrix(NA_real_, nboot, ncol(G))
  seGam <- dat$seGamma[idx]; seG <- dat$seG[idx, , drop = FALSE]
  for (b in seq_len(nboot)) {
    Gb <- G + matrix(stats::rnorm(length(G), 0, seG), nrow = nrow(G))
    Gamb <- Gam + stats::rnorm(length(Gam), 0, seGam)
    bb[b, ] <- tryCatch(lad(Gb, Gamb, w0), error = function(e) rep(NA, ncol(G)))
  }
  se <- apply(bb, 2L, stats::sd, na.rm = TRUE); zq <- stats::qnorm(1 - alpha/2)
  list(beta = beta, se = se, ci = cbind(beta - zq*se, beta + zq*se))
}

## ---- MV-PRESSO: heterogeneity-based outlier removal, then MV-IVW ----
## Iteratively removes the instrument with the largest leave-one-out reduction
## in the weighted residual sum of squares while a global heterogeneity test
## rejects, then refits MV-IVW on the retained set. Distinct from TSHT, which
## thresholds residual votes at the Bonferroni level.
mv_presso <- function(dat, idx = NULL, alpha = 0.05, het_alpha = 0.05) {
  if (is.null(idx)) idx <- screen_instruments(dat)
  keep <- idx
  for (it in seq_along(idx)) {
    fit <- mv_ivw(dat, keep, alpha)
    st  <- residual_stats(fit$beta, dat, keep)
    rss <- sum(st$z^2); dfree <- length(keep) - ncol(dat$G)
    if (dfree <= 1) break
    pglobal <- stats::pchisq(rss, dfree, lower.tail = FALSE)
    if (pglobal > het_alpha) break                    # no residual heterogeneity
    j <- which.max(st$z^2)                             # largest outlier contribution
    if (st$z[j]^2 < stats::qchisq(1 - het_alpha / length(keep), 1)) break
    keep <- keep[-j]
    if (length(keep) < ncol(dat$G) + 1) { keep <- idx; break }
  }
  fit <- mv_ivw(dat, keep, alpha)
  c(fit, list(valid = keep))
}

## ---- MV-Lasso: L1-penalised pleiotropy (sisVIVE/CIIV style), select, refit ----
## Augmented model Gamma = G beta + pi + e with beta unpenalised and pi
## L1-penalised; solved by weighted block coordinate descent over a lambda grid,
## lambda chosen by BIC. The selected valid set {pi_j = 0} is refit by MV-IVW.
mv_lasso <- function(dat, idx = NULL, alpha = 0.05) {
  if (is.null(idx)) idx <- screen_instruments(dat)
  G <- dat$G[idx, , drop = FALSE]; Gam <- dat$Gamma[idx]
  w <- 1 / dat$seGamma[idx]^2; p <- length(idx); q <- ncol(G)
  soft <- function(x, t) sign(x) * pmax(abs(x) - t, 0)
  fit_lambda <- function(lambda) {
    beta <- as.numeric(solve(crossprod(G, G * w), crossprod(G, Gam * w)))
    pi <- numeric(p)
    for (i in 1:200) {
      pi_new <- soft(Gam - as.numeric(G %*% beta), lambda / w)
      beta_new <- as.numeric(solve(crossprod(G, G * w), crossprod(G, (Gam - pi_new) * w)))
      if (max(abs(beta_new - beta)) < 1e-7 && max(abs(pi_new - pi)) < 1e-7) {
        beta <- beta_new; pi <- pi_new; break }
      beta <- beta_new; pi <- pi_new
    }
    list(beta = beta, pi = pi)
  }
  res <- abs(Gam - as.numeric(G %*% solve(crossprod(G, G*w), crossprod(G, Gam*w))))
  lam_max <- max(w * res); lambdas <- lam_max * exp(seq(0, -4, length.out = 30))
  best <- NULL; bestbic <- Inf
  for (lam in lambdas) {
    f <- fit_lambda(lam)
    rss <- sum(w * (Gam - as.numeric(G %*% f$beta) - f$pi)^2)
    k <- sum(abs(f$pi) > 1e-8) + q
    bic <- rss + log(p) * k
    if (bic < bestbic) { bestbic <- bic; best <- f }
  }
  valid <- idx[abs(best$pi) <= 1e-8]
  if (length(valid) < q + 1) valid <- idx
  fit <- mv_ivw(dat, valid, alpha)
  c(fit, list(valid = valid))
}

## ============================================================
##  First-stage diagnostics (Section 5 of the paper)
## ============================================================

## ---- Sanderson-Windmeijer conditional F-statistic, summary form ----
## Conditional F for exposure k given the others (Sanderson, Windmeijer and
## Bowden 2019; Sanderson, Spiller and Bowden 2021), computed from summary
## statistics. For exposure k, delta minimizes the conditional Q-statistic
##   Q_k(delta) = sum_j (G_jk - G_{j,-k} delta)^2 / s2_j(delta),
##   s2_j(delta) = Sigma_j[k,k] - 2 delta' Sigma_j[-k,k] + delta' Sigma_j[-k,-k] delta,
## where Sigma_j is the per-SNP qxq covariance of the SNP-exposure estimates
## (built from dat$gencov or dat$Sigma_g; diagonal seG_jk^2 when absent). The
## reported statistic is F_k = Q_k(deltahat) / (p - q + 1), the standard MVMR
## conditional F that reduces to the mean squared marginal z-statistic when
## q = 1. With only marginal SEs and correlated exposures this overstates
## conditional strength, so dat$gencov should be supplied (see build_gencov).
conditional_F <- function(dat, idx = NULL) {
  if (is.null(idx)) idx <- screen_instruments(dat)
  G <- dat$G[idx, , drop = FALSE]; seG <- dat$seG[idx, , drop = FALSE]
  p <- nrow(G); q <- ncol(G)
  if (q == 1L) return(mean((G[, 1] / seG[, 1])^2))
  cov_j <- .snp_cov_list(dat, idx)
  if (is.null(cov_j)) cov_j <- lapply(seq_len(p), function(j) diag(seG[j, ]^2, q))
  Fk <- numeric(q)
  for (k in seq_len(q)) {
    ok <- setdiff(seq_len(q), k)
    Xo <- G[, ok, drop = FALSE]; y <- G[, k]
    s2 <- function(delta) vapply(seq_len(p), function(j) {
      S <- cov_j[[j]]
      as.numeric(S[k, k] - 2 * crossprod(delta, S[ok, k]) +
                   crossprod(delta, S[ok, ok] %*% delta)) }, numeric(1))
    delta <- as.numeric(stats::coef(stats::lm(y ~ Xo - 1)))   # OLS start
    for (it in seq_len(100)) {
      v <- s2(delta); v[v <= 0] <- min(v[v > 0], na.rm = TRUE) / 10
      w <- 1 / v
      dn <- as.numeric(solve(crossprod(Xo, Xo * w), crossprod(Xo, y * w)))
      if (max(abs(dn - delta)) < 1e-8) { delta <- dn; break }
      delta <- dn
    }
    v <- s2(delta); v[v <= 0] <- min(v[v > 0], na.rm = TRUE) / 10
    Qk <- sum((y - as.numeric(Xo %*% delta))^2 / v)
    Fk[k] <- Qk / (p - q + 1)
  }
  Fk
}

## ---- build a per-SNP SNP-exposure covariance object from a phenotypic
## correlation matrix of the exposures (the standard same-sample approximation).
## `R_pheno` is the q x q observational correlation among the exposures, from
## individual-level data or a reference; the per-SNP covariance is then
##   Sigma_{gamma,j} = diag(seG_j) R_pheno diag(seG_j).
## Returns a value to assign to dat$gencov. With R_pheno = I this is the
## independence (marginal-SE) default. Supplying it makes conditional_F and the
## residual variance account for correlated exposure estimates, as required when
## the exposure GWAS share a sample.
build_gencov <- function(R_pheno) {
  R_pheno <- as.matrix(R_pheno)
  stopifnot(nrow(R_pheno) == ncol(R_pheno),
            isTRUE(all.equal(diag(R_pheno), rep(1, nrow(R_pheno)))))
  R_pheno
}

## variance inflation factor per exposure (collinearity of the first stage)
vif_exposure <- function(dat, idx = NULL) {
  if (is.null(idx)) idx <- screen_instruments(dat)
  G <- dat$G[idx, , drop = FALSE]; q <- ncol(G)
  vif <- numeric(q)
  for (k in seq_len(q)) {
    if (q == 1) { vif[k] <- 1; next }
    r2 <- summary(lm(G[, k] ~ G[, -k, drop = FALSE]))$r.squared
    vif[k] <- 1 / (1 - r2)
  }
  vif
}

## partial condition number kappa_k controlling projection length (Cor. 1)
## kappa_k = sqrt( [ (G_V^T V_{G,VV}^{-1} G_V)^{-1} ]_kk ), valid set V (or Shat)
partial_condition_number <- function(dat, idx = NULL) {
  if (is.null(idx)) idx <- screen_instruments(dat)
  G <- dat$G[idx, , drop = FALSE]; seG <- dat$seG[idx, , drop = FALSE]
  W <- 1 / rowMeans(seG^2)                       # diagonal first-stage precision
  info <- crossprod(G, G * W)
  sqrt(diag(solve(info)))
}

## count of instruments individually strongly relevant to exposure k
strong_count <- function(dat, idx = NULL, alpha = 0.05) {
  if (is.null(idx)) idx <- screen_instruments(dat)
  za <- z_alpha(length(idx), alpha)
  G <- dat$G[idx, , drop = FALSE]; seG <- dat$seG[idx, , drop = FALSE]
  colSums(abs(G / seG) > za)
}

## ---- Oracle test-inversion benchmark (supplement sanity check) ----
## Knows the true valid set `dat$valid` and inverts the residual test on it
## (no Wald, no thresholding). Returns whether beta* is accepted; averaging
## over replications gives the nominal-coverage check in the Supplementary
## Material. Contrast with mv_oracle(), which is the Wald (MV-IVW) oracle.
oracle_test_inversion <- function(dat, beta = dat$beta_star, alpha = 0.05) {
  V  <- dat$valid
  rs <- residual_stats(beta, dat, idx = V)
  za <- z_alpha(length(V), alpha)
  list(accept = max(rs$z) <= za, za = za, max_z = max(rs$z))
}

## Coverage of the oracle test-inversion benchmark over `n_rep` draws of a
## population generator `gen` (e.g. SETTINGS[[\"S2\"]]$gen).
oracle_test_inversion_cover <- function(gen, n_rep = 300, alpha = 0.05) {
  acc <- logical(n_rep)
  for (r in seq_len(n_rep)) {
    dat <- draw_sample(gen())
    acc[r] <- oracle_test_inversion(dat, alpha = alpha)$accept
  }
  mean(acc)
}

## ---- Scalar (q=1) searching baseline, applied one exposure at a time ----
## Inverts the univariate residual majority test on a single exposure k,
## ignoring the other exposures (the method of Guo 2023 applied marginally).
## Demonstrates in simulation why the multivariable construction is needed:
## under cross-loading the region is typically empty and never covers beta*_k.
univariate_searching <- function(dat, k, alpha = 0.05, ng = 1201, Rmax = 3) {
  g <- dat$G[, k]; seg <- dat$seG[, k]; rel <- abs(g) / seg
  idx <- which(rel >= sqrt(log(dat$n)))
  if (length(idx) < 3)
    idx <- order(rel, decreasing = TRUE)[seq_len(max(3, floor(length(g) / 2)))]
  nS <- length(idx); za <- z_alpha(nS, alpha); Tmaj <- floor(nS / 2)
  gk <- dat$G[idx, k]; segk <- dat$seG[idx, k]
  Gam <- dat$Gamma[idx]; seGam <- dat$seGamma[idx]
  grid <- seq(-Rmax, Rmax, length.out = ng)
  acc <- vapply(grid, function(b) {
    r <- Gam - b * gk; se <- sqrt(seGam^2 + (segk * b)^2)
    sum(abs(r) / se > za) <= (Tmaj - 1)
  }, logical(1))
  b0 <- dat$beta_star[k]
  if (!any(acc)) return(list(cover = FALSE, empty = TRUE, lower = NA, upper = NA))
  bs <- grid[acc]
  list(cover = b0 >= min(bs) && b0 <= max(bs), empty = FALSE,
       lower = min(bs), upper = max(bs))
}

## Marginal (univariate) IVW point estimate for coordinate k; its bias for
## the direct effect quantifies the blend of direct and indirect paths.
marginal_ivw <- function(dat, k) {
  g <- dat$G[, k]; w <- 1 / dat$seGamma^2
  sum(w * g * dat$Gamma) / sum(w * g^2)
}

## Scalar-baseline summary over a generator: coverage of beta*_k, empty-region
## fraction, and mean marginal-IVW bias (Table "Necessity of the
## multivariable construction").
scalar_baseline <- function(gen, k = 1, n_rep = 300, alpha = 0.05) {
  cov <- emp <- bias <- numeric(n_rep)
  for (r in seq_len(n_rep)) {
    dat <- draw_sample(gen())
    u <- univariate_searching(dat, k, alpha = alpha)
    cov[r] <- as.numeric(u$cover); emp[r] <- as.numeric(u$empty)
    bias[r] <- marginal_ivw(dat, k) - dat$beta_star[k]
  }
  c(coverage = mean(cov), empty_fraction = mean(emp), marg_ivw_bias = mean(bias))
}
