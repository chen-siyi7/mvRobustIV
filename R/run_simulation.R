## ============================================================
##  run_simulation.R
##
##  Runs the S1-S6 study (seed 2024, 500 replications by default) and saves
##  per-replication metrics plus aggregates to results/sim_results.rds.
##  These feed make_tables.R (Tables for coverage-by-object and the
##  method comparison) and plot_simulations.R (the coverage figure).
##
##  Usage:
##    source("R/run_simulation.R")
##    res <- run_all_settings(n_rep = 500, seed = 2024)     # full
##    res <- run_all_settings(n_rep = 50,  quick = TRUE)    # fast check
## ============================================================

source(file.path("R", "mvrobustiv_methods.R"))
source(file.path("R", "simulation_settings.R"))

## per-coordinate CI coverage and joint coverage for a competitor fit
.cov_from_ci <- function(ci, beta_star) {
  cov_k <- (beta_star >= ci[, 1]) & (beta_star <= ci[, 2])
  list(joint = all(cov_k), cov1 = unname(cov_k[1]),
       len1 = unname(ci[1, 2] - ci[1, 1]),
       vol  = unname(prod(ci[, 2] - ci[, 1])))
}

## one replication of one setting -> named numeric vector of metrics
run_once <- function(pop_gen, ng_q2 = 121, ng_q3 = 41) {
  pop <- pop_gen()
  dat <- draw_sample(pop)
  bstar <- dat$beta_star; q <- length(bstar)
  idx <- screen_instruments(dat)
  if (length(idx) < q + 1) return(NULL)         # degenerate draw, skip

  ## proposed searching region (projection is reference-free; see methods)
  sr <- searching_region(dat, idx = idx)
  b0 <- sr$beta0
  sec1 <- sr$section[1, ]; pr1 <- sr$projection[1, ]
  sec_cov1  <- !is.na(sec1["lower"]) && bstar[1] >= sec1["lower"] && bstar[1] <= sec1["upper"]
  proj_cov1 <- !is.na(pr1["lower"])  && bstar[1] >= pr1["lower"]  && bstar[1] <= pr1["upper"]
  joint_cov <- sr$in_joint(bstar)
  sec_len1  <- unname(sec1["upper"] - sec1["lower"])
  proj_len1 <- unname(pr1["upper"]  - pr1["lower"])

  ## competitors (all nine methods of the full comparison)
  f_ivw <- mv_ivw(dat, idx);            ivw <- .cov_from_ci(f_ivw$ci, bstar)
  f_egg <- mv_egger(dat, idx);          egg <- .cov_from_ci(f_egg$ci, bstar)
  f_med <- mv_median(dat, idx = idx);   med <- .cov_from_ci(f_med$ci, bstar)
  f_prs <- mv_presso(dat, idx = idx);   prs <- .cov_from_ci(f_prs$ci, bstar)
  f_las <- mv_lasso(dat, idx = idx);    las <- .cov_from_ci(f_las$ci, bstar)
  f_tsh <- mv_tsht(dat, idx = idx);     tsh <- .cov_from_ci(f_tsh$ci, bstar)
  f_orc <- mv_oracle(dat, dat$valid);   orc <- .cov_from_ci(f_orc$ci, bstar)
  smp   <- sampling_region(dat, idx = idx, beta_eval = bstar)
  smp_len1 <- unname(smp$projection[1, "upper"] - smp$projection[1, "lower"])
  smp_mid1 <- unname((smp$projection[1, "lower"] + smp$projection[1, "upper"]) / 2)

  c(
    ## searching region, by object (Table: coverage by object)
    sear_sec_cov = as.numeric(sec_cov1),
    sear_proj_cov = as.numeric(proj_cov1),
    sear_joint_cov = as.numeric(joint_cov),
    sear_sec_len1 = sec_len1,
    sear_proj_len1 = proj_len1,
    sear_proj_unbd = unname(pr1["unbounded"]),
    sear_vol = sr$volume,
    sear_bias1 = sr$center[1] - bstar[1],
    ## method comparison cells (cov = joint coverage; len1 = first coordinate)
    ivw_cov = as.numeric(ivw$joint), ivw_len1 = ivw$len1, ivw_vol = ivw$vol, ivw_bias1 = f_ivw$beta[1]-bstar[1],
    egg_cov = as.numeric(egg$joint), egg_len1 = egg$len1, egg_vol = egg$vol, egg_bias1 = f_egg$beta[1]-bstar[1],
    med_cov = as.numeric(med$joint), med_len1 = med$len1, med_vol = med$vol, med_bias1 = f_med$beta[1]-bstar[1],
    prs_cov = as.numeric(prs$joint), prs_len1 = prs$len1, prs_vol = prs$vol, prs_bias1 = f_prs$beta[1]-bstar[1],
    las_cov = as.numeric(las$joint), las_len1 = las$len1, las_vol = las$vol, las_bias1 = f_las$beta[1]-bstar[1],
    tsh_cov = as.numeric(tsh$joint), tsh_len1 = tsh$len1, tsh_vol = tsh$vol, tsh_bias1 = f_tsh$beta[1]-bstar[1],
    orc_cov = as.numeric(orc$joint), orc_len1 = orc$len1, orc_vol = orc$vol, orc_bias1 = f_orc$beta[1]-bstar[1],
    smp_cov = as.numeric(smp$cover), smp_len1 = smp_len1, smp_vol = smp$volume, smp_bias1 = smp_mid1-bstar[1]
  )
}

## run all replications for one setting
run_setting <- function(pop_gen, n_rep = 500, ...) {
  out <- vector("list", n_rep)
  for (r in seq_len(n_rep)) out[[r]] <- run_once(pop_gen, ...)
  do.call(rbind, out[!vapply(out, is.null, logical(1))])
}

## Monte Carlo standard error of a coverage estimate
mc_se <- function(p, m) sqrt(p * (1 - p) / m)

## aggregate one setting's replication matrix into reported summaries
summarise_setting <- function(M) {
  m <- nrow(M)
  agg <- function(col, f = mean) f(M[, col], na.rm = TRUE)
  list(
    n_rep = m,
    ## coverage-by-object table
    section_cov   = agg("sear_sec_cov"),
    projection_cov= agg("sear_proj_cov"),
    joint_cov     = agg("sear_joint_cov"),
    proj_len1     = agg("sear_proj_len1", stats::median),
    proj_unbd     = agg("sear_proj_unbd"),
    cov_se        = mc_se(agg("sear_joint_cov"), m),
    ## method-comparison table: all nine methods (coverage, volume, bias, length)
    methods = local({
      keys <- c(ivw="MV-IVW", egg="MV-Egger", med="MV-Median", prs="MV-PRESSO",
                las="MV-Lasso", tsh="TSHT", smp="Sampling", sear="Searching", orc="Oracle")
      cov <- function(k) if (k == "sear") agg("sear_sec_cov") else agg(paste0(k, "_cov"))
      len <- function(k) if (k == "sear") agg("sear_sec_len1", stats::median)
                         else agg(paste0(k, "_len1"), stats::median)
      data.frame(
        method = unname(keys),
        cov    = vapply(names(keys), cov, numeric(1)),
        vol    = vapply(names(keys), function(k) agg(paste0(k, "_vol"), stats::median), numeric(1)),
        bias1  = vapply(names(keys), function(k) agg(paste0(k, "_bias1")), numeric(1)),
        len1   = vapply(names(keys), len, numeric(1)),
        row.names = NULL, stringsAsFactors = FALSE
      )
    })
  )
}

run_all_settings <- function(n_rep = 500, seed = 2024, quick = FALSE,
                             ng_q2 = 121, ng_q3 = 41, save = TRUE) {
  if (quick) { n_rep <- min(n_rep, 50); ng_q2 <- 81; ng_q3 <- 31 }
  set.seed(seed)
  raw <- list(); summ <- list()
  for (nm in names(SETTINGS)) {
    message(sprintf("Running %s (%d reps) ...", nm, n_rep))
    M <- run_setting(SETTINGS[[nm]]$gen, n_rep = n_rep, ng_q2 = ng_q2, ng_q3 = ng_q3)
    raw[[nm]]  <- M
    summ[[nm]] <- summarise_setting(M)
  }
  res <- list(raw = raw, summary = summ, labels = lapply(SETTINGS, `[[`, "label"),
              n_rep = n_rep, seed = seed)
  if (save) { dir.create("results", showWarnings = FALSE)
    saveRDS(res, file.path("results", "sim_results.rds")) }
  res
}

if (sys.nframe() == 0L) {
  res <- run_all_settings(n_rep = 500, seed = 2024)
  print(lapply(res$summary, function(s) s$methods))
}
