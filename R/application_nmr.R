## ============================================================
##  application_nmr.R   (all data from OpenGWAS)
##
##  Real-data analysis for the manuscript application: highly correlated
##  Nightingale NMR lipoprotein subfractions on coronary artery disease
##  (CAD). Exposures and outcome are both pulled from OpenGWAS, harmonized,
##  clumped, screened, diagnosed, and run through the estimators and the
##  searching/sampling regions. Produces the three LaTeX tables and a
##  paste-ready summary block for the manuscript placeholders (\PH).
##
##  It reuses the reference method functions in mvrobustiv_methods.R;
##  nothing here re-derives a result.
##
##  ONE-TIME SETUP
##    install.packages("remotes")
##    remotes::install_github("MRCIEU/ieugwasr")
##    remotes::install_github("MRCIEU/TwoSampleMR")
##  Get a free token at https://api.opengwas.io/ and paste it just below.
##
##  TO RUN (from the package root, the directory that contains R/):
##    source("R/application_nmr.R")
##  The run block at the bottom executes automatically and prints the
##  summary; `res` holds the full results.
## ============================================================

## ============================================================
##  >>>>>>>>>>>>>>>>  PUT YOUR OPENGWAS TOKEN HERE  <<<<<<<<<<<<<<<<
##  Paste the JWT string from https://api.opengwas.io/ between the quotes.
##  Treat it like a password: do not share it or commit it to git.
## ============================================================
OPENGWAS_JWT <- "PASTE_YOUR_OPENGWAS_TOKEN_HERE"

## ---- what to analyze (the manuscript's final, separable triple) ----
EXPOSURE_IDS    <- c("met-d-IDL_C", "met-d-LDL_C", "met-d-HDL_C")
EXPOSURE_LABELS <- c("IDL-C", "LDL-C", "HDL-C")
OUTCOME_ID      <- "ieu-a-7"   # CARDIoGRAMplusC4D 2015 (Nikpay), CAD

## ---- exposure-estimate covariance for same-sample exposures (REVIEW-CRITICAL) ----
## The Nightingale NMR subfractions are measured on the SAME UK Biobank
## individuals, so the per-SNP SNP-exposure estimates gamma_hat_{j.} are
## correlated and the off-diagonals of Sigma_{gamma,j} are NOT zero. Supply the
## q x q phenotypic correlation of the exposures (rows/cols in EXPOSURE_IDS
## order), from UK Biobank individual-level data or an LDSC genetic-correlation
## proxy. The per-SNP covariance is then
##   Sigma_{gamma,j} = diag(seG_j) %*% R_pheno %*% diag(seG_j),
## which the conditional F, the residual-contrast variance, and the region all
## use. Leaving this NULL invokes the working-independence approximation
## (diagonal Sigma); with correlated exposures that OVERSTATES conditional
## strength, and the reported conditional F and region are then explicitly
## conditional on independence (stated as such in the output and the paper).
EXPOSURE_PHENO_COR <- NULL
## example (replace with the UK Biobank values for IDL-C, LDL-C, HDL-C):
## EXPOSURE_PHENO_COR <- matrix(c(1.00, 0.92, 0.10,
##                                0.92, 1.00, 0.05,
##                                0.10, 0.05, 1.00), 3, 3, byrow = TRUE)

## ---- data provenance recorded with the run (printed and returned) ----
## Edit if you change ancestry, build, or clumping. Palindromic SNPs are
## resolved by TwoSampleMR harmonisation at strictness 2 (ambiguous palindromes
## with intermediate allele frequencies are dropped). Two-sample independence is
## assumed; see overlap_sensitivity() for the partial-overlap stress check.
PROVENANCE <- list(
  exposure_source = "UK Biobank Nightingale NMR, OpenGWAS 'met-d' batch (Julkunen et al. 2023)",
  outcome_source  = "CARDIoGRAMplusC4D 2015 (Nikpay et al.), OpenGWAS ieu-a-7",
  ancestry        = "European",
  genome_build    = "GRCh37/hg19 (OpenGWAS harmonised)",
  clump_r2        = 0.001, clump_kb = 10000,
  harmonise_strictness = 2,
  sample_overlap  = "exposure and outcome samples not confirmed disjoint; two-sample independence assumed, overlap stress-tested")

source(file.path("R", "mvrobustiv_methods.R"))

## register the token for ieugwasr / TwoSampleMR (this session)
## strip accidental surrounding whitespace, newlines, or a "Bearer " prefix
OPENGWAS_JWT <- trimws(gsub("[\r\n]", "", OPENGWAS_JWT))
OPENGWAS_JWT <- sub("^[Bb]earer[[:space:]]+", "", OPENGWAS_JWT)
if (OPENGWAS_JWT != "PASTE_YOUR_OPENGWAS_TOKEN_HERE" && nchar(OPENGWAS_JWT) > 0)
  Sys.setenv(OPENGWAS_JWT = OPENGWAS_JWT)

## ------------------------------------------------------------
##  OpenGWAS helpers: token check, id verification, sample sizes
## ------------------------------------------------------------
check_opengwas <- function() {
  if (!requireNamespace("ieugwasr", quietly = TRUE))
    stop("install ieugwasr (remotes::install_github('MRCIEU/ieugwasr')).")
  tok <- ieugwasr::get_opengwas_jwt()
  if (nchar(tok) == 0)
    stop("No OpenGWAS token. Paste your JWT into OPENGWAS_JWT at the top of ",
         "this script (or set OPENGWAS_JWT in ~/.Renviron), then re-source.")
  tryCatch(invisible(ieugwasr::user()),     # 401 here => token bad/expired
    error = function(e) stop(
      "OpenGWAS rejected the token (", conditionMessage(e), ").\n",
      "  Most likely it has expired or was copied with extra/missing characters.\n",
      "  Generate a fresh token at https://api.opengwas.io/profile/, paste the\n",
      "  whole string (no spaces) into OPENGWAS_JWT, and re-source. The token\n",
      "  should be a few hundred characters with exactly two '.' separators.",
      call. = FALSE))
}

## confirm ids exist and report their reported sample sizes (ncase+ncontrol
## for a binary outcome, sample_size for quantitative exposures)
opengwas_n <- function(ids) {
  gi <- ieugwasr::gwasinfo(ids)
  if (!nrow(gi)) stop("none of these ids were found: ", paste(ids, collapse = ", "))
  n <- gi$sample_size
  if (is.null(n) || all(is.na(n))) {
    nc <- suppressWarnings(as.numeric(gi$ncase)); nk <- suppressWarnings(as.numeric(gi$ncontrol))
    n <- ifelse(is.na(nc), NA, nc) + ifelse(is.na(nk), 0, nk)
  }
  stats::setNames(n, gi$id)
}

## ------------------------------------------------------------
##  OpenGWAS acquisition via TwoSampleMR
##  Verify the exact met-d ids first with
##    ao <- TwoSampleMR::available_outcomes()
##    ao[grepl("^met-d", ao$id) & grepl("LDL|IDL", ao$trait),
##       c("id","trait","sample_size")]
##  Subfraction-cholesterol ids look like "met-d-L_LDL_C","met-d-S_LDL_C",
##  "met-d-M_LDL_C","met-d-IDL_C","met-d-L_VLDL_C", etc.
## ------------------------------------------------------------
fetch_opengwas <- function(exposure_ids, outcome_id = "ieu-a-7",
                           clump_r2 = 0.001, clump_kb = 10000) {
  if (!requireNamespace("TwoSampleMR", quietly = TRUE))
    stop("install TwoSampleMR (remotes::install_github('MRCIEU/TwoSampleMR')).")
  check_opengwas()
  ## instruments associated with ANY exposure, LD-clumped jointly: the
  ## standard MVMR instrument pool. Returns the SNP-exposure matrix already
  ## clumped at the given thresholds.
  exp_dat <- TwoSampleMR::mv_extract_exposures(
    exposure_ids, clump_r2 = clump_r2, clump_kb = clump_kb)
  out_dat <- TwoSampleMR::extract_outcome_data(
    snps = unique(exp_dat$SNP), outcomes = outcome_id)
  mvh <- TwoSampleMR::mv_harmonise_data(exp_dat, out_dat, harmonise_strictness = 2)
  ## mvh$exposure_beta : p x q (columns named by exposure id)
  ## mvh$exposure_se   : p x q
  ## mvh$outcome_beta  : length p (log-odds for a binary outcome)
  ## mvh$outcome_se    : length p
  list(G = mvh$exposure_beta, seG = mvh$exposure_se,
       Gamma = as.numeric(mvh$outcome_beta), seGamma = as.numeric(mvh$outcome_se),
       id_order = colnames(mvh$exposure_beta))
}

## assemble the summary-statistic `dat` object the methods expect.
## `pheno_cor`, when given, is the q x q phenotypic correlation of the
## exposures in `cols` order; it is attached as dat$gencov so the conditional
## F, residual variance, and region account for same-sample exposure
## correlation. NULL keeps the working-independence (diagonal) covariance.
build_dat_nmr <- function(raw, cols, labels, n_exp, n_out, pheno_cor = NULL) {
  miss <- setdiff(cols, colnames(raw$G))
  if (length(miss))
    stop("requested exposures not found in the harmonized data: ",
         paste(miss, collapse = ", "),
         "\n  available columns: ", paste(colnames(raw$G), collapse = ", "))
  G   <- raw$G[, cols, drop = FALSE]
  seG <- raw$seG[, cols, drop = FALSE]
  colnames(G) <- colnames(seG) <- labels
  keep <- stats::complete.cases(G, seG, raw$Gamma, raw$seGamma)
  dat <- list(Gamma = raw$Gamma[keep], seGamma = raw$seGamma[keep],
              G = G[keep, , drop = FALSE], seG = seG[keep, , drop = FALSE],
              n = round((n_exp + n_out) / 2))
  if (!is.null(pheno_cor)) {
    R <- as.matrix(pheno_cor)
    if (nrow(R) != length(cols) || ncol(R) != length(cols))
      stop("pheno_cor must be ", length(cols), " x ", length(cols),
           " in the order of `cols`.")
    dat$gencov <- build_gencov(R)        # validated q x q correlation
  }
  dat
}

## standardized first-stage condition number (conditioning of G^T W G)
condition_number <- function(dat, idx = NULL) {
  if (is.null(idx)) idx <- screen_instruments(dat)
  G <- dat$G[idx, , drop = FALSE]; seG <- dat$seG[idx, , drop = FALSE]
  W <- 1 / rowMeans(seG^2)
  info <- crossprod(G, G * W)
  ev <- eigen(info, symmetric = TRUE, only.values = TRUE)$values
  sqrt(max(ev) / min(ev))
}

## ---- formatting helpers (match the manuscript table style) ----
.fmt  <- function(x, d = 2) formatC(x, format = "f", digits = d)
.intv <- function(lo, hi, unb = 0, d = 2) {
  if (isTRUE(unb == 1) || is.na(lo)) "unbounded"
  else sprintf("$(%s,\\,%s)$", .fmt(lo, d), .fmt(hi, d))
}
.pt   <- function(est, lo, hi, d = 2)         # "0.39 (0.28,0.50)"
  sprintf("$%s\\ (%s,%s)$", .fmt(est, d), .fmt(lo, d), .fmt(hi, d))

## ---- one estimator row (point estimate + Wald CI per coordinate) ----
.est_row <- function(name, fit, q)
  sprintf("%s & %s \\\\", name,
    paste(vapply(seq_len(q),
      function(k) .pt(fit$beta[k], fit$ci[k, 1], fit$ci[k, 2]),
      character(1)), collapse = " & "))

## ---- first-stage diagnostics table + returned values ----
diagnostics_nmr <- function(dat, labels, path = "tables/tab_app_diag.tex",
                            label = "tab:app_diag", model = "three-subfraction") {
  idx  <- screen_instruments(dat)
  Fk   <- conditional_F(dat, idx)
  vif  <- vif_exposure(dat, idx)
  sc   <- strong_count(dat, idx)
  kap  <- condition_number(dat, idx)
  Tmaj <- floor(length(idx) / 2)
  bnd  <- ifelse(sc >= Tmaj, "yes", "no")
  basis <- if (is.null(dat$gencov)) "under working independence across exposures"
           else "with the supplied exposure covariance"
  rows <- vapply(seq_along(labels), function(k)
    sprintf("%s & $%s$ & $%s$ & $%d$ & %s \\\\", labels[k],
            .fmt(Fk[k], 1), .fmt(vif[k], 2), sc[k], bnd[k]), character(1))
  tex <- c("\\begin{table}[ht]", "\\centering",
    sprintf("\\caption{First-stage diagnostics for the %s model. Conditional $F$ is the Sanderson--Windmeijer multivariable conditional $F$, computed from summary statistics %s. ``Strong'' counts instruments whose marginal association with the subfraction exceeds the Bonferroni cutoff $z_{\\alpha}$; a coordinate's searching section is bounded only if this count reaches the majority threshold $\\lfloor|\\hat{\\mathcal{S}}|/2\\rfloor=%d$.}", model, basis, Tmaj),
    sprintf("\\label{%s}", label), "\\small", "\\begin{tabular}{lrrrc}", "\\toprule",
    "Subfraction & Conditional $F$ & VIF & Strong & Section bounded \\\\", "\\midrule",
    rows, "\\bottomrule", "\\end{tabular}", "\\end{table}")
  writeLines(tex, path); message("wrote ", path)
  list(idx = idx, nS = length(idx), Tmaj = Tmaj, za = z_alpha(length(idx)),
       F = Fk, vif = vif, strong = sc, kappa = kap, bounded = bnd,
       cov_basis = if (is.null(dat$gencov)) "working-independence" else "supplied-covariance")
}

## ---- full estimate table (six estimators + searching section/projection) ----
estimate_nmr <- function(dat, labels, caption, label, path) {
  idx <- screen_instruments(dat); q <- length(labels)
  fits <- list(
    "MV-IVW"    = mv_ivw(dat, idx),
    "MV-Egger"  = mv_egger(dat, idx),
    "MV-Median" = mv_median(dat, idx = idx),
    "MV-PRESSO" = mv_presso(dat, idx = idx),
    "MV-Lasso"  = mv_lasso(dat, idx = idx),
    "TSHT"      = mv_tsht(dat, idx = idx))
  sr <- searching_region(dat, idx = idx, certify = TRUE)
  sec_row <- sprintf("Searching (section) & %s \\\\",
    paste(vapply(seq_len(q), function(k)
      .intv(sr$section[k, "lower"], sr$section[k, "upper"], sr$section[k, "unbounded"]),
      character(1)), collapse = " & "))
  pr_row <- sprintf("Searching (projection) & %s \\\\",
    paste(vapply(seq_len(q), function(k)
      .intv(sr$projection[k, "lower"], sr$projection[k, "upper"], sr$projection[k, "unbounded"]),
      character(1)), collapse = " & "))
  est_rows <- vapply(names(fits), function(nm) .est_row(nm, fits[[nm]], q), character(1))
  tex <- c("\\begin{table}[ht]", "\\centering",
    sprintf("\\caption{%s}", caption), sprintf("\\label{%s}", label), "\\small",
    sprintf("\\begin{tabular}{l%s}", paste(rep("c", q), collapse = "")), "\\toprule",
    paste0("Method & ", paste(labels, collapse = " & "), " \\\\"), "\\midrule",
    est_rows, "\\midrule", sec_row, pr_row, "\\bottomrule",
    "\\end{tabular}", "\\end{table}")
  writeLines(tex, path); message("wrote ", path)
  list(fits = fits, sr = sr, unbounded = which(sr$projection[, "unbounded"] == 1),
       certified = sr$certified, Rmax = sr$Rmax, proj_doubled = sr$proj_doubled)
}

## ---- driver (OpenGWAS only) ----
run_nmr_application <- function(exposure_ids, exposure_labels,
                                outcome_id = "ieu-a-7",
                                n_exp = NULL, n_out = NULL,
                                pheno_cor = EXPOSURE_PHENO_COR) {
  stopifnot(length(exposure_ids) == 3, length(exposure_labels) == 3)
  if (!is.null(pheno_cor)) {
    pheno_cor <- as.matrix(pheno_cor)
    stopifnot(nrow(pheno_cor) == 3, ncol(pheno_cor) == 3)
  }
  dir.create("tables", showWarnings = FALSE)
  raw <- fetch_opengwas(exposure_ids, outcome_id)
  if (is.null(n_exp)) n_exp <- round(stats::median(opengwas_n(exposure_ids), na.rm = TRUE))
  if (is.null(n_out)) n_out <- round(unname(opengwas_n(outcome_id)))
  cols <- exposure_ids

  ## ----- full three-subfraction model -----
  dat3  <- build_dat_nmr(raw, cols, exposure_labels, n_exp, n_out, pheno_cor = pheno_cor)
  diag3 <- diagnostics_nmr(dat3, exposure_labels)
  est3  <- estimate_nmr(dat3, exposure_labels,
    "Three-subfraction model. Point estimators report estimate (interval); the searching region reports intervals only, with ``unbounded'' marking a coordinate it does not bound. Effects are log-odds ratios per standard-deviation increment.",
    "tab:app_est", "tables/tab_app_est.tex")

  ## ----- reduced two-subfraction model: drop the unbounded coordinate -----
  drop_k <- if (length(est3$unbounded)) est3$unbounded[1] else which.min(diag3$strong)
  keep2  <- setdiff(seq_along(exposure_labels), drop_k)
  pc2    <- if (is.null(pheno_cor)) NULL else pheno_cor[keep2, keep2, drop = FALSE]
  dat2   <- build_dat_nmr(raw, cols[keep2], exposure_labels[keep2], n_exp, n_out, pheno_cor = pc2)
  diag2  <- diagnostics_nmr(dat2, exposure_labels[keep2],
    path = "tables/tab_app_two_diag.tex", label = "tab:app_two_diag",
    model = "reduced two-subfraction")
  est2   <- estimate_nmr(dat2, exposure_labels[keep2],
    "Reduced two-subfraction model. Layout as in Table~\\ref{tab:app_est}; ``unbounded'' marks a coordinate the searching region does not bound. Effects are log-odds ratios per standard-deviation increment.",
    "tab:app_two", "tables/tab_app_two.tex")

  ## ----- paste-ready summary for the LaTeX \PH cells -----
  cat("\n================  PASTE-READY SUMMARY  ================\n")
  cat(sprintf("n_exp = %d ; n_out = %d ; |S-hat| = %d ; z_alpha = %.2f ; Tmaj = %d\n",
              n_exp, n_out, diag3$nS, diag3$za, diag3$Tmaj))
  cat(sprintf("standardized first-stage condition number = %.2f\n", diag3$kappa))
  cat(sprintf("exposure-covariance basis: %s%s\n", diag3$cov_basis,
      if (identical(diag3$cov_basis, "working-independence"))
        "  (set EXPOSURE_PHENO_COR for same-sample exposures)" else ""))
  cat(sprintf("provenance: %s ; outcome %s ; ancestry %s ; clump r2=%.3f kb=%d ; %s\n",
      PROVENANCE$exposure_source, PROVENANCE$outcome_source, PROVENANCE$ancestry,
      PROVENANCE$clump_r2, PROVENANCE$clump_kb, PROVENANCE$sample_overlap))
  cat("\n-- three-subfraction diagnostics (Table tab:app_diag) --\n")
  for (k in seq_along(exposure_labels))
    cat(sprintf("  %-14s  F=%.1f  VIF=%.2f  strong=%d  section_bounded=%s\n",
        exposure_labels[k], diag3$F[k], diag3$vif[k], diag3$strong[k], diag3$bounded[k]))
  cat(sprintf("\nunbounded coordinate in full model: %s\n",
      if (length(est3$unbounded)) paste(exposure_labels[est3$unbounded], collapse = ", ")
      else "none (all projections finite)"))
  cat(sprintf("projection certification (full model): search radius Rmax=%.0f, re-checked at 2*Rmax with halved cell;\n  per-coordinate stable = %s\n",
      est3$Rmax, paste(sprintf("%s:%s", exposure_labels,
        ifelse(est3$certified, "yes", "NO")), collapse = ", ")))
  cat(sprintf("reduced model retains: %s\n", paste(exposure_labels[keep2], collapse = ", ")))
  cat(sprintf("\n-- reduced two-subfraction diagnostics (Table tab:app_two_diag) --\n"))
  cat(sprintf("   |S-hat|_2 = %d ; Tmaj_2 = %d ; condition number = %.2f\n",
              diag2$nS, diag2$Tmaj, diag2$kappa))
  for (k in seq_along(keep2))
    cat(sprintf("  %-14s  F=%.1f  VIF=%.2f  strong=%d  section_bounded=%s\n",
        exposure_labels[keep2][k], diag2$F[k], diag2$vif[k], diag2$strong[k], diag2$bounded[k]))
  cat(sprintf("reduced projection certification: Rmax=%.0f ; per-coordinate stable = %s\n",
      est2$Rmax, paste(sprintf("%s:%s", exposure_labels[keep2],
        ifelse(est2$certified, "yes", "NO")), collapse = ", ")))
  cat("\nFull tables written to tables/tab_app_diag.tex, tab_app_est.tex, tab_app_two.tex\n")
  cat("Supplement counts: per-subfraction strong = (",
      paste(diag3$strong, collapse = ", "), "), |S-hat| = ", diag3$nS,
      ", threshold = ", diag3$Tmaj, "\n", sep = "")
  cat("======================================================\n")

  invisible(list(dat3 = dat3, dat2 = dat2, diag3 = diag3, diag2 = diag2,
                 est3 = est3, est2 = est2, provenance = PROVENANCE,
                 pheno_cor = pheno_cor))
}

## ---- sample-overlap sensitivity (two-sample independence is maintained) ----
## Two-sample independence is the theoretical assumption. This stress-tests
## partial overlap: for each assumed overlap correlation rho between the
## SNP-outcome and SNP-exposure estimates, the per-SNP cross-covariance
## cov(hatGamma_j, hatgamma_jk) = rho * seGamma_j * seG_jk is injected (via
## dat$cross) and the searching projection recomputed. A stable bounded/
## unbounded pattern across rho shows the conclusions do not hinge on exact
## non-overlap. Returns a matrix of unbounded flags (1 = unbounded).
overlap_sensitivity <- function(dat, labels, rhos = c(0, 0.1, 0.2), alpha = 0.05) {
  q <- ncol(dat$G)
  rows <- lapply(rhos, function(rho) {
    d <- dat
    d$cross <- if (rho == 0) NULL
               else (dat$seGamma %o% rep(1, q)) * dat$seG * rho
    searching_region(d, alpha = alpha)$projection[, "unbounded"]
  })
  M <- do.call(rbind, rows)
  rownames(M) <- sprintf("rho=%.2f", rhos); colnames(M) <- labels
  cat("\nSample-overlap sensitivity (1 = unbounded projection):\n"); print(M)
  invisible(M)
}

## ---- quick screen: conditioning diagnostics for a candidate set ----
## Cheap check before a full run: fetches the candidate exposures + outcome
## and prints |S-hat|, conditional F, VIF, and the condition number, so you
## can pick a set that is correlated but still separable (aim for conditional
## F well above ~10 and VIF elevated but not in the dozens).
screen_exposures <- function(exposure_ids, exposure_labels = exposure_ids,
                             outcome_id = "ieu-a-7") {
  raw   <- fetch_opengwas(exposure_ids, outcome_id)
  n_exp <- round(stats::median(opengwas_n(exposure_ids), na.rm = TRUE))
  n_out <- round(unname(opengwas_n(outcome_id)))
  dat   <- build_dat_nmr(raw, exposure_ids, exposure_labels, n_exp, n_out)
  idx   <- screen_instruments(dat)
  Fk <- conditional_F(dat, idx); vif <- vif_exposure(dat, idx); sc <- strong_count(dat, idx)
  cat(sprintf("|S-hat| = %d ; condition number = %.1f ; Tmaj = %d\n",
              length(idx), condition_number(dat, idx), floor(length(idx) / 2)))
  for (k in seq_along(exposure_labels))
    cat(sprintf("  %-16s  F=%.1f  VIF=%.1f  strong=%d\n",
        exposure_labels[k], Fk[k], vif[k], sc[k]))
  invisible(dat)
}

## ============================================================
##  RUN  (executes when you source this file)
##  To load the functions WITHOUT running (e.g. to screen candidate
##  exposure sets cheaply), set  SKIP_RUN <- TRUE  before sourcing, then call
##  screen_exposures(c("met-d-IDL_C","met-d-LDL_C","met-d-HDL_C"),
##                   c("IDL-C","LDL-C","HDL-C"))
## ============================================================
if (!exists("SKIP_RUN") || !isTRUE(SKIP_RUN)) {
  if (OPENGWAS_JWT == "PASTE_YOUR_OPENGWAS_TOKEN_HERE" &&
      nchar(Sys.getenv("OPENGWAS_JWT")) == 0)
    stop("Paste your OpenGWAS token into OPENGWAS_JWT at the top of this script.")

  res <- run_nmr_application(
    exposure_ids    = EXPOSURE_IDS,
    exposure_labels = EXPOSURE_LABELS,
    outcome_id      = OUTCOME_ID)
}
