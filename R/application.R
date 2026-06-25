## ============================================================
##  application.R
##
##  Regenerates the lipid / coronary-heart-disease application: the
##  first-stage diagnostics table, the three-trait estimate table, and the
##  two-trait (LDL-C, HDL-C) table. Reads a harmonised, pre-clumped
##  two-sample summary-statistics file; the published numbers come from your
##  GLGC (exposure) and CARDIoGRAMplusC4D (outcome) extract.
##
##  Expected input CSV (one row per clumped instrument), columns:
##    SNP,
##    beta_LDL, se_LDL, beta_HDL, se_HDL, beta_TG, se_TG,   (SNP-exposure)
##    beta_CHD, se_CHD                                       (SNP-outcome, log-odds)
##  with n_exp and n_out passed to run_application().
##  Clumping (r^2<0.001, 10 Mb) is assumed done upstream with an LD panel.
##
##  Usage:
##    source("R/application.R")
##    run_application("data/lipids_chd.csv", n_exp = 188577, n_out = 184305)
## ============================================================

source(file.path("R", "mvrobustiv_methods.R"))

load_application_data <- function(path) {
  d <- utils::read.csv(path, stringsAsFactors = FALSE)
  req <- c("beta_LDL","se_LDL","beta_HDL","se_HDL","beta_TG","se_TG","beta_CHD","se_CHD")
  miss <- setdiff(req, names(d))
  if (length(miss)) stop("input is missing columns: ", paste(miss, collapse = ", "))
  d
}

## build the summary-statistic `dat` object for a chosen set of exposures
build_dat <- function(d, exposures, n_exp, n_out) {
  G   <- as.matrix(d[, paste0("beta_", exposures), drop = FALSE])
  seG <- as.matrix(d[, paste0("se_",   exposures), drop = FALSE])
  colnames(G) <- colnames(seG) <- exposures
  list(Gamma = d$beta_CHD, seGamma = d$se_CHD, G = G, seG = seG,
       n = round((n_exp + n_out) / 2))
}

.fmt  <- function(x, d = 2) formatC(x, format = "f", digits = d)
.intv <- function(lo, hi, unb, d = 2)
  if (isTRUE(unb == 1) || is.na(lo)) "unbounded" else sprintf("(%s,\\,%s)", .fmt(lo,d), .fmt(hi,d))

## ---- first-stage diagnostics table ----
diagnostics_table <- function(dat, exposures, path = "tables/tab_app_diag.tex") {
  idx <- screen_instruments(dat)
  Fk  <- conditional_F(dat, idx)
  vif <- vif_exposure(dat, idx)
  sc  <- strong_count(dat, idx)
  Tmaj <- floor(length(idx) / 2)
  bounded <- ifelse(sc >= Tmaj, "yes", "no")
  rows <- vapply(seq_along(exposures), function(k)
    sprintf("%s & $%s$ & $%s$ & $%d$ & %s \\\\", exposures[k],
            .fmt(Fk[k], 1), .fmt(vif[k], 2), sc[k], bounded[k]), character(1))
  tex <- c("\\begin{table}[ht]","\\centering",
    sprintf("\\caption{First-stage diagnostics. ``Strong'' counts instruments whose marginal association exceeds the Bonferroni cutoff; a coordinate's section is bounded only if this count reaches the majority threshold $%d$.}", Tmaj),
    "\\label{tab:app_diag}","\\small","\\begin{tabular}{lrrrc}","\\toprule",
    "Trait & Conditional $F$ & VIF & Strong & Section bounded \\\\","\\midrule",
    rows, "\\bottomrule","\\end{tabular}","\\end{table}")
  writeLines(tex, path); message("wrote ", path)
  invisible(list(idx = idx, F = Fk, vif = vif, strong = sc, Tmaj = Tmaj))
}

## ---- estimate table for a given model ----
estimate_table <- function(dat, exposures, label, path) {
  idx <- screen_instruments(dat)
  ivw <- mv_ivw(dat, idx); tsh <- mv_tsht(dat, idx = idx)
  sr  <- searching_region(dat, idx = idx)
  q <- length(exposures)
  ci_row <- function(name, ci) sprintf("%s & %s \\\\", name,
    paste(vapply(seq_len(q), function(k) .intv(ci[k,1], ci[k,2], 0), character(1)),
          collapse = " & "))
  sec_row <- sprintf("Searching (section) & %s \\\\",
    paste(vapply(seq_len(q), function(k) .intv(sr$section[k,"lower"], sr$section[k,"upper"], sr$section[k,"unbounded"]), character(1)), collapse = " & "))
  pr_row  <- sprintf("Searching (projection) & %s \\\\",
    paste(vapply(seq_len(q), function(k) .intv(sr$projection[k,"lower"], sr$projection[k,"upper"], sr$projection[k,"unbounded"]), character(1)), collapse = " & "))
  tex <- c("\\begin{table}[ht]","\\centering",
    sprintf("\\caption{%s. Intervals are %d\\%% confidence regions; the searching region reports section and projection, with ``unbounded'' where the projection is not finite.}", label, 95),
    sprintf("\\label{tab:app_%s}", gsub("[^a-z]", "", tolower(label))),
    "\\small", sprintf("\\begin{tabular}{l%s}", paste(rep("c", q), collapse = "")),
    "\\toprule",
    paste0("Method & ", paste(exposures, collapse = " & "), " \\\\"), "\\midrule",
    ci_row("MV-IVW", ivw$ci), ci_row("TSHT", tsh$ci), sec_row, pr_row,
    "\\bottomrule","\\end{tabular}","\\end{table}")
  writeLines(tex, path); message("wrote ", path)
}

run_application <- function(path, n_exp = 188577, n_out = 184305) {
  dir.create("tables", showWarnings = FALSE)
  d <- load_application_data(path)
  ## three-trait model
  dat3 <- build_dat(d, c("LDL","HDL","TG"), n_exp, n_out)
  diagnostics_table(dat3, c("LDL-C","HDL-C","Triglycerides"))
  estimate_table(dat3, c("LDL-C","HDL-C","Triglycerides"),
                 "Three-trait model", "tables/tab_app_est.tex")
  ## two-trait model in LDL-C and HDL-C
  dat2 <- build_dat(d, c("LDL","HDL"), n_exp, n_out)
  estimate_table(dat2, c("LDL-C","HDL-C"),
                 "Two-trait model in LDL-C and HDL-C", "tables/tab_app_two.tex")
  message("Application tables written to tables/.")
}

if (sys.nframe() == 0L) {
  message("Provide a harmonised CSV: run_application('data/lipids_chd.csv').")
}
