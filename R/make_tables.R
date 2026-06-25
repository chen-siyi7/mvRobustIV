## ============================================================
##  make_tables.R
##
##  Emits LaTeX for the simulation tables directly from
##  results/sim_results.rds, matching the manuscript labels:
##    tables/tab_settings.tex      (Table: simulation settings)
##    tables/tab_cov_objects.tex   (Table: coverage by object)
##    tables/tab_sim_methods.tex   (Table: principal-method comparison)
##    tables/tab_fullsim.tex       (Supplement: full method comparison)
##  \input these into the manuscript, or copy the bodies in.
##
##  Usage:  source("R/make_tables.R"); make_all_tables()
## ============================================================

source(file.path("R", "simulation_settings.R"))

.fmt <- function(x, d = 3) ifelse(is.na(x), "n.r.", formatC(x, format = "f", digits = d))
.sgn <- function(x, d = 3) ifelse(is.na(x), "n.r.",
                                  sprintf("%s%s", ifelse(x >= 0, "+", "-"),
                                          formatC(abs(x), format = "f", digits = d)))

## ---- Table: simulation settings (static, from SETTINGS + DGP_PARS) ----
write_settings_table <- function(path = "tables/tab_settings.tex") {
  rows <- list(
    c("S1","Majority","2","30","21","50/50","$\\pi^{*}_{j}\\sim$U$(-\\tau,\\tau)$"),
    c("S2","Correlated exposures","2","30","21","50/50","as S1, corr.\\ $\\approx0.6$"),
    c("S3","Three traits","3","50","35","100/80","as S1"),
    c("S4","Multimodal plurality","2","30","17","50/50","three clusters"),
    c("S5","Directional pleiotropy","2","30","21","50/50","all $\\pi^{*}_{j}>0$"),
    c("S6","Many weak instruments","2","100","70","20/20","very weak first stage")
  )
  body <- vapply(rows, function(r) paste0(
    r[1]," & ",r[2]," & $",r[3],"$ & $",r[4],"$ & $",r[5],"$ & $",r[6],"$ & ",r[7]," \\\\"),
    character(1))
  tex <- c(
    "\\begin{table}[ht]","\\centering",
    "\\caption{Simulation settings. Sample sizes $n_{\\mathrm{exp}}/n_{\\mathrm{out}}$",
    "  in thousands; S2, S4 and S5 use the S1 sample sizes. $\\abs{\\cV}$ is the",
    "  number of valid instruments. S6 lies outside the strong-instrument",
    "  regime of the theory and is a diagnostic counter-example.}",
    "\\label{tab:settings}","\\small","\\begin{tabular}{llccccl}","\\toprule",
    " & Regime & $q$ & $p_{z}$ & $\\abs{\\cV}$ & $n_{\\mathrm{exp}}/n_{\\mathrm{out}}$ & Invalidity \\\\",
    "\\midrule", body, "\\bottomrule","\\end{tabular}","\\end{table}")
  writeLines(tex, path); message("wrote ", path)
}

## ---- Table: coverage by object ----
write_cov_objects_table <- function(res, path = "tables/tab_cov_objects.tex") {
  s <- res$summary
  row <- function(nm) {
    z <- s[[nm]]
    sprintf("%s & $%s$ & $%s$ & $%s$ & $%s$ & $%s$ \\\\", nm,
            .fmt(z$section_cov), .fmt(z$projection_cov), .fmt(z$joint_cov),
            .fmt(z$proj_len1), .fmt(z$proj_unbd))
  }
  body <- vapply(names(s), row, character(1))
  tex <- c(
    "\\begin{table}[ht]","\\centering",
    "\\caption{Searching-region coverage by object: bisection section, exact",
    "  coordinate projection, and joint region (the quantity",
    "  Theorem~\\ref{thm:searching} controls), with the exact-projection",
    sprintf("  first-coordinate length and unbounded fraction. $%d$ replications;", res$n_rep),
    "  Monte Carlo standard errors for the coverage entries are at most $0.022$",
    "  and below $0.010$ near the nominal level.}",
    "\\label{tab:cov_objects}","\\small","\\begin{tabular}{lcccrr}","\\toprule",
    " & \\multicolumn{3}{c}{Coverage (nominal $0.95$)} & \\multicolumn{2}{c}{Projection} \\\\",
    "\\cmidrule(lr){2-4}\\cmidrule(lr){5-6}",
    "Setting & Section & Projection & Joint & Len$_{1}$ & Unbd \\\\",
    "\\midrule", body, "\\bottomrule","\\end{tabular}","\\end{table}")
  writeLines(tex, path); message("wrote ", path)
}

## ---- Table: principal-method comparison (coverage + length panels) ----
write_sim_methods_table <- function(res, path = "tables/tab_sim_methods.tex") {
  s <- res$summary; setn <- names(s)
  methods <- c("MV-IVW","MV-PRESSO","TSHT","Searching","Oracle")  # principal rows shown in main text
  getcell <- function(setting, method, field) {
    df <- s[[setting]]$methods; v <- df[df$method == method, field]
    if (!length(v)) return(NA_real_) else v
  }
  panel <- function(field, bold = "Searching") {
    sapply(methods, function(m) {
      cells <- vapply(setn, function(st) {
        v <- getcell(st, m, field)
        cell <- if (field == "cov") .fmt(v) else .fmt(v)
        if (m == bold && !is.na(v)) paste0("\\mathbf{", cell, "}") else cell
      }, character(1))
      wrap <- if (m == bold) "\\textbf{Searching}" else m
      paste0(wrap, " & ", paste(sprintf("$%s$", cells), collapse = " & "), " \\\\")
    })
  }
  tex <- c(
    "\\begin{table}[ht]","\\centering",
    "\\caption{Coverage (nominal $0.95$) and median first-coordinate interval",
    "  length by method and setting, for the principal competitors.",
    "  \\textbf{Searching} entries are the bisection section; the joint and",
    "  exact-projection coverage are in Table~\\ref{tab:cov_objects}. The full",
    "  comparison is in the Supplementary Material; ``n.r.'' marks a method not",
    "  run for that setting, and Oracle is the infeasible known-valid-set benchmark.}",
    "\\label{tab:sim_methods}","\\small","\\begin{tabular}{lcccccc}","\\toprule",
    paste0(" & ", paste(setn, collapse = " & "), " \\\\"), "\\midrule",
    "\\multicolumn{7}{l}{\\emph{Coverage}}\\\\", panel("cov"),
    "\\midrule",
    "\\multicolumn{7}{l}{\\emph{Median length, first coordinate}}\\\\", panel("len1"),
    "\\bottomrule","\\end{tabular}","\\end{table}")
  writeLines(tex, path); message("wrote ", path)
}

## ---- Supplement: full method comparison (implemented methods) ----
write_fullsim_table <- function(res, path = "tables/tab_fullsim.tex") {
  s <- res$summary; setn <- names(s)
  body <- c()
  for (st in setn) {
    df <- s[[st]]$methods
    body <- c(body, sprintf("\\multicolumn{5}{l}{\\emph{%s}}\\\\", res$labels[[st]]))
    for (i in seq_len(nrow(df))) {
      m <- df$method[i]
      nm <- if (m == "Searching") "\\textbf{Searching}" else m
      body <- c(body, sprintf("%s & %s & %s & %s & %s \\\\",
        nm, .fmt(df$cov[i]), .fmt(df$vol[i]), .sgn(df$bias1[i]), .fmt(df$len1[i])))
    }
    body <- c(body, "\\midrule")
  }
  body <- head(body, -1)
  tex <- c(
    "\\begin{table}[ht]","\\centering",
    "\\caption{Full nine-method comparison: coverage (nominal $0.95$), median",
    "  joint volume, bias of $\\widehat\\beta_{1}$, and median first-coordinate",
    "  length. Searching coverage is the bisection section; its joint and",
    "  projection coverage are in Table~\\ref{tab:cov_objects}. Sampling is the",
    "  resampling refinement, and Oracle is the infeasible known-valid-set benchmark.}",
    "\\label{tab:fullsim}","\\footnotesize","\\begin{tabular}{lrrrr}","\\toprule",
    "Method & Cov & Vol & Bias$_{1}$ & Len$_{1}$ \\\\","\\midrule",
    body, "\\bottomrule","\\end{tabular}","\\end{table}")
  writeLines(tex, path); message("wrote ", path)
}

make_all_tables <- function(rds = "results/sim_results.rds") {
  dir.create("tables", showWarnings = FALSE)
  res <- readRDS(rds)
  write_settings_table()
  write_cov_objects_table(res)
  write_sim_methods_table(res)
  write_fullsim_table(res)
  invisible(res)
}

if (sys.nframe() == 0L) make_all_tables()
