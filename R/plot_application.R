## ============================================================
##  plot_application.R
##
##  Application figures: per-exposure interval ("forest") plots of the
##  estimate and confidence interval for each method, for the three-trait
##  and two-trait lipid/CHD models. Reads the same harmonised CSV as
##  application.R. Base graphics, grayscale-safe.
##    make_application_plots("data/lipids_chd.csv")
##      -> figures/app_forest_three.pdf, figures/app_forest_two.pdf
## ============================================================

if (!exists("searching_region")) source(file.path("R", "mvrobustiv_methods.R"))
if (!exists("build_dat"))        source(file.path("R", "application.R"))

## intervals (q x 2) for each method on a fitted data object
.app_intervals <- function(dat, idx) {
  list(
    "MV-IVW"    = mv_ivw(dat, idx)$ci,
    "MV-Egger"  = mv_egger(dat, idx)$ci,
    "MV-Median" = mv_median(dat, idx = idx)$ci,
    "MV-PRESSO" = mv_presso(dat, idx = idx)$ci,
    "MV-Lasso"  = mv_lasso(dat, idx = idx)$ci,
    "TSHT"      = mv_tsht(dat, idx = idx)$ci,
    "Searching" = searching_region(dat, idx = idx)$projection[, 1:2, drop = FALSE]
  )
}

.forest_one <- function(intervals, k, label, xlim = NULL) {
  meths <- names(intervals); nm <- length(meths)
  lo <- vapply(intervals, function(ci) ci[k, 1], numeric(1))
  hi <- vapply(intervals, function(ci) ci[k, 2], numeric(1))
  est <- (lo + hi) / 2
  fin <- is.finite(lo) & is.finite(hi)
  if (is.null(xlim)) {
    rng <- range(c(lo[fin], hi[fin], 0), na.rm = TRUE)
    xlim <- rng + c(-1, 1) * 0.08 * diff(rng)
  }
  plot(NA, xlim = xlim, ylim = c(0.5, nm + 0.5), yaxt = "n",
       xlab = "Effect on CHD (log odds)", ylab = "", main = label)
  axis(2, at = nm:1, labels = meths, las = 1, cex.axis = 0.85)
  abline(v = 0, col = "grey60", lty = 3)
  for (i in seq_len(nm)) {
    y <- nm - i + 1
    if (!fin[i]) {                       # unbounded: draw an arrowed bar
      arrows(xlim[1], y, xlim[2], y, length = 0.06, code = 3, angle = 30,
             col = "grey50", lty = 2)
      text(mean(xlim), y + 0.25, "unbounded", cex = 0.7, col = "grey40")
    } else {
      segments(lo[i], y, hi[i], y, lwd = 1.6)
      points(est[i], y, pch = if (meths[i] == "Searching") 19 else 16, cex = 1.1)
    }
  }
}

make_application_plots <- function(path, n_exp = 188577, n_out = 184305) {
  dir.create("figures", showWarnings = FALSE)
  d <- load_application_data(path)

  dat3 <- build_dat(d, c("LDL", "HDL", "TG"), n_exp, n_out)
  iv3  <- .app_intervals(dat3, screen_instruments(dat3))
  pdf("figures/app_forest_three.pdf", width = 9, height = 3.4)
  op <- par(mfrow = c(1, 3), mar = c(4, 6.5, 1.6, 0.6), cex = 0.9)
  .forest_one(iv3, 1, "LDL-C"); .forest_one(iv3, 2, "HDL-C"); .forest_one(iv3, 3, "Triglycerides")
  par(op); dev.off(); message("wrote figures/app_forest_three.pdf")

  dat2 <- build_dat(d, c("LDL", "HDL"), n_exp, n_out)
  iv2  <- .app_intervals(dat2, screen_instruments(dat2))
  pdf("figures/app_forest_two.pdf", width = 7, height = 3.4)
  op <- par(mfrow = c(1, 2), mar = c(4, 6.5, 1.6, 0.6), cex = 0.9)
  .forest_one(iv2, 1, "LDL-C"); .forest_one(iv2, 2, "HDL-C")
  par(op); dev.off(); message("wrote figures/app_forest_two.pdf")
}

if (sys.nframe() == 0L)
  message("Provide a harmonised CSV: make_application_plots('data/lipids_chd.csv').")
