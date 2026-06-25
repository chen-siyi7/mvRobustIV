## ============================================================
##  plot_simulations.R
##
##  Clean simulation figures, base graphics only (grayscale-safe).
##  Design choices that keep the panels readable:
##    - methods are DODGED within each setting (no connecting lines, which
##      would falsely imply continuity across categorical settings);
##    - a single legend sits BELOW the panels, never over the data;
##    - the proposed region is drawn as a larger filled marker.
##    make_coverage_plot()    -> figures/sim_coverage.pdf
##    make_trajectory_plot()  -> figures/sim_trajectory.pdf
## ============================================================

if (!exists("searching_region")) source(file.path("R", "mvrobustiv_methods.R"))
if (!exists("SETTINGS"))        source(file.path("R", "simulation_settings.R"))

## five principal methods, matching the main-text table, with markers
.METH <- c("MV-IVW", "MV-PRESSO", "TSHT", "Searching", "Oracle")
.PCH  <- c("MV-IVW" = 1, "MV-PRESSO" = 0, "TSHT" = 2, "Searching" = 19, "Oracle" = 5)
.CEX  <- c("MV-IVW" = 1.1, "MV-PRESSO" = 1.1, "TSHT" = 1.1, "Searching" = 1.5, "Oracle" = 1.1)

.cov_of <- function(z, m) if (m == "Searching") z$joint_cov else {
  v <- z$methods[z$methods$method == m, "cov"]; if (length(v)) v else NA_real_ }
.len_of <- function(z, m) {
  v <- z$methods[z$methods$method == m, "len1"]; if (length(v)) v else NA_real_ }

make_coverage_plot <- function(rds = "results/sim_results.rds",
                               out = "figures/sim_coverage.pdf") {
  dir.create(dirname(out), showWarnings = FALSE)
  s <- readRDS(rds)$summary; setn <- names(s)
  cov <- sapply(setn, function(st) vapply(.METH, function(m) .cov_of(s[[st]], m), numeric(1)))
  len <- sapply(setn, function(st) vapply(.METH, function(m) .len_of(s[[st]], m), numeric(1)))
  rownames(cov) <- rownames(len) <- .METH
  ns <- length(setn); nm <- length(.METH)
  off <- seq(-0.28, 0.28, length.out = nm)

  pdf(out, width = 9.6, height = 4.6)
  op <- par(mfrow = c(1, 2), oma = c(4.2, 0.5, 0.5, 0.5),
            mar = c(3.0, 4.2, 2.2, 1.0), cex = 0.95, las = 1, tcl = -0.3,
            mgp = c(2.4, 0.6, 0))

  ## ---- left: coverage by setting, methods dodged ----
  plot(NA, xlim = c(0.5, ns + 0.5), ylim = c(0, 1.06), xaxt = "n",
       xlab = "", ylab = "Empirical coverage", main = "Coverage by setting")
  abline(h = seq(0, 1, 0.2), col = "grey93")
  abline(v = seq(1.5, ns - 0.5), col = "grey93")
  abline(h = 0.95, lty = 2, col = "grey45")
  axis(1, at = seq_len(ns), labels = setn, tick = FALSE)
  for (j in seq_len(nm)) {
    m <- .METH[j]
    points(seq_len(ns) + off[j], cov[m, ], pch = .PCH[m], cex = .CEX[m],
           lwd = 1.5, col = "grey15")
  }
  mtext("Setting", side = 1, line = 1.8, cex = 0.95)

  ## ---- right: coverage-width trade-off ----
  cv <- as.vector(cov); ln <- as.vector(len)
  mp <- rep(.METH, times = ns); ok <- is.finite(cv) & is.finite(ln) & ln > 0
  plot(ln[ok], cv[ok], log = "x", pch = .PCH[mp[ok]], cex = .CEX[mp[ok]],
       lwd = 1.5, col = "grey15", ylim = c(0, 1.06),
       xlab = "Median first-coordinate length (log)",
       ylab = "Empirical coverage", main = "Coverage vs. width")
  abline(h = 0.95, lty = 2, col = "grey45")

  ## ---- shared legend, below both panels ----
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0),
      new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  legend("bottom", legend = .METH, pch = .PCH[.METH], pt.cex = .CEX[.METH],
         pt.lwd = 1.5, col = "grey15", horiz = TRUE, bty = "n",
         x.intersp = 0.7, text.width = NA, cex = 0.92)
  par(op); dev.off(); message("wrote ", out)
}

## ---- coverage trajectory in S2 as the sample size grows ----
make_trajectory_plot <- function(out = "figures/sim_trajectory.pdf",
                                 n_grid = c(5e4, 1e5, 1.5e5, 2e5),
                                 n_rep = 100, seed = 2024, M_samp = 60) {
  dir.create(dirname(out), showWarnings = FALSE)
  set.seed(seed)
  cov <- matrix(NA_real_, length(n_grid), 4,
                dimnames = list(NULL, c("section", "projection", "joint", "sampling")))
  for (i in seq_along(n_grid)) {
    n <- n_grid[i]; acc <- matrix(NA_real_, n_rep, 4)
    for (r in seq_len(n_rep)) {
      pop <- pop_S2(); pop$n_exp <- n; pop$n_out <- n; pop$n <- n
      dat <- draw_sample(pop); idx <- screen_instruments(dat)
      if (length(idx) < ncol(dat$G) + 1) next
      b <- dat$beta_star; sr <- searching_region(dat, idx = idx)
      s1 <- sr$section[1, ]; p1 <- sr$projection[1, ]
      acc[r, 1] <- !is.na(s1["lower"]) && b[1] >= s1["lower"] && b[1] <= s1["upper"]
      acc[r, 2] <- !is.na(p1["lower"]) && b[1] >= p1["lower"] && b[1] <= p1["upper"]
      acc[r, 3] <- sr$in_joint(b)
      acc[r, 4] <- sampling_region(dat, idx = idx, M = M_samp, beta_eval = b)$cover
    }
    cov[i, ] <- colMeans(acc, na.rm = TRUE)
  }
  obj <- colnames(cov); pch <- c(1, 2, 19, 5); lty <- c(3, 3, 1, 3)
  pdf(out, width = 6.6, height = 4.4)
  op <- par(oma = c(3.2, 0.5, 0.5, 0.5), mar = c(3.4, 4.2, 2.2, 1.0),
            cex = 0.95, las = 1, tcl = -0.3, mgp = c(2.4, 0.6, 0))
  ylo <- min(0.9, min(cov, na.rm = TRUE) - 0.02)
  plot(NA, xlim = range(n_grid), ylim = c(ylo, 1.01), xaxt = "n",
       xlab = "", ylab = "Empirical coverage",
       main = "S2 coverage as the sample size grows")
  abline(h = 0.95, lty = 2, col = "grey45")
  axis(1, at = n_grid, labels = formatC(n_grid, format = "d", big.mark = ","))
  for (k in seq_along(obj)) {
    lines(n_grid, cov[, k], lty = lty[k], col = "grey35", lwd = 1.3)
    points(n_grid, cov[, k], pch = pch[k], cex = 1.3, lwd = 1.5, col = "grey15")
  }
  mtext(expression(n[exp] == n[out]), side = 1, line = 2.0, cex = 0.95)
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  legend("bottom", legend = obj, pch = pch, lty = lty, col = "grey20",
         horiz = TRUE, bty = "n", cex = 0.92, x.intersp = 0.7)
  par(op); dev.off(); message("wrote ", out)
}

if (sys.nframe() == 0L) { make_coverage_plot(); make_trajectory_plot() }
