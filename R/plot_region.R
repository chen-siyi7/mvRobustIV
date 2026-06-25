## ============================================================
##  plot_region.R
##
##  Clean illustration of the searching region for one q = 2 data set:
##  the region drawn as a solid filled raster (not scattered squares),
##  its projection hull, the section through the centre, the truth, and
##  competitor estimates. Legend sits in the empty corner. Grayscale-safe.
##    make_region_plot() -> figures/region_geometry.pdf
## ============================================================

if (!exists("searching_region")) source(file.path("R", "mvrobustiv_methods.R"))
if (!exists("SETTINGS"))        source(file.path("R", "simulation_settings.R"))

make_region_plot <- function(out = "figures/region_geometry.pdf",
                             setting = "S1", seed = 11, npix = 260) {
  dir.create(dirname(out), showWarnings = FALSE)
  set.seed(seed)
  pop <- SETTINGS[[setting]]$gen(); dat <- draw_sample(pop)
  if (ncol(dat$G) != 2) stop("region plot is for q = 2 settings")
  bstar <- dat$beta_star; idx <- screen_instruments(dat)
  za <- z_alpha(length(idx)); Tmaj <- floor(length(idx) / 2)
  sr <- searching_region(dat, idx = idx); pr <- sr$projection

  padx <- 0.18 * (pr[1, 2] - pr[1, 1] + 0.1); pady <- 0.18 * (pr[2, 2] - pr[2, 1] + 0.1)
  ax <- seq(pr[1, 1] - padx, pr[1, 2] + padx, length.out = npix)
  ay <- seq(pr[2, 1] - pady, pr[2, 2] + pady, length.out = npix)
  G  <- as.matrix(expand.grid(b1 = ax, b2 = ay))
  acc <- accept_many(G, dat, idx, za, Tmaj)
  zmat <- matrix(as.numeric(acc), nrow = length(ax))   # rows = ax, cols = ay

  ivw <- mv_ivw(dat, idx)$beta
  tsh <- mv_tsht(dat, idx = idx)$beta
  med <- mv_median(dat, idx = idx)$beta

  pdf(out, width = 5.8, height = 5.4)
  op <- par(mar = c(4.0, 4.0, 2.0, 1.0), cex = 0.95, las = 1, tcl = -0.3,
            mgp = c(2.3, 0.6, 0))
  ## solid region as a raster, drawn first
  image(ax, ay, zmat, col = c("white", "grey86"), useRaster = TRUE,
        xlab = expression(beta[1]), ylab = expression(beta[2]),
        main = paste0("Searching region (", setting, ")"))
  box()
  rect(pr[1, 1], pr[2, 1], pr[1, 2], pr[2, 2], border = "grey30", lty = 2, lwd = 1.4)
  abline(v = sr$center[1], col = "grey65", lty = 3)
  abline(h = sr$center[2], col = "grey65", lty = 3)
  points(ivw[1], ivw[2], pch = 1, cex = 1.5, lwd = 1.6)
  points(tsh[1], tsh[2], pch = 2, cex = 1.5, lwd = 1.6)
  points(med[1], med[2], pch = 6, cex = 1.5, lwd = 1.6)
  points(bstar[1], bstar[2], pch = 8, cex = 1.9, lwd = 2.2)

  ## legend in whichever corner the region avoids (region descends L->R here)
  legend("topright",
         legend = c("region", "projection hull", "section", "truth",
                    "MV-IVW", "TSHT", "MV-Median"),
         pch = c(15, NA, NA, 8, 1, 2, 6),
         lty = c(NA, 2, 3, NA, NA, NA, NA),
         col = c("grey86", "grey30", "grey65", "black", "black", "black", "black"),
         pt.cex = c(1.6, 1, 1, 1.3, 1.2, 1.2, 1.2), lwd = 1.4,
         bty = "n", cex = 0.82, y.intersp = 1.05)
  par(op); dev.off(); message("wrote ", out)
}

if (sys.nframe() == 0L) make_region_plot()
