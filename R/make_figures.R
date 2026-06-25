## ============================================================
##  make_figures.R  -- regenerate every figure
##
##  source("R/make_figures.R")
##  make_all_figures()                       # coverage + region (fast)
##  make_all_figures(trajectory = TRUE)      # also the S2 trajectory (slower)
##  make_application_plots("data/lipids_chd.csv")   # application forests (needs data)
## ============================================================

source(file.path("R", "plot_simulations.R"))
source(file.path("R", "plot_region.R"))
source(file.path("R", "plot_application.R"))

make_all_figures <- function(rds = "results/sim_results.rds", trajectory = FALSE) {
  dir.create("figures", showWarnings = FALSE)
  make_coverage_plot(rds)          # figures/sim_coverage.pdf
  make_region_plot()               # figures/region_geometry.pdf
  if (trajectory) make_trajectory_plot()   # figures/sim_trajectory.pdf (runs its own draws)
  invisible(TRUE)
}

if (sys.nframe() == 0L) make_all_figures(trajectory = TRUE)
