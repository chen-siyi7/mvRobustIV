#' @export
print.searching_region <- function(x, ...) {
  cat(sprintf("<searching_region> (alpha = %.3f, %d selected instruments)\n",
              x$alpha, length(x$S)))
  tab <- data.frame(
    exposure = x$exposure,
    lower = ifelse(is.infinite(x$sections[, 1]), "-Inf", formatC(x$sections[, 1], digits = 3, format = "f")),
    upper = ifelse(is.infinite(x$sections[, 2]),  "Inf", formatC(x$sections[, 2], digits = 3, format = "f")),
    bounded = x$bounded,
    stringsAsFactors = FALSE
  )
  print(tab, row.names = FALSE)
  if (!all(x$bounded)) {
    cat("note: unbounded coordinates are not identified under unknown invalidity.\n")
  }
  invisible(x)
}

#' Summarize a searching region
#'
#' @param object A `"searching_region"` from [searching_region()].
#' @param ... Unused.
#'
#' @return A data frame with the coordinate sections and their bounded status
#'   (invisibly), printed to the console.
#' @export
summary.searching_region <- function(object, ...) {
  cat("Searching confidence region\n")
  cat(sprintf("  level: %.1f%%   selected instruments: %d\n",
              100 * (1 - object$alpha), length(object$S)))
  cat(sprintf("  bounded coordinates: %d of %d\n",
              sum(object$bounded), object$q))
  tab <- data.frame(exposure = object$exposure,
                    lower = object$sections[, 1],
                    upper = object$sections[, 2],
                    bounded = object$bounded,
                    row.names = NULL)
  print(tab, row.names = FALSE)
  invisible(tab)
}

#' Plot a two-exposure searching region
#'
#' Draws the bisection-section rectangle of a two-exposure region with the
#' reference center marked. For \eqn{q \ne 2} a message is shown instead.
#'
#' @param x A `"searching_region"` from [searching_region()].
#' @param ... Passed to [graphics::plot()].
#'
#' @return `x`, invisibly.
#' @export
plot.searching_region <- function(x, ...) {
  if (x$q != 2) {
    message("Plotting is implemented for q = 2 only.")
    return(invisible(x))
  }
  s <- x$sections
  finite_box <- function(v, c0, pad = 1) {
    lo <- if (is.finite(v[1])) v[1] else c0 - pad
    hi <- if (is.finite(v[2])) v[2] else c0 + pad
    c(lo, hi)
  }
  xr <- finite_box(s[1, ], x$center[1])
  yr <- finite_box(s[2, ], x$center[2])
  graphics::plot(NA, xlim = xr, ylim = yr,
                 xlab = x$exposure[1], ylab = x$exposure[2],
                 main = "Searching region (coordinate sections)", ...)
  graphics::rect(s[1, 1], s[2, 1], s[1, 2], s[2, 2],
                 border = "#2c7fb8", lwd = 2,
                 col = grDevices::adjustcolor("#2c7fb8", alpha.f = 0.12))
  graphics::points(x$center[1], x$center[2], pch = 19, col = "#d95f02")
  graphics::abline(h = 0, v = 0, lty = 3, col = "grey60")
  invisible(x)
}
