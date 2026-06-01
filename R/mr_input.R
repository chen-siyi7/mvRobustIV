#' Build a two-sample summary-data object
#'
#' Assembles and validates the summary statistics used by every other function
#' in the package: SNP-outcome associations and a SNP-exposure association
#' matrix with their standard errors.
#'
#' @param beta_Y Numeric vector of length \eqn{p}: SNP-outcome associations.
#' @param se_Y Numeric vector of length \eqn{p}: their standard errors
#'   (strictly positive).
#' @param beta_X Numeric matrix \eqn{p \times q}: SNP-exposure associations, one
#'   column per exposure. A vector is treated as a single exposure.
#' @param se_X Numeric matrix \eqn{p \times q} of standard errors for `beta_X`
#'   (strictly positive), or a vector when \eqn{q = 1}.
#' @param n_X,n_Y Optional exposure- and outcome-sample sizes. `n_X` is used by
#'   [relevance_screen()] for the \eqn{q \log n} cutoff.
#' @param snp Optional character vector of SNP identifiers of length \eqn{p}.
#' @param exposure Optional character vector of exposure names of length
#'   \eqn{q}.
#'
#' @return An object of class `"mr_input"`: a list with entries `beta_Y`,
#'   `se_Y`, `beta_X`, `se_X`, `p`, `q`, `n_X`, `n_Y`, `snp`, `exposure`.
#'
#' @examples
#' sim <- simulate_mvmr(q = 2, p = 30, seed = 1)
#' sim
#'
#' @export
mr_input <- function(beta_Y, se_Y, beta_X, se_X,
                     n_X = NULL, n_Y = NULL, snp = NULL, exposure = NULL) {
  beta_X <- as.matrix(beta_X)
  se_X <- as.matrix(se_X)
  beta_Y <- as.numeric(beta_Y)
  se_Y <- as.numeric(se_Y)

  p <- length(beta_Y)
  q <- ncol(beta_X)

  if (nrow(beta_X) != p || length(se_Y) != p) {
    stop("beta_Y, se_Y and the rows of beta_X must all have length p.", call. = FALSE)
  }
  if (!all(dim(se_X) == dim(beta_X))) {
    stop("se_X must have the same dimensions as beta_X.", call. = FALSE)
  }
  if (anyNA(beta_Y) || anyNA(se_Y) || anyNA(beta_X) || anyNA(se_X)) {
    stop("Missing values are not allowed in the summary statistics.", call. = FALSE)
  }
  if (any(se_Y <= 0) || any(se_X <= 0)) {
    stop("Standard errors must be strictly positive.", call. = FALSE)
  }
  if (p <= q) {
    stop("The number of instruments p must exceed the number of exposures q.", call. = FALSE)
  }
  if (is.null(exposure)) exposure <- paste0("exposure", seq_len(q))
  if (length(exposure) != q) stop("exposure must have length q.", call. = FALSE)
  if (!is.null(snp) && length(snp) != p) stop("snp must have length p.", call. = FALSE)

  colnames(beta_X) <- exposure
  colnames(se_X) <- exposure

  structure(
    list(beta_Y = beta_Y, se_Y = se_Y, beta_X = beta_X, se_X = se_X,
         p = p, q = q, n_X = n_X, n_Y = n_Y, snp = snp, exposure = exposure),
    class = "mr_input"
  )
}

#' @export
print.mr_input <- function(x, ...) {
  cat("<mr_input>\n")
  cat(sprintf("  instruments p = %d, exposures q = %d\n", x$p, x$q))
  cat(sprintf("  exposures: %s\n", paste(x$exposure, collapse = ", ")))
  if (!is.null(x$n_X) || !is.null(x$n_Y)) {
    cat(sprintf("  n_X = %s, n_Y = %s\n",
                format(x$n_X %||% NA), format(x$n_Y %||% NA)))
  }
  invisible(x)
}
