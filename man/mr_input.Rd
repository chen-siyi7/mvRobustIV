\name{mr_input}
\alias{mr_input}
\alias{print.mr_input}
\title{Build a two-sample summary-data object}
\usage{
mr_input(beta_Y, se_Y, beta_X, se_X, n_X = NULL, n_Y = NULL,
         snp = NULL, exposure = NULL)
}
\arguments{
\item{beta_Y}{Numeric vector of length p: SNP-outcome associations.}
\item{se_Y}{Positive numeric vector of length p: their standard errors.}
\item{beta_X}{Numeric p-by-q matrix of SNP-exposure associations.}
\item{se_X}{Positive numeric p-by-q matrix of standard errors for \code{beta_X}.}
\item{n_X, n_Y}{Optional exposure- and outcome-sample sizes.}
\item{snp}{Optional character vector of SNP identifiers (length p).}
\item{exposure}{Optional character vector of exposure names (length q).}
}
\value{An object of class \code{"mr_input"}.}
\description{
Assembles and validates the summary statistics used throughout the package.
}
\examples{
sim <- simulate_mvmr(q = 2, p = 30, seed = 1)
sim
}
