\name{simulate_mvmr}
\alias{simulate_mvmr}
\title{Simulate two-sample multivariable summary data with invalid instruments}
\usage{
simulate_mvmr(q = 2, p = 30, n_valid = 21, beta = NULL,
              strength = c(0.02, 0.06), n_X = 5e4, n_Y = 5e4, rho = 0,
              invalidity = c("symmetric", "plurality", "directional"),
              tau = 0.2, seed = NULL)
}
\arguments{
\item{q}{Number of exposures.}
\item{p}{Number of instruments.}
\item{n_valid}{Number of valid instruments.}
\item{beta}{True causal vector of length q; a default is used when NULL.}
\item{strength}{Range \code{c(low, high)} for first-stage magnitudes.}
\item{n_X, n_Y}{Exposure- and outcome-sample sizes.}
\item{rho}{Correlation injected between exposures' first-stage vectors.}
\item{invalidity}{Direct-effect pattern.}
\item{tau}{Direct-effect scale.}
\item{seed}{Optional random seed.}
}
\value{An \code{\link{mr_input}} object with attributes \code{true_beta} and
\code{valid_set}.}
\description{
Generates SNP-exposure and SNP-outcome summary statistics under a linear model
with a known valid set, for examples, tests and reproduction.
}
\examples{
sim <- simulate_mvmr(q = 2, p = 30, n_valid = 21, seed = 1)
attr(sim, "true_beta")
}
