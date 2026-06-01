\name{searching_region}
\alias{searching_region}
\title{Searching confidence region}
\usage{
searching_region(object, alpha = 0.05, center = NULL, S = NULL,
                 max_abs = 50, tol = 1e-3)
}
\arguments{
\item{object}{An \code{\link{mr_input}} object.}
\item{alpha}{Nominal level; coverage is at least 1 - alpha.}
\item{center}{Feasible reference point; defaults to \code{\link{mv_ivw}}.}
\item{S}{Optional selected set.}
\item{max_abs}{Outer bracket for the bisection.}
\item{tol}{Bisection tolerance.}
}
\value{An object of class \code{"searching_region"} with \code{sections},
\code{bounded}, \code{center}, \code{alpha} and \code{S}.}
\description{
Computes the searching confidence region for the causal vector: the
per-coordinate bisection sections, the boundedness diagnostic and the reference
center. The region has coverage at least 1 - alpha, uniformly over the class of
data-generating processes in the reference.
}
\examples{
searching_region(simulate_mvmr(seed = 1))
}
