\name{searching_section}
\alias{searching_section}
\title{Coordinate section of the searching region}
\usage{
searching_section(object, k, alpha = 0.05, center = NULL, S = NULL,
                  max_abs = 50, tol = 1e-3)
}
\arguments{
\item{object}{An \code{\link{mr_input}} object.}
\item{k}{Index of the exposure whose section is computed.}
\item{alpha}{Nominal level.}
\item{center}{Feasible reference point (length q); defaults to \code{\link{mv_ivw}}.}
\item{S}{Optional selected set.}
\item{max_abs}{Outer bracket for the bisection on each side.}
\item{tol}{Bisection tolerance.}
}
\value{Numeric vector \code{c(lower, upper)}; infinite endpoints mark an
unbounded coordinate.}
\description{
One-dimensional slice of the searching region along exposure \code{k}, computed
by bisection with the other coordinates held at \code{center}.
}
\examples{
searching_section(simulate_mvmr(seed = 1), k = 1)
}
