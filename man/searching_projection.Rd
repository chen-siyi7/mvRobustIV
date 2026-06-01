\name{searching_projection}
\alias{searching_projection}
\title{Approximate coordinate projection of the searching region}
\usage{
searching_projection(object, k, alpha = 0.05, center = NULL, S = NULL,
                     grid = 41, span = 5, max_abs = 50, tol = 1e-3)
}
\arguments{
\item{object}{An \code{\link{mr_input}} object.}
\item{k}{Index of the target exposure.}
\item{alpha}{Nominal level.}
\item{center}{Feasible reference point.}
\item{S}{Optional selected set.}
\item{grid}{Number of grid points per non-target coordinate.}
\item{span}{Half-width of the grid around the center.}
\item{max_abs}{Outer bracket for the bisection.}
\item{tol}{Bisection tolerance.}
}
\value{Numeric vector \code{c(lower, upper)}; infinite endpoints mark an
unbounded projection.}
\description{
Grid approximation to the exact coordinate projection (the range of beta_k over
the region). Wider than the section and feasible only for small q.
}
\examples{
searching_projection(simulate_mvmr(seed = 1), k = 1)
}
