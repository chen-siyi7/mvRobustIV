\name{sampling_region}
\alias{sampling_region}
\alias{print.sampling_region}
\title{Resampling refinement of the searching region}
\usage{
sampling_region(object, alpha = 0.05, S = NULL, M = 500, level = 0.95,
                seed = NULL)
}
\arguments{
\item{object}{An \code{\link{mr_input}} object.}
\item{alpha}{Nominal level.}
\item{S}{Optional selected set.}
\item{M}{Number of resamples.}
\item{level}{Envelope level per coordinate.}
\item{seed}{Optional random seed.}
}
\value{An object of class \code{"sampling_region"} with \code{intervals},
\code{centers}, \code{M} and \code{level}.}
\description{
A coverage-oriented sensitivity refinement secondary to
\code{\link{searching_region}}, which carries the uniform guarantee.
}
\examples{
sampling_region(simulate_mvmr(seed = 1), M = 100, seed = 2)$intervals
}
