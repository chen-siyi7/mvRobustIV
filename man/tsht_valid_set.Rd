\name{tsht_valid_set}
\alias{tsht_valid_set}
\title{Iterated hard-thresholding valid-set estimate}
\usage{
tsht_valid_set(object, alpha = 0.05, S = NULL, max_iter = 10)
}
\arguments{
\item{object}{An \code{\link{mr_input}} object.}
\item{alpha}{Nominal level.}
\item{S}{Optional selected set.}
\item{max_iter}{Maximum number of refitting iterations.}
}
\value{A list with \code{valid}, \code{estimate}, \code{surplus} and \code{S}.}
\description{
A practical two-stage hard-thresholding valid-set estimate used for the surplus
diagnostic and as a comparator. The uniform searching region does not use this
selection.
}
\examples{
tsht_valid_set(simulate_mvmr(seed = 1))$surplus
}
