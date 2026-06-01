\name{summary.searching_region}
\alias{summary.searching_region}
\title{Summarize a searching region}
\usage{
\method{summary}{searching_region}(object, ...)
}
\arguments{
\item{object}{A \code{"searching_region"} from \code{\link{searching_region}}.}
\item{...}{Unused.}
}
\value{A data frame of coordinate sections and bounded status, invisibly.}
\description{Prints a short summary of a searching region.}
\examples{
summary(searching_region(simulate_mvmr(seed = 1)))
}
