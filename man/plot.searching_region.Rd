\name{plot.searching_region}
\alias{plot.searching_region}
\title{Plot a two-exposure searching region}
\usage{
\method{plot}{searching_region}(x, ...)
}
\arguments{
\item{x}{A \code{"searching_region"} from \code{\link{searching_region}}.}
\item{...}{Passed to \code{\link[graphics]{plot}}.}
}
\value{\code{x}, invisibly.}
\description{
Draws the bisection-section rectangle of a two-exposure region with the
reference center marked. For q != 2 a message is shown.
}
\examples{
plot(searching_region(simulate_mvmr(seed = 1)))
}
