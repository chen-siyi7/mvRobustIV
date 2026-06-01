\name{mv_ivw}
\alias{mv_ivw}
\alias{print.mv_ivw}
\title{Multivariable inverse-variance-weighted estimate}
\usage{
mv_ivw(object, S = NULL, alpha = 0.05)
}
\arguments{
\item{object}{An \code{\link{mr_input}} object.}
\item{S}{Optional selected set.}
\item{alpha}{Level for the Wald intervals.}
}
\value{An object of class \code{"mv_ivw"} with \code{estimate}, \code{se} and
\code{ci}.}
\description{
Inverse-variance-weighted multivariable estimate over the selected set, used as
a reference center and a conventional comparator. Assumes all selected
instruments are valid.
}
\examples{
mv_ivw(simulate_mvmr(seed = 1))
}
