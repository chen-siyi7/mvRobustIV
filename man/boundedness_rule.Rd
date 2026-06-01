\name{boundedness_rule}
\alias{boundedness_rule}
\title{Per-coordinate boundedness rule}
\usage{
boundedness_rule(object, alpha = 0.05, S = NULL)
}
\arguments{
\item{object}{An \code{\link{mr_input}} object.}
\item{alpha}{Nominal level.}
\item{S}{Optional selected set from \code{\link{relevance_screen}}.}
}
\value{A data frame with \code{exposure}, \code{strong}, \code{threshold} and
\code{bounded}.}
\description{
A coordinate section is bounded iff at least floor(|S|/2) selected instruments
are individually strongly associated with that exposure.
}
\examples{
boundedness_rule(simulate_mvmr(seed = 1))
}
