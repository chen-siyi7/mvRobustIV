\name{searching_membership}
\alias{searching_membership}
\title{Membership test for the searching confidence region}
\usage{
searching_membership(object, beta, alpha = 0.05, S = NULL)
}
\arguments{
\item{object}{An \code{\link{mr_input}} object.}
\item{beta}{Numeric vector of length q.}
\item{alpha}{Nominal level; coverage is at least 1 - alpha.}
\item{S}{Optional selected set from \code{\link{relevance_screen}}.}
}
\value{A list with \code{in_region}, \code{n_flagged}, \code{max_flagged} and
\code{flagged}.}
\description{
Tests whether a candidate causal vector lies in the searching region by
counting flagged instruments against the majority threshold.
}
\examples{
sim <- simulate_mvmr(seed = 1)
searching_membership(sim, attr(sim, "true_beta"))$in_region
}
