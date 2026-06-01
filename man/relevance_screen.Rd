\name{relevance_screen}
\alias{relevance_screen}
\title{First-stage relevance screen}
\usage{
relevance_screen(object, n = NULL, chisq_level = 0.999)
}
\arguments{
\item{object}{An \code{\link{mr_input}} object.}
\item{n}{Sample size for the q*log(n) cutoff; defaults to \code{object$n_X}.}
\item{chisq_level}{Tail level for the fallback chi-squared(q) cutoff.}
}
\value{Integer vector of selected instrument indices.}
\description{
Selects the instrument set by the joint first-stage rule
n * gamma_j^T inverse(V) gamma_j >= q log(n).
}
\examples{
length(relevance_screen(simulate_mvmr(seed = 1)))
}
