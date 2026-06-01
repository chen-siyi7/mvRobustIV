\name{iv_diagnostics}
\alias{iv_diagnostics}
\alias{print.iv_diagnostics}
\title{First-stage diagnostics}
\usage{
iv_diagnostics(object, alpha = 0.05, S = NULL)
}
\arguments{
\item{object}{An \code{\link{mr_input}} object.}
\item{alpha}{Nominal level.}
\item{S}{Optional selected set.}
}
\value{An object of class \code{"iv_diagnostics"} with \code{n_instruments},
\code{per_exposure}, \code{majority_threshold}, \code{surplus} and
\code{joint_bounded}.}
\description{
Collects the diagnostics governing whether the searching region is bounded and
well conditioned: marginal-relevance counts with the boundedness rule, an
approximate conditional F, variance inflation factors, partial condition
numbers and the hard-thresholding surplus. The conditional F is read in one
direction only.
}
\examples{
iv_diagnostics(simulate_mvmr(q = 3, p = 50, seed = 2))
}
