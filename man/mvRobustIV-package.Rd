\name{mvRobustIV-package}
\alias{mvRobustIV-package}
\alias{mvRobustIV}
\docType{package}
\title{Uniform Inference for Multiple Causal Effects with Invalid Instruments}
\description{
Confidence regions for a vector of causal effects in multivariable
instrumental-variable and Mendelian randomization analyses when some
instruments may be invalid and the valid set is unknown. The package implements
the searching confidence region obtained by inverting scalar residual-validity
tests under a majority rule, with relevance screening, boundedness and strength
diagnostics, an iterated hard-thresholding valid-set estimate, and an optional
resampling refinement.
}
\details{
Main entry points: \code{\link{mr_input}} builds a summary-data object,
\code{\link{searching_region}} computes the region and its coordinate sections,
\code{\link{iv_diagnostics}} reports relevance, boundedness and strength
diagnostics, and \code{\link{simulate_mvmr}} generates example data.
}
\references{
Chen, S. (2026). Uniform inference for multiple causal effects with invalid
instruments.
}
\author{Siyi Chen}
\keyword{package}
