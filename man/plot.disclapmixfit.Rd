\name{plot.disclapmixfit}
\alias{plot.disclapmixfit}
\title{Plot a disclapmixfit}
\description{Is able to plot a wide of different things for a \code{\link{disclapmixfit}} object.}
\usage{
  \method{plot}{disclapmixfit}(x, which = 1, \dots)
}
\arguments{
  \item{x}{a \code{\link{disclapmixfit}} object}
  \item{which}{a vector of which plots to do}
  \item{...}{not used}
}
\details{
  Which can be one or more of the following:
  1: Prior probabilities of belonging to a center (tau)
  2: v gain (see \code{\link{disclapmix}})
  3: marginal log likelihood
  4: full log likelihood
}
\author{
  Mikkel Meyer Andersen <mikl@math.aau.dk> and Poul Svante Eriksen
}

\seealso{
\code{\link{disclapmix}}
\code{\link{disclapmixfit}}
\code{\link{predict.disclapmixfit}}
\code{\link{extractMarginalBIC}}
\code{\link{extractFullBIC}}
}
\keyword{disclap}
\keyword{disclapmix}
\keyword{disclapmixfit}
\keyword{discrete Laplace}
\keyword{EM algorithm}

