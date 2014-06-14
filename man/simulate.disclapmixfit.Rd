\name{simulate.disclapmixfit}
\alias{simulate.disclapmixfit}
\title{Simulate from a disclapmixfit}
\description{Simulate from a \code{\link{disclapmixfit}} object.}
\usage{
  \method{simulate}{disclapmixfit}(object, nsim = 1L, seed = NULL, \dots)
}
\arguments{
  \item{object}{a \code{\link{disclapmixfit}} object, usually from a result of a call to \code{disclapmix}.}
  \item{nsim}{number of haplotypes to generate.}
  \item{seed}{not used}
  \item{...}{not used}
}
\value{
  A matrix where the rows correspond to the simulated haplotypes.
}
\author{
  Mikkel Meyer Andersen <mikl@math.aau.dk> and Poul Svante Eriksen
}

\seealso{
\code{\link{disclapmix}}
\code{\link{disclapmixfit}}
\code{\link{predict.disclapmixfit}}
\code{\link{print.disclapmixfit}}
\code{\link{plot.disclapmixfit}}
\code{\link{summary.disclapmixfit}}
\code{\link{haplotype_diversity}}
\code{\link{clusterdist}}
}
\keyword{print}

