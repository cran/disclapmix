\name{haplotype_diversity}
\alias{haplotype_diversity}
\title{Calculate haplotype diversity from a disclapmixfit}
\description{Calculate haplotype diversity from a \code{\link{disclapmixfit}} object. The method is based on simulating a huge database that approximates the population.}
\usage{
  haplotype_diversity(object, nsim = 1e4L)
}
\arguments{
  \item{object}{a \code{\link{disclapmixfit}} object, usually from a result of a call to \code{disclapmix}.}
  \item{nsim}{number of haplotypes to generate for calculating the haplotype diversity.}
}
\value{
  The calculated haplotype diversity.
}
\author{
  Mikkel Meyer Andersen <mikl@math.aau.dk> and Poul Svante Eriksen
}

\seealso{
\code{\link{disclapmix}}
\code{\link{disclapmixfit}}
\code{\link{predict.disclapmixfit}}
\code{\link{print.disclapmixfit}}
\code{\link{summary.disclapmixfit}}
\code{\link{simulate.disclapmixfit}}
\code{\link{haplotype_diversity}}
\code{\link{clusterdist}}
}
\keyword{print}

