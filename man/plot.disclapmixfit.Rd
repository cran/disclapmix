\name{plot.disclapmixfit}
\alias{plot.disclapmixfit}
\title{Plot a disclapmixfit}
\description{Plot a \code{\link{disclapmixfit}} object.}
\usage{
  \method{plot}{disclapmixfit}(x, which = 1L, \dots)  
}
\arguments{
  \item{x}{a \code{\link{disclapmixfit}} object, usually from a result of a call to \code{disclapmix}.}
  \item{which}{What plot to make. 1L = clusters and their distances.}
  \item{...}{not used}
}
\value{  
  A data frame with discrete Laplace distributions for each cluster and locus. Side effect: A plot.
}
\author{
  Mikkel Meyer Andersen <mikl@math.aau.dk> and Poul Svante Eriksen
}

\seealso{
\code{\link{disclapmix}}
\code{\link{disclapmixfit}}
\code{\link{predict.disclapmixfit}}
\code{\link{print.disclapmixfit}}
\code{\link{simulate.disclapmixfit}}
\code{\link{summary.disclapmixfit}}
\code{\link{haplotype_diversity}}
\code{\link{clusterdist}}
}
\examples{
data(danes)
db <- as.matrix(danes[rep(1:nrow(danes), danes$n), 1:(ncol(danes)-1)])
fit <- disclapmix(db, clusters = 4L)
plot(fit)
}
\keyword{plot}

