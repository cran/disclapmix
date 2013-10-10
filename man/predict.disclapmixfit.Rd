\name{predict.disclapmixfit}
\alias{predict.disclapmixfit}
\title{Predict from a disclapmixfit}
\description{Is able to predict haplotype frequencies using a \code{\link{disclapmixfit}} object.}
\usage{
  \method{predict}{disclapmixfit}(object, newdata, \dots)
}
\arguments{
  \item{object}{a \code{\link{disclapmixfit}} object}
  \item{newdata}{the haplotypes in matrix format to estimate haplotype probabilities for}
  \item{...}{not used}
}
\author{
  Mikkel Meyer Andersen <mikl@math.aau.dk> and Poul Svante Eriksen
}

\seealso{
\code{\link{disclapmix}}
\code{\link{disclapmixfit}}
\code{\link{print.disclapmixfit}}
\code{\link{summary.disclapmixfit}}
\code{\link{clusterdist}}
\code{\link{clusterprob}}
}
\keyword{predict}

