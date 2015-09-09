\name{clusterdist}
\alias{clusterdist}
\title{Calculate distance between clusters}
\description{\code{clusterdist} calculates the distance between each pair of clusters. The distance measure is based on a symmetric Kullback-Leibler divergence.}
\usage{
  clusterdist(fit, \dots)
}
\arguments{
  \item{fit}{A \code{\link{disclapmixfit}} object.}
  \item{...}{Not used}
}
\value{
  A distance matrix
}
\author{
  Mikkel Meyer Andersen <mikl@math.aau.dk> and Poul Svante Eriksen
}

\seealso{
\code{\link{disclapmix-package}}
\code{\link{disclapmix}}
\code{\link{disclapmixfit}}
\code{\link{clusterprob}}
\code{\link{predict.disclapmixfit}}
\code{\link{print.disclapmixfit}}
\code{\link{summary.disclapmixfit}}
\code{\link{simulate.disclapmixfit}}
%\code{\link{haplotype_diversity}}
\code{\link[disclap:disclap-package]{disclap}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{distance}
\keyword{clusters}

