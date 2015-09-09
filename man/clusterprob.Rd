\name{clusterprob}
\alias{clusterprob}
\title{Cluster origin probabilities for haplotypes}
\description{\code{clusterprob} calculates the cluster origin probabilities for haplotypes.}
\usage{
clusterprob(fit, newdata, \dots)
}
\arguments{
  \item{fit}{A \code{\link{disclapmixfit}} object.}
  \item{newdata}{The haplotypes to predict the cluster origin probabilities for.}
  \item{...}{Not used}
}
\value{
  A matrix where the rows correspond to the rows in \code{newdata} and the sum of each row is 1.
}
\author{
  Mikkel Meyer Andersen <mikl@math.aau.dk> and Poul Svante Eriksen
}

\seealso{
\code{\link{disclapmix-package}}
\code{\link{disclapmix}}
\code{\link{disclapmixfit}}
\code{\link{clusterdist}}
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

