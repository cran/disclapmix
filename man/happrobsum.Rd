\name{happrobsum}
\alias{happrobsum}
\title{Haplotypes probability sum}
\description{Calculate sum of haplotype probabilities where haplotypes are defined by allele ranges.}
\usage{
  happrobsum(object, alleles, \dots)
}
\arguments{
  \item{object}{a \code{\link{disclapmixfit}} object}
  \item{alleles}{matrix defining the allele ranges: a column per locus where first row is minimum allele and second row is maximum allele}
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
\code{\link{simulate.disclapmixfit}}
\code{\link{plot.disclapmixfit}}
\code{\link{predict.disclapmixfit}}
\code{\link{haplotype_diversity}}
\code{\link{clusterdist}}
\code{\link{clusterprob}}
}
\keyword{happrobsum}

