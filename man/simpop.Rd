\name{simpop}
\alias{simpop}
\docType{data}
\title{Simulated Y-STR haplotypes}
\description{
  500 simulated Y-STR 7 loci haplotypes. 
  
  First a population has been simulated with \code{\link[fwsim:fwsim]{fwsim}} using the following parameters: initial population size of 10000, 500 generations, mutation rate of 0.001, and a population growth such that expected population size was 20,000,000 after the 500 generations. The actual end population sized was 19397385 consisting of 34180 different haplotypes.
  
  Afterwards a database of 500 has been drawn randomly from this population.
}
\usage{data(simpop)}
\format{
  A data frame with 500 observations: on the following 7 loci (\code{n} is the number of times each haplotype has been observed)
  \describe{
    \item{\code{L1}}{Locus 1}
    \item{\code{L2}}{Locus 2}
    \item{\code{L3}}{Locus 3}
    \item{\code{L4}}{Locus 4}
    \item{\code{L5}}{Locus 5}
    \item{\code{L6}}{Locus 6}
    \item{\code{L7}}{Locus 7}
    \item{\code{Freq}}{Frequency in the population}
    \item{\code{n}}{Number of times the haplotype has been observed in the database}
  }
}
\seealso{
\code{\link[fwsim:fwsim]{fwsim}}
}
\keyword{datasets}
