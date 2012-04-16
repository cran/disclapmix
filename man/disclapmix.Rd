\name{disclapmix}
\alias{disclapmix}
\title{Discrete Laplace mixture inference using the EM algorithm}
\description{\code{disclapmix} makes inference in a mixture of Discrete Laplace distributions using the EM algorithm. After the EM algorithm has converged, the centers are moved if the marginal likelihood increases by doing so. And then the EM algorithm is run again. This continues until the centers are not moved.}
\usage{
disclapmix(x, centers = 1:5, use.parallel = FALSE, iterations = 25, eps = 0.001, calculate.logLs = FALSE, plots.prefix = NULL, verbose = 1)
}
\arguments{
  \item{x}{Dataset}
  \item{centers}{A vector of center counts to fit the model for}
  \item{use.parallel}{Whether to do the computations in parallel (requires \code{verbose = 0})}
  \item{iterations}{Maximum number of iterations in the EM-algorithm}
  \item{eps}{Convergence stop criteria in the EM algorithm which is compared to \eqn{\frac{\max \{ v_{new} - v_{old} \}}{\max \{ v_{old} \}}}{| max (v\_new - v\_old) | / max(v\_old)}, where \code{v} is a matrix of each observation's probability of belonging to a certain center}
  \item{calculate.logLs}{Whether to calculate the model's log likelihood at each iterations (might be slightly slower than not doing it)}
  \item{plots.prefix}{Where to put plots (NULL to avoid plotting)}
  \item{verbose}{from 0 to 3 (both including): 0 for silent, 3 for extra verbose}
}
\value{
  \describe{
    \item{\code{best.fit}}{The \code{\link{disclapmixfit}} object with the lowest marginal BIC score}
    \item{\code{best.fit.index}}{The index such that \code{best.fit = fits[[best.fit.index]]}}
    \item{\code{fits}}{A list of the \code{\link{disclapmixfit}} objects corresponding to the specified \code{centers} vector}
  }
}
\author{
Mikkel Meyer Andersen <mikl@math.aau.dk> and Poul Svante Eriksen
}

\seealso{
\code{\link{disclapmix-package}}
\code{\link{disclapmix}}
\code{\link{disclapmixfit}}
\code{\link{plot.disclapmixfit}}
\code{\link{predict.disclapmixfit}}
\code{\link{extractMarginalBIC}}
\code{\link{extractFullBIC}}
\code{\link[disclap:disclap-package]{disclap}}
}
\examples{
\dontrun{
demo(danes)
demo(simpop)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{disclapmix}
\keyword{centers}
\keyword{use.parallel}
\keyword{eps}

