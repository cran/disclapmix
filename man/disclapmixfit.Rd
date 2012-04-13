\name{disclapmixfit}
\alias{disclapmixfit}
\title{Discrete Laplace mixture fit by the EM algorithm}
\description{A \code{disclapmixfit} object contains various information about the Discrete Laplace mixture fit by the EM algorithm obtained using the \code{\link{disclapmix}} function.}
\value{
  \describe{
    \item{\code{changed.center}}{A vector with the iteration numbers where the centers have changed}
    \item{\code{disclapdata}}{Structure containing the database \code{x}, the centers \code{y}, the matrix \code{v} of each observation's probability of belonging to a certain center, the data frame with weights used when fitting the model with \code{\link{glm}} \code{dat}, the prior probabilities of belonging to a center \code{tau}, \code{n = nrow(x)}, \code{c = nrow(y)}, and \code{r = ncol(x)}}
    \item{\code{fit}}{The \code{\link{glm}} fit of the final model}
    \item{\code{iterations}}{Number of iterations performed in total (including moving centers and reestimating using the EM algorithm)}
    \item{\code{pred.ps}}{A list where each entry is a vector of locus parameters, each entry in the list corresponds to a center}
    \item{\code{logL.full}}{Full log likelihood of the final model}
    \item{\code{logL.marginal}}{Marginal log likelihood of the final model}
    \item{\code{logL.full.iterations}}{Full log likelihood of the models during the iterations}
    \item{\code{logL.marginal.iterations}}{Marginal log likelihood of the models during the iterations}
    \item{\code{BIC.full.iterations}}{BIC based on full log likelihood of the models during the iterations}
    \item{\code{BIC.marginal.iterations}}{BIC based on marginal log likelihood of the models during the iterations}
    \item{\code{v.gain.iterations}}{The gain \eqn{\frac{\max \{ v_{new} - v_{old} \}}{\max \{ v_{old} \}}}{| max (v\_new - v\_old) | / max(v\_old)}, where \code{v} is the matrix mentioned above under \code{disclapdata}, during the iterations}
    \item{\code{tau.iterations}}{The prior probability of belonging to the centers during the iterations}
  }
}
\author{
Mikkel Meyer Andersen <mikl@math.aau.dk> and Poul Svante Eriksen
}

\seealso{
\code{\link{disclapmix-package}}
\code{\link{disclapmix}}
\code{\link{plot.disclapmixfit}}
\code{\link{predict.disclapmixfit}}
\code{\link{extractMarginalBIC}}
\code{\link{extractFullBIC}}
\code{\link[disclap:disclap-package]{disclap}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{disclap}
\keyword{disclapmix}
\keyword{disclapmixfit}
\keyword{discrete Laplace}
\keyword{EM algorithm}

