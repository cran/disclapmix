\name{disclapmix}
\alias{disclapmix}
\title{Discrete Laplace mixture inference using the EM algorithm}
\description{\code{disclapmix} makes inference in a mixture of Discrete Laplace distributions using the EM algorithm. After the EM algorithm has converged, the centers are moved if the marginal likelihood increases by doing so. And then the EM algorithm is run again. This continues until the centers are not moved.}
\usage{
disclapmix(x, clusters, init_y = NULL, iterations = 100L, eps = 0.001, 
  verbose = 0L, glm_method = "internal_coef", glm_control_maxit = 50L, 
  glm_control_eps = 1e-6, init_y_method = "pam", \dots)
}
\arguments{
  \item{x}{Dataset.}
  \item{clusters}{The number of clusters/components to fit the model for.}
  \item{init_y}{Initial central haplotypes, if NULL, these will be estimated as described under the \code{init_y_method} argument.}
  \item{iterations}{Maximum number of iterations in the EM-algorithm.}
  \item{eps}{Convergence stop criteria in the EM algorithm which is compared to \eqn{\frac{\max \{ v_{new} - v_{old} \}}{\max \{ v_{old} \}}}{| max (v\_new - v\_old) | / max(v\_old)}, where \code{v} is a matrix of each observation's probability of belonging to a certain center.}
  \item{verbose}{from 0 to 2 (both including): 0 for silent, 2 for extra verbose.}
  \item{glm_method}{\code{internal_coef}, \code{internal_dev} or \code{glm.fit}. Please see details.}
  \item{glm_control_maxit}{Integer giving the maximal number of IWLS iterations.}
%  \item{glm_control_eps}{Positive convergence tolerance epsilon; the iterations converge when \code{|dev - dev_{old}|/(|dev| + 0.1) < epsilon}. For \code{internal_coef}, \code{max(abs(beta_correction - beta_correction_old)/(|beta_correction| + 0.1))} is used instead (that is, the maximum of the changes in the parameter vector).}
  \item{glm_control_eps}{Positive convergence tolerance epsilon; the iterations converge when \code{|x - x_{old}|/(|x| + 0.1) < epsilon}, where \code{x = beta_correction} for \code{internal_coef} and \code{x = deviance} otherwise.}
  \item{init_y_method}{Which cluster method to use for finding initial central haplotypes, y: \code{pam} (recommended) or \code{clara}. Ignored if \code{init_y} is supplied.}
  \item{...}{Used to detect obsolete usage (when using parameters \code{centers}, \code{use.parallel}, \code{calculate.logLs} or \code{plots.prefix}).}
}
\details{
  \code{glm_method}: 
  \code{internal_coef} is the fastest as it uses the relative changes in the coefficients as a stopping criterium, hence it does not need to compute the deviance until the very end. In normal situations, it would not be a problem to use this method. \code{internal_dev} is the reasonably fast method that uses the deviance as a stopping criterium (like \code{glm.fit}). \code{glm.fit} to use the traditional \code{glm.fit} IWLS implementation and is slow compared to the other two methods.
  
  \code{init_y_method}: 
  For \code{init_y_method = 'clara'}, the sampling parameters are: \code{samples = 100}, \code{sampsize = min(ceiling(nrow(x)/2), 100 + 2*clusters)} and the random number generator in R is used.
}
\value{
  A \code{\link{disclapmixfit}} object with estimated parameters.
}
\author{
Mikkel Meyer Andersen <mikl@math.aau.dk> and Poul Svante Eriksen
}

\seealso{
\code{\link{disclapmix-package}}
\code{\link{disclapmix}}
\code{\link{disclapmixfit}}
\code{\link{predict.disclapmixfit}}
\code{\link{print.disclapmixfit}}
\code{\link{summary.disclapmixfit}}
\code{\link{simulate.disclapmixfit}}
\code{\link{haplotype_diversity}}
\code{\link{clusterdist}}
\code{\link{clusterprob}}
\code{\link{glm.fit}}
\code{\link[disclap]{disclap}}
\code{\link[cluster]{pam}}
\code{\link[cluster]{clara}}
}
\examples{
# Generate sample database
db <- matrix(disclap::rdisclap(1000, 0.3), nrow = 250, ncol = 4)

# Add location parameters
db <- sapply(1:ncol(db), function(i) as.integer(db[, i]+13+i))

head(db)

fit1 <- disclapmix(db, clusters = 1L, verbose = 1L, glm_method = "glm.fit")
fit1$disclap_parameters
fit1$y

# Generate another type of database
db2 <- matrix(disclap::rdisclap(2000, 0.1), nrow = 500, ncol = 4)
db2 <- sapply(1:ncol(db2), function(i) as.integer(db2[, i]+14+i))
fit2 <- disclapmix(rbind(db, db2), clusters = 2L, verbose = 1L)
fit2$disclap_parameters
fit2$y
fit2$tau
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{disclapmix}
\keyword{clusters}
\keyword{eps}

