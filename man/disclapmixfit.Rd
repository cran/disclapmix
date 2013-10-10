\name{disclapmixfit}
\alias{disclapmixfit}
\title{Discrete Laplace mixture fit by the EM algorithm}
\description{A \code{disclapmixfit} object contains various information about the Discrete Laplace mixture fit by the EM algorithm obtained using the \code{\link{disclapmix}} function.}
\value{
  \describe{
    \item{\code{glm_method}}{The supplied GLM method.}
    \item{\code{init_y}}{The supplied initial central haplotypes, \code{init_y}.}
    \item{\code{init_y_method}}{The supplied method for choosing initial central haplotypes (only used if \code{init_y} is \code{NULL}).}

    \item{\code{converged}}{Whether the estimation converged or not.}
       
    \item{\code{y}}{The central haplotypes, \code{y}.}
    \item{\code{tau}}{The prior probabilities of belonging to a cluster, \code{tau}.}
    \item{\code{v_matrix}}{The matrix \code{v} of each observation's probability of belonging to a certain cluster. The rows are in the same order as the observations in \code{x} used to generate this fit.}
    \item{\code{disclap_parameters}}{A matrix with the estimated dicrete Laplace parameters.}
    \item{\code{glm_coef}}{The coefficients from the last GLM fit (used to calculate \code{disclap_parameters}).}
    
    \item{\code{model_observations}}{Number of observations.}
    \item{\code{model_parameters}}{Number of parameters in the model.}
    \item{\code{iterations}}{Number of iterations performed in total (including moving centers and re-estimating using the EM algorithm).}

    \item{\code{logL_full}}{Full log likelihood of the final model.}
    \item{\code{logL_marginal}}{Marginal log likelihood of the final model.}
    \item{\code{BIC_full}}{BIC based on the full log likelihood of the final model.}
    \item{\code{BIC_marginal}}{BIC based on the marginal log likelihood of the final model.}

    \item{\code{v_gain_iterations}}{The gain \eqn{\frac{\max \{ v_{new} - v_{old} \}}{\max \{ v_{old} \}}}{| max (v\_new - v\_old) | / max(v\_old)}, where \code{v} is \code{vic_matrix} mentioned above, during the iterations.}
    \item{\code{tau_iterations}}{The prior probability of belonging to the centers during the iterations.}

    \item{\code{changed_center}}{A vector with the iteration numbers where the centers have changed.}
    \item{\code{centers_iterations}}{The centers before the changes in \code{changed_center}.}
    
    \item{\code{logL_full_iterations}}{Full log likelihood of the models during the iterations (only calculated when \code{verbose = 2L}).}
    \item{\code{logL_marginal_iterations}}{Marginal log likelihood of the models during the iterations (only calculated when \code{verbose = 2L}).}
    \item{\code{BIC_full_iterations}}{BIC based on full log likelihood of the models during the iterations (only calculated when \code{verbose = 2L}).}
    \item{\code{BIC_marginal_iterations}}{BIC based on marginal log likelihood of the models during the iterations (only calculated when \code{verbose = 2L}).}
  }
}
\author{
Mikkel Meyer Andersen <mikl@math.aau.dk> and Poul Svante Eriksen
}

\seealso{
\code{\link{disclapmix-package}}
\code{\link{disclapmix}}
\code{\link{clusterdist}}
\code{\link{clusterprob}}
\code{\link{predict.disclapmixfit}}
\code{\link{print.disclapmixfit}}
\code{\link{summary.disclapmixfit}}
\code{\link[disclap:disclap-package]{disclap}}
}

\keyword{disclapmixfit}

