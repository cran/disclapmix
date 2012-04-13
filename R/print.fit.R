print.fit <-
function(fit, disclapdata, digits = 3) {
  logL.full <- get.loglikelihood.full(fit, disclapdata$tau, disclapdata$r)
  
  cat("IWLS itererations:       ", fit$iter, "\n", sep = "")
  cat("IWLS converged:          ", fit$converged, "\n", sep = "")
  cat("IWLS boundary:           ", fit$boundary, "\n", sep = "")
  cat("full log likelihood:     ", round(logL.full, digits), "\n", sep = "")
  cat("marginal log likelihood: ", round(get.loglikelihood.marginal(disclapdata = disclapdata, 
    p = get.p.parameters(disclapdata, fit)), digits), "\n", sep = "")
  cat("Null Deviance:           ", round(fit$null.deviance, digits), 
    " (on ", fit$df.null, " degrees of freedom)\n", sep = "")
  cat("Residual Deviance:       ", round(fit$deviance, digits), 
    " (on ", fit$df.residual, " degrees of freedom)\n", sep = "")
  cat("BIC:                     ", round(get.full.BIC(logL.full, disclapdata), digits), 
    "\n", sep = "")
}

