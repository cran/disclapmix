summary.disclapmixfit <-
function(object, ...) {
  
  cat("Number of observations        = ", 
    formatC(object$disclapdata$n), "\n", sep = "")
  cat("Number of loci                = ", 
    formatC(object$disclapdata$r), "\n", sep = "")
  cat("Number of centers             = ", 
    formatC(object$disclapdata$c), "\n", sep = "")

  cat("Full log likelihood           = ", 
    formatC(object$logL.full), "\n", sep = "")
  cat("Marginal log likelihood       = ", 
    formatC(object$logL.marginal), "\n", sep = "")
  cat("BIC (full log likelihood)     = ", 
    formatC(extractFullBIC(object)), "\n", sep = "")
  cat("BIC (marginal log likelihood) = ", 
    formatC(extractMarginalBIC(object)), "\n", sep = "")

  cat("Center information:\n")
  
  for (j in 1:object$disclapdata$c) {
    distdb <- as.vector(apply(object$disclapdata$x, 1, function(h) sum(abs(h - object$disclapdata$y[1, ]))))
  
    cat("  y[", j, "] = (", 
      paste(object$disclapdata$y[j, ], collapse = ", "), 
      ") with tau[", j, "] = ", 
      formatC(object$disclapdata$tau[j]), " and L1 norm to nearest haplotype in database = ", min(distdb), "\n", sep = "")
  }

  cat("Centers changed               = ", 
    length(object$changed.center), " times", sep = "")

  if (length(object$changed.center) > 0) {
    cat(" after the following number of total iterations: ", 
      paste(object$changed.center, collapse = ", "), sep = "")
  }  
  cat("\n")
  
  return(invisible(object))
}

