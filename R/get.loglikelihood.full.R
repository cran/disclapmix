get.loglikelihood.full <-
function(fit, tau, r) { 
  theta <- fit$linear.predictors
  p <- exp(theta)
  d <- fit$data$dist
  v <- fit$data$v
  den <- ddisclap(d, p)
  
  if (abs(sum(tau) - 1) > 1e-15) {
    print(tau)
    print(sum(tau))
    stop("sum(tau) != 1")
  }
  
  zi <- fit$data$center
  
  tau <- tau^(1/r)  
  tau <- tau[zi]
    
  logL <- sum(v*log(tau*den))
  
  return(logL)
}

