get.w <-
function(fit, disclapdata) {
  new.dat <- cbind(disclapdata$dat, 
    theta = fit$linear.predictors, 
    p = exp(fit$linear.predictors))

  new.w <- aggregate(ddisclap(dist, p) ~ individual + center, new.dat, prod)
  
  if (ncol(new.w) != 3) stop("ncol(new.w) != 3")
  if (nrow(new.w) != disclapdata$n*disclapdata$c) {
    stop("nrow(new.w) != disclapdata$n*disclapdata$c")
  }
  
  new.w.mat <- matrix(NA, nrow = disclapdata$n, ncol = disclapdata$c)
  
  for (j in 1:disclapdata$c) {
    new.w.mat[, j] <- disclapdata$tau[j] * subset(new.w, new.w$center == j)$"ddisclap(dist, p)"
  }
  
  return(new.w.mat)
}

