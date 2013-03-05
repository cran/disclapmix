disclapmix <-
function(x, centers = 1:5, use.parallel = FALSE, iterations = 25, eps = 0.001, calculate.logLs = FALSE, plots.prefix = NULL, verbose = 1) {
  if (use.parallel == TRUE && verbose != 0) {
    stop("You cannot get verbose output when running in parallel")
  }
  
  em.res.s <- NULL
  
  if (length(centers) > 1 && use.parallel == TRUE) {
    no.cores <- NA
    try(no.cores <- detectCores(), silent = TRUE)   
    if (is.na(no.cores)) {
      no.cores <- 4
    }
    
    cluster.arguments <- list(x = x, iterations = iterations, eps = eps, 
        calculate.logLs = calculate.logLs, plots.prefix = plots.prefix, 
        verbose = verbose)
    
    cl <- makeCluster(no.cores)
    
    em.res.s <- clusterApplyLB(cl = cl, x = centers, fun = function(centercount, cluster.arguments) {
      x <- cluster.arguments$x
      iterations <- cluster.arguments$iterations
      eps <- cluster.arguments$eps
      calculate.logLs <- cluster.arguments$calculate.logLs
      plots.prefix <- cluster.arguments$plots.prefix
      verbose <- cluster.arguments$verbose
      
      em.res <- do.em(x = x, centers = centercount, 
        iterations = iterations, eps = eps, calculate.logLs = calculate.logLs, 
        plots.prefix = plots.prefix, verbose = verbose)
      
      return(em.res)
    }, cluster.arguments = cluster.arguments)
    
    stopCluster(cl)
  } else {
    em.res.s <- lapply(centers, function(centercount) {
      em.res <- do.em(x = x, centers = centercount, iterations = iterations, eps = eps, calculate.logLs = calculate.logLs, plots.prefix = plots.prefix, verbose = verbose)

      if (verbose >= 1) {
        cat("disclapmix for ", centercount, " centers done\n", sep = "")
      }
      
      return(em.res)
    })
  }
  
  marginal.logL <- sapply(em.res.s, function(em.res) em.res$logL.marginal)
  marginal.BIC <- sapply(1:length(em.res.s), function(em.res.i) get.marginal.BIC(logL.marginal = marginal.logL[em.res.i], disclapdata = em.res.s[[em.res.i]]$disclapdata))
  
  best.model.index <- which.min(marginal.BIC)
  em.res <- em.res.s[[best.model.index]]
  
  if (length(centers) > 1 && em.res$disclapdata$c %in% range(centers)) {
    warning("The best number of centers was on the boundary of the specified allowed number of centers; consider expanding the allowed number of centers.")
  }
  
  res <- list(best.fit = em.res, best.fit.index = best.model.index, fits = em.res.s)
  
  return(res)
}

