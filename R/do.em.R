do.em <-
function(x, centers, iterations, eps, calculate.logLs = FALSE, plots.prefix = NULL, verbose = 0) {
  converged <- FALSE
  warning.messages <- NULL
  
  ##############################################################################

  if (length(centers) != 1) stop("length(centers) != 1")

  ##############################################################################
  
  if (!is.null(plots.prefix)) {
    plots.prefix <- paste(plots.prefix, "c", centers, "-", sep = "")
  }
  
  ##############################################################################
  
  y <- NULL
  
  if (centers == 1) {
    y <- matrix(as.vector(apply(x, 2, median)), nrow = 1)
  } else {
    centers.result <- infer.y(x = x, k = centers)
    
    if (!is.null(plots.prefix)) {
      plot.centers(centers.result, plots.prefix)
    }
    
    y <- centers.result$medoids
  }
  
  ##############################################################################
  
  disclapdata <- create.disclapdata.with.weights(x = x, y = y)

  if (verbose >= 1) {
    cat(nrow(y), " center", ifelse(nrow(y) == 1, "", "s"), ":\n", sep = "")
    print(y)
    cat("\n")
    cat("object.size(disclapdata):\n")
    print(object.size(disclapdata), units = "auto")
    cat("\n")
  }

  if (disclapdata$c == 1) {
    # Only one center, only one iteration in the EM algorithm is 
    # necessary (the weights are all 1 and no reestimation is needed)
    iterations <- 1
  }

  iterations.total <- 0
  logL.full.iterations <- NULL
  logL.marginal.iterations <- NULL
  BIC.full.iterations <- NULL
  BIC.marginal.iterations <- NULL
  v.gain.iterations <- NULL
  tau.iterations <- disclapdata$tau
  changed.center <- NULL
  
  fit <- NULL
  #contr <- glm.control(epsilon = 1e-12, maxit = 100, trace = verbose >= 2)
  contr <- glm.control(trace = verbose >= 2)
 
  repeat {   
    v.max.last <- 0
    last.logL.full <- -Inf
    last.logL.marginal <- -Inf 

    for (iter in 1:iterations) {
      iterations.total <- iterations.total + 1
      
      #################################
      # M-step ########################
      #################################
      if (disclapdata$c == 1) {
        fit <- glm(dist ~ locus - 1, 
          weights = disclapdata$dat$v,
          data = disclapdata$dat, 
          family = DiscreteLaplace(), 
          control = contr)
      } else {
        fit <- glm(dist ~ center + locus - 1, 
          weights = disclapdata$dat$v, 
          data = disclapdata$dat, 
          family = DiscreteLaplace(), 
          control = contr)
      }
      
      if (fit$converged == FALSE) {
        msg <- paste("glm fit did not converge in iteration ", 
          iterations.total, sep = "")
        warning(msg)
        
        if (FALSE) {
          if (verbose >= 1) {
            cat("\nWarning: ", msg, "\n", sep = "")
          } else {          
            warning.messages <- c(warning.messages, msg)
          }
        }
      }
      
      #################################
      # E-step ########################
      #################################
      new.w <- get.w(fit, disclapdata)
      new.v <- get.v(new.w)     
      new.disclapdata <- update.disclapdata.with.new.v(disclapdata, new.v)
      tau.iterations <- rbind(tau.iterations, new.disclapdata$tau)
      
      #################################
      # v gain ########################
      #################################
      v.max.diff <- max(abs(new.v - disclapdata$v))
      v.gain <- v.max.diff / v.max.last
      v.gain.iterations <- c(v.gain.iterations, v.gain)
      
      #################################
      # Updating ######################
      #################################
      disclapdata <- new.disclapdata
      
      if (calculate.logLs == TRUE || verbose >= 1) {
        logL.full <- get.loglikelihood.full(fit, disclapdata$tau, disclapdata$r)      
        logL.full.iterations <- c(logL.full.iterations, logL.full)

        BIC.full <- get.full.BIC(logL.full, disclapdata)
        BIC.full.iterations <- c(BIC.full.iterations, BIC.full)
                
        #->
        pred.ps <- get.p.parameters(disclapdata, fit)
        logL.marginal <- get.loglikelihood.marginal(disclapdata = disclapdata, p = pred.ps)
        logL.marginal.iterations <- c(logL.marginal.iterations, logL.marginal)
        
        BIC.marginal <- get.marginal.BIC(logL.marginal, disclapdata)
        BIC.marginal.iterations <- c(BIC.marginal.iterations, BIC.marginal)
        # <-

        if (verbose >= 1) {
          cat("Iteration ", iter, ":\n", sep = "")
          print.fit(fit, disclapdata)
          cat("max|vij - vijk_old| / max(vijk_old) = ", 
            v.max.diff, " / ", v.max.last, " = ", v.gain, "\n", sep = "")
        }
      }

      if (!is.na(v.gain) && v.gain < eps) {
        converged <- TRUE

        if (verbose >= 1) {
          cat("\nStopping after ", iterations.total, 
            " iterations due to convergence, ", eps, " > ", 
            v.gain, "\n\n", sep = "")
        }

        break
      }
      
      v.max.last <- max(abs(new.v))
      
      if (verbose >= 1) {
        cat("\n")
      }
      
      if (converged) {
        break
      }
    }

    ##############################################################################
    # Move centers
    ##############################################################################
    new.y <- move.centers(disclapdata)
    dist.new.y <- sum(abs(disclapdata$y - new.y))
    
    if (dist.new.y == 0) {
      if (verbose >= 1) {
        cat("Centers were optimal, no need to more EM iterations\n")
      }
      
      break
    } else {
      if (verbose >= 1) {      
        cat("Current centers not optimal:\n")
        print(disclapdata$y)
        cat("With score = ", centers.score(disclapdata$x, disclapdata$y, disclapdata$v), "\n", sep = "")
        cat("New centers:\n")
        print(new.y)
        cat("Differences:\n")
        print(new.y - disclapdata$y)
        cat("With score = ", centers.score(disclapdata$x, new.y, disclapdata$v), "\n", sep = "")
        cat("Number of stepwise mutations between center configurations = ", dist.new.y, "\n", sep = "")
        cat("Doing EM again with the new centers...\n")
      }
      
      converged <- FALSE
      changed.center <- c(changed.center, iterations.total)
      disclapdata <- update.disclapdata.with.new.centers(disclapdata, new.y) 
    }
  }
  
  if (verbose == 0) {
    logL.full <- get.loglikelihood.full(fit, disclapdata$tau, disclapdata$r)      

    pred.ps <- get.p.parameters(disclapdata, fit)
    logL.marginal <- get.loglikelihood.marginal(disclapdata = disclapdata, p = pred.ps)
  }

  ##############################################################################
  #dev.off()
  ##############################################################################
  
  if (disclapdata$c > 1 && !converged) {
    msg <- paste("EM did not converge according to the specified eps = ", 
      eps, " (only reached ", v.gain, ")", sep = "")
    warning(msg)
    
    if (FALSE) {
      if (verbose >= 1) {
        cat("\nWarning: ", msg, "\n", sep = "")
      } else {          
        warning.messages <- c(warning.messages, msg)
      }
    }
  }
 
  ##############################################################################
  # Ready for prediction
  ##############################################################################
  pred.ps <- get.p.parameters(disclapdata, fit)
  ##############################################################################
  if (!is.null(plots.prefix) && !is.null(logL.full.iterations)) {
    pdf(paste(plots.prefix, "logL-full.pdf", sep = ""))
    plot(logL.full.iterations, type = "b", 
      xlab = "Iteration", ylab = "full log likelihood")
    for (iter in changed.center) abline(v = iter + 0.5, lty = 2)
    dev.off()
  }
  
  if (!is.null(plots.prefix) && !is.null(logL.marginal.iterations)) {
    pdf(paste(plots.prefix, "logL-marginal.pdf", sep = ""))
    plot(logL.marginal.iterations, type = "b", 
      xlab = "Iteration", ylab = "sum log match probability")
    for (iter in changed.center) abline(v = iter + 0.5, lty = 2)
    dev.off()
  }
  
  if (!is.null(plots.prefix) && !is.null(BIC.full.iterations)) {
    pdf(paste(plots.prefix, "BIC-full.pdf", sep = ""))
    plot(BIC.full.iterations, type = "b", xlab = "Iteration", ylab = "BIC")
    for (iter in changed.center) abline(v = iter + 0.5, lty = 2)
    dev.off()
  }
  
  if (!is.null(plots.prefix) && !is.null(BIC.marginal.iterations)) {
    pdf(paste(plots.prefix, "BIC-marginal.pdf", sep = ""))
    plot(BIC.marginal.iterations, type = "b", xlab = "Iteration", ylab = "BIC")
    for (iter in changed.center) abline(v = iter + 0.5, lty = 2)
    dev.off()
  }

  if (!is.null(plots.prefix) && disclapdata$c >= 2) {
    pdf(paste(plots.prefix, "tau.pdf", sep = ""))
    plot(tau.iterations[, 1], type = "l", 
      xlab = "Iteration", ylab = expression(tau), 
      ylim = range(tau.iterations), col = 1)  
    for (tau.i in 2:disclapdata$c) lines(tau.iterations[, tau.i], type = "l", col = tau.i)
    for (iter in changed.center) abline(v = iter + 0.5, lty = 2)
    dev.off()
  }
  
  if (iterations.total >= 2) {
    is <- floor(iterations.total / 2):iterations.total

    if (!is.null(plots.prefix) && !is.null(logL.full.iterations)) {
      pdf(paste(plots.prefix, "logL-full-last-half.pdf", sep = ""))
      plot(is, logL.full.iterations[is], type = "b", xlab = "Iteration", ylab = "full log likelihood")
      for (iter in changed.center) abline(v = iter + 0.5, lty = 2)
      dev.off()
      
      logL.full.iterations.sign <- sign(sapply(2:iterations.total, function(i) logL.full.iterations[i] - logL.full.iterations[i-1]))
      pdf(paste(plots.prefix, "logL-full-sign.pdf", sep = ""))
      plot(logL.full.iterations.sign, type = "b", main = "full log L sgn", 
        ylab = expression(sgn(logLik[i] - logLik[i-1])))
      dev.off()
    }

    if (!is.null(plots.prefix) && !is.null(logL.marginal.iterations)) {
      pdf(paste(plots.prefix, "logL-marginal-last-half.pdf", sep = ""))
      plot(is, logL.marginal.iterations[is], type = "b", 
        xlab = "Iteration", ylab = "marginal log likelihood")
      for (iter in changed.center) abline(v = iter + 0.5, lty = 2)
      dev.off()
      
      logL.marginal.iterations.sign <- sign(sapply(2:iterations.total, 
        function(i) logL.marginal.iterations[i] - logL.marginal.iterations[i-1]))
      pdf(paste(plots.prefix, "logL-marginal-sign.pdf", sep = ""))
      plot(logL.marginal.iterations.sign, type = "b", 
        main = "marginal log L sgn", ylab = expression(sgn(logLik[i] - logLik[i-1])))
      dev.off()
    }
        
    if (!is.null(plots.prefix) && any(!is.na(v.gain.iterations))) {
      pdf(paste(plots.prefix, "v-gain.pdf", sep = ""))
      plot(v.gain.iterations, type = "b", xlab = "Iteration", ylab = "max|v - v_old| / max v_old")
      for (iter in changed.center) abline(v = iter + 0.5, lty = 2)
      dev.off()
    }
  }
  
  if (disclapdata$c == 1) {
    converged <- TRUE
  }

  #warning.messages = warning.messages,  
  
  ans <- list(
    converged = converged,
    changed.center = changed.center,
    disclapdata = disclapdata, 
    fit = fit, 
    iterations = iterations.total,
    pred.ps = pred.ps,
    logL.full = logL.full,
    logL.marginal = logL.marginal,  
    logL.full.iterations = logL.full.iterations,
    logL.marginal.iterations = logL.marginal.iterations,
    BIC.full.iterations = BIC.full.iterations, 
    BIC.marginal.iterations = BIC.marginal.iterations, 
    v.gain.iterations = v.gain.iterations,
    tau.iterations = tau.iterations)
  
  class(ans) <- "disclapmixfit"

  return(ans)
}

