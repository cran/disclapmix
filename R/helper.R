check_x <- function(x) {
  if (!is.matrix(x) || !is.integer(x)) {
    stop("x is not an integer matrix, please use the str function to ensure that the object is a matrix and only as integer columns (and not numeric)")
  } else if (nrow(x) <= 1) {
    stop("x must have more than one row")
  }
}

check_y <- function(y, clusters, loci) {
  if (!is.matrix(y) || !is.integer(y)) {
    stop("y is not an integer matrix, please use the str function to ensure that the object is a matrix and only as integer columns (and not numeric)")
  } 
  #else if (nrow(y) <= 1) {
  #  stop("y must have more than one row")
  #} 
  else if (nrow(y) != clusters) {
    stop("y must have exactly as many rows as the specified number of clusters")
  } else if (ncol(y) != loci) {
    stop("y must have the same number of columns as x")
  }
}

create_model_matrix <- function(x, clusters, model.matrix.function) {
  mat <- rcpp_create_design_matrix(x, clusters)
  colnames(mat) <- c("cluster", "locus")
  
  mat <- as.data.frame(mat)
  mat$cluster <- factor(mat$cluster)
  mat$locus <- factor(mat$locus)
  
  if (clusters == 1L && ncol(x) == 1L) {
    mat$cluster <- NULL
    mat$locus <- NULL
    mat <- model.matrix.function(~ 1, mat)
  } else if (clusters == 1L) {
    mat$cluster <- NULL
    mat <- model.matrix.function(~ locus - 1, mat, contrasts = list(locus = "contr.treatment"))
  } else if (ncol(x) == 1L) {
    mat$locus <- NULL
    mat <- model.matrix.function(~ cluster - 1, mat, contrasts = list(cluster = "contr.treatment"))
  } else {
    mat <- model.matrix.function(~ cluster + locus - 1, mat, contrasts = list(cluster = "contr.treatment", locus = "contr.treatment"))
  }
  
  return(mat)
}

create_response_vector <- function(x, y) {
  response <- rcpp_create_response_vector(x, y)
  return(response)
}

create_initial_y <- function(x, clusters, init_y_method = "pam") {  
  if (is.null(init_y_method) | !is.character(init_y_method) | length(init_y_method) != 1) {
    stop("Invalid init_y_method argument.")
  }
  
  y <- NULL
  
  if (clusters == 1L) {
    y <- matrix(as.integer(round(apply(x, 2L, median))), nrow = 1L)
  } else {
    centers.result <- NULL
        
    if (init_y_method == "pam") {      
      centers.result <- cluster::pam(x, k = clusters, diss = FALSE, stand = FALSE, metric = "manhattan", do.swap = TRUE)
      y <- centers.result$medoids
    } else if (init_y_method == "clara") {
      centers.result <- cluster::clara(x, k = clusters, metric = "manhattan", 
        stand = FALSE, samples = 100, 
        sampsize = min(ceiling(nrow(x)/2), 100 + 2*clusters),
        medoids.x = TRUE, keep.data = FALSE,
        rngR = TRUE, pamLike = TRUE)
      y <- centers.result$medoids
    } else {
      stop("Unsupported init_y_method chosen")
    }

    y <- as.matrix(apply(y, 2L, function(r) as.integer(round(r))))
  }

  return(y)
}

convert_coef_to_disclap_parameters <- function(beta, clusters) {
  if (clusters == 1L) {
    pcl <- matrix(exp(beta), nrow = 1)
    return(pcl)  
  } else {
    theta_cs <- beta[1L:clusters]
    theta_ls <- c(0L, beta[-(1L:clusters)])
    theta <- outer(theta_cs, theta_ls, "+")
    pcl <- exp(theta)
    return(pcl)
  }
}

convert_coef_to_disclap_parameters_internal <- function(beta, clusters) {
  loci <- length(beta) - clusters
  theta_ls <- beta[1L:loci]
  theta_cs <- beta[-(1L:loci)]
  theta <- outer(theta_cs, theta_ls, "+")
  pcl <- exp(theta)
  return(pcl)
}

move_centers <- function(x, y, v_matrix) {
  y_candidates <- apply(x, 2L, range)
  
  new_y <- sapply(1L:ncol(x), function(l) {
    ycls <- y_candidates[1L, l]:y_candidates[2L, l]
    res <- sapply(ycls, function(ycl) {
      sapply(1L:nrow(y), function(cluster) {
        sum(v_matrix[, cluster]*abs(x[, l] - ycl))
      })
    })
    
    if (nrow(y) == 1L) {
       indices <- which.min(res)
    } else {
       indices <- apply(res, 1L, which.min)
    }

    return(ycls[indices])
  })
  
  if (nrow(y) == 1L) {
    new_y <- matrix(new_y, nrow = 1L)
  }
  
  colnames(new_y) <- colnames(y)
  return(new_y)
}

get_loglikelihood_full <- function(fit, clusters, response_vector, weight_vector, tau_norm) { 
  theta <- fit$linear.predictors
  p <- exp(theta)
  #print(head(theta))
  #print(head(p))
  #print(theta[p < 0 | p >= 1])
  #print(p[p < 0 | p >= 1])
  
  den <- ddisclap(response_vector, p)
  zi <- rep(1L:clusters, each = length(theta)/clusters)
  tau_norm <- tau_norm[zi]
  logL <- sum(weight_vector*log(tau_norm*den))
  
  return(logL)
}

get_loglikelihood_marginal <- function(x, y, disclap_parameters, tau_vector) {
  xprobs <- rcpp_calculate_haplotype_probabilities(x, y, disclap_parameters, tau_vector)
  logL <- sum(log(xprobs))
  return(logL)
}

get_AIC <- function(logL, individuals, clusters, loci) {
        #coord  #glm         # tau (sums to 1)
  k <- (clusters * loci) + (loci + clusters - 1) + (clusters - 1)

  # Our observations is individuals' haplotype
  # We have a model for the haplotype
  return(2*k - 2*logL)
}

get_BIC <- function(logL, individuals, clusters, loci) {
        #coord  #glm         # tau (sums to 1)
  k <- (clusters * loci) + (loci + clusters - 1) + (clusters - 1)

  # Our observations is individuals' haplotype
  # We have a model for the haplotype
  return(log(individuals)*k - 2*logL)
}

predict.disclapmixfit <-
function(object, newdata, ...) {
  if (!is(object, "disclapmixfit")) stop("object must be a disclapmixfit")
  probs <- rcpp_calculate_haplotype_probabilities(newdata, object$y, object$disclap_parameters, object$tau)
  return(probs)
}

print.disclapmixfit <-
function(x, ...) {
  if (!is(x, "disclapmixfit")) stop("x must be a disclapmixfit")
  
  cat("disclapmixfit from ", 
    formatC(nrow(x$v_matrix)), " observations on ", 
    formatC(ncol(x$y)), " loci with ", 
    formatC(nrow(x$y)), " clusters.\n", sep = "")
  
  return(invisible(x))
}

summary.disclapmixfit <-
function(object, ...) {
  if (!is(object, "disclapmixfit")) stop("object must be a disclapmixfit")
  
  print(object)
  
  cat("\n")
  cat("EM converged:                                                       ", object$converged, "\n", sep = "")
  cat("Number of central haplotype changes:                                ", length(object$changed_center), "\n", sep = "")
  cat("Total number of EM iterations:                                      ", formatC(object$iterations), "\n", sep = "")
  cat("Model observations (n*loci*clusters):                               ", formatC(object$model_observations), "\n", sep = "")
  cat("Model parameters ((clusters*loci)+(loci+clusters-1)+(clusters-1)):  ", formatC(object$model_parameters), "\n", sep = "")
  cat("GLM method:                                                         ", formatC(object$glm_method), "\n", sep = "")
  cat("Initial central haplotypes supplied:                                ", !is.null(object$init_y), "\n", sep = "")
  
  if (is.null(object$init_y)) {
  cat("Method to find initial central haplotypes:                          ", object$init_y_method, "\n", sep = "")
  }

  return(invisible(object))
}

simulate.disclapmixfit <- function(object, nsim = 1L, seed = NULL, ...) {
  if (!is(object, "disclapmixfit")) stop("object must be a disclapmixfit")
  
  if (!is.null(seed)) {
    stop("seed must be null or else ")
  }
  
  if (is.null(nsim) || length(nsim) != 1L || !is.integer(nsim) || nsim <= 0L) {
    stop("nsim must be >= 1L (note the L postfix for integer)")
  }
  
  tau_cumsum <- cumsum(object$tau)
  res <- rcpp_simulate(nsim, object$y, tau_cumsum, object$disclap_parameters)

  return(res)
}

haplotype_diversity <- function(object, nsim = 1e4L) {
  if (!is(object, "disclapmixfit")) stop("object must be a disclapmixfit")
  
  if (is.null(nsim) || length(nsim) != 1L || !is.integer(nsim) || nsim <= 0L) {
    stop("nsim must be >= 1L (note the L postfix for integer)")
  }
  
  get_db_counts <- function(x) { 
    order_x <- do.call(order, as.data.frame(x))
    equal.to.previous <- rowSums(x[tail(order_x, -1),] != x[head(order_x, -1),]) == 0
    indices <- split(order_x, cumsum(c(TRUE, !equal.to.previous)))
    Ns <- unlist(lapply(indices, length))
    return(Ns)
  }
  
  db <- simulate.disclapmixfit(object, nsim = nsim)
  Ns <- get_db_counts(db)
  freqs <- Ns / nsim
  #D <- 1 - sum(freqs^2)
  D <- (nsim / (nsim - 1)) * (1 - sum(freqs^2))
  
  return(D)
}

INTERNAL_glmfit <- function(loci, clusters, individuals, response_vector, apriori_probs, weight_vector, vmat, 
  verbose = FALSE, stop_by_deviance = TRUE, epsilon = 1e-4, maxit = 25L) {
  
  converged <- FALSE
  devold <- Inf
  dev <- 0

  ks <- 1:loci
  js <- 1:clusters
  
  mi <- min(loci, clusters)
  di <- max(loci, clusters) - mi
    
  ############################################################
  # Initialisation
  ############################################################
  tmp_d <- array(response_vector, c(loci, clusters, individuals))
 
  RkSj <- matrix(0, nrow = loci, ncol = clusters)
  for (k in ks) {
    for (j in js) {
      RkSj[k, j] <- RkSj[k, j] + sum(vmat[, j]*tmp_d[k, j, ])
    }
  }
  
  Rk <- apply(RkSj, 1, sum)
  Sj <- apply(RkSj, 2, sum)
  ###########################################################
  
  d_vec <- as.numeric(response_vector)
  
  beta <- c(rep(0.5, loci), rep(-1, clusters))
  beta_correction <- beta
  beta_correction_old <- beta

  lin.pred <- rep(beta[ks], clusters * individuals)
  lin.pred <- lin.pred + rep(rep(beta[-(ks)], each = loci), individuals)

  for (iter in 1L:maxit) {
    if (stop_by_deviance) {
      mu_m <- disclapglm_linkinv(lin.pred)
      
      # Deviance:
      dev <- disclapglm_deviance(d_vec, mu_m, weight_vector)
      
      if (verbose == TRUE) {
        cat("  IRLS iteration ", iter, ", deviance = ", dev, "\n", sep = "")
      }
      
      # Stop if devience is NaN
      if (is.null(dev) || is.na(dev) || is.nan(dev) || is.infinite(dev)) {
        converged <- FALSE
        break
      }
      
      if (abs(dev - devold)/(0.1 + abs(dev)) < epsilon) {
        converged <- TRUE
        break
      }
      
      devold <- dev
    } else if (iter > 1) {
      beta_change <- max(abs(beta_correction - beta_correction_old)/(0.1 + abs(beta_correction)))
      
      if (verbose == TRUE) {
        cat("  IRLS iteration ", iter, ", max relative coefficient change = ", beta_change, "\n", sep = "")
      }
      
      if (beta_change < epsilon) {
        converged <- TRUE
        break
      }

      beta_correction_old <- beta_correction
    }
    ############################################################################
    
    lin_pred_generic <- outer(beta[ks], beta[(loci+1):(loci+clusters)], "+")
    
    mu_generic <- apply(lin_pred_generic, 2, disclapglm_linkinv)
    if (!is.matrix(mu_generic)) { # == 1 cluster
      mu_generic <- matrix(mu_generic, nrow = loci, ncol = clusters) 
    }
    
    var_generic <- apply(mu_generic, 2, disclapglm_varfunc)    
    if (!is.matrix(var_generic)) { # == 1 cluster
      var_generic <- matrix(var_generic, nrow = loci, ncol = clusters)
    }
    
    apriori_probs_indv <- apriori_probs * individuals

    offD <- matrix(0, nrow = loci, ncol = clusters)
    mu_generic_apriori_probs <- matrix(0, nrow = loci, ncol = clusters)
    for (k in ks) {
      offD[k, ] <- apriori_probs_indv * var_generic[k, ]
      mu_generic_apriori_probs[k, ] <- apriori_probs_indv * mu_generic[k, ]
    }
    
    D1 <- apply(offD, 1, sum)
    D2 <- apply(offD, 2, sum)
    y1 <- Rk - apply(mu_generic_apriori_probs, 1, sum)
    y2 <- Sj - apply(mu_generic_apriori_probs, 2, sum)    
    
    H <- D1^-.5*offD
    H <- t(D2^-.5*t(H))
    
    dec <- svd(H, loci, clusters)

    OD <- 1-dec$d^2
    OD[1] <- 0.25
    
    if (mi >= 2L) { # > 1 cluster
      OD[2:mi] <- 1/OD[2:mi]
    }
    
    if (loci == mi) {
      dr <- OD
    } else {
      dr <- c(OD, rep(1, di))
    }
    
    if (clusters == mi) {
      ds <- OD
    } else {
      ds <- c(OD, rep(1, di))
    }
    
    dia <- (dec$d-2*dec$d*OD)/(1+dec$d^2)
    if (length(dia) > 1) { # > 1 cluster
      dia <- diag(dia)
    }

    K <- matrix(0, loci, clusters)
    
    if (loci <= clusters) {
      K[ks, ks] <- dia
    } else {
      K[1:clusters, 1:clusters] <- dia
    }
    
    if (length(dr) > 1) {
      dr <- diag(dr)
    }
    
    if (length(ds) > 1) {
      ds <- diag(ds)
    }

    mI <- matrix(0, loci+clusters, loci+clusters)
    U <- D1^-.5*dec$u
    V <- D2^-.5*dec$v

    mI[ks, ks] <- U %*% dr %*% t(U)
    mI[ks, (loci+1):(loci+clusters)] <- U %*% K %*% t(V)
    mI[(loci+1):(loci+clusters), ks] <- t(mI[ks,(loci+1):(loci+clusters)])
    mI[(loci+1):(loci+clusters), (loci+1):(loci+clusters)] <- V %*% ds %*% t(V)  
    beta_correction <- mI %*% c(y1, y2)
    lin_pred_correction <- rep(beta_correction[ks, 1], clusters * individuals)
    lin_pred_correction <- lin_pred_correction + rep(rep(beta_correction[(loci+1):(loci+clusters), 1], each = loci), individuals)
    ##################    
    beta <- beta + beta_correction
    lin.pred <- lin.pred + lin_pred_correction
    ##################
  }

  coefficients <- as.numeric(beta)

  if (stop_by_deviance == FALSE) {
    mu_m <- disclapglm_linkinv(lin.pred)
    dev <- disclapglm_deviance(d_vec, mu_m, weight_vector)
  }
  
  ans <- list(
    coefficients = coefficients,
    converged = converged,
    deviance = dev,
    linear.predictors = lin.pred
  )

  return(ans)
}

