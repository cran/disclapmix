create.model.data.matrix.with.weights <-
function(x, y, v) {
  if (any(abs(apply(v, 1, sum) - 1) > 1e-15)) {
    stop("v_{i+} != 1")
  }

  n <- nrow(x) # individuals
  r <- ncol(x) # loci
  c <- nrow(y) # centers
  
  dat <- matrix(NA, nrow = n*c*r, ncol = 7)
  colnames(dat) <- c("individual", "center", "locus", "v", "x", "y", "dist")
  
  dat[, 1] <- as.vector(sapply(1:n, rep.int, times = r*c))
  dat[, 2] <- rep(as.vector(sapply(1:c, rep.int, times = r)), n)
  dat[, 3] <- rep(1:r, n*c)
  
  cols <- t(apply(dat, 1, function(row) {
    i <- row[1]
    j <- row[2]
    k <- row[3]
    return(c(v[i, j], x[i, k], y[j, k]))
  }))
  
  dat[, 4] <- cols[, 1]
  dat[, 5] <- cols[, 2]
  dat[, 6] <- cols[, 3]
  dat[, 7] <- abs(dat[, 5] - dat[, 6])
  
  dat <- data.frame(dat)
  dat$individual <- factor(dat$individual)
  
  if (c > 1) {
    dat$center <- factor(dat$center)
  }
  
  if (r > 1) {
    dat$locus <- factor(dat$locus)
  }

  return(dat)
}

