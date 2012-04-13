update.model.data.matrix.new.v <-
function(dat, v) {
  if (any(abs(apply(v, 1, sum) - 1) > 1e-15)) {
    stop("v_{i+} != 1")
  }
  
  dat[, 4] <- v[matrix(c(dat[, 1], dat[, 2]), ncol = 2)]
  return(dat)
}

