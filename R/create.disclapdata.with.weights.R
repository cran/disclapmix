create.disclapdata.with.weights <-
function(x, y) {
  v <- matrix(1/nrow(y), nrow = nrow(x), ncol = nrow(y))  
  tau <- get.tau(v)
  dat <- create.model.data.matrix.with.weights(x, y, v)
  disclapdata <- list(n = nrow(x), r = ncol(x), c = nrow(y), x = x, y = y, dat = dat, v = v, tau = tau)
  class(disclapdata) <- c("disclapdataweighted", class(disclapdata)) 
  return(disclapdata)
}

