get.tau <-
function(v) {
  if (ncol(v) == 1) {
    return(1)
  }
  
  tau <- apply(v, 2, mean)
  
  if (abs(sum(tau) - 1) > 1e-15) {
    print(tau)
    print(sum(tau))
    stop("sum(tau) != 1")
  }
  
  return(tau)
}

