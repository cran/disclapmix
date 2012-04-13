get.marginal.BIC <-
function(logL.marginal, disclapdata) {
  n <- disclapdata$n
  c <- disclapdata$c
  r <- disclapdata$r
  
  k <- c*(r + 2) # r + 1 + 1 for each center: r coords, 1 tau and 1 alpha_j
  
  return(log(n)*k - 2*logL.marginal)
}

