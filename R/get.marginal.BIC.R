get.marginal.BIC <-
function(logL.marginal, disclapdata) {
  c <- disclapdata$c
  r <- disclapdata$r
  
        #coord  #glm         # tau (sums to 1)
  k <- (c*r) + (r + c - 1) + (c-1)
    
  return(log(disclapdata$n)*k - 2*logL.marginal)
}

