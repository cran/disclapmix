get.loglikelihood.marginal <-
function(disclapdata, p) {
  match.probs <- sapply(1:disclapdata$n, 
    function(i) {
      estimate.haplotype.frequency(h = disclapdata$x[i, ], 
        disclapdata = disclapdata, p = p)
    })
  
  logL <- sum(log(match.probs))
  return(logL)
}

