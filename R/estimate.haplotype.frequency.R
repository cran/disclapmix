estimate.haplotype.frequency <-
function(h, disclapdata, p) {
  h <- as.numeric(h)

  if (disclapdata$c == 1) {
    pred.prob <- ddisclap(abs(h - disclapdata$y[1, ]), p)
    pred <- prod(pred.prob)
    return(pred)  
  }

  tau <- disclapdata$tau
  
  pred.prob <- lapply(1:disclapdata$c, function(j) { 
      ddisclap(abs(h - disclapdata$y[j, ]), p[[j]])
    })
  pred.prod <- sapply(pred.prob, prod)
  pred <- sum(disclapdata$tau * pred.prod)

  return(pred)
}

