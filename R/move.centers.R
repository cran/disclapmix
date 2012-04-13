move.centers <-
function(disclapdata) {
  x <- disclapdata$x
  y.candidates <- apply(x, 2, range)
  v <- disclapdata$v
  
  new.y <- sapply(1:disclapdata$r, function(k) {
    yjks <- y.candidates[1, k]:y.candidates[2, k]
    res <- sapply(yjks, function(yjk) {
      sapply(1:disclapdata$c, function(j) {
        sum(v[, j]*abs(x[, k] - yjk))
      })
    })
    
    if (disclapdata$c == 1) {
       indices <- which.min(res)
    } else {
       indices <- apply(res, 1, which.min)
    }

    return(yjks[indices])
  })
  
  if (disclapdata$c == 1) {
    new.y <- matrix(new.y, nrow = 1)
    colnames(new.y) <- colnames(disclapdata$y)
    return(new.y)
  }
  
  colnames(new.y) <- colnames(disclapdata$y)
  return(new.y)
}

