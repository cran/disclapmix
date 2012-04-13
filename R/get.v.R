get.v <-
function(w) {
  if (ncol(w) == 1) {
    new.v <- matrix(rep(1, nrow(w), ncol = 1, nrow = nrow(w)))
    return(new.v)
  }
  
  new.v <- t(apply(w, 1, function(wij) wij / sum(wij)))
  return(new.v)
}

