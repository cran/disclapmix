centers.score <-
function(x, y, v) {
  score <- 0

  for (k in 1:ncol(x)) {
    for (j in 1:nrow(y)) {
      score <- score + sum(v[, j]*abs(x[, k] - y[j, k]))
    }
  }

  return(score)
}

