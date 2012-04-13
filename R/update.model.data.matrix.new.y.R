update.model.data.matrix.new.y <-
function(dat, y) {
  dat[, 6] <- y[matrix(c(dat[, 2], dat[, 3]), ncol = 2)]
  dat[, 7] <- abs(dat[, 5] - dat[, 6])
  return(dat)
}

