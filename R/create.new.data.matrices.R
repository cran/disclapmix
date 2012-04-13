create.new.data.matrices <-
function(disclapdata) {
  r <- disclapdata$r # loci
  y <- disclapdata$y
  c <- disclapdata$c # centers

  dat <- matrix(NA, nrow = r, ncol = 3)
  dat <- data.frame(dat)
  colnames(dat) <- c("center", "locus", "v")
  dat[, 2] <- factor(1:r, levels(disclapdata$dat$locus))
  dat[, 3] <- 1
  
  if (c == 1) {
    dat[, 1] <- 1
    return(dat)
  } 
  
  c.lvls <- levels(disclapdata$dat$center)
  
  res <- lapply(1:c, function(j) {
    d <- dat
    d[, 1] <- factor(j, c.lvls)
    return(d)
  })  

  return(res)
}

