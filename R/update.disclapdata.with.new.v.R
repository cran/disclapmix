update.disclapdata.with.new.v <-
function(disclapdata, new.v) {
  new.disclapdata <- disclapdata
  new.disclapdata$v <- new.v
  new.disclapdata$tau <- get.tau(new.v)
  new.disclapdata$dat <- update.model.data.matrix.new.v(disclapdata$dat, new.v)
  
  return(new.disclapdata)
}

