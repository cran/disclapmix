update.disclapdata.with.new.centers <-
function(disclapdata, new.y) {
  new.disclapdata <- disclapdata
  new.disclapdata$y <- new.y
  new.disclapdata$dat <- update.model.data.matrix.new.y(disclapdata$dat, new.y)
  
  return(new.disclapdata)
}

