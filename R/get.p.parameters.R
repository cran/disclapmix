get.p.parameters <-
function(disclapdata, fit) {
  new.data <- create.new.data.matrices(disclapdata)

  if (disclapdata$c == 1) {
    pred.mus <- predict(fit, newdata = new.data, type = "response")
    pred.ps <- sapply(pred.mus, get.p.from.mu)
  } else {
    pred.mus <- lapply(new.data, function(nd) { 
        predict(fit, newdata = nd, type = "response") 
      })
    pred.ps <- lapply(pred.mus, get.p.from.mu)
  }

  return(pred.ps)
}

