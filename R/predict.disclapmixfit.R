predict.disclapmixfit <-
function(object, newdata, ...) {
  if (!is(object, "disclapmixfit")) stop("object must be a disclapmixfit")

  disclapdata <- object$disclapdata
  p <- object$pred.ps
  
  disclap.estimates <- apply(newdata, 1, function(haplotype) estimate.haplotype.frequency(haplotype, disclapdata, p))

  return(disclap.estimates)
}

