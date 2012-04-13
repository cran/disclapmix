extractFullBIC <- 
function(fit) {
  if (!is(fit, "disclapmixfit")) stop("fit must be a disclapmixfit")
  return(get.full.BIC(logL.full = fit$logL.full, disclapdata = fit$disclapdata))
}

