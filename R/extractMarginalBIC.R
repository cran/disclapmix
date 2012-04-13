extractMarginalBIC <- 
function(fit) {
  if (!is(fit, "disclapmixfit")) stop("fit must be a disclapmixfit")
  return(get.marginal.BIC(logL.marginal = fit$logL.marginal, disclapdata = fit$disclapdata))
}

