infer.y <-
function(x, k) { 
  res <- pam(x, k, diss = FALSE, stand = FALSE, metric = "manhattan")
  return(res)
}

