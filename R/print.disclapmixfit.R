print.disclapmixfit <-
function(x, ...) {
  cat("disclapmixfit from ", 
    formatC(x$disclapdata$n), " observations on ", 
    formatC(x$disclapdata$r), " loci with ", 
    formatC(x$disclapdata$c), " centers.\n", sep = "")
  
  return(invisible(x))
}

