plot.disclapmixfit <-
function(x, which = 1, ...) {
  if (!is(x, "disclapmixfit")) stop("x must be a disclapmixfit")
  
  plot.centers.changed <- function(ys = c(0, 1), legend.pos = "topright") {
    if (length(x$changed.center) > 0) {
      for (cc in x$changed.center) {
        lines(c(cc, cc), ys, col = 1, lty = 2, )
      }
      
      legend(legend.pos, legend = "Centers changed", col = 1, lty = 2)
    }
  }
  
  which <- match.call(expand.dots=TRUE)$which

  if (1 %in% which) {
    plot(x$tau.iterations[, 1], col = 1, type = "l", 
      ylim = c(0, 1), xlab = "Iteration", ylab = expression(tau[j]^(i)))
    
    if (x$disclapdata$c > 1) {
      for (i in 2:x$disclapdata$c) {
        lines(x$tau.iterations[, i], col = i)
      }
    }
    
    plot.centers.changed()
  }

  if (2 %in% which) {
    plot(x$v.gain.iterations, col = 1, type = "b", xlab = "Iteration", ylab = "v gain")
    plot.centers.changed()
  }

  if (3 %in% which) {
    plot(x$logL.marginal.iterations, col = 1, type = "b", xlab = "Iteration", ylab = "marginal log likelihood")
    plot.centers.changed(ys = range(x$logL.marginal.iterations), legend.pos = "bottomright")
  }

  if (4 %in% which) {
    plot(x$logL.full.iterations, col = 1, type = "b", xlab = "Iteration", ylab = "full log likelihood")
    plot.centers.changed(ys = range(x$logL.full.iterations), legend.pos = "bottomright")
  }

  
  return(invisible(NULL))
}

