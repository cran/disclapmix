plot.centers <-
function(centers.result, plots.prefix) {
  if (is.null(plots.prefix)) {
    return(FALSE)
  }

  center.plot <- plot
  
  # Silhouette plot only for >= 2 centers 
  if (nrow(centers.result$medoids) == 1) {
    center.plot <- clusplot
  }
  
  pdf(paste(plots.prefix, "centers-partition.pdf", sep = ""))
  center.plot(centers.result)
  dev.off()
  
  return(TRUE)
}

