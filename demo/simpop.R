data(simpop)

db <- simpop[rep(1:nrow(simpop), simpop$n), 1:7]
res <- disclapmix(db, centers = 1:5, use.parallel = TRUE, calculate.logLs = TRUE, verbose = 0)

par(ask = TRUE)

################################################################################
# Best fit: 
#   Basic information
################################################################################
res$best.fit
summary(res$best.fit)

################################################################################
# Plot the BIC (based on marginal likelihood) for the different number of centers
################################################################################
marginalBICs <- sapply(res$fits, extractMarginalBIC)
plot(marginalBICs, xlab = "Number of centers", ylab = "Marginal BIC")

################################################################################
# Best fit: 
#   Predict haplotype probabilities
################################################################################
disclap.estimates <- predict(res$best.fit, newdata = simpop[, 1:7])

# Robbins' (1968) estimate of unobserved probability mass
sum(simpop$n == 1) / (sum(simpop$n) + 1)

# disclapmix' estimate
1 - sum(disclap.estimates)

plot(simpop$Freq, disclap.estimates, log = "xy", 
  xlab = "Frequency in the population of 20,000,000 from which the database was sampled",
  ylab = "Predicted frequency using the discrete Laplace mixture model")

par(ask = FALSE)

