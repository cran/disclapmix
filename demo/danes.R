data(danes)

db <- danes[rep(1:nrow(danes), danes$n), 1:(ncol(danes)-1)]
res <- disclapmix(db, centers = 1:5, use.parallel = TRUE, calculate.logLs = TRUE, verbose = 0)

par(ask = TRUE)

################################################################################
# Best fit: 
#   Basic information
################################################################################
res$best.fit
summary(res$best.fit)

################################################################################
# Best fit: 
#   Plot how prior probabily (tau) of belonging to the centers has changed
################################################################################
plot(res$best.fit, which = 1)

################################################################################
# Best fit: 
#   Plot how marginal log likelihood has changed
################################################################################
plot(res$best.fit, which = 3)

################################################################################
# Plot the BIC (based on marginal likelihood) for the different number of centers
################################################################################
marginalBICs <- sapply(res$fits, extractMarginalBIC)
plot(marginalBICs, xlab = "Number of centers", ylab = "Marginal BIC")

################################################################################
# Best fit: 
#   Predict haplotype probabilities
################################################################################
disclap.estimates <- predict(res$best.fit, newdata = danes[, 1:(ncol(danes) - 1)])

# Robbins' (1968) estimate of unobserved probability mass
sum(danes$n == 1) / (sum(danes$n) + 1)

# disclapmix' estimate
1 - sum(disclap.estimates)

plot(danes$n, disclap.estimates, 
  xlab = "Number of times the haplotype has been observed",
  ylab = "Predicted frequency using the discrete Laplace mixture model")

par(ask = FALSE)

