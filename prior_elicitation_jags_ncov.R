###########################################
### Prior elicitation for re-implementationn of C. Althaus' CFR estimation using JAGS
### Johannes Bracher, johannes.bracher@uzh.ch
###########################################

# try to come up with a reasonable prior reproducing the CIs in Linton et al:

# assume that shape is fixed, choose prior for rate such that quantiles
# coincide with CIs for median
# (motivation: quantiles of gamma distribution are linearly related to the inverse rate;
# but not sure this is the best way to go)

# point estimates from Linton et al:
death_shape <- 5.049793
death_rate <- 0.3294643

# set up a grid for the rate parameter:
grid_rate <- seq(from = 0.2, to = 0.6, by = 0.01)
# evaluate the implied median duration to death:
vals_med <- sapply(grid_rate, qgamma, p = 0.5, shape = death_shape)

# limits of a CI for rate implied by the CI for the median:
upper_bound_rate <- grid_rate[which.min(vals_med > 12.2)]
lower_bound_rate <- grid_rate[which.min(vals_med > 16.7)]

# visualize:
plot(grid_rate, vals_med, type = "l", xlab = "rate", ylab = "median duration to death")
abline(h = c(12.2, 14.3, 16.7))
abline(v = c(lower_bound_rate, death_rate, upper_bound_rate))

# slightly skewed... a possible transformation to induce symmetry is x^(-2)
plot(grid_rate^-2, vals_med, type = "l", xlab = "rate^(-2)", ylab = "median")
abline(h = c(12.2, 14.3, 16.7))
abline(v = c(lower_bound_rate, death_rate, upper_bound_rate)^-2)

# compute standard deviation on the one-over-square scale:
- (upper_bound_rate^(-2) - death_rate^(-2))/1.96
(lower_bound_rate^(-2) - death_rate^(-2))/1.96
# choose 1.355 as standard deviation of death_rate^(-2)

# check the limits this prior implies for the median of the duration to death:
qgamma(0.5, shape = death_shape, rate = upper_bound_rate) # 12.2; should be 4.8
qgamma(0.5, shape = death_shape, rate = lower_bound_rate) # 16.3; should be 16.7
# this works.

# check the limits this prior implies for the 5% quantile of the duration to death:
qgamma(0.05, shape = death_shape, rate = upper_bound_rate) # 5.1; should be 4.8
qgamma(0.05, shape = death_shape, rate = lower_bound_rate) # 6.9; should be 7.9

# and the limits for the 95% quantile of the duration to death:
qgamma(0.95, shape = death_shape, rate = upper_bound_rate) # 23.6; should be 21.8
qgamma(0.95, shape = death_shape, rate = lower_bound_rate) # 31.8; should be 34.7

# somewhat too narrow, but roughly okay

# What does this prior look like?
prior_samples_rate <- rnorm(100000, mean = death_rate^(-2),
                            1.35)^(-0.5)
plot(density(prior_samples_rate), xlab = "death_rate")
abline(v = death_rate)

# plot distribution of time to death for 0.025, 0.5 and 0.975 prior quantiles of death_rate:
png("figures/prior_choice_death_rate.png", height = 350, width = 700)
curve(dgamma(x, rate = upper_bound_rate, shape = death_shape), 0, 40,
      col = "darkgreen", xlab = "Time from onset to death (days)", ylab = "Probability density",
      frame = FALSE, main = "Plausible range of distributions of time from onset\n to death under chosen prior")
curve(dgamma(x, rate = lower_bound_rate, shape = death_shape), col = "red", add = TRUE)
curve(dgamma(x, rate = death_rate, shape = death_shape), col = "steelblue", add = TRUE)

abline(v = c(qgamma(0.5, shape = death_shape, rate = upper_bound_rate),
             14.3,
             qgamma(0.5, shape = death_shape, rate = lower_bound_rate)),
       col = c("darkgreen", "steelblue", "red"), lty = 2)

legend("topright", legend = c(paste0("rate = ",
                                     c(upper_bound_rate, round(death_rate, 2), lower_bound_rate)),
                              "medians"),
       col = c("darkgreen", "steelblue", "red", "black"), lty = c(rep(1, 3), 2), bty = "n")
dev.off()
