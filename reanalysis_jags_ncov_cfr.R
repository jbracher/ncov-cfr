###########################################
### Re-implementationn of C. Althaus' CFR estimation using JAGS
### Johannes Bracher, johannes.bracher@uzh.ch
###########################################

# Short summary: Bayesian and likelihood inference yield almost identical results
# and the effect of carrying forward the uncertainty about the parameters
# of the gamma distribution from Linton et al is vanishingly small.
# Caveat: Eliciting a prior from Linton et al (see separate file) is not straightforward
# (but all more or less reasonable choices should yield quite similar results)

#-------------------------------------------------------------------------#
# re-running code by C. Althaus for likelihood-ased inference (for comparison):

# Load libraries
library(lubridate)
library(bbmle)
library(plotrix)

# Load 2019-nCoV cases (n=130) identified outside of China
# Source: WHO Novel Coronavirus(2019-nCoV) Situation Report - 19, and media reports for the death on the Philippines
exports <- read.csv("data/ncov_cases_20200215.csv")
begin <- ymd(exports$date[1])
cases <- exports$cases
deaths <- exports$deaths
n_cases <- sum(cases)
n_deaths <- sum(deaths)
n_days <- length(cases)
days <- 1:n_days
interval <- seq(1, n_days + 7, 7)
onset <- rep(days, cases)

# Estimating gamma distribution of onset of death from Linton et al. (http://dx.doi.org/10.1101/2020.01.26.20018754)
linton <- function(shape, rate) {
  ssr <- sum((qgamma(c(0.05, 0.5, 0.95), shape = exp(shape), rate = exp(rate)) - c(6.1, 14.3, 28.0))^2)
  return(ssr)
}
free_linton <- c(shape = log(4), rate = log(4/14.3))
fit_linton <- mle2(linton, start = as.list(free_linton), method = "Nelder-Mead")

# Likelihood and expected mortality function
nll <- function(cfr, death_shape, death_rate) {
  cfr <- plogis(cfr)
  expected <- numeric(n_days)
  for(i in days) {
    for(j in 1:n_cases) {
      d <- i - onset[j]
      if(d >= 0) {
        expected[i] <- expected[i] + cfr*diff(pgamma(c(d - 0.5, d + 0.5), shape = death_shape, rate = death_rate))
      }
    }
  }
  ll <- sum(dpois(deaths, expected, log = TRUE)) # Daily incidence of deaths
  #ll <- dpois(1, sum(expected), log = TRUE) # Cumulative incidence of deaths
  return(-ll)
}

# Fit the model
free <- c(cfr = 0)
fixed <- c(death_shape = exp(coef(fit_linton)[[1]]), death_rate = exp(coef(fit_linton)[[2]]))
fit <- mle2(nll, start = as.list(free), fixed = as.list(fixed), method = "Brent", lower = -100, upper = 100)
summary(fit)
plogis(coef(fit)[1])
confidence <- confint(fit)
plogis(confidence)

#-------------------------------------------------------------------------#
# New code for Bayesian inference using JAGS:

library(R2jags)
library(rjags)

# extract point estimates of shape and scale parameter:
death_shape = exp(coef(fit_linton)[[1]])
death_rate = exp(coef(fit_linton)[[2]])

#------------------------------------------#
# the JAGS model:
jags.mod <- function(){
  cfr ~ dunif(0, 1) # "uninformative" prior

  # prior on rate parameter of gamma distribution describing time to death
  # see below for elicitation
  # use a transformed parameter for which a symetrical prior is reasonable
  death_rate_trafo ~ dnorm(9.2, 1.355^(-2))
  death_rate = death_rate_trafo^(-0.5)

  # manually implement weights resulting from gamma distribution
  # all weights are shifted by n_days in order to assign zero weights
  # in case of negative time differences
  # note that the integrals over the gamma distribution are approximated
  # by the gamma density at the center of the respective interval
  for(i in 1:(n_days - 1)){
    wgts[i] = 0
  }
  for(i in n_days:(2*n_days)){
    wgts[i] = death_rate^death_shape/
      ((exp(loggam(death_shape))))*
      (i - n_days)^(death_shape - 1)*exp(-(i - n_days)*death_rate)
  }

  # actual model definition:
  for(i in 1:n_days){
    for(p in 1:n_cases){
      contribution[i, p] = wgts[i - onset[p] + n_days]
    }
    expected[i] = cfr*sum(contribution[i, ])
    deaths[i] ~ dpois(expected[i])
  }
}

#------------------------------------------#
# Fit no. 1: with distribution of onset of death fixed (for comparison with ML estimation)

# data to be passed to JAGS:
death_rate_trafo <- death_rate^(-2)
dat.jags <- list("deaths", "n_days", "n_cases", "onset", "death_shape", "death_rate_trafo")
# also passing death_rate_trafo here

# parameters for which to store posterior:
jags.mod.params <- c("cfr", "death_rate")

# initial values for three chains:
jags.mod.inits <- list(list("cfr" = 0.01), list("cfr" = 0.05), list("cfr" = 0.1))

# run:
set.seed(111) # set seed to make reproducible
jags.mod.fit <- jags(data = dat.jags, inits = jags.mod.inits,
                      parameters.to.save = jags.mod.params, n.chains = 3,
                     n.iter = 1000000, n.burnin = 10000, model.file = jags.mod)

jags.mod.fit

# smooth posterior density:
post.cfr <- density(jags.mod.fit$BUGSoutput$sims.list$cfr)
plot(post.cfr)
# get different point estimates:
(post.mod.cfr <- post.cfr$x[which.max(post.cfr$y)])
(post.mean.cfr <- jags.mod.fit$BUGSoutput$mean$cfr)
(post.median.cfr <- jags.mod.fit$BUGSoutput$median$cfr)
# differ as posterior is somewhat skewed

# compare posterior mode and ML estimator:
post.mod.cfr; quantile(jags.mod.fit$BUGSoutput$sims.list$cfr, c(0.05, 0.95))
plogis(coef(fit)[1]); plogis(confidence)

#------------------------------------------#
# Fit no. 2: with priors assigned to parameters of
# the distribution of onset of death

# Re-run the model without providing death_rate parameter in the data, i.e using the
# elicited prior:
dat.jags2 <- list("deaths", "n_days", "n_cases", "onset", "death_shape")
set.seed(111) # set seed to make reproducible
jags.mod.fit2 <- jags(data = dat.jags2, inits = jags.mod.inits,
                     parameters.to.save = jags.mod.params, n.chains = 3,
                     n.iter = 100000, n.burnin = 10000, model.file = jags.mod)
jags.mod.fit2
jags.mod.fit
# virually no difference exceeding the usual Monte Carlo error.

# smooth posterior and compute mode:
post.cfr2 <- density(jags.mod.fit2$BUGSoutput$sims.list$cfr)
(post.mod.cfr2 <- post.cfr2$x[which.max(post.cfr2$y)])

#------------------------------------------#
# compare results with and without prior on death_rate as well as likelihood-based:
png("figures/cfr_comparison_ML_vs_JAGS.png", height = 350, width = 700)
plot(post.cfr, main = "", xlab = "Time from onset to death (days)")
abline(v = post.mod.cfr, col = "black", lty = 2)
lines(post.cfr2, col = "blue")
abline(v = post.mod.cfr2, col = "blue", lty = 2)

# add likelihood curve for comparison:
vals_cfr <- 0:150/1000
lik_cfr <- exp(-sapply(qlogis(vals_cfr), FUN = nll,
                       death_shape = death_shape, death_rate = death_rate))
lik_cfr_scaled <- lik_cfr/sum(lik_cfr)*1000 # re-scale likelihood
lines(vals_cfr, lik_cfr_scaled, col = "red")
abline(v = vals_cfr[which.max(lik_cfr_scaled)], col = "red", lty = 2)

legend("topright", legend = c("posterior density, death_rate fixed",
                              "posterior density, with prior on death_rate",
                              "scaled likelihood function",
                              "posterior mode / ML estimator"),
       col = c("black", "blue", "red"), lty = c(1, 1, 1, 2), bty = "n")
dev.off()
# pretty much identical results.