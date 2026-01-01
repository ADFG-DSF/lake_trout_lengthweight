## (1)
## This script first fits a Von B growth model to a large Lake Trout
## length-age dataset I found online

## https://figshare.com/articles/dataset/Lake_Trout_growth_data_csv_file_/6957620

## I was curious how sample quantiles might compare to estimated L_inf.  In this
## dataset, the 90th percentile length gave a good estimate of L_inf.


## (2)
## Then, the same model is fit to a smaller dataset (Chandler lake)
## Again, approximately 90% of the lengths were smaller than the estimated L_inf

## Oddly, both models estimate L_inf at the 89th percentile!


## (3)
## Then, a nifty model was constructed to see if mortality (exponential decay)
## and selectivity (logistic) combined would be consistent with the age
## distribution observed in the large dataset.  It worked surprisingly well!

## Note: the same model was used with the Chandler Lake dataset.  It didn't
## work at all, which was likely a result of the much smaller dataset.


## (4)
## Finally, a simulation was constructed using the mortality parameter estimated
## in (3) with the vonB parameters estimated in (1) to see if the results were
## still consistent with the 90th percentile method.

## They weren't. Still not sure what to do with that.


library(tidyverse)
library(jagsUI)
library(jagshelper)

## (1) reading large dataset and defining lvb model

LML <- read_csv(file="flat_data/LML_LT_Growth.csv")

par(mfrow=c(1,1))
with(LML, plot(LENGTH~FINAL_AGE))
with(LML, plot(WEIGHT~FINAL_AGE, log="xy"))

# {\displaystyle L(a)=L_{\infty }(1-\exp(-k(a-t_{0})))}

# specify model, which is written to a temporary file
vb_jags <- tempfile()
cat('model {
  for(i in 1:n) {
    y[i] ~ dnorm(mu[i], tau)
    ypp[i] ~ dnorm(mu[i], tau)
    mu[i] <- Linf*(1-exp(-k*(x[i]-t0)))
  }

  tau <- pow(sig, -2)
  sig ~ dunif(0, 200)
  sig_prior ~ dunif(0, 200)
  Linf ~ dnorm(700, 0.0001)
  Linf_prior ~ dnorm(700, 0.0001)
  t0 ~ dnorm(0, 0.01)
  t0_prior ~ dnorm(0, 0.01)
  # k ~ dnorm(0, 0.01)
  # k_prior ~ dnorm(0, 0.01)
  k ~ dexp(0.1)
  k_prior ~ dexp(0.1)
}', file=vb_jags)


# bundle data to pass into JAGS
vb_data <- list(x = LML$FINAL_AGE[!is.na(LML$FINAL_AGE) & !is.na(LML$LENGTH)],
                y = LML$LENGTH[!is.na(LML$FINAL_AGE) & !is.na(LML$LENGTH)],
                n = sum(!is.na(LML$FINAL_AGE) & !is.na(LML$LENGTH)))

# JAGS controls
niter <- 100000
# ncores <- 3
ncores <- min(10, parallel::detectCores()-1)

{
  tstart <- Sys.time()
  print(tstart)
  vb_jags_out <- jagsUI::jags(model.file=vb_jags, data=vb_data,
                              parameters.to.save=c("sig","Linf","t0","k","mu","ypp",
                                                   "sig_prior","Linf_prior","t0_prior","k_prior"),
                              n.chains=ncores, parallel=T, n.iter=niter,
                              n.burnin=niter/2, n.thin=niter/2000)
  print(Sys.time() - tstart)
}

# nbyname(vb_jags_out)
plotRhats(vb_jags_out)
traceworstRhat(vb_jags_out, parmfrow = c(3, 3))
comparepriors(vb_jags_out, parmfrow = c(3, 3))

par(mfrow=c(2,2))
qq_postpred(vb_jags_out, p="ypp", y=vb_data$y)
# ts_postpred(vb_jags_out, p="ypp", y=vb_data$y)
ts_postpred(vb_jags_out, p="ypp", y=vb_data$y, x=vb_data$x)


envelope_vec <- function(x, ...) {  # making this as a raw vector first
  xdf <- data.frame(x,x,x)
  envelope(xdf, x=par("usr")[c(1,2,2)], add=TRUE, ...=...)
  ## ok found a bug in envelope: does not seem to work with a 2-column df!
}


plot(vb_data)
envelope(vb_jags_out, p="mu", x=vb_data$x, add=TRUE)
envelope_vec(x=vb_jags_out$sims.list$Linf, col=2, dark=.1)
plot(vb_data)
envelope(vb_jags_out, p="ypp", x=vb_data$x, add=TRUE)
envelope_vec(x=vb_jags_out$sims.list$Linf, col=2, dark=.1)


## ok moment of truth, how do the length quantiles compare to modeled Linf??
mean(vb_data$y < vb_jags_out$q50$Linf)

Linf_post <- vb_jags_out$sims.list$Linf
quantile_post <- colMeans(outer(vb_data$y, Linf_post, "<"))
plotdens(quantile_post)


## (2) Fitting the same model with Chandler lake dataset

## that was fun, let's do it again
chandler <- read_csv(file="flat_data/Copy of ChandlerLake_LakeTroutOtoliths_Tyers.csv")
chan_data <- list(x = chandler$Age[!is.na(chandler$Age)],
                  y = chandler$Length.mmFL[!is.na(chandler$Age)],
                  n = sum(!is.na(chandler$Age)))
{
  tstart <- Sys.time()
  print(tstart)
  chan_jags_out <- jagsUI::jags(model.file=vb_jags, data=chan_data,
                              parameters.to.save=c("sig","Linf","t0","k","mu","ypp",
                                                   "sig_prior","Linf_prior","t0_prior","k_prior"),
                              n.chains=ncores, parallel=T, n.iter=niter,
                              n.burnin=niter/2, n.thin=niter/2000)
  print(Sys.time() - tstart)
}

# nbyname(vb_jags_out)
plotRhats(chan_jags_out)
traceworstRhat(chan_jags_out, parmfrow = c(3, 3))
comparepriors(chan_jags_out, parmfrow = c(3, 3))

par(mfrow=c(2,2))
qq_postpred(chan_jags_out, p="ypp", y=chan_data$y)
ts_postpred(chan_jags_out, p="ypp", y=chan_data$y)

plot(chan_data)
envelope(chan_jags_out, p="mu", x=chan_data$x, add=TRUE)
envelope_vec(x=chan_jags_out$sims.list$Linf, col=2, dark=.1)
plot(chan_data)
envelope(chan_jags_out, p="ypp", x=chan_data$x, add=TRUE)
envelope_vec(x=chan_jags_out$sims.list$Linf, col=2, dark=.1)


## ok moment of truth, how do the length quantiles compare to modeled Linf??
mean(chan_data$y < chan_jags_out$q50$Linf)

Linf_post <- chan_jags_out$sims.list$Linf
quantile_post <- colMeans(outer(chan_data$y, Linf_post, "<"))
plotdens(quantile_post)



## quick & dirty sim to look into how universal this 90% thing is
## interestingly, the 90% thing worked great for Chandler and not for LML

## actually it's pretty weird that it worked because the age sample should
## not be uniform!

# t0 <- vb_jags_out$q50$t0
# k <- vb_jags_out$q50$k
# Linf <- vb_jags_out$q50$Linf
# sig <- vb_jags_out$q50$sig

t0 <- chan_jags_out$q50$t0
k <- chan_jags_out$q50$k
Linf <- chan_jags_out$q50$Linf
sig <- chan_jags_out$q50$sig

nsim <- 1000
xsim <- runif(nsim, 0, 50)
# xsim <- rexp(nsim, .01)
musim <- Linf*(1-exp(-k*(xsim-t0)))
ysim <- rnorm(n=nsim, mean=musim, sd=sig)
plot(xsim, musim)
plot(xsim, ysim)
mean(ysim < Linf)



## (3) Defining a model that fits natural mortality (exponential decay) and
## selectivity (logistic)

# specify model, which is written to a temporary file
mortrate_jags <- tempfile()
cat('model {
  for(i in 1:n) {
    y[i] ~ dbin(p[i]*pi*(lambda^i), 10000)
    ypp[i] ~ dbin(p[i]*pi*(lambda^i), 10000)
    logit(p[i]) <- b0 + b1*age[i]
  }

  # N[1] <- 1000
  # for(j in 2:n) {
  #   N[j] ~ dbin(lambda, N[j-1])
  # }

  lambda ~ dbeta(5,5)
  pi ~ dbeta(5,5)
  b0 ~ dnorm(0, 0.001)
  b1 ~ dnorm(0, 0.001)
}', file=mortrate_jags)



rawtable <- table(vb_data$x) #%>% plot
agetable <- rep(0, max(as.numeric(names(rawtable))))
agetable[as.numeric(names(rawtable))] <- rawtable
plot(agetable)


mortrate_data <- list(y = agetable,
                      age = 1:length(agetable),
                      n = length(agetable))

# JAGS controls
niter <- 100000
# ncores <- 3
ncores <- min(10, parallel::detectCores()-1)

{
  tstart <- Sys.time()
  print(tstart)
  mortrate_jags_out <- jagsUI::jags(model.file=mortrate_jags, data=mortrate_data,
                                    parameters.to.save=c("lambda","pi","b0","b1","N","ypp","p"),
                                    n.chains=ncores, parallel=T, n.iter=niter,
                                    n.burnin=niter/2, n.thin=niter/2000)
  print(Sys.time() - tstart)
}

# nbyname(mortrate_jags_out)
plotRhats(mortrate_jags_out)
traceworstRhat(mortrate_jags_out, parmfrow = c(3, 3))
qq_postpred(mortrate_jags_out, p="ypp", y=mortrate_data$y)
ts_postpred(mortrate_jags_out, p="ypp", y=mortrate_data$y)
par(mfrow=c(3,2))
envelope(mortrate_jags_out, "ypp")
points(mortrate_data$y)
envelope(mortrate_jags_out, p="p")
envelope(outer(mortrate_jags_out$sims.list$lambda, mortrate_data$age, "^"),
         main="lambda effect")

envelope(outer(mortrate_jags_out$sims.list$lambda, mortrate_data$age, "^") *
           mortrate_jags_out$sims.list$p,
         main="lambda x p")
envelope(outer(mortrate_jags_out$sims.list$lambda, mortrate_data$age, "^") *
           mortrate_jags_out$sims.list$p *
           mortrate_jags_out$sims.list$pi,
         main="lambda x p x pi")
plotcor_jags(mortrate_jags_out)




## (4) Simulating data from this and see if it's consistent with the vb model & output

# defining selectivity for a given sample: hopefully this will be self-normalized later
selectivity_mat <- outer(mortrate_jags_out$sims.list$lambda, mortrate_data$age, "^") *
  mortrate_jags_out$sims.list$p *
  mortrate_jags_out$sims.list$pi

# defining sample of ages (x 10000 mcmc rows)
nsample <- 500  # number of fish in each simulated sample
sample_mat <- matrix(nrow=nrow(selectivity_mat), ncol=ncol(selectivity_mat))
for(i in 1:nrow(sample_mat)) {
  sample_mat[i, ] <- rmultinom(1, nsample, selectivity_mat[i,])
  # sample_mat[i, ] <- rmultinom(1, nsample, colMeans(selectivity_mat))
}

# turning this into raw age data (mcmc rows x nsample columns)
age_sim <- matrix(nrow=nrow(sample_mat), ncol=nsample)
for(i in 1:nrow(sample_mat)) {
  mcmc_vec <- NULL
  for(j in 1:ncol(sample_mat)) {  # this is a hack
    mcmc_vec <- c(mcmc_vec, rep(j, sample_mat[i,j]))
  }
  age_sim[i,] <- mcmc_vec
}

# now simulating!
# mu[i] <- Linf*(1-exp(-k*(x[i]-t0)))
length_sim <- NA*age_sim
for(j in 1:ncol(age_sim)) { # probably a better way to do this, oh well
  length_sim[, j] <- rnorm(n = nrow(age_sim),
                           mean = vb_jags_out$sims.list$Linf *
                             (1 - exp(-vb_jags_out$sims.list$k *
                                        (age_sim[, j] - vb_jags_out$sims.list$t0))),
                           sd = vb_jags_out$sims.list$sig)
}
plot(vb_data)
points(age_sim[2,], length_sim[2,], pch=16, col=2)
boxplot(rowMeans(length_sim < vb_jags_out$q50$Linf),
        rowMeans(length_sim < vb_jags_out$sims.list$Linf))

length_sim <- NA*age_sim
for(j in 1:ncol(age_sim)) { # probably a better way to do this, oh well
  length_sim[, j] <- rnorm(n = nrow(age_sim),
                           mean = vb_jags_out$q50$Linf *
                             (1 - exp(-vb_jags_out$q50$k *
                                        (age_sim[, j] - vb_jags_out$q50$t0))),
                           sd = vb_jags_out$q50$sig)
}
plot(age_sim[2,], length_sim[2,])
boxplot(rowMeans(length_sim < vb_jags_out$q50$Linf),
        rowMeans(length_sim < vb_jags_out$sims.list$Linf))
mean(length_sim < vb_jags_out$q50$Linf)

# Lester did "the mean fork length of the largest 10% in our sample,
# after removing fish smaller than 300 mm."
length_lester <- function(x) {
  x1 <- x[x >= 300]
  x2 <- x1[x1 >= quantile(x1, 0.9, na.rm=TRUE)]
  return(mean(x2, na.rm=TRUE))
}
mean(apply(length_sim, 1, length_lester))
mean(length_sim[length_sim > quantile(length_sim, 0.9)])



### Tried it again with the Chandler data, it didn't work

# rawtable <- table(chan_data$x) #%>% plot
# agetable <- rep(0, max(as.numeric(names(rawtable))))
# agetable[as.numeric(names(rawtable))] <- rawtable
# plot(agetable)
#
#
# mortrate_data <- list(y = agetable,
#                       age = 1:length(agetable),
#                       n = length(agetable))
#
# # JAGS controls
# niter <- 1000000
# # ncores <- 3
# ncores <- min(10, parallel::detectCores()-1)
#
# {
#   tstart <- Sys.time()
#   print(tstart)
#   mortrate_jags_out <- jagsUI::jags(model.file=mortrate_jags, data=mortrate_data,
#                                     parameters.to.save=c("lambda","pi","b0","b1","N","ypp","p"),
#                                     n.chains=ncores, parallel=T, n.iter=niter,
#                                     n.burnin=niter/2, n.thin=niter/2000)
#   print(Sys.time() - tstart)
# }
#
# # nbyname(mortrate_jags_out)
# plotRhats(mortrate_jags_out)
# traceworstRhat(mortrate_jags_out, parmfrow = c(3, 3))
# qq_postpred(mortrate_jags_out, p="ypp", y=mortrate_data$y)
# ts_postpred(mortrate_jags_out, p="ypp", y=mortrate_data$y)
# par(mfrow=c(3,2))
# envelope(mortrate_jags_out, "ypp")
# points(mortrate_data$y)
# envelope(mortrate_jags_out, p="p")
# envelope(outer(mortrate_jags_out$sims.list$lambda, mortrate_data$age, "^"),
#          main="lambda effect")
#
# envelope(outer(mortrate_jags_out$sims.list$lambda, mortrate_data$age, "^") *
#            mortrate_jags_out$sims.list$p,
#          main="lambda x p")
# envelope(outer(mortrate_jags_out$sims.list$lambda, mortrate_data$age, "^") *
#            mortrate_jags_out$sims.list$p *
#            mortrate_jags_out$sims.list$pi,
#          main="lambda x p x pi")
# plotcor_jags(mortrate_jags_out)




############ delete all this, it didn't work

# specify model, which is written to a temporary file
regtest1 <- tempfile()
cat('model {
  for(i in 1:n) {
    mu[i] <- b0 + b1*x[i]
    y[i] ~ dnorm(mu[i], tau)
  }

  tau <- pow(sig, -2)
  sig ~ dunif(0,10)
  b0 ~ dnorm(0, 0.001)
  b1 ~ dnorm(0, 0.001)
}', file=regtest1)

regtest2 <- tempfile()
cat('model {
  for(i in 1:n) {
    mu[i] <- b1*(b0 + x[i])
    y[i] ~ dnorm(mu[i], tau)
  }

  tau <- pow(sig, -2)
  sig ~ dunif(0,10)
  b0 ~ dnorm(0, 0.001)
  b1 ~ dnorm(0, 0.001)
}', file=regtest2)

regtest_data <- list(x = rnorm(100, mean=100, sd=10),
                     n = 100)
regtest_data$y <- rnorm(100, mean=regtest_data$x+100, sd=3)
plot(regtest_data)

# JAGS controls
niter <- 1000
# ncores <- 3
ncores <- min(10, parallel::detectCores()-1)

{
  tstart <- Sys.time()
  print(tstart)
  regtest1_jags_out <- jagsUI::jags(model.file=regtest1, data=regtest_data,
                                    parameters.to.save=c("b0", "b1","sig"),
                                    n.chains=ncores, parallel=T, n.iter=niter,
                                    n.burnin=niter/2, n.thin=niter/1000)
  print(Sys.time() - tstart)
}
plotRhats(regtest1_jags_out)
traceworstRhat(regtest1_jags_out, parmfrow = c(3, 3))
plotcor_jags(regtest1_jags_out)

{
  tstart <- Sys.time()
  print(tstart)
  regtest2_jags_out <- jagsUI::jags(model.file=regtest2, data=regtest_data,
                                    parameters.to.save=c("b0", "b1","sig"),
                                    n.chains=ncores, parallel=T, n.iter=niter,
                                    n.burnin=niter/2, n.thin=niter/1000)
  print(Sys.time() - tstart)
}
plotRhats(regtest2_jags_out)
traceworstRhat(regtest2_jags_out, parmfrow = c(3, 3))
plotcor_jags(regtest2_jags_out)
