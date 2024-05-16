## This script just fits a Von B growth model to a large Lake Trout
## length-age dataset I found online

## https://figshare.com/articles/dataset/Lake_Trout_growth_data_csv_file_/6957620

## I was curious how sample quantiles might compare to estimated L_inf

library(tidyverse)
library(jagsUI)
library(jagshelper)

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
ts_postpred(vb_jags_out, p="ypp", y=vb_data$y)

plot(vb_data)
envelope(vb_jags_out, p="mu", x=vb_data$x, add=TRUE)
plot(vb_data)
envelope(vb_jags_out, p="ypp", x=vb_data$x, add=TRUE)


## ok moment of truth, how do the length quantiles compare to modeled Linf??
mean(vb_data$y < vb_jags_out$q50$Linf)

Linf_post <- vb_jags_out$sims.list$Linf
quantile_post <- colMeans(outer(vb_data$y, Linf_post, "<"))
plotdens(quantile_post)



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
plot(chan_data)
envelope(chan_jags_out, p="ypp", x=chan_data$x, add=TRUE)


## ok moment of truth, how do the length quantiles compare to modeled Linf??
mean(chan_data$y < chan_jags_out$q50$Linf)

Linf_post <- chan_jags_out$sims.list$Linf
quantile_post <- colMeans(outer(chan_data$y, Linf_post, "<"))
plotdens(quantile_post)



## quick & dirty sim to look into how universal this 90% thing is

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




rawtable <- table(chan_data$x) #%>% plot
agetable <- rep(0, max(as.numeric(names(rawtable))))
agetable[as.numeric(names(rawtable))] <- rawtable
plot(agetable)


mortrate_data <- list(y = agetable,
                      age = 1:length(agetable),
                      n = length(agetable))

# JAGS controls
niter <- 1000000
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
