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
  sig ~ dunif(0, 100)
  sig_prior ~ dunif(0, 100)
  Linf ~ dnorm(700, 0.0001)
  Linf_prior ~ dnorm(700, 0.0001)
  t0 ~ dnorm(0, 0.01)
  t0_prior ~ dnorm(0, 0.01)
  k ~ dnorm(0, 0.01)
  k_prior ~ dnorm(0, 0.01)
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
