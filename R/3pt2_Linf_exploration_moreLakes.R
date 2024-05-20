## After the master lake trout dataset was updated, there are now several lakes
## with viable length & age measurements!!

## This script fits Von Bertalanffy growth models to the set of lakes with the
## most measurements.  Currently the cutoff is n=150, with a total of 4 lakes.

## Sample quantiles were then compared to the modeled L_inf.
## Interestingly, for these lakes, L_inf was at the 93-98th percentile values.

## Finally, a comparison of the 80th, 90th, and 95th percentile values, as well
## as Lester's method, for all lakes.  From this, it looks like there is much
## more difference between lakes than difference between method, so it might
## not matter super much which percentile value is chosen!


library(tidyverse)
library(jagsUI)
library(jagshelper)

source("R/1_laketrout_lwdata.R")



## filtering to a dataset only including the lakes with a minimum number of
## viable age x length samples (currently 150)

lake_tbl <- table(laketrout_all$LakeName, !is.na(laketrout_all$Age) & !is.na(laketrout_all$ForkLength_mm))

min_sample <- 150 # let's play with this
lakeswithmin <- rownames(lake_tbl)[lake_tbl[,2] >= min_sample]


## plotting

laketrout_all %>%
  filter(LakeName %in% lakeswithmin) %>%
  filter(Age < 50) %>%
  ggplot(aes(x=Age, y=ForkLength_mm, col=LakeName)) +
  geom_point() +
  theme_bw()



## defining VB length model

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

# JAGS controls
niter <- 200000   # biggerize this
# ncores <- 3
ncores <- min(10, parallel::detectCores()-1)

vb_outs <- list()
vb_datas <- list()



## first constructing a list of data objects (to be used later), one for each lake
## then fitting the VB model for each lake

for(i_lake in 1:length(lakeswithmin)) {
  lake_i <- filter(laketrout_all, LakeName == lakeswithmin[i_lake]) %>%
  filter(Age < 50)
  vb_datas[[i_lake]] <- list(x = lake_i$Age[!is.na(lake_i$Age) & !is.na(lake_i$ForkLength_mm)],
                  y = lake_i$ForkLength_mm[!is.na(lake_i$Age) & !is.na(lake_i$ForkLength_mm)],
                  n = sum(!is.na(lake_i$Age) & !is.na(lake_i$ForkLength_mm)))
}
for(i_lake in 1:length(lakeswithmin)) {
  {
    tstart <- Sys.time()
    print(tstart)
    vb_outs[[i_lake]] <- jagsUI::jags(model.file=vb_jags, data=vb_datas[[i_lake]],
                                parameters.to.save=c("sig","Linf","t0","k","mu","ypp",
                                                     "sig_prior","Linf_prior","t0_prior","k_prior"),
                                n.chains=ncores, parallel=T, n.iter=niter,
                                n.burnin=niter/2, n.thin=niter/2000)
    print(Sys.time() - tstart)
  }
}

## plotting some fit diagnostics & model output

par(mfrow=c(2,3))
for(i_lake in 1:length(vb_outs)) {
  plotRhats(vb_outs[[i_lake]])
}

for(i_lake in 1:length(vb_outs)) {
  par(mfrow=c(4,3))
  traceworstRhat(vb_outs[[i_lake]])
}

for(i_lake in 1:length(vb_outs)) {
  par(mfrow=c(2,2))
  comparepriors(vb_outs[[i_lake]])
}

par(mfrow=c(2,2))
comparecat(vb_outs, p="Linf")
comparecat(vb_outs, p="k", ylim=0:1)
comparecat(vb_outs, p="t0")
comparecat(vb_outs, p="sig")


## comparing inferences between lakes.  The proportion of the length sample
## greater than median L_inf is printed to the console for each lake.

maxx <- max(sapply(vb_datas, function(x) max(x$x)))
maxy <- max(sapply(vb_datas, function(x) max(x$y)))
par(mfrow=c(2,2))
for(i_lake in 1:length(vb_outs)) {
  plot(vb_datas[[i_lake]]$x, vb_datas[[i_lake]]$y, xlim=c(0,maxx), ylim=c(0,maxy))
  envelope(vb_outs[[i_lake]], p="ypp", add=TRUE, x=vb_datas[[i_lake]]$x)
  abline(h=vb_outs[[i_lake]]$q50$Linf, lty=3)
  print(mean(vb_datas[[i_lake]]$y < vb_outs[[i_lake]]$q50$Linf))
}



## finally, visualizing the differences between 80th, 90th, and 95th percentile
## values for each lake, as well as Lester's method.

## note that Lester's method seems to be sensitive to outliers in the positive direction!

lake_tbl2 <- table(laketrout_all$LakeName, !is.na(laketrout_all$ForkLength_mm))
minn_quantile <- 50  # play with this
lakeswithmin2 <- rownames(lake_tbl2)[lake_tbl2[,2]>=minn_quantile]

lakeswithlengths <- filter(laketrout_all, LakeName %in% lakeswithmin2)
plot(tapply(lakeswithlengths$ForkLength_mm, lakeswithlengths$LakeName, quantile, 0.95, na.rm=TRUE))
points(tapply(lakeswithlengths$ForkLength_mm, lakeswithlengths$LakeName, quantile, 0.90, na.rm=TRUE))
points(tapply(lakeswithlengths$ForkLength_mm, lakeswithlengths$LakeName, quantile, 0.80, na.rm=TRUE))
segments(x0=seq_along(lakeswithmin2),
         y0=tapply(lakeswithlengths$ForkLength_mm, lakeswithlengths$LakeName, quantile, 0.95, na.rm=TRUE),
          y1=tapply(lakeswithlengths$ForkLength_mm, lakeswithlengths$LakeName, quantile, 0.80, na.rm=TRUE))
# Lester did "the mean fork length of the largest 10% in our sample,
# after removing fish smaller than 300 mm."
length_lester <- function(x) {
  x1 <- x[x >= 300]
  x2 <- x1[x1 >= quantile(x1, 0.9, na.rm=TRUE)]
  return(mean(x2, na.rm=TRUE))
}
points(tapply(lakeswithlengths$ForkLength_mm, lakeswithlengths$LakeName, length_lester), pch=16, col=2)
