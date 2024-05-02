# source("R/1_laketrout_lwdata.R")


# bundle data to pass into JAGS
laketrout_lw <- filter(laketrout_lw, !is.na(SurfaceArea_h))

logarea <- log(with(laketrout_lw, tapply(SurfaceArea_h, LakeName, median, na.rm=T)))
lat <- with(laketrout_lw, tapply(Latitude_WGS84, LakeName, median))

lt_data <- list(x=log(laketrout_lw$ForkLength_mm) - mean(log(laketrout_lw$ForkLength_mm)),
                y=log(laketrout_lw$Weight_g/1000),
                n=nrow(laketrout_lw),
                lake=as.numeric(as.factor(laketrout_lw$LakeName)),
                nlake=length(unique(laketrout_lw$LakeName)),
                meanx = mean(log(laketrout_lw$ForkLength_mm)),
                lat=lat - mean(lat),
                area=logarea - mean(logarea))
lt_data$xpred <- seq(from=min(lt_data$x), to=max(lt_data$x), length.out=30)
lt_data$npred <- length(lt_data$xpred)

## Ultimately, what I want to compare is the following:

## - (1) hierarchical intercept, hierarchical slope
## - (2) hierarchical intercept, trending slope (lat)
## - (3) hierarchical intercept, trending slope (logarea)
## - (4) hierarchical intercept, trending slope (lat & logarea)
## - (5) trending intercept (lat), hierarchical slope
## - (6) trending intercept (logarea), hierarchical slope
## - (7) trending intercept (lat & logarea), hierarchical slope
## - (8) trending intercept (lat), trending slope (lat)
## - (9) trending intercept (lat), trending slope (logarea)
## - (10) trending intercept (lat), trending slope (lat & logarea)
## - (11) trending intercept (logarea), trending slope (lat)
## - (12) trending intercept (logarea), trending slope (logarea)
## - (13) trending intercept (logarea), trending slope (lat & logarea)
## - (14) trending intercept (lat & logarea), trending slope (lat)
## - (15) trending intercept (lat & logarea), trending slope (logarea)
## - (16) trending intercept (lat & logarea), trending slope (lat & logarea)

## write this as a matrix of ones and zeroes
## actually this is way easier
controlmat <- expand.grid(0:1,0:1,0:1,0:1)

intmodel <- slopemodel <- rep(NA, nrow(controlmat))
intmodel[controlmat[,1]==0 & controlmat[,2]==0] <- "hierarchical"
intmodel[controlmat[,1]==1 & controlmat[,2]==0] <- "trends with lat"
intmodel[controlmat[,1]==0 & controlmat[,2]==1] <- "trends with area"
intmodel[controlmat[,1]==1 & controlmat[,2]==1] <- "trends with lat & area"
slopemodel[controlmat[,3]==0 & controlmat[,4]==0] <- "hierarchical"
slopemodel[controlmat[,3]==1 & controlmat[,4]==0] <- "trends with lat"
slopemodel[controlmat[,3]==0 & controlmat[,4]==1] <- "trends with area"
slopemodel[controlmat[,3]==1 & controlmat[,4]==1] <- "trends with lat & area"
modeldescription <- paste0("intercept ", intmodel,", ", "slope ", slopemodel)
modeldescription


## write a model that can be extended to all possible configurations
## - supply data that will turn on/off trending terms as needed

lt_jags <- tempfile()
cat('model {
  for(i in 1:n) {
    y[i] ~ dnorm(mu[i], tau)
    ypp[i] ~ dnorm(mu[i], tau)
    mu[i] <- b0[lake[i]] + b1[lake[i]]*x[i]
  }

  for(j in 1:nlake) {
    b0[j] ~ dnorm(mu_b0[j], tau_b0)
    mu_b0[j] <- b0_int
                + b0_lat*lat[j]   * intlat
                + b0_area*area[j] * intarea

    b1[j] ~ dnorm(mu_b1[j], tau_b1)
    mu_b1[j] <- b1_int
                + b1_lat*lat[j]   * slopelat
                + b1_area*area[j] * slopearea
  }

  sig_b0 ~ dunif(0, 10)
  tau_b0 <- pow(sig_b0, -2)

  sig_b1 ~ dunif(0, 10)
  tau_b1 <- pow(sig_b1, -2)

  b0_int ~ dnorm(0, 0.001)
  b0_lat ~ dnorm(0, 0.001)
  b0_area ~ dnorm(0, 0.001)

  b1_int ~ dnorm(0, 0.001)
  b1_lat ~ dnorm(0, 0.001)
  b1_area ~ dnorm(0, 0.001)

  tau <- pow(sig, -2)
  sig ~ dunif(0, 10)

}', file=lt_jags)


lt_jags_commonint <- tempfile()
cat('model {
  for(i in 1:n) {
    y[i] ~ dnorm(mu[i], tau)
    ypp[i] ~ dnorm(mu[i], tau)
    mu[i] <- b0[lake[i]] + b1[lake[i]]*x[i]
  }

  for(j in 1:nlake) {
    b0[j] <- b0_int
    # b0[j] ~ dnorm(mu_b0[j], tau_b0)
    # mu_b0[j] <- b0_int
    #             + b0_lat*lat[j]   * intlat
    #             + b0_area*area[j] * intarea

    b1[j] ~ dnorm(mu_b1[j], tau_b1)
    mu_b1[j] <- b1_int
                + b1_lat*lat[j]   * slopelat
                + b1_area*area[j] * slopearea
  }

  sig_b0 ~ dunif(0, 10)
  tau_b0 <- pow(sig_b0, -2)

  sig_b1 ~ dunif(0, 10)
  tau_b1 <- pow(sig_b1, -2)

  b0_int ~ dnorm(0, 0.001)
  b0_lat ~ dnorm(0, 0.001)
  b0_area ~ dnorm(0, 0.001)

  b1_int ~ dnorm(0, 0.001)
  b1_lat ~ dnorm(0, 0.001)
  b1_area ~ dnorm(0, 0.001)

  tau <- pow(sig, -2)
  sig ~ dunif(0, 10)

}', file=lt_jags_commonint)


lt_jags_commonslope <- tempfile()
cat('model {
  for(i in 1:n) {
    y[i] ~ dnorm(mu[i], tau)
    ypp[i] ~ dnorm(mu[i], tau)
    mu[i] <- b0[lake[i]] + b1[lake[i]]*x[i]
  }

  for(j in 1:nlake) {
    b0[j] ~ dnorm(mu_b0[j], tau_b0)
    mu_b0[j] <- b0_int
                + b0_lat*lat[j]   * intlat
                + b0_area*area[j] * intarea

    b1[j] <- b1_int
    # b1[j] ~ dnorm(mu_b1[j], tau_b1)
    # mu_b1[j] <- b1_int
    #             + b1_lat*lat[j]   * slopelat
    #             + b1_area*area[j] * slopearea
  }

  sig_b0 ~ dunif(0, 10)
  tau_b0 <- pow(sig_b0, -2)

  sig_b1 ~ dunif(0, 10)
  tau_b1 <- pow(sig_b1, -2)

  b0_int ~ dnorm(0, 0.001)
  b0_lat ~ dnorm(0, 0.001)
  b0_area ~ dnorm(0, 0.001)

  b1_int ~ dnorm(0, 0.001)
  b1_lat ~ dnorm(0, 0.001)
  b1_area ~ dnorm(0, 0.001)

  tau <- pow(sig, -2)
  sig ~ dunif(0, 10)

}', file=lt_jags_commonslope)


lt_jags_common <- tempfile()
cat('model {
  for(i in 1:n) {
    y[i] ~ dnorm(mu[i], tau)
    ypp[i] ~ dnorm(mu[i], tau)
    mu[i] <- b0[lake[i]] + b1[lake[i]]*x[i]
  }

  for(j in 1:nlake) {
    b0[j] <- b0_int
    # b0[j] ~ dnorm(mu_b0[j], tau_b0)
    # mu_b0[j] <- b0_int
    #             + b0_lat*lat[j]   * intlat
    #             + b0_area*area[j] * intarea

    b1[j] <- b1_int
    # b1[j] ~ dnorm(mu_b1[j], tau_b1)
    # mu_b1[j] <- b1_int
    #             + b1_lat*lat[j]   * slopelat
    #             + b1_area*area[j] * slopearea
  }

  sig_b0 ~ dunif(0, 10)
  tau_b0 <- pow(sig_b0, -2)

  sig_b1 ~ dunif(0, 10)
  tau_b1 <- pow(sig_b1, -2)

  b0_int ~ dnorm(0, 0.001)
  b0_lat ~ dnorm(0, 0.001)
  b0_area ~ dnorm(0, 0.001)

  b1_int ~ dnorm(0, 0.001)
  b1_lat ~ dnorm(0, 0.001)
  b1_area ~ dnorm(0, 0.001)

  tau <- pow(sig, -2)
  sig ~ dunif(0, 10)

}', file=lt_jags_common)


## will compare
## - DIC
## - sig
## - b0_interp, b0, b1


# JAGS controls
niter <- 50000
ncores <- min(10, parallel::detectCores()-1)
par(mfrow=c(4,4))


### then build a master model and take stuff out

alloutputs <- list()

for(imodel in 1:nrow(controlmat)) {
  lt_data$intlat <- controlmat[imodel, 1]
  lt_data$intarea <- controlmat[imodel, 2]
  lt_data$slopelat <- controlmat[imodel, 3]
  lt_data$slopearea <- controlmat[imodel, 4]

  {
    tstart <- Sys.time()
    print(tstart)
    lt_jags_out <- jagsUI::jags(model.file=lt_jags, data=lt_data,
                                parameters.to.save=c("b0","b1","sig",
                                                     "mu_b0","sig_b0","mu_b1","sig_b1",
                                                     "fit","pred","ypp","b0_interp",
                                                     "b0_int","b1_int",
                                                     "b0_lat","b0_area",
                                                     "b1_lat","b1_area"),
                                n.chains=ncores, parallel=T, n.iter=niter,
                                n.burnin=niter/2, n.thin=niter/2000)
    print(Sys.time() - tstart)
    plotRhats(lt_jags_out)
    alloutputs[[imodel]] <- lt_jags_out
  }
}

## appending models with common intercept
controlmat <- rbind(controlmat,
                    expand.grid(0, 0, 0:1, 0:1))
for(imodel in (imodel+1):nrow(controlmat)) {
  lt_data$intlat <- controlmat[imodel, 1]
  lt_data$intarea <- controlmat[imodel, 2]
  lt_data$slopelat <- controlmat[imodel, 3]
  lt_data$slopearea <- controlmat[imodel, 4]

  {
    tstart <- Sys.time()
    print(tstart)
    lt_jags_out <- jagsUI::jags(model.file=lt_jags_commonint, data=lt_data,
                                parameters.to.save=c("b0","b1","sig",
                                                     "mu_b0","sig_b0","mu_b1","sig_b1",
                                                     "fit","pred","ypp","b0_interp",
                                                     "b0_int","b1_int",
                                                     "b0_lat","b0_area",
                                                     "b1_lat","b1_area"),
                                n.chains=ncores, parallel=T, n.iter=niter,
                                n.burnin=niter/2, n.thin=niter/2000)
    print(Sys.time() - tstart)
    plotRhats(lt_jags_out)
    alloutputs[[imodel]] <- lt_jags_out
  }
}

## appending models with common slope
controlmat <- rbind(controlmat,
                    expand.grid(0:1, 0:1, 0, 0))
for(imodel in (imodel+1):nrow(controlmat)) {
  lt_data$intlat <- controlmat[imodel, 1]
  lt_data$intarea <- controlmat[imodel, 2]
  lt_data$slopelat <- controlmat[imodel, 3]
  lt_data$slopearea <- controlmat[imodel, 4]

  {
    tstart <- Sys.time()
    print(tstart)
    lt_jags_out <- jagsUI::jags(model.file=lt_jags_commonslope, data=lt_data,
                                parameters.to.save=c("b0","b1","sig",
                                                     "mu_b0","sig_b0","mu_b1","sig_b1",
                                                     "fit","pred","ypp","b0_interp",
                                                     "b0_int","b1_int",
                                                     "b0_lat","b0_area",
                                                     "b1_lat","b1_area"),
                                n.chains=ncores, parallel=T, n.iter=niter,
                                n.burnin=niter/2, n.thin=niter/2000)
    print(Sys.time() - tstart)
    plotRhats(lt_jags_out)
    alloutputs[[imodel]] <- lt_jags_out
  }
}

## appending model with common all
controlmat <- rbind(controlmat,
                    expand.grid(0, 0, 0, 0))
for(imodel in (imodel+1):nrow(controlmat)) {
  lt_data$intlat <- controlmat[imodel, 1]
  lt_data$intarea <- controlmat[imodel, 2]
  lt_data$slopelat <- controlmat[imodel, 3]
  lt_data$slopearea <- controlmat[imodel, 4]

  {
    tstart <- Sys.time()
    print(tstart)
    lt_jags_out <- jagsUI::jags(model.file=lt_jags_common, data=lt_data,
                                parameters.to.save=c("b0","b1","sig",
                                                     "mu_b0","sig_b0","mu_b1","sig_b1",
                                                     "fit","pred","ypp","b0_interp",
                                                     "b0_int","b1_int",
                                                     "b0_lat","b0_area",
                                                     "b1_lat","b1_area"),
                                n.chains=ncores, parallel=T, n.iter=niter,
                                n.burnin=niter/2, n.thin=niter/2000)
    print(Sys.time() - tstart)
    plotRhats(lt_jags_out)
    alloutputs[[imodel]] <- lt_jags_out
  }
}

theDIC <- sapply(alloutputs, function(x) x$DIC)[1:16]
theDIC
theDIC - min(theDIC)
# [1] 2.25661808 1.89560947 0.26957872 0.83358178 1.98750118 3.20484979 0.00000000 0.02992762
# [9] 1.61076528 1.17305266 1.22793261 0.86946090 1.01976078 1.64548233 0.74567197 0.46453610

modeldescription[which.min(theDIC)]

par(mfrow=c(1, 1))
plot(x=theDIC - min(theDIC), y=rank(theDIC),
     xlim=c(0, 10))
text(x=theDIC - min(theDIC), y=rank(theDIC),
     labels=modeldescription, pos=4)

# par(mfrow=c(1,1))
# comparecat(alloutputs, p="sig")
# thesigs <- sapply(alloutputs, function(x) jags_df(x, p="sig", exact=T))
# thesigs %>% as.data.frame %>% caterpillar

par(mfrow=c(2,2))
comparecat(alloutputs[controlmat[,1]==1], p="b0_lat")
abline(h=0, lty=2)
comparecat(alloutputs[controlmat[,2]==1], p="b0_area")
abline(h=0, lty=2)
comparecat(alloutputs[controlmat[,3]==1], p="b1_lat")
abline(h=0, lty=2)
comparecat(alloutputs[controlmat[,4]==1], p="b1_area")
abline(h=0, lty=2)



par(mfrow=c(2,2))
caterpillar(alloutputs[[1]], p="b0", x=lat, xlab="latitude")
caterpillar(alloutputs[[1]], p="b1", x=lat, xlab="latitude")
caterpillar(alloutputs[[1]], p="b0", x=logarea, xlab="log area")
caterpillar(alloutputs[[1]], p="b1", x=logarea, xlab="log area")

caterpillar(alloutputs[[7]], p="b0", x=lat, xlab="latitude")
caterpillar(alloutputs[[7]], p="b1", x=lat, xlab="latitude")
envelope(alloutputs[[7]], p="mu_b1", x=lat, add=TRUE)
caterpillar(alloutputs[[7]], p="b0", x=logarea, xlab="log area")
envelope(alloutputs[[7]], p="mu_b0", x=logarea, add=TRUE)
caterpillar(alloutputs[[7]], p="b1", x=logarea, xlab="log area")

caterpillar(alloutputs[[16]], p="b0", x=lat, xlab="latitude")
caterpillar(alloutputs[[16]], p="b1", x=lat, xlab="latitude")
caterpillar(alloutputs[[16]], p="b0", x=logarea, xlab="log area")
caterpillar(alloutputs[[16]], p="b1", x=logarea, xlab="log area")

comparecat(alloutputs[c(1,7,16)], p="b0", ylim=c(-0.4,0.3))
comparecat(alloutputs[c(1,7,16)], p="b1", ylim=c(2.4,4.2))

for(imodel in c(1,7,16)) {
  print(paste("model",imodel))
  print(mean(alloutputs[[imodel]]$sims.list$b0_lat > 0))
  print(mean(alloutputs[[imodel]]$sims.list$b0_area > 0))
  print(mean(alloutputs[[imodel]]$sims.list$b1_lat > 0))
  print(mean(alloutputs[[imodel]]$sims.list$b1_area > 0))
}



data.frame(modeldescription,
           deltaDIC=theDIC - min(theDIC))
