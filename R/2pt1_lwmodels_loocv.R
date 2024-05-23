## loading all data
source("R/1_laketrout_lwdata.R")

## load full model runs for comparison
load(file="alloutputs.Rdata")
theDICs <- sapply(alloutputs, function(x) x$DIC)

## loop over lakes in base dataset
  ## take out each lake in turn and redefine data bundle
  ## loop over best models
    ## store intercept & slope prediction


# define global data
laketrout_lw <- filter(laketrout_lw, !is.na(SurfaceArea_h))

logarea_all <- log(with(laketrout_lw, tapply(SurfaceArea_h, LakeName, median, na.rm=T)))
lat_all <- with(laketrout_lw, tapply(Latitude_WGS84, LakeName, median))
lakenames_all <- names(lat_all)



## defining model as a matrix of ones and zeroes
## actually this is way easier!
controlmat <- expand.grid(0:1,0:1,0:1,0:1)
bestmodels <- which(theDICs %in% tapply(theDICs, rowSums(controlmat), min))


## create a vector describing the model from the control matrix of ones and zeroes
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
# [1] "intercept hierarchical, slope hierarchical"
# [2] "intercept trends with lat, slope hierarchical"
# [3] "intercept trends with area, slope hierarchical"
# [4] "intercept trends with lat & area, slope hierarchical"
# [5] "intercept hierarchical, slope trends with lat"
# [6] "intercept trends with lat, slope trends with lat"
# [7] "intercept trends with area, slope trends with lat"
# [8] "intercept trends with lat & area, slope trends with lat"
# [9] "intercept hierarchical, slope trends with area"
# [10] "intercept trends with lat, slope trends with area"
# [11] "intercept trends with area, slope trends with area"
# [12] "intercept trends with lat & area, slope trends with area"
# [13] "intercept hierarchical, slope trends with lat & area"
# [14] "intercept trends with lat, slope trends with lat & area"
# [15] "intercept trends with area, slope trends with lat & area"
# [16] "intercept trends with lat & area, slope trends with lat & area"


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
    b0_interp[j] <- b0[j] - b1[j]*meanx
    mu_b0[j] <- b0_int
                + b0_lat*lat[j]   * intlat
                + b0_area*area[j] * intarea

    b1[j] ~ dnorm(mu_b1[j], tau_b1)
    mu_b1[j] <- b1_int
                + b1_lat*lat[j]   * slopelat
                + b1_area*area[j] * slopearea
  }

  new_b0 ~ dnorm(new_mu_b0, tau_b0)
  new_b0_interp <- new_b0 - new_b1*meanx
  new_mu_b0 <- b0_int
                + b0_lat*newlat   * intlat
                + b0_area*newarea * intarea
  new_b1 ~ dnorm(new_mu_b1, tau_b1)
  new_mu_b1 <- b1_int
                + b1_lat*newlat   * slopelat
                + b1_area*newarea * slopearea

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

## initialize place to put the output
b0_interp_arr <- b0_arr <- b1_arr <- mu_b0_arr <- mu_b1_arr <- array(dim=c(nrow(alloutputs[[1]]$sims.list$b0),
                                                          length(lakenames_all),
                                                          length(bestmodels)))

# JAGS controls
niter <- 10*1000#20*1000
ncores <- min(10, parallel::detectCores()-1)
par(mfrow=c(4,4))

for(ilake in seq_along(lakenames_all)) {

  laketrout_loocv <- filter(laketrout_lw, LakeName != lakenames_all[ilake])
  logarea_loocv <- log(with(laketrout_loocv, tapply(SurfaceArea_h, LakeName, median, na.rm=T)))
  lat_loocv <- with(laketrout_loocv, tapply(Latitude_WGS84, LakeName, median))


lt_data <- list(x=log(laketrout_loocv$ForkLength_mm) - mean(log(laketrout_lw$ForkLength_mm)),
                # x=log(laketrout_loocv$ForkLength_mm) - mean(log(laketrout_loocv$ForkLength_mm)),
                y=log(laketrout_loocv$Weight_g/1000),
                n=nrow(laketrout_loocv),
                lake=as.numeric(as.factor(laketrout_loocv$LakeName)),
                nlake=length(unique(laketrout_loocv$LakeName)),
                meanx = mean(log(laketrout_lw$ForkLength_mm)),
                # meanx = mean(log(laketrout_loocv$ForkLength_mm)),
                lat=lat_loocv - mean(lat_loocv),
                area=logarea_loocv - mean(logarea_loocv),
                newlat = lat_all[ilake] - mean(lat_loocv),
                newarea = logarea_all[ilake] - mean(logarea_loocv))






### defining a list to put all model outputs in
# alloutputs <- list()

for(imodel in seq_along(bestmodels)) {
  lt_data$intlat <- controlmat[bestmodels[imodel], 1]
  lt_data$intarea <- controlmat[bestmodels[imodel], 2]
  lt_data$slopelat <- controlmat[bestmodels[imodel], 3]
  lt_data$slopearea <- controlmat[bestmodels[imodel], 4]

  {
    tstart <- Sys.time()
    print(tstart)
    print(ilake)
    lt_jags_out <- jagsUI::jags(model.file=lt_jags, data=lt_data,
                                parameters.to.save=c("b0","b1","sig",
                                                     "mu_b0","sig_b0","mu_b1","sig_b1",
                                                     "fit","pred","b0_interp",
                                                     "b0_int","b1_int",
                                                     "b0_lat","b0_area",
                                                     "b1_lat","b1_area",
                                                     "new_b0","new_b1","new_b0_interp",
                                                     "new_mu_b0","new_mu_b1"),#,"ypp"
                                n.chains=ncores, parallel=T, n.iter=niter,
                                n.burnin=niter/2, n.thin=niter/2000)
    print(Sys.time() - tstart)
    plotRhats(lt_jags_out)
    # alloutputs[[imodel]] <- lt_jags_out

    b0_arr[,ilake,imodel] <- lt_jags_out$sims.list$new_b0
    b0_interp_arr[,ilake,imodel] <- lt_jags_out$sims.list$new_b0_interp
    b1_arr[,ilake,imodel] <- lt_jags_out$sims.list$new_b1
    mu_b0_arr[,ilake,imodel] <- lt_jags_out$sims.list$new_mu_b0
    mu_b1_arr[,ilake,imodel] <- lt_jags_out$sims.list$new_mu_b1
  }
}
}

########################################################
# save(b0_arr, b0_interp_arr, b1_arr, mu_b0_arr, mu_b1_arr, file="lw_loocv.Rdata")

par(mfrow=c(1,1))
comparecat(list(as.data.frame(b0_arr[,,1]),
                as.data.frame(b0_arr[,,2]),
                as.data.frame(b0_arr[,,3]),
                as.data.frame(b0_arr[,,4]),
                as.data.frame(b0_arr[,,5])), col=1:5)
for(i in seq_along(bestmodels)) lines(alloutputs[[bestmodels[i]]]$q50$b0)

comparecat(list(as.data.frame(b0_interp_arr[,,1]),
                as.data.frame(b0_interp_arr[,,2]),
                as.data.frame(b0_interp_arr[,,3]),
                as.data.frame(b0_interp_arr[,,4]),
                as.data.frame(b0_interp_arr[,,5])), col=1:5)
for(i in seq_along(bestmodels)) lines(alloutputs[[bestmodels[i]]]$q50$b0_interp)

comparecat(list(as.data.frame(b1_arr[,,1]),
                as.data.frame(b1_arr[,,2]),
                as.data.frame(b1_arr[,,3]),
                as.data.frame(b1_arr[,,4]),
                as.data.frame(b1_arr[,,5])), col=1:5)
for(i in seq_along(bestmodels)) lines(alloutputs[[bestmodels[i]]]$q50$b1)

comparecat(list(as.data.frame(mu_b0_arr[,,1]),
                as.data.frame(mu_b0_arr[,,2]),
                as.data.frame(mu_b0_arr[,,3]),
                as.data.frame(mu_b0_arr[,,4]),
                as.data.frame(mu_b0_arr[,,5])), col=1:5)
for(i in seq_along(bestmodels)) lines(alloutputs[[bestmodels[i]]]$q50$b0)
for(i in seq_along(bestmodels)) lines(alloutputs[[bestmodels[i]]]$q50$mu_b0, lty=2)

comparecat(list(as.data.frame(mu_b1_arr[,,1]),
                as.data.frame(mu_b1_arr[,,2]),
                as.data.frame(mu_b1_arr[,,3]),
                as.data.frame(mu_b1_arr[,,4]),
                as.data.frame(mu_b1_arr[,,5])), col=1:5)
for(i in seq_along(bestmodels)) lines(alloutputs[[bestmodels[i]]]$q50$b1)
for(i in seq_along(bestmodels)) lines(alloutputs[[bestmodels[i]]]$q50$mu_b1, lty=2)


## plot pp quantiles for new arrays vs truth (figure out what truth is!)
# let's say truth is naive (independently  hierarchical) model
bestmodels
# [1]  1  3  7 12 16
truthmodel <- 1

b0_truth <- alloutputs[[truthmodel]]$q50$b0
b0_qq <- (b0_arr < array(matrix(b0_truth, nrow=dim(b0_arr)[1], ncol=length(b0_truth),
                                byrow=TRUE), dim=dim(b0_arr))) %>% apply(c(2,3), mean)
par(mfrow=c(2,3))
for(i in 1:5) {
  plot(NA, xlim=c(1,22), ylim=0:1)
  points(sort(b0_qq[,i]))
  abline(0,1/22)
}

b0_interp_truth <- alloutputs[[truthmodel]]$q50$b0_interp
b0_interp_qq <- (b0_interp_arr < array(matrix(b0_interp_truth, nrow=dim(b0_interp_arr)[1], ncol=length(b0_interp_truth),
                                byrow=TRUE), dim=dim(b0_interp_arr))) %>% apply(c(2,3), mean)
par(mfrow=c(2,3))
for(i in 1:5) {
  plot(NA, xlim=c(1,22), ylim=0:1)
  points(sort(b0_interp_qq[,i]))
  abline(0,1/22)
}

b1_truth <- alloutputs[[truthmodel]]$q50$b1
b1_qq <- (b1_arr < array(matrix(b1_truth, nrow=dim(b1_arr)[1], ncol=length(b1_truth),
                                byrow=TRUE), dim=dim(b1_arr))) %>% apply(c(2,3), mean)
par(mfrow=c(2,3))
for(i in 1:5) {
  plot(NA, xlim=c(1,22), ylim=0:1)
  points(sort(b1_qq[,i]))
  abline(0,1/22)
}

mu_b0_truth <- alloutputs[[truthmodel]]$q50$mu_b0
mu_b0_qq <- (mu_b0_arr < array(matrix(mu_b0_truth, nrow=dim(mu_b0_arr)[1], ncol=length(mu_b0_truth),
                                      byrow=TRUE), dim=dim(mu_b0_arr))) %>% apply(c(2,3), mean)
par(mfrow=c(2,3))
for(i in 1:5) {
  plot(NA, xlim=c(1,22), ylim=0:1)
  points(sort(mu_b0_qq[,i]))
  abline(0,1/22)
}

mu_b1_truth <- alloutputs[[truthmodel]]$q50$mu_b1
mu_b1_qq <- (mu_b1_arr < array(matrix(mu_b1_truth, nrow=dim(mu_b1_arr)[1], ncol=length(mu_b1_truth),
                                      byrow=TRUE), dim=dim(mu_b1_arr))) %>% apply(c(2,3), mean)

mu_b1_qq <- matrix(nrow=length(mu_b1_truth), ncol=length(bestmodels))
for(ilake in 1:length(mu_b1_truth)) {
  for(imodel in 1:length(bestmodels)) {
    mu_b1_qq[ilake, imodel] <- mean(mu_b1_arr[, ilake, imodel] < mu_b1_truth[ilake])
  }
}

par(mfrow=c(2,3))
for(i in 1:5) {
  plot(NA, xlim=c(1,22), ylim=0:1)
  points(sort(mu_b1_qq[,i]))
  abline(0,1/22)
}


## plot predictive accuracy of new array (medians) vs truth
par(mfrow=c(3,3))
b0_medians <- apply(b0_arr, 2:3, median)
boxplot(b0_medians-as.vector(b0_truth))
apply(b0_medians-as.vector(b0_truth), 2, function(x) sqrt(mean(x^2))) %>% plot
apply(b0_medians-as.vector(b0_truth), 2, function(x) sqrt(mean(abs(x)))) %>% plot

b0_interp_medians <- apply(b0_interp_arr, 2:3, median)
boxplot(b0_interp_medians-as.vector(b0_interp_truth))
apply(b0_interp_medians-as.vector(b0_interp_truth), 2, function(x) sqrt(mean(x^2))) %>% plot
apply(b0_interp_medians-as.vector(b0_interp_truth), 2, function(x) sqrt(mean(abs(x)))) %>% plot

b1_medians <- apply(b1_arr, 2:3, median)
boxplot(b1_medians-as.vector(b1_truth))
apply(b1_medians-as.vector(b1_truth), 2, function(x) sqrt(mean(x^2))) %>% plot
apply(b1_medians-as.vector(b1_truth), 2, function(x) sqrt(mean(abs(x)))) %>% plot

plot(b0_truth[order(b0_truth)], b0_medians[,1][order(b0_truth)], type='l',
     ylim=range(b0_medians), lwd=2)
for(i in 2:5) lines(b0_truth[order(b0_truth)], b0_medians[,i][order(b0_truth)], col=i, lwd=2)
abline(0,1,lty=2)
plot(b0_interp_truth[order(b0_interp_truth)], b0_interp_medians[,1][order(b0_interp_truth)], type='l',
     ylim=range(b0_interp_medians), lwd=2)
for(i in 2:5) lines(b0_interp_truth[order(b0_interp_truth)], b0_interp_medians[,i][order(b0_interp_truth)], col=i, lwd=2)
abline(0,1,lty=2)
plot(b1_truth[order(b1_truth)], b1_medians[,1][order(b1_truth)], type='l',
     ylim=range(b1_medians), lwd=2)
for(i in 2:5) lines(b1_truth[order(b1_truth)], b1_medians[,i][order(b1_truth)], col=i, lwd=2)
abline(0,1,lty=2)


# b0_truths <- sapply(alloutputs[bestmodels], function(x) x$q50$b0)
# b0_interp_truths <- sapply(alloutputs[bestmodels], function(x) x$q50$b0_interp)
# b1_truths <- sapply(alloutputs[bestmodels], function(x) x$q50$b1)
b0_truths <- sapply(alloutputs, function(x) x$q50$b0)
b0_interp_truths <- sapply(alloutputs, function(x) x$q50$b0_interp)
b1_truths <- sapply(alloutputs, function(x) x$q50$b1)

# all models, all truths
b0_rmse <- b0_interp_rmse <- b1_rmse <- matrix(nrow=ncol(b0_truths), ncol=5)
for(itruth in 1:ncol(b0_truths)) {
  for(imodel in 1:5) {
    b0_rmse[itruth, imodel] <- sqrt(mean((b0_truths[,itruth]-b0_medians[,imodel])^2))
    b0_interp_rmse[itruth, imodel] <- sqrt(mean((b0_interp_truths[,itruth]-b0_interp_medians[,imodel])^2))
    b1_rmse[itruth, imodel] <- sqrt(mean((b1_truths[,itruth]-b1_medians[,imodel])^2))
  }
}
boxplot(b0_rmse)
boxplot(b0_interp_rmse)
boxplot(b1_rmse)

plot(b0_rmse[1,], type="l", ylim=range(b0_rmse))
for(i in 1:ncol(b0_truths)) lines(b0_rmse[i,], lwd=2, col=i)
plot(b0_interp_rmse[1,], type="l", ylim=range(b0_interp_rmse))
for(i in 1:ncol(b0_truths)) lines(b0_interp_rmse[i,], lwd=2, col=i)
plot(b1_rmse[1,], type="l", ylim=range(b1_rmse))
for(i in 1:ncol(b0_truths)) lines(b1_rmse[i,], lwd=2, col=i)


## better still would be residual sd of DATA!
pred_r <- pred_sd <- pred_rmse <- matrix(nrow=length(bestmodels), ncol=length(lakenames_all))
for(imodel in seq_along(bestmodels)) {
  for(ilake in seq_along(lakenames_all)) {
    # ypred <- b0_interp_medians[ilake, imodel] +
    #   b1_medians[ilake, imodel]*log(laketrout_lw$ForkLength_mm[laketrout_lw$LakeName==lakenames_all[ilake]])
    ypred <- b0_medians[ilake, imodel] +
      b1_medians[ilake, imodel]*
      (log(laketrout_lw$ForkLength_mm[laketrout_lw$LakeName==lakenames_all[ilake]]) -
         mean(log(laketrout_lw$ForkLength_mm)))
    ytrue <- log(laketrout_lw$Weight_g[laketrout_lw$LakeName==lakenames_all[ilake]]/1000)
    pred_sd[imodel, ilake] <- sd(ypred-ytrue)  # should this be rmse instead??
    pred_rmse[imodel, ilake] <- sqrt(mean(((ypred-ytrue)^2)))
    pred_r[imodel, ilake] <- cor(ypred, ytrue)
  }
}
boxplot(pred_sd)
boxplot(t(pred_sd))
boxplot(pred_rmse)
boxplot(t(pred_rmse))
boxplot(pred_r)
boxplot(t(pred_r))

apply(pred_rmse, 1, median) %>% plot
apply(pred_rmse, 1, mean) %>% plot
apply(pred_r, 1, median) %>% plot
apply(pred_r, 1, mean) %>% plot

plot(pred_rmse[,1], ylim=range(pred_rmse), type='l')
for(i in 1:ncol(pred_rmse)) lines(pred_rmse[,i], col=i)

nperlake <- table(laketrout_lw$LakeName)
apply(pred_rmse, 1, weighted.mean, w=nperlake) %>% plot
apply(pred_rmse, 1, weighted.mean, w=sqrt(nperlake)) %>% plot

