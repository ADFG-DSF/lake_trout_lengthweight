## loading all data
source("R/1_laketrout_lwdata.R")

## prep a clean dataset with lengths
laketrout_linf_filter1 <- laketrout_all %>%
  # mutate(Weight_g = ifelse(Weight_g > 100000, Weight_g/1000, Weight_g)) %>%
  mutate(Weight_g = ifelse(Weight_g < 50, Weight_g*1000, Weight_g)) %>%
  filter(is.na(Weight_g) | Weight_g < 100000) %>%  # might not even be needed for linf
  filter(ForkLength_mm > 100) %>%
  filter(!ProjectTitle %in% c("Mark-Recapture Event 1 - (September - 2002)",
                              "Mark-Recapture Event 1 - (September - 2003)",
                              "Mark-Recapture Event 1 - (September - 2004)",
                              "Mark-Recapture Event 2 - (May - 2003)"))

plot(laketrout_linf_filter1$Weight_g ~ laketrout_linf_filter1$ForkLength_mm, log="xy")

laketrout_lm_for_linf_filter <- lm(log(laketrout_linf_filter1$Weight_g) ~
                                     log(laketrout_linf_filter1$ForkLength_mm))
par(mfrow=c(2,2))
plot(laketrout_lm_for_linf_filter)
plot(laketrout_lm_for_linf_filter$residuals)

all_resids <- log(laketrout_linf_filter1$Weight_g) -
  predict(laketrout_lm_for_linf_filter, newdata=laketrout_linf_filter1)
plot(all_resids)
sd_resids <- sd(all_resids, na.rm=TRUE)

laketrout_linf_filter2 <- laketrout_linf_filter1 %>%
  filter(is.na(all_resids) |
           (all_resids > -4*sd_resids & all_resids < 4*sd_resids))

plot(laketrout_linf_filter2$Weight_g ~ laketrout_linf_filter2$ForkLength_mm, log="xy")
# %>% lm %>% plot

par(mfrow=c(2,2))
lm(log(laketrout_linf_filter2$Weight_g) ~ log(laketrout_linf_filter2$ForkLength_mm)) %>% plot

laketrout_linf_filter2$ForkLength_mm %>% plot

plot(laketrout_linf_filter2$ForkLength_mm,
     col=rcolors(200)[as.numeric(as.factor(laketrout_linf_filter2$LakeName))])
plot(laketrout_linf_filter2$ForkLength_mm,
     col=rcolors(200)[as.numeric(as.factor(laketrout_linf_filter2$ProjectTitle))])
abline(h=c(400,450,500,600))
# boxplot(laketrout_linf_filter2$ForkLength_mm ~
#      as.factor(laketrout_linf_filter2$ProjectTitle))

laketrout_linf <- filter(laketrout_linf_filter2, ForkLength_mm < 1400)


nlengths <- table(laketrout_linf$LakeName[!is.na(laketrout_linf$ForkLength_mm) &
                                            !is.na(laketrout_linf$SurfaceArea_h)])
nweights <- table(laketrout_linf$LakeName[!is.na(laketrout_linf$Weight_g) &
                    !is.na(laketrout_linf$SurfaceArea_h)])
plot(ecdf(nlengths), xlim=c(0,100))
abline(v=20)

n_greaterthan <- function(x, y) {
  colSums(outer(x, y, ">"))
}
p_greaterthan <- function(x, y) {
  colMeans(outer(x, y, ">"))
}
n_greaterthan(x=nlengths, y=1:100)  %>% plot(ylim=c(0,60))
n_greaterthan(x=nweights, y=1:100)  %>% points


## define input data for Lester lake area relationship

minn <- 25  ## minimum number of length samples to accept per lake (investigate this!)
laketrout_lengtharea <- filter(laketrout_linf,
                               LakeName %in% names(nlengths[nlengths >= minn]))

# Lester did "the mean fork length of the largest 10% in our sample,
# after removing fish smaller than 300 mm."
length_lester <- function(x) {
  x1 <- x[x >= 300]
  x2 <- x1[x1 >= quantile(x1, 0.9, na.rm=TRUE)]
  return(mean(x2, na.rm=TRUE))
}
linf1 <- with(laketrout_lengtharea,
             tapply(ForkLength_mm, LakeName, length_lester))

# also looking at 95th percentile value (can even investigate other percentages)
linf2 <- with(laketrout_lengtharea,
              tapply(ForkLength_mm, LakeName, quantile, p=0.90))

larea <- with(laketrout_lengtharea,
              tapply(SurfaceArea_h, LakeName, median, na.rm=TRUE))

# comparing Linf approximation methods
lester_curve <- function(x) 957*(1-exp(-0.14*(1+log(x))))
plot(larea, linf1, log="x")
curve(lester_curve(x), lty=3, add=TRUE)  # overlay Lester curve
plot(larea, linf2, log="x")
curve(lester_curve(x), lty=3, add=TRUE)  # overlay Lester curve

## first look at min number of length samples
minns <- seq(15, 60, by=5)  # can even go finer
rmse1 <- rmse80 <- rmse90 <- rmse95 <- rmse99 <- NA*minns
for(i_minn in seq_along(minns)) {
  laketrout_lengtharea_test <- filter(laketrout_linf,
                                 LakeName %in% names(nlengths[nlengths >= minns[i_minn]]))
  larea_test <- with(laketrout_lengtharea_test,
                tapply(SurfaceArea_h, LakeName, median, na.rm=TRUE))
  linf1_test <- with(laketrout_lengtharea_test,
                tapply(ForkLength_mm, LakeName, length_lester))
  linf80_test <- with(laketrout_lengtharea_test,
                      tapply(ForkLength_mm, LakeName, quantile, p=0.8))
  linf90_test <- with(laketrout_lengtharea_test,
                      tapply(ForkLength_mm, LakeName, quantile, p=0.9))
  linf95_test <- with(laketrout_lengtharea_test,
                      tapply(ForkLength_mm, LakeName, quantile, p=0.95))
  linf99_test <- with(laketrout_lengtharea_test,
                      tapply(ForkLength_mm, LakeName, quantile, p=0.99))
  rmse1[i_minn] <- sqrt(mean((linf1_test-lester_curve(larea_test))^2, na.rm=T))
  rmse80[i_minn] <- sqrt(mean((linf80_test-lester_curve(larea_test))^2, na.rm=T))
  rmse90[i_minn] <- sqrt(mean((linf90_test-lester_curve(larea_test))^2, na.rm=T))
  rmse95[i_minn] <- sqrt(mean((linf95_test-lester_curve(larea_test))^2, na.rm=T))
  rmse99[i_minn] <- sqrt(mean((linf99_test-lester_curve(larea_test))^2, na.rm=T))
}
plot(minns, rmse1, ylim=range(rmse1, rmse99, rmse95, rmse90, rmse80))
lines(minns, rmse80, col=2)
lines(minns, rmse90, col=3)
lines(minns, rmse95, col=4)
lines(minns, rmse99, col=5)



## actually it might make more sense to fit our own curve, and see if there is
## a data source that is most consistent (enough data for fit, lower amount of low-data noise)


## trying to re-fit Lester's L_inf from lake area model
# specify model, which is written to a temporary file
linf_jags <- tempfile()
cat('model {
  for(i in 1:n) {
    y[i] ~ dnorm(mu[i], tau)
    ypp[i] ~ dnorm(mu[i], tau)
    mu[i] <- b0 * (1 - exp(-b1 * (1 + log(x[i]))))
  }

  tau <- pow(sig, -2)
  sig ~ dunif(0, 300)
  b0 ~ dnorm(957, pow(0.3*957, -2))
  b0_prior ~ dnorm(957, pow(0.3*957, -2))
  b1 ~ dnorm(0.14, pow(1*0.14, -2))
  b1_prior ~ dnorm(0.14, pow(1*0.14, -2))
}', file=linf_jags)


# bundle data to pass into JAGS
# xx <- area_ln_q[!is.na(area_ln_q) & !is.na(ln_q)]
# yy <- ln_q[!is.na(area_ln_q) & !is.na(ln_q)]
linf_data <- list(x=larea[!is.na(larea) & !is.na(linf1)],
                  y=linf1[!is.na(larea) & !is.na(linf1)],
                  n=sum(!is.na(larea) & !is.na(linf1)))

# JAGS controls
niter <- 10000
# ncores <- 3
ncores <- min(10, parallel::detectCores()-1)

{
  tstart <- Sys.time()
  print(tstart)
  linf_jags_out <- jagsUI::jags(model.file=linf_jags, data=linf_data,
                                parameters.to.save=c("b0","b1",
                                                     "b0_prior","b1_prior",
                                                     "sig", "mu","ypp"),
                                n.chains=ncores, parallel=T, n.iter=niter,
                                n.burnin=niter/2, n.thin=niter/2000)
  print(Sys.time() - tstart)
}
# nbyname(linf_jags_out)
plotRhats(linf_jags_out)
traceworstRhat(linf_jags_out, parmfrow=c(3,3))

par(mfrow=c(3,2))
comparepriors(linf_jags_out)
plot(linf_data$x, linf_data$y, log="x",
     xlab="Lake Area (ha)", ylab="L_inf (mm)",
     main="Trend")
envelope(linf_jags_out, p="mu", x=linf_data$x, add=TRUE)
curve(957*(1-exp(-0.14*(1+log(x)))), lty=3, add=TRUE)

envelope(linf_jags_out, p="ypp", x=linf_data$x, log="x",
         xlab="Lake Area (ha)", ylab="L_inf (mm)",
         main="Post Predictive")
points(linf_data$x, linf_data$y)
curve(957*(1-exp(-0.14*(1+log(x)))), lty=3, add=TRUE)

qq_postpred(linf_jags_out, p="ypp", y=linf_data$y)
ts_postpred(linf_jags_out, p="ypp", y=linf_data$y, x=linf_data$x, log="x")
cor(linf_jags_out$q50$mu, linf_data$y)^2


## CLEAN THIS ALL UP!!!

# try different minns again
# try different quantiles again
# this time Bayesian fit, and record R^2 and median sig
# heck, record rmse as well

rmse <- function(x, y) sqrt(mean((x - y)^2, na.rm=TRUE))

minns <- seq(15, 60, by=5)  # can even go finer
rmse_test <- sig_test <- matrix(nrow=length(minns), ncol=5)
for(i_minn in seq_along(minns)) {
  laketrout_lengtharea_test <- filter(laketrout_linf,
                                      LakeName %in% names(nlengths[nlengths >= minns[i_minn]]))
  larea_test <- with(laketrout_lengtharea_test,
                     tapply(SurfaceArea_h, LakeName, median, na.rm=TRUE))

  linfs <- list(with(laketrout_lengtharea_test,
                     tapply(ForkLength_mm, LakeName, length_lester)),
                with(laketrout_lengtharea_test,
                                    tapply(ForkLength_mm, LakeName, quantile, p=0.8)),
                with(laketrout_lengtharea_test,
                                    tapply(ForkLength_mm, LakeName, quantile, p=0.9)),
                with(laketrout_lengtharea_test,
                                    tapply(ForkLength_mm, LakeName, quantile, p=0.95)),
                with(laketrout_lengtharea_test,
                                    tapply(ForkLength_mm, LakeName, quantile, p=0.99)))
  for(i_data in 1:5) {
    xx <- larea_test
    yy <- linfs[[i_data]]
    linf_data_test <- list(x = xx[!is.na(xx) & !is.na(yy)],
                           y = yy[!is.na(xx) & !is.na(yy)],
                           n = sum(!is.na(xx) & !is.na(yy)))
    tstart <- Sys.time()
    print(tstart)
    linf_jags_out_test <- jagsUI::jags(model.file=linf_jags, data=linf_data_test,
                                  parameters.to.save=c("b0","b1",
                                                       "b0_prior","b1_prior",
                                                       "sig", "mu","ypp"),
                                  n.chains=ncores, parallel=T, n.iter=niter,
                                  n.burnin=niter/2, n.thin=niter/2000)
    print(Sys.time() - tstart)

    rmse_test[i_minn, i_data] <- rmse(x=linf_data_test$y,
                                      y=linf_jags_out_test$q50$mu)
    sig_test[i_minn, i_data] <- linf_jags_out_test$q50$sig
  }
}
which.min(apply(rmse_test, 1, min))
which.min(apply(rmse_test, 2, min))
which.min(apply(sig_test, 1, min))
which.min(apply(sig_test, 2, min))

### NOT SURE IF THIS IS MEANINGFUL
# yes, we have the best fit, but does that say anything about L_inf ITSELF??

i_minn <- 8   ############ hard coded this
laketrout_lengtharea_test <- filter(laketrout_linf,
                                    LakeName %in% names(nlengths[nlengths >= minns[i_minn]]))
larea_test <- with(laketrout_lengtharea_test,
                   tapply(SurfaceArea_h, LakeName, median, na.rm=TRUE))

linfs <- list(with(laketrout_lengtharea_test,
                   tapply(ForkLength_mm, LakeName, length_lester)),
              with(laketrout_lengtharea_test,
                   tapply(ForkLength_mm, LakeName, quantile, p=0.8)),
              with(laketrout_lengtharea_test,
                   tapply(ForkLength_mm, LakeName, quantile, p=0.9)),
              with(laketrout_lengtharea_test,
                   tapply(ForkLength_mm, LakeName, quantile, p=0.95)),
              with(laketrout_lengtharea_test,
                   tapply(ForkLength_mm, LakeName, quantile, p=0.99)))
for(i_data in 2) {   ############ hard coded this
  xx <- larea_test
  yy <- linfs[[i_data]]
  linf_data_test <- list(x = xx[!is.na(xx) & !is.na(yy)],
                         y = yy[!is.na(xx) & !is.na(yy)],
                         n = sum(!is.na(xx) & !is.na(yy)))
  tstart <- Sys.time()
  print(tstart)
  linf_jags_out_test <- jagsUI::jags(model.file=linf_jags, data=linf_data_test,
                                     parameters.to.save=c("b0","b1",
                                                          "b0_prior","b1_prior",
                                                          "sig", "mu","ypp"),
                                     n.chains=ncores, parallel=T, n.iter=niter,
                                     n.burnin=niter/2, n.thin=niter/2000)
  print(Sys.time() - tstart)

  # rmse_test[i_minn, i_data] <- rmse(x=linf_data_test$y,
  #                                   y=linf_jags_out_test$q50$mu)
  # sig_test[i_minn, i_data] <- linf_jags_out_test$q50$sig
}
plotRhats(linf_jags_out_test)
traceworstRhat(linf_jags_out_test, parmfrow=c(3,3))

par(mfrow=c(3,2))
comparepriors(linf_jags_out_test)
plot(linf_data_test$x, linf_data_test$y, log="x",
     xlab="Lake Area (ha)", ylab="L_inf (mm)",
     main="Trend")
envelope(linf_jags_out_test, p="mu", x=linf_data_test$x, add=TRUE)
curve(957*(1-exp(-0.14*(1+log(x)))), lty=3, add=TRUE)

envelope(linf_jags_out_test, p="ypp", x=linf_data_test$x, log="x",
         xlab="Lake Area (ha)", ylab="L_inf (mm)",
         main="Post Predictive")
points(linf_data_test$x, linf_data_test$y)
curve(957*(1-exp(-0.14*(1+log(x)))), lty=3, add=TRUE)

qq_postpred(linf_jags_out_test, p="ypp", y=linf_data_test$y)
ts_postpred(linf_jags_out_test, p="ypp", y=linf_data_test$y, x=linf_data_test$x, log="x")
cor(linf_jags_out_test$q50$mu, linf_data_test$y)^2
