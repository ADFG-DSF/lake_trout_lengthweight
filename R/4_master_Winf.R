library(tidyverse)
library(jagsUI)
library(jagshelper)

## load data, filter bad observations BUT NOT empty lakes
## - bad = obviously bad W~L, obviously bad L~age, outlying W, L, age

# morphometry <- read_csv("flat_data/lake_morphometry.csv", skip=2)
morphometry <- read_csv("flat_data/lake_morphometry2.csv", skip=2)

# # is lake name unique?  YES
# sum(!is.na(morphometry$LakeName))
# length(unique(morphometry$LakeName))

# laketrout_all <- read_csv("flat_data/length_weight.csv", skip=2) %>%
#   left_join(morphometry)
laketrout_all <- read_csv("flat_data/length_weight2.csv", skip=2) %>%
  left_join(morphometry)
nrow(laketrout_all)
sapply(laketrout_all, function(x) sum(is.na(x)))

laketrout1 <- laketrout_all %>%
  filter(!ProjectTitle %in% c("Mark-Recapture Event 1 - (September - 2002)",
                              "Mark-Recapture Event 1 - (September - 2003)",
                              "Mark-Recapture Event 1 - (September - 2004)",
                              "Mark-Recapture Event 2 - (May - 2003)"))
  # filter(is.na(ForkLength_mm) | ForkLength_mm )
lm1 <- with(laketrout1, lm(log(Weight_g) ~ log(ForkLength_mm)))
resids1 <- log(laketrout1$Weight_g) - predict(lm1, newdata=laketrout1)
laketrout2 <- filter(laketrout1, is.na(resids1) | abs(resids1) < 4*sd(resids1, na.rm=TRUE))

laketrout <- laketrout2 %>%
  filter(is.na(Age) | Age < 50) %>%
  mutate(LakeNum = as.numeric(as.factor(LakeName)))

nrow(laketrout)

### plotting weight ~ length
laketrout_all %>%
  # ggplot(aes(x=ForkLength_mm, y=Weight_g, color=LakeName)) +
  ggplot(aes(x=ForkLength_mm, y=Weight_g)) +
  geom_point() +
  scale_y_log10() +
  scale_x_log10() +
  theme_bw()

laketrout %>%
  # ggplot(aes(x=ForkLength_mm, y=Weight_g, color=LakeName)) +
  ggplot(aes(x=ForkLength_mm, y=Weight_g)) +
  geom_point() +
  scale_y_log10() +
  scale_x_log10() +
  theme_bw()


## plotting length ~ age
laketrout_all %>%
  # ggplot(aes(y=ForkLength_mm, x=Age, color=LakeName)) +
  ggplot(aes(y=ForkLength_mm, x=Age)) +
  geom_point() +
  theme_bw()

laketrout %>%
  # ggplot(aes(y=ForkLength_mm, x=Age, color=LakeName)) +
  ggplot(aes(y=ForkLength_mm, x=Age)) +
  geom_point() +
  theme_bw()

sapply(laketrout, function(x) sum(is.na(x)))
table(is.na(laketrout$Latitude_WGS84),
      is.na(laketrout$SurfaceArea_h),
      is.na(laketrout$Weight_g))




## run W~L model, using post predictive when no data exists
WL_jags <- tempfile()
cat('model {

  # loop over all WL data
  for(i in whichdata) {
    y[i] ~ dnorm(mu[i], tau)
    # ypp[i] ~ dnorm(mu[i], tau)    # might take this out if model seems ok
    mu[i] <- b0[lake[i]] + b1[lake[i]]*x[i]
  }

  # loop over lakes     -- need to somehow predict lakes without lat & area???
  for(j in whichlakes) {
    b0[j] ~ dnorm(mu_b0[j], tau_b0)
    b0_interp[j] <- b0[j] - b1[j]*meanx
    mu_b0[j] <- b0_int
                + b0_area*area[j]

    b1[j] ~ dnorm(mu_b1[j], tau_b1)
    mu_b1[j] <- b1_int
                + b1_lat*lat[j]
  }

  for(jc in whichlakesc) {
    b0[jc] ~ dnorm(mu_b0[jc], tau_b0)
    b0_interp[jc] <- b0[jc] - b1[jc]*meanx
    mu_b0[jc] <- b0_int
                + b0_area*areasim[jc]

    b1[jc] ~ dnorm(mu_b1[jc], tau_b1)
    mu_b1[jc] <- b1_int
                + b1_lat*latsim[jc]
    areasim[jc] ~ dnorm(mean(area[whichlakes]), pow(sd(area[whichlakes]), -2))
    latsim[jc] ~ dnorm(mean(lat[whichlakes]), pow(sd(lat[whichlakes]), -2))
  }

  sig_b0 ~ dunif(0, 10)
  tau_b0 <- pow(sig_b0, -2)

  sig_b1 ~ dunif(0, 10)
  tau_b1 <- pow(sig_b1, -2)

  b0_int ~ dnorm(0, 0.001)
  b0_area ~ dnorm(0, 0.001)

  b1_int ~ dnorm(0, 0.001)
  b1_lat ~ dnorm(0, 0.001)

  tau <- pow(sig, -2)
  sig ~ dunif(0, 10)

}', file=WL_jags)

WL_data <- list(y = log(laketrout$Weight_g/1000),
                x = log(laketrout$ForkLength_mm) - mean(log(laketrout$ForkLength_mm), na.rm=TRUE),
                meanx = mean(log(laketrout$ForkLength_mm), na.rm=TRUE),
                whichdata = which(!is.na(laketrout$ForkLength_mm) &
                                    !is.na(laketrout$SurfaceArea_h) &
                                    !is.na(laketrout$Latitude_WGS84)),
                whichlakes = sort(unique(laketrout$LakeNum[!is.na(laketrout$Latitude_WGS84) & !is.na(laketrout$SurfaceArea_h)])),
                whichlakesc = sort(unique(laketrout$LakeNum[is.na(laketrout$Latitude_WGS84) | is.na(laketrout$SurfaceArea_h)])),
                lat = tapply(laketrout$Latitude_WGS84, laketrout$LakeNum, median) -
                  mean(tapply(laketrout$Latitude_WGS84, laketrout$LakeNum, median), na.rm=TRUE),
                area = log(tapply(laketrout$SurfaceArea_h, laketrout$LakeNum, median)) -
                  mean(log(tapply(laketrout$SurfaceArea_h, laketrout$LakeNum, median)), na.rm=TRUE),
                lake = laketrout$LakeNum)


# JAGS controls
niter <- 10*1000#20*1000
ncores <- min(10, parallel::detectCores()-1)

{
  tstart <- Sys.time()
  print(tstart)
  WL_jags_out <- jagsUI::jags(model.file=WL_jags, data=WL_data,
                              parameters.to.save=c("b0","b1","sig",
                                                   "mu_b0","sig_b0","mu_b1","sig_b1",
                                                   "fit","pred","b0_interp",
                                                   "b0_int","b1_int",
                                                   "b0_lat","b0_area",
                                                   "b1_lat","b1_area"),#,"ypp"
                              n.chains=ncores, parallel=T, n.iter=niter,
                              n.burnin=niter/2, n.thin=niter/2000)
  print(Sys.time() - tstart)
}

par(mfrow=c(1,1))
nbyname(WL_jags_out)
plotRhats(WL_jags_out)
traceworstRhat(WL_jags_out, parmfrow=c(3,3))

# # qq_postpred(WL_jags_out, p="ypp", y=WL_data$y)
# # ts_postpred(WL_jags_out, p="ypp", y=WL_data$y[1:32509])
# # ts_postpred(WL_jags_out, p="ypp", y=WL_data$y, x=WL_data$x)
# qq_postpred(WL_jags_out$sims.list$ypp[,1:20000], y=WL_data$y[1:20000])
# # qq_postpred(WL_jags_out$sims.list$ypp[,(1:20000)[!is.na(WL_data$y[1:20000])]],
# #             y=WL_data$y[1:20000][!is.na(WL_data$y[1:20000])])

caterpillar(WL_jags_out, p="b1")
caterpillar(WL_jags_out, p="b0")
caterpillar(WL_jags_out, p="b0_interp")

# plot regression bands for each lake, overlay data and lester line
xpredict <- seq(from = min(log(laketrout$ForkLength_mm[!is.na(laketrout$Weight_g)]), na.rm=TRUE),
                  to = max(log(laketrout$ForkLength_mm[!is.na(laketrout$Weight_g)]), na.rm=TRUE),
                  length.out=50)
par(mfrow=c(3,3))
for(ilake in 1:max(laketrout$LakeNum)) {
  logweight_predict <- WL_jags_out$sims.list$b0_interp[,ilake] +
    outer(WL_jags_out$sims.list$b1[,ilake], xpredict)
  plot(NA, xlim=range(xpredict), ylim=range(log(laketrout$Weight_g/1000), na.rm=TRUE))
  if(all(!is.na(logweight_predict))) envelope(logweight_predict, x=xpredict, add=TRUE)
  points(x = log(laketrout$ForkLength_mm[laketrout$LakeNum==ilake]),
         y = log(laketrout$Weight_g[laketrout$LakeNum==ilake]/1000))
  abline(a=-19.56, b=3.2, lty=3)
}

par(mfrow=c(1,1))
caterpillar(WL_jags_out$sims.list$b0_interp + WL_jags_out$sims.list$b1*6.5, transform="exp")



### start building master dataframe to report
getn <- function(x) sum(!is.na(x))
laketrout_Winf <- data.frame(n_Weight = tapply(laketrout$Weight_g, laketrout$LakeName, getn),
                             n_Length = tapply(laketrout$ForkLength_mm, laketrout$LakeName, getn),
                             n_Age = tapply(laketrout$Age, laketrout$LakeName, getn))
laketrout_Winf$Wt_quantile <- laketrout_Winf$n_Weight > 200
laketrout_Winf$Ln_Age <- laketrout_Winf$n_Age > 150
laketrout_Winf$Ln_quantile <- laketrout_Winf$n_Length > 200

laketrout_Winf




## run L~Age model for lakes with both (also include the big bonus lake??)

LA_jags <- tempfile()
cat('model {
  for(i in whichdata) {
    y[i] ~ dnorm(mu[i], tau)
    ypp[i] ~ dnorm(mu[i], tau)
    mu[i] <- Linf[lake[i]]*(1-exp(-k[lake[i]]*(x[i]-t0[lake[i]])))
  }

  for(j in whichlakes) {
    # Linf[j] ~ dnorm(700, 0.0001)
    # Linf_prior[j] ~ dnorm(700, 0.0001)
    # t0[j] ~ dnorm(0, 0.1)
    # t0_prior[j] ~ dnorm(0, 0.1)
    # k[j] ~ dexp(0.1)
    # k_prior[j] ~ dexp(0.1)
    Linf[j] ~ dnorm(mu_Linf, tau_Linf)
    Linf_prior[j] ~ dnorm(mu_Linf, tau_Linf)
    t0[j] ~ dnorm(mu_t0, tau_t0)
    t0_prior[j] ~ dnorm(mu_t0, tau_t0)
    k[j] ~ dlnorm(mu_k, tau_k)
    k_prior[j] ~ dlnorm(mu_k, tau_k)

    for(ifit in 1:nfit) {
      yfit[ifit, j] <- Linf[j]*(1-exp(-k[j]*(xfit[ifit]-t0[j])))
    }
  }

  mu_Linf ~ dnorm(700, 0.0001)
  sig_Linf ~ dunif(0, 300)
  tau_Linf <- pow(sig_Linf, -2)

  mu_t0 ~ dnorm(0, 0.1)
  sig_t0 ~ dunif(0, 1)
  tau_t0 <- pow(sig_t0, -2)

  mu_k ~ dnorm(0, 0.1)
  sig_k ~ dunif(0, 2)
  tau_k <- pow(sig_k, -2)

  tau <- pow(sig, -2)
  sig ~ dunif(0, 200)
  sig_prior ~ dunif(0, 200)

}', file=LA_jags)


# bundle data to pass into JAGS
LA_data <- list(x = laketrout$Age,
                y = laketrout$ForkLength_mm,
                lake = laketrout$LakeNum,
                whichlakes = which(laketrout_Winf$n_Age > 0),  # was > 150
                xfit = 1:max(laketrout$Age, na.rm=TRUE),
                nfit = max(laketrout$Age, na.rm=TRUE))
LA_data$whichdata <- which((laketrout$LakeNum %in% LA_data$whichlakes) &
                            (!is.na(laketrout$Age)))

# JAGS controls
niter <- 20*1000#20*1000
ncores <- min(10, parallel::detectCores()-1)

{
  tstart <- Sys.time()
  print(tstart)
  LA_jags_out <- jagsUI::jags(model.file=LA_jags, data=LA_data,
                              parameters.to.save=c("Linf","t0","k",
                                                   "Linf_prior","t0_prior","k_prior",
                                                   "sig","sig_prior","ypp",
                                                   "mu_Linf","sig_Linf",
                                                   "mu_k", "sig_k",
                                                   "mu_t0","sig_t0",
                                                   "yfit"),#
                              n.chains=ncores, parallel=T, n.iter=niter,
                              n.burnin=niter/2, n.thin=niter/2000)
  print(Sys.time() - tstart)
}

par(mfrow=c(1,1))
# nbyname(LA_jags_out)
plotRhats(LA_jags_out)
traceworstRhat(LA_jags_out, parmfrow=c(4,3))

par(mfrow=c(1,2))
caterpillar(LA_jags_out, p="Linf")
caterpillar(LA_jags_out, p="Linf_prior")
caterpillar(LA_jags_out, p="k")
caterpillar(LA_jags_out, p="k_prior")
caterpillar(LA_jags_out, p="t0")
caterpillar(LA_jags_out, p="t0_prior")

par(mfrow=c(3,3))
for(ilake in LA_data$whichlakes) {
  envelope(LA_jags_out$sims.list$yfit[,,ilake],
           ylim=c(0, max(LA_data$y[!is.na(LA_data$x)], na.rm=TRUE)))
  points(x=LA_data$x[LA_data$lake==ilake], y=LA_data$y[LA_data$lake==ilake])
}


## run Linf~Area model, using Lester priors



## .....do the big gnarly predictive thingy from mcmc vectors
