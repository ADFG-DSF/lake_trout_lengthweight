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
niter <- 20*1000#20*1000
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

cols <- 3 - 1*(sort(unique(laketrout$LakeNum)) %in% WL_data$whichlakesc)
caterpillar(WL_jags_out, p="b1", col=cols)
caterpillar(WL_jags_out, p="b0", col=cols)
caterpillar(WL_jags_out, p="b0_interp", col=cols)

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

    # L_quantile[j] <- mean(x[whichdata] <= Linf[j])   #########
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
                                                   "yfit",
                                                   "L_quantile"),#
                              n.chains=ncores, parallel=T, n.iter=niter,
                              n.burnin=niter/2, n.thin=niter/2000)
  print(Sys.time() - tstart)
}

par(mfrow=c(1,1))
# nbyname(LA_jags_out)
plotRhats(LA_jags_out)
traceworstRhat(LA_jags_out, parmfrow=c(3,3))

par(mfrow=c(1,2))
caterpillar(LA_jags_out, p="Linf")
caterpillar(LA_jags_out, p="Linf_prior")
caterpillar(LA_jags_out, p="k")
caterpillar(LA_jags_out, p="k_prior")
caterpillar(LA_jags_out, p="t0")
caterpillar(LA_jags_out, p="t0_prior")

## calculate a vector of Linf quantile
L_quantile <- matrix(nrow=nrow(LA_jags_out$sims.list$Linf),
                     ncol=ncol(LA_jags_out$sims.list$Linf))
# for(i in 1:nrow(L_quantile)) {
#   for(j in 1:ncol(L_quantile)) {   # might be able to make this more efficient
#     L_quantile[i,j] <- mean(LA_data$y[LA_data$lake==j] < LA_jags_out$sims.list$Linf[i,j], na.rm=TRUE)
#   }
# }
for(j in 1:ncol(L_quantile)) {
  L_quantile[,j] <- colMeans(outer(LA_data$y[LA_data$lake==j], LA_jags_out$sims.list$Linf[,j], "<"), na.rm=TRUE)
  # str(LA_jags_out$sims.list$Linf[,j])
}
caterpillar(L_quantile, col=3-(laketrout_Winf$n_Age > 100))

## imputing for data-poor lakes
L_quantile_imputed <- matrix(nrow=nrow(L_quantile), ncol=nrow(laketrout_Winf))
for(j in which(laketrout_Winf$n_Age < 100)) {
  L_quantile_imputed[,j] <- sample(L_quantile[, which(laketrout_Winf$n_Age > 100)], nrow(L_quantile))
}
for(j in which(laketrout_Winf$n_Age > 100)) {
  L_quantile_imputed[,j] <- L_quantile[,j]
}
caterpillar(L_quantile_imputed)


par(mfrow=c(3,3))
for(ilake in LA_data$whichlakes) {
  envelope(LA_jags_out$sims.list$yfit[,,ilake],
           ylim=c(0, max(LA_data$y[!is.na(LA_data$x)], na.rm=TRUE)))
  points(x=LA_data$x[LA_data$lake==ilake], y=LA_data$y[LA_data$lake==ilake])
}


## run Linf~Area model, using Lester priors

## trying to re-fit Lester's L_inf from lake area model
# specify model, which is written to a temporary file
LinfArea_jags <- tempfile()
cat('model {
  for(i in whichlakes) {
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
}', file=LinfArea_jags)

# bundle data to pass into JAGS
larea <- tapply(laketrout$SurfaceArea_h, laketrout$LakeNum, median, na.rm=TRUE)
quantile_meds <- apply(L_quantile_imputed, 2, median, na.rm=TRUE)
Linf_med <- rep(NA, nrow(laketrout_Winf))
for(ilake in 1:length(Linf_med)) {
  Linf_med[ilake] <- quantile(laketrout$ForkLength_mm[laketrout$LakeNum==ilake],
                              quantile_meds[ilake], na.rm=TRUE)
}
LinfArea_data <- list(x = larea,
                  y = Linf_med,
                  whichlakes = which(!is.na(larea)))

# JAGS controls
niter <- 100000
# ncores <- 3
ncores <- min(10, parallel::detectCores()-1)

{
  tstart <- Sys.time()
  print(tstart)
  LinfArea_jags_out <- jagsUI::jags(model.file=LinfArea_jags, data=LinfArea_data,
                                parameters.to.save=c("b0","b1",
                                                     "b0_prior","b1_prior",
                                                     "sig", "mu","ypp"),
                                n.chains=ncores, parallel=T, n.iter=niter,
                                n.burnin=niter/2, n.thin=niter/2000)
  print(Sys.time() - tstart)
}
# nbyname(LinfArea_jags_out)
plotRhats(LinfArea_jags_out)
traceworstRhat(LinfArea_jags_out, parmfrow=c(3,3))

par(mfrow=c(3,2))
comparepriors(LinfArea_jags_out)
plot(LinfArea_data$x, LinfArea_data$y, log="x",
     xlab="Lake Area (ha)", ylab="L_inf (mm)",
     main="Trend")
envelope(LinfArea_jags_out, p="mu", x=LinfArea_data$x, add=TRUE)
curve(957*(1-exp(-0.14*(1+log(x)))), lty=3, add=TRUE)

plot(LinfArea_data$x, LinfArea_data$y, log="x",
     xlab="Lake Area (ha)", ylab="L_inf (mm)",
     main="Trend")
envelope(LinfArea_jags_out, p="ypp", x=LinfArea_data$x, log="x",
         xlab="Lake Area (ha)", ylab="L_inf (mm)",
         main="Post Predictive", add=TRUE)
points(LinfArea_data$x, LinfArea_data$y)
curve(957*(1-exp(-0.14*(1+log(x)))), lty=3, add=TRUE)

qq_postpred(LinfArea_jags_out, p="ypp", y=LinfArea_data$y)
# ts_postpred(LinfArea_jags_out, p="ypp", y=LinfArea_data$y, x=LinfArea_data$x, log="x")




### all Winf methods:
nlake <- nrow(laketrout_Winf)
Winf_all <- array(dim=c(10000, nlake, 4))  # n mcmc, n lakes, n methods

# Winf from quantile
# need sufficient weights sample + L_quantile_imputed
for(ilake in 1:nlake) {
  print(ilake)
  for(imcmc in 1:10000) {
    Winf_all[imcmc, ilake, 1] <- quantile(laketrout$Weight_g[laketrout$LakeNum==ilake]/1000,
                                          L_quantile_imputed[imcmc, ilake],
                                          na.rm=TRUE)
  }
}
caterpillar(Winf_all[,,1])
# for(ilake in 1:nlake) {
#   Winf_all[, ilake, 1] <- outer(laketrout$Weight_g[laketrout$LakeNum==ilake]/1000,
#             L_quantile_imputed[, ilake],
#             quantile, na.rm=TRUE)
# }


# Winf from Linf
# need Linf + WL relationship
Linf_mat <- matrix(nrow=10000, ncol=nlake)
Linf_mat[, 1:ncol(LA_jags_out$sims.list$Linf)] <- LA_jags_out$sims.list$Linf

# nbyname(WL_jags_out)
b0 <- WL_jags_out$sims.list$b0_interp
b1 <- WL_jags_out$sims.list$b1


Winf_all[,,2] <- exp(b0 + b1*log(Linf_mat))
caterpillar(Winf_all[,,2])



# Winf from (Linf from quantile)
# need sufficient lengths sample + L_quantile_imputed + WL relationship
Linf_placeholder <- NA*Winf_all[,,3]
for(ilake in 1:nlake) {
  print(ilake)
  for(imcmc in 1:10000) {
    Linf_placeholder[imcmc, ilake] <- quantile(laketrout$ForkLength_mm[laketrout$LakeNum==ilake],
                                          L_quantile_imputed[imcmc, ilake],
                                          na.rm=TRUE)
  }
}
Winf_all[,,3] <- exp(b0 + b1*log(Linf_placeholder))
caterpillar(Winf_all[,,3])



# Winf from (Linf from area)
# need ypp Linf  + WL relationship
ypp_mat <- LinfArea_jags_out$sims.list$ypp
Winf_all[,,4] <- exp(b0 + b1*log(ypp_mat))
caterpillar(Winf_all[,,4])


# par(mfrow=c(1,1))
# comparecat(list(data.frame(Winf_all[,,1]),
#                 data.frame(Winf_all[,,2]),
#                 data.frame(Winf_all[,,3]),
#                 data.frame(Winf_all[,,4])), ylim=c(0,15))
caterpillar(Winf_all[,,1], ylim=c(0,15), col=1)
caterpillar(Winf_all[,,2], add=TRUE, col=2, x=1:nlake+.2)
caterpillar(Winf_all[,,3], add=TRUE, col=3, x=1:nlake+.4)
caterpillar(Winf_all[,,4], add=TRUE, col=4, x=1:nlake+.6)


## should only include those with adequate data to support them
Winf_all_all <- Winf_all
# Winf_all <- Winf_all_all

# Winf from quantile
# need sufficient weights sample + L_quantile_imputed
Winf_all[, laketrout_Winf$n_Weight < 150, 1] <- NA
Winf_all[, laketrout_Winf$n_Weight >= 150, 2:4] <- NA

# Winf from Linf
# need Linf + WL relationship
Winf_all[, laketrout_Winf$n_Age < 50, 2] <- NA
Winf_all[, laketrout_Winf$n_Age >= 50, 3:4] <- NA

# Winf from (Linf from quantile)
# need sufficient lengths sample + L_quantile_imputed + WL relationship
Winf_all[, laketrout_Winf$n_Length < 150, 3] <- NA
Winf_all[, laketrout_Winf$n_Length >= 150, 4] <- NA

# Winf from (Linf from area)
# need ypp Linf  + WL relationship


caterpillar(Winf_all[,,1], ylim=c(0,15), col=1, x=2*(1:nlake))
caterpillar(Winf_all[,,2], add=TRUE, col=2, x=2*(1:nlake)+.3)
caterpillar(Winf_all[,,3], add=TRUE, col=3, x=2*(1:nlake)+.6)
caterpillar(Winf_all[,,4], add=TRUE, col=4, x=2*(1:nlake)+.9)

## .....do the big gnarly predictive thingy from mcmc vectors
