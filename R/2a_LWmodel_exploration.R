library(jagsUI)
library(jagshelper)
library(tidyverse)


load(file="Rdata/laketrout_sampling_formodel.Rdata")


##### starting a summarized data.frame
getn <- function(x) sum(!is.na(x))
lakenums <- factor(laketrout$LakeNum, levels=1:nrow(morphometry))
laketrout$LakeNum <- factor(laketrout$LakeNum, levels=1:nrow(morphometry))
laketrout_Winf <- data.frame(LakeName = morphometry$LakeName,
                             LakeNum = morphometry$LakeNum,
                             n_Weight = tapply(laketrout$Weight_g, lakenums, getn),
                             n_Length = tapply(laketrout$ForkLength_mm, lakenums, getn),
                             n_Age = tapply(laketrout$Age, lakenums, getn),
                             Latitude = morphometry$Latitude_WGS84,
                             Area_ha = morphometry$SurfaceArea_h,
                             Elevation = morphometry$`Elevation (m)`,
                             Temperature = morphometry$`Temp (C)`,
                             use_fish = morphometry$use_fish,
                             make_estimates = morphometry$make_estimates) %>%
  mutate(n_Weight = ifelse(is.na(n_Weight), 0, n_Weight)) %>%
  mutate(n_Length = ifelse(is.na(n_Length), 0, n_Length)) %>%
  mutate(n_Age = ifelse(is.na(n_Age), 0, n_Age))
lakenames <- morphometry$LakeName




int_Winf_jags <- tempfile()
cat('model {

  ## ---- Length ~ Age Component ---- ##
  for(i in whichdata_Lt) {   # fish-level data from lakes with paired lengths & ages
    # L[i] ~ dnorm(mu_Lt[i], tau_Lt)
    # Lpp[i] ~ dnorm(mu_Lt[i], tau_Lt)
    # mu_Lt[i] <- Linf[lake[i]]*(1-exp(-k[lake[i]]*(Age[i]-t0[lake[i]])))
    L[i] ~ dlnorm(logmu_Lt[i], tau_Lt)
    Lpp[i] ~ dlnorm(logmu_Lt[i], tau_Lt)
    logmu_Lt[i] <- log(Linf[lake[i]]*(1-exp(-k[lake[i]]*(Age[i]-t0[lake[i]]))))
  }

  tau_Lt <- pow(sig_Lt, -2)
  # sig_Lt ~ dunif(0, 100)
  # sig_Lt_prior ~ dunif(0, 100)
  sig_Lt ~ dunif(0, 3)
  sig_Lt_prior ~ dunif(0, 3)


  for(j in whichlakes_Lt) {   # lake level
    t0[j] ~ dnorm(mu_t0, tau_t0)T(,1)
    # t0_prior[j] ~ dnorm(mu_t0, tau_t0)
    k[j] ~ dlnorm(mu_k, tau_k)
    # k_prior[j] ~ dlnorm(mu_k, tau_k)

    # making plottable LVB curve envelopes for each lake
    for(ifit in 1:nfit) {
      Lfit[ifit, j] <- Linf[j]*(1-exp(-k[j]*(Agefit[ifit]-t0[j])))
    }
  }

  mu_t0 ~ dnorm(0, 1)   # was 1
  mu_t0_prior ~ dnorm(0, 1)
  # sig_t0 ~ dunif(0, 0.2)
  sig_t0 ~ dunif(0, 10)
  sig_t0_prior ~ dunif(0, 10)
  # sig_t0 ~ dexp(10)
  tau_t0 <- pow(sig_t0, -2)

  mu_k ~ dnorm(0, 0.1)
  mu_k_prior ~ dnorm(0, 0.1)
  sig_k ~ dunif(0, 3)
  sig_k_prior ~ dunif(0, 3)
  tau_k <- pow(sig_k, -2)



  ## ---- Linf ~ Area component ---- ##
  for(j in whichlakes_LA) {  # lake-level data where there are areas
    # Linf[j] ~ dnorm(mu_LA[j], tau_LA)T(1,)
    # mu_LA[j] <- gam * (1 - exp(-lam * (1 + log(Area[j]))))
    Linf[j] ~ dlnorm(logmu_LA[j], tau_LA)
    logmu_LA[j] <- log(gam * (1 - exp(-lam * (1 + log(Area[j])))))
  }
  for(j in whichlakes_LA_c) {  # lake-level priors where there are no areas
    Linf[j] ~ dnorm(600, 0.00001)T(1,)
    # Linf[j] ~ dnorm(mu_Linf_noarea, tau_Linf_noarea)
  }
  # mu_Linf_noarea ~ dnorm(600, 0.0001)
  # tau_Linf_noarea <- pow(sig_Linf_noarea, -2)
  # sig_Linf_noarea ~ dunif(0, 500)
  # sig_Linf_noarea_prior ~ dunif(0, 500)

  tau_LA <- pow(sig_LA, -2)
  # sig_LA ~ dunif(0, 300)
  # sig_LA_prior ~ dunif(0, 300)
  sig_LA ~ dunif(0, 3)
  sig_LA_prior ~ dunif(0, 3)
  # sig_LA ~ dexp(0.05)
  # sig_LA_prior ~ dexp(0.05)

  gam ~ dnorm(gam_lester, pow(cv_gam_lester*gam_lester, -2))T(0.01,)
  gam_prior ~ dnorm(gam_lester, pow(cv_gam_lester*gam_lester, -2))T(0.01,)
  lam ~ dnorm(lam_lester, pow(cv_lam_lester*lam_lester, -2))T(0.01,)
  lam_prior ~ dnorm(lam_lester, pow(cv_lam_lester*lam_lester, -2))T(0.01,)



  ## ---- Length quantile ~ Linf component ---- ##
  for(j in whichlakes_L) {
    qL[j] ~ dnorm(Linf[j], tau_L[j])
    qLpp[j] ~ dnorm(Linf[j], tau_L[j])
    tau_L[j] <- pow(sig_L[j], -2)
    sig_L[j] <- pow((eta_L^2) + ((zeta_L^2)/nL[j]), 0.5)
  }

  eta_L ~ dunif(0, eta_L_cap)
  eta_L_prior ~ dunif(0, eta_L_cap)
  zeta_L ~ dunif(0, zeta_L_cap)
  zeta_L_prior ~ dunif(0, zeta_L_cap)

  # eta_L ~ dexp(0.1)
  # eta_L_prior ~ dexp(0.1)
  # zeta_L ~ dexp(0.1)
  # zeta_L_prior ~ dexp(0.1)



  ## ---- log Weight ~ log Length component ---- ##
  # loop over all WL data
  for(i in whichdata_WL) {
    logW[i] ~ dnorm(mu_WL[i], tau_WL)
    # logWpp[i] ~ dnorm(mu_WL[i], tau_WL)    # might take this out if model seems ok
    mu_WL[i] <- b0[lake[i]] + b1[lake[i]]*logLc[i]
  }

  # loop over lakes with lat & area data
  for(j in whichlakes_WL) {
    b0[j] ~ dnorm(mu_b0[j], tau_b0)
    b0_interp[j] <- b0[j] - b1[j]*meanlogLc
    mu_b0[j] <- b0_int
                # + b0_area*logareac[j]

    b1[j] ~ dnorm(mu_b1[j], tau_b1)
    mu_b1[j] <- b1_int
                # + b1_lat*latc[j]
  }

  # # loop over lakes without area & lat data (c for complement)
  # # draw log(area) and lat from Normal distributions (mean & sd of actual data)
  # # corresponding b0, b0_interp, and b1 will effectively be from the post pred
  # for(jc in whichlakes_WL_c) {
  #   b0[jc] ~ dnorm(mu_b0[jc], tau_b0)
  #   b0_interp[jc] <- b0[jc] - b1[jc]*meanlogLc
  #   mu_b0[jc] <- b0_int
  #               # + b0_area*areasim[jc]
  #
  #   b1[jc] ~ dnorm(mu_b1[jc], tau_b1)
  #   mu_b1[jc] <- b1_int
  #               # + b1_lat*latsim[jc]
  #   areasim[jc] ~ dnorm(mean(logareac[whichlakes_WL]), pow(sd(logareac[whichlakes_WL]), -2))
  #   latsim[jc] ~ dnorm(mean(latc[whichlakes_WL]), pow(sd(latc[whichlakes_WL]), -2))
  # }

  # global priors
  sig_b0 ~ dunif(0, 10)
  sig_b0_prior ~ dunif(0, 10)
  tau_b0 <- pow(sig_b0, -2)

  sig_b1 ~ dunif(0, 10)
  sig_b1_prior ~ dunif(0, 10)
  tau_b1 <- pow(sig_b1, -2)

  b0_int ~ dnorm(0, 0.001)
  b0_area ~ dnorm(0, 0.001)
  b0_int_prior ~ dnorm(0, 0.001)
  b0_area_prior ~ dnorm(0, 0.001)

  b1_int ~ dnorm(0, 0.001)
  b1_lat ~ dnorm(0, 0.001)
  b1_int_prior ~ dnorm(0, 0.001)
  b1_lat_prior ~ dnorm(0, 0.001)

  tau_WL <- pow(sig_WL, -2)
  sig_WL ~ dunif(0, 10)
  sig_WL_prior ~ dunif(0, 10)



  ## ---- Weight quantile ~ Winf component ---- ##
  for(j in whichlakes_W) {
    qW[j] ~ dnorm(Winf[j], tau_W[j])
    qWpp[j] ~ dnorm(Winf[j], tau_W[j])
    tau_W[j] <- pow(sig_W[j], -2)
    sig_W[j] <- pow((eta_W^2) + ((zeta_W^2)/nW[j]), 0.5)
  }

  eta_W ~ dunif(0, eta_W_cap)
  eta_W_prior ~ dunif(0, eta_W_cap)
  zeta_W ~ dunif(0, zeta_W_cap)
  zeta_W_prior ~ dunif(0, zeta_W_cap)

  # eta_W ~ dexp(0.1)
  # eta_W_prior ~ dexp(0.1)
  # zeta_W ~ dexp(0.1)
  # zeta_W_prior ~ dexp(0.1)



  ## ---- Winf from Linf! ---- ##
  for(j in whichlakes_WL) {  # alllakes
    Winf[j] <- exp(b0_interp[j] + b1[j]*log(Linf[j]))
  }


}', file=int_Winf_jags)



# bundle data to pass into JAGS
q_input <- 0.95   # assumed quantile value for asymptotic size
int_Winf_data <- list(

  # fish-level data
  Age = laketrout$Age,
  L = laketrout$ForkLength_mm,
  logW = log(laketrout$Weight_g/1000),
  logLc = log(laketrout$ForkLength_mm) - mean(log(laketrout$ForkLength_mm), na.rm=TRUE),
  lake = as.numeric(as.character(laketrout$LakeNum)),
  whichdata_WL = which(laketrout$use_fish &
                         !is.na(laketrout$ForkLength_mm) &
                         !is.na(laketrout$SurfaceArea_h) &
                         !is.na(laketrout$Latitude_WGS84) &
                         !is.na(laketrout$Weight_g)),  # note: only 14 use_fish are missing surface area
  # whichdata_WL = which(laketrout$use_fish &
  #                        !is.na(laketrout$ForkLength_mm) &
  #                        !is.na(laketrout$SurfaceArea_h) &
  #                        !is.na(laketrout$Latitude_WGS84)),  # note: only 14 use_fish are missing surface area
  # # ^^^ this may be vestigial
  # # whichdata_WL = which(laketrout$use_fish),  # note: only 14 use_fish are missing surface area

  # lake-level data
  Area = laketrout_Winf$Area_ha,
  latc = morphometry$Latitude_WGS84 - mean(morphometry$Latitude_WGS84, na.rm=TRUE),
  logareac = log(morphometry$SurfaceArea_h) - mean(log(morphometry$SurfaceArea_h), na.rm=TRUE),

  whichlakes_L = which(morphometry$use_fish & laketrout_Winf$n_Length > 10),
  qL = tapply(laketrout$ForkLength_mm, laketrout$LakeNum, quantile, p=q_input, na.rm=TRUE),
  nL = laketrout_Winf$n_Length,
  whichlakes_W = which(morphometry$use_fish & laketrout_Winf$n_Weight > 10),
  qW = tapply(laketrout$Weight_g/1000, laketrout$LakeNum, quantile, p=q_input, na.rm=TRUE),
  nW = laketrout_Winf$n_Weight,

  whichlakes_Lt = which(morphometry$use_fish & laketrout_Winf$n_Age > 0 & laketrout_Winf$n_Length > 0),
  whichlakes_LA = which(!is.na(laketrout_Winf$Area_ha)),# & laketrout_Winf$n_Length > 0),
  whichlakes_LA_c = which(!(!is.na(laketrout_Winf$Area_ha))),# & laketrout_Winf$n_Length > 0)),
  whichlakes_WL = which(!is.na(laketrout_Winf$Lat) & !is.na(laketrout_Winf$Area_ha)),
  whichlakes_WL_c = which(is.na(laketrout_Winf$Lat) | is.na(laketrout_Winf$Area_ha)),
  alllakes = 1:nrow(laketrout_Winf),

  # inputs for fitted vals of L~Age
  Agefit = 1:max(laketrout$Age, na.rm=TRUE),
  nfit = max(laketrout$Age, na.rm=TRUE),

  # priors from Lester
  gam_lester = 957,
  cv_gam_lester = 0.2*2.5,
  lam_lester = 0.14,
  cv_lam_lester = 1*2.5,

  # global values
  eta_L_cap = 500,
  zeta_L_cap = 2000,
  eta_W_cap = 5,
  zeta_W_cap = 20,
  meanlogLc = mean(log(laketrout$ForkLength_mm), na.rm=TRUE)
)

# can i make this within the above?
int_Winf_data$whichdata_Lt <- which((laketrout$LakeNum %in% int_Winf_data$whichlakes_Lt) &
                                      (!is.na(laketrout$Age)) & (!is.na(laketrout$ForkLength_mm)))


parameters <- c("sig_Lt", "sig_Lt_prior",
                "t0", "t0_prior",
                "k", "k_prior",
                "mu_t0","mu_t0_prior",
                "sig_t0","sig_t0_prior",
                "mu_k","mu_k_prior",
                "sig_k","sig_k_prior",
                "Lfit",
                "Linf", "Winf",
                "mu_LA",
                "sig_LA", "sig_LA_prior",
                "gam", "gam_prior",
                "lam", "lam_prior",
                "eta_L", "eta_L_prior",
                "zeta_L", "zeta_L_prior",
                "sig_L",
                "Lpp","qLpp",
                "mu_Linf_noarea", "sig_Linf_noarea", "sig_Linf_noarea_prior",
                "b0","b1","sig_WL","sig_WL_prior",
                "mu_b0","mu_b1",
                "sig_b0","sig_b0_prior",
                "sig_b1","sig_b1_prior",
                "b0_interp",
                "b0_int","b1_int","b0_int_prior","b1_int_prior",
                "b0_lat","b0_area","b0_lat_prior","b0_area_prior",
                "b1_lat","b1_area","b1_lat_prior","b1_area_prior",
                "eta_W", "eta_W_prior",
                "zeta_W", "zeta_W_prior",
                "sig_W",
                "Wpp","qWpp")




# JAGS controls
# niter <- 2*1000
# niter <- 20*1000
niter <- 50*1000      # 50k in 15 minutes
# niter <- 100*1000
# niter <- 200*1000     # 1.1 hr
# niter <- 500*1000  # 3 hr

# ncores <- 3
ncores <- 8
# ncores <- min(10, parallel::detectCores()-1)


###### RUNNING THE MODEL #####
{
  tstart <- Sys.time()
  print(tstart)
  int_Winf_jags_out <- jagsUI::jags(model.file=int_Winf_jags, data=int_Winf_data,
                                    parameters.to.save=parameters,
                                    n.chains=ncores, parallel=T, n.iter=niter,
                                    n.burnin=niter/2, n.thin=niter/2000)
  print(Sys.time() - tstart)
}
int_Winf_jags_out$DIC


wl_lakes <- subset(laketrout, !is.na(Weight_g) & !is.na(ForkLength_mm))$LakeNum %>%
  unique %>%
  sort %>%
  as.character %>%
  as.numeric

b0 <- int_Winf_jags_out$q50$b0[wl_lakes]
b1 <- int_Winf_jags_out$q50$b1[wl_lakes]
elev <- morphometry$`Elevation (m)`[wl_lakes]
temp <- morphometry$`Temp (C)`[wl_lakes]
larea <- log(morphometry$SurfaceArea_h[wl_lakes])
lat <- morphometry$Latitude_WGS84[wl_lakes]

par(mfrow=c(2,2))
plot(b0 ~ elev)
plot(b0 ~ temp)
plot(b0 ~ larea)
plot(b0 ~ lat)

plot(b1 ~ elev)
plot(b1 ~ temp)
plot(b1 ~ larea)
plot(b1 ~ lat)

lm(b0 ~ elev) %>% summary
lm(b0 ~ temp) %>% summary
lm(b0 ~ larea) %>% summary  # this is the only one!
lm(b0 ~ lat) %>% summary

lm(b1 ~ elev) %>% summary
lm(b1 ~ temp) %>% summary
lm(b1 ~ larea) %>% summary
lm(b1 ~ lat) %>% summary  # this one just a little

vars <- c("elev", "temp", "larea", "lat")
sapply(1:4,
       \(i) apply(combn(1:4, i), 2,
                  \(x) paste(vars[x], collapse="+"))) %>%
  unlist %>%
  c(1, .) %>%
  paste("b0 ~", .) %>%
  lapply(lm) %>%
  sapply(AIC) %>%
  (\(x) x - min(x)) %>%
  plot
abline(h=2)
sapply(1:4,
       \(i) apply(combn(1:4, i), 2,
                  \(x) paste(vars[x], collapse="+"))) %>%
  unlist %>%
  c(1, .) %>%
  paste("b1 ~", .) %>%
  lapply(lm) %>%
  sapply(AIC) %>%
  (\(x) x - min(x)) %>%
  plot
abline(h=2)

modelbuilder <- function(x) {
  sapply(1:4,
         \(i) apply(combn(1:4, i), 2,
                    \(x) paste(vars[x], collapse="+"))) %>%
    unlist %>%
    c(1, .) %>%
    paste(x, "~", .)
}
elev0 <- elev
temp0 <- temp
larea0 <- larea
lat0 <- lat
b00 <- b0
b10 <- b1
rmse <- function(x,y) sqrt(mean((x-y)^2, na.rm=TRUE))

# rmses <- rep(NA, length(elev))
b0_preds <- b1_preds <- matrix(nrow=length(elev0), ncol=16)
for(i in 1:length(elev0)) {
  elev <- elev0[-i]
  temp <- temp0[-i]
  larea <- larea0[-i]
  lat <- lat0[-i]
  b0 <- b00[-i]
  b1 <- b10[-i]
  thelms_b0 <- lapply(modelbuilder("b0"), lm)
  thelms_b1 <- lapply(modelbuilder("b1"), lm)
  b0_preds[i,] <- sapply(thelms_b0, predict, newdata=data.frame(elev=elev0[i], temp=temp0[i], larea=larea0[i], lat=lat0[i]))
  b1_preds[i,] <- sapply(thelms_b1, predict, newdata=data.frame(elev=elev0[i], temp=temp0[i], larea=larea0[i], lat=lat0[i]))
}
b0_rmses <- sqrt(colMeans((as.vector(b00)-b0_preds)^2))
plot(b0_rmses)
b1_rmses <- sqrt(colMeans((as.vector(b10)-b1_preds)^2))
plot(b1_rmses)
modelbuilder("b0")[which.min(b0_rmses)]
modelbuilder("b1")[which.min(b1_rmses)]





cat('model {

  ## ---- Length ~ Age Component ---- ##
  for(i in whichdata_Lt) {   # fish-level data from lakes with paired lengths & ages
    # L[i] ~ dnorm(mu_Lt[i], tau_Lt)
    # Lpp[i] ~ dnorm(mu_Lt[i], tau_Lt)
    # mu_Lt[i] <- Linf[lake[i]]*(1-exp(-k[lake[i]]*(Age[i]-t0[lake[i]])))
    L[i] ~ dlnorm(logmu_Lt[i], tau_Lt)
    Lpp[i] ~ dlnorm(logmu_Lt[i], tau_Lt)
    logmu_Lt[i] <- log(Linf[lake[i]]*(1-exp(-k[lake[i]]*(Age[i]-t0[lake[i]]))))
  }

  tau_Lt <- pow(sig_Lt, -2)
  # sig_Lt ~ dunif(0, 100)
  # sig_Lt_prior ~ dunif(0, 100)
  sig_Lt ~ dunif(0, 3)
  sig_Lt_prior ~ dunif(0, 3)


  for(j in whichlakes_Lt) {   # lake level
    t0[j] ~ dnorm(mu_t0, tau_t0)T(,1)
    # t0_prior[j] ~ dnorm(mu_t0, tau_t0)
    k[j] ~ dlnorm(mu_k, tau_k)
    # k_prior[j] ~ dlnorm(mu_k, tau_k)

    # making plottable LVB curve envelopes for each lake
    for(ifit in 1:nfit) {
      Lfit[ifit, j] <- Linf[j]*(1-exp(-k[j]*(Agefit[ifit]-t0[j])))
    }
  }

  mu_t0 ~ dnorm(0, 1)   # was 1
  mu_t0_prior ~ dnorm(0, 1)
  # sig_t0 ~ dunif(0, 0.2)
  sig_t0 ~ dunif(0, 10)
  sig_t0_prior ~ dunif(0, 10)
  # sig_t0 ~ dexp(10)
  tau_t0 <- pow(sig_t0, -2)

  mu_k ~ dnorm(0, 0.1)
  mu_k_prior ~ dnorm(0, 0.1)
  sig_k ~ dunif(0, 3)
  sig_k_prior ~ dunif(0, 3)
  tau_k <- pow(sig_k, -2)



  ## ---- Linf ~ Area component ---- ##
  for(j in whichlakes_LA) {  # lake-level data where there are areas
    # Linf[j] ~ dnorm(mu_LA[j], tau_LA)T(1,)
    # mu_LA[j] <- gam * (1 - exp(-lam * (1 + log(Area[j]))))
    Linf[j] ~ dlnorm(logmu_LA[j], tau_LA)
    logmu_LA[j] <- log(gam * (1 - exp(-lam * (1 + log(Area[j])))))
  }
  for(j in whichlakes_LA_c) {  # lake-level priors where there are no areas
    Linf[j] ~ dnorm(600, 0.00001)T(1,)
    # Linf[j] ~ dnorm(mu_Linf_noarea, tau_Linf_noarea)
  }
  # mu_Linf_noarea ~ dnorm(600, 0.0001)
  # tau_Linf_noarea <- pow(sig_Linf_noarea, -2)
  # sig_Linf_noarea ~ dunif(0, 500)
  # sig_Linf_noarea_prior ~ dunif(0, 500)

  tau_LA <- pow(sig_LA, -2)
  # sig_LA ~ dunif(0, 300)
  # sig_LA_prior ~ dunif(0, 300)
  sig_LA ~ dunif(0, 3)
  sig_LA_prior ~ dunif(0, 3)
  # sig_LA ~ dexp(0.05)
  # sig_LA_prior ~ dexp(0.05)

  gam ~ dnorm(gam_lester, pow(cv_gam_lester*gam_lester, -2))T(0.01,)
  gam_prior ~ dnorm(gam_lester, pow(cv_gam_lester*gam_lester, -2))T(0.01,)
  lam ~ dnorm(lam_lester, pow(cv_lam_lester*lam_lester, -2))T(0.01,)
  lam_prior ~ dnorm(lam_lester, pow(cv_lam_lester*lam_lester, -2))T(0.01,)



  ## ---- Length quantile ~ Linf component ---- ##
  for(j in whichlakes_L) {
    qL[j] ~ dnorm(Linf[j], tau_L[j])
    qLpp[j] ~ dnorm(Linf[j], tau_L[j])
    tau_L[j] <- pow(sig_L[j], -2)
    sig_L[j] <- pow((eta_L^2) + ((zeta_L^2)/nL[j]), 0.5)
  }

  eta_L ~ dunif(0, eta_L_cap)
  eta_L_prior ~ dunif(0, eta_L_cap)
  zeta_L ~ dunif(0, zeta_L_cap)
  zeta_L_prior ~ dunif(0, zeta_L_cap)

  # eta_L ~ dexp(0.1)
  # eta_L_prior ~ dexp(0.1)
  # zeta_L ~ dexp(0.1)
  # zeta_L_prior ~ dexp(0.1)



  ## ---- log Weight ~ log Length component ---- ##
  # loop over all WL data
  for(i in whichdata_WL) {
    logW[i] ~ dnorm(mu_WL[i], tau_WL)
    # logWpp[i] ~ dnorm(mu_WL[i], tau_WL)    # might take this out if model seems ok
    mu_WL[i] <- b0[lake[i]] + b1[lake[i]]*logLc[i]
  }

  # loop over lakes with lat & area data
  for(j in whichlakes_WL) {
    b0[j] ~ dnorm(mu_b0[j], tau_b0)
    b0_interp[j] <- b0[j] - b1[j]*meanlogLc
    mu_b0[j] <- b0_int
                + Z[1]*b0_elev*elevc[j]
                + Z[2]*b0_temp*tempc[j]
                + Z[3]*b0_area*logareac[j]
                + Z[4]*b0_lat*latc[j]

    b1[j] ~ dnorm(mu_b1[j], tau_b1)
    mu_b1[j] <- b1_int
                + Z[5]*b1_elev*elevc[j]
                + Z[6]*b1_temp*tempc[j]
                + Z[7]*b1_area*logareac[j]
                + Z[8]*b1_lat*latc[j]
  }

  # # loop over lakes without area & lat data (c for complement)
  # # draw log(area) and lat from Normal distributions (mean & sd of actual data)
  # # corresponding b0, b0_interp, and b1 will effectively be from the post pred
  # for(jc in whichlakes_WL_c) {
  #   b0[jc] ~ dnorm(mu_b0[jc], tau_b0)
  #   b0_interp[jc] <- b0[jc] - b1[jc]*meanlogLc
  #   mu_b0[jc] <- b0_int
  #               # + b0_area*areasim[jc]
  #
  #   b1[jc] ~ dnorm(mu_b1[jc], tau_b1)
  #   mu_b1[jc] <- b1_int
  #               # + b1_lat*latsim[jc]
  #   areasim[jc] ~ dnorm(mean(logareac[whichlakes_WL]), pow(sd(logareac[whichlakes_WL]), -2))
  #   latsim[jc] ~ dnorm(mean(latc[whichlakes_WL]), pow(sd(latc[whichlakes_WL]), -2))
  # }

  # global priors
  sig_b0 ~ dunif(0, 10)
  sig_b0_prior ~ dunif(0, 10)
  tau_b0 <- pow(sig_b0, -2)

  sig_b1 ~ dunif(0, 10)
  sig_b1_prior ~ dunif(0, 10)
  tau_b1 <- pow(sig_b1, -2)

  b0_int ~ dnorm(0, 0.001)
  b0_elev ~ dnorm(0, 0.001)
  b0_temp ~ dnorm(0, 0.001)
  b0_area ~ dnorm(0, 0.001)
  b0_lat ~ dnorm(0, 0.001)
  b0_int_prior ~ dnorm(0, 0.001)
  b0_area_prior ~ dnorm(0, 0.001)

  b1_int ~ dnorm(0, 0.001)
  b1_elev ~ dnorm(0, 0.001)
  b1_temp ~ dnorm(0, 0.001)
  b1_area ~ dnorm(0, 0.001)
  b1_lat ~ dnorm(0, 0.001)
  b1_int_prior ~ dnorm(0, 0.001)
  b1_lat_prior ~ dnorm(0, 0.001)

  tau_WL <- pow(sig_WL, -2)
  sig_WL ~ dunif(0, 10)
  sig_WL_prior ~ dunif(0, 10)



  ## ---- Weight quantile ~ Winf component ---- ##
  for(j in whichlakes_W) {
    qW[j] ~ dnorm(Winf[j], tau_W[j])
    qWpp[j] ~ dnorm(Winf[j], tau_W[j])
    tau_W[j] <- pow(sig_W[j], -2)
    sig_W[j] <- pow((eta_W^2) + ((zeta_W^2)/nW[j]), 0.5)
  }

  eta_W ~ dunif(0, eta_W_cap)
  eta_W_prior ~ dunif(0, eta_W_cap)
  zeta_W ~ dunif(0, zeta_W_cap)
  zeta_W_prior ~ dunif(0, zeta_W_cap)

  # eta_W ~ dexp(0.1)
  # eta_W_prior ~ dexp(0.1)
  # zeta_W ~ dexp(0.1)
  # zeta_W_prior ~ dexp(0.1)



  ## ---- Winf from Linf! ---- ##
  for(j in whichlakes_WL) {  # alllakes
    Winf[j] <- exp(b0_interp[j] + b1[j]*log(Linf[j]))
  }


}', file=int_Winf_jags)




# defining a matrix to express all possible 1-variable scenarios
thegrid <- as.matrix(expand.grid(0:4, c(0,5:8)))
themat <- matrix(0, nrow=nrow(thegrid), ncol=8)
for(i in 1:nrow(thegrid)) {
  themat[i, thegrid[i,]] <- 1
}

# REDEFINING THIS, using the models we're actually considering
toconsider <- c(1,4,14,24)
thegrid <- as.matrix(expand.grid(0:4, c(0,5:8)))[toconsider,]
themat <- matrix(0, nrow=nrow(thegrid), ncol=8)
for(i in 1:nrow(thegrid)) {
  themat[i, thegrid[i,]] <- 1
}

# adding elevation and temp to data object
int_Winf_data$tempc <- morphometry$`Temp (C)` - mean(morphometry$`Temp (C)`, na.rm=TRUE)
int_Winf_data$elevc <- morphometry$`Elevation (m)` - mean(morphometry$`Elevation (m)`, na.rm=TRUE)

# JAGS controls
# niter <- 2*1000
# niter <- 20*1000
# niter <- 50*1000      # 50k in 15 minutes per
niter <- 100*1000    # about 13 hours total, for all 25 models
# niter <- 200*1000
# niter <- 500*1000

# ncores <- 3
ncores <- 8
# ncores <- min(10, parallel::detectCores()-1)

outs <- list()

for(i in 1:nrow(thegrid)) {
  int_Winf_data$Z <- themat[i,]




  ###### RUNNING THE MODEL #####
  {
    print(i)
    tstart <- Sys.time()
    print(tstart)
    outs[[i]] <- jagsUI::jags(model.file=int_Winf_jags, data=int_Winf_data,
                                      parameters.to.save=parameters,
                                      n.chains=ncores, parallel=T, n.iter=niter,
                                      n.burnin=niter/2, n.thin=niter/2000)
    print(Sys.time() - tstart)
  }
}

thedics <- sapply(outs, \(x) x$DIC)
plot(thedics)
thedics - min(thedics)

aa <- c("none", "elev", "temp", "area", "lat")
expand.grid(aa,aa)[which.min(thedics),]



## ------  trying this again, using k-fold cross validation

# JAGS controls
niter <- 2*1000
# niter <- 20*1000
# niter <- 50*1000      # 50k in 15 minutes per
# niter <- 100*1000    # about 13 hours total, for all 25 models
# niter <- 200*1000
# niter <- 500*1000

# ncores <- 3
ncores <- 8
# ncores <- min(10, parallel::detectCores()-1)

outs <- list()

for(i in 1:nrow(thegrid)) {
  int_Winf_data$Z <- themat[i,]




  ###### RUNNING THE MODEL #####
  {
    print(i)
    tstart <- Sys.time()
    print(tstart)
    outs[[i]] <- kfold(p="logW",
                       k=5,
                       model.file=int_Winf_jags, data=int_Winf_data,
                              # parameters.to.save=parameters,
                              n.chains=ncores, parallel=T, n.iter=niter,
                              n.burnin=niter/2, n.thin=niter/2000)
    print(Sys.time() - tstart)
  }
}
