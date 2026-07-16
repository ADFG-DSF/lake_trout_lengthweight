library(jagsUI)
library(jagshelper)
library(tidyverse)


load(file="Rdata/laketrout_sampling_formodel.Rdata")

load(file="interim_posts/int_Winf_modelrun.Rdata")


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
                             make_estimates = morphometry$make_estimates,
                             censor = morphometry$censor) %>%
  mutate(n_Weight = ifelse(is.na(n_Weight), 0, n_Weight)) %>%
  mutate(n_Length = ifelse(is.na(n_Length), 0, n_Length)) %>%
  mutate(n_Age = ifelse(is.na(n_Age), 0, n_Age))
lakenames <- morphometry$LakeName


normalcurve <- function(p) {
  curve(dnorm(x,
              mean=int_Winf_jags_out$q50[[p]],
              sd=int_Winf_jags_out$sd[[p]]),
        add=TRUE)
}

## make priors from full model post predictive
# global
# par(mfrow=c(2,2))
# sig_Lt
plotdens(int_Winf_jags_out, "sig_Lt", exact=TRUE)
normalcurve("sig_Lt")
sig_Lt_mu <- int_Winf_jags_out$q50$sig_Lt
sig_Lt_tau <- int_Winf_jags_out$sd$sig_Lt^(-2)

# gam
plotdens(int_Winf_jags_out, "gam", exact=TRUE)
normalcurve("gam")
gam_mu <- int_Winf_jags_out$q50$gam
gam_tau <- int_Winf_jags_out$sd$gam^(-2)

# lam
plotdens(int_Winf_jags_out, "lam", exact=TRUE)
normalcurve("lam")
lam_mu <- int_Winf_jags_out$q50$lam
lam_tau <- int_Winf_jags_out$sd$lam^(-2)

# sig_WL
plotdens(int_Winf_jags_out, "sig_WL", exact=TRUE)
normalcurve("sig_WL")
sig_WL_mu <- int_Winf_jags_out$q50$sig_WL
sig_WL_tau <- int_Winf_jags_out$sd$sig_WL^(-2)

#
# lake level
# t0
nmc <- length(int_Winf_jags_out$sims.list$mu_t0)
normalcurve2 <- function(y) {
  curve(dnorm(x, mean=median(y), sd=sd(y)), add=TRUE)
}
lognormalcurve2 <- function(y) {
  curve(dlnorm(x, meanlog=median(log(y)), sdlog=sd(log(y))), add=TRUE)
}
t0_pp <- rnorm(n = nmc,
               mean = int_Winf_jags_out$sims.list$mu_t0,
               sd = int_Winf_jags_out$sims.list$sig_t0)
plotdens(t0_pp)
normalcurve2(t0_pp)
t0_mu <- median(t0_pp)     ## not totally thrilled with these
t0_tau <- sd(t0_pp)^(-2)

# k
k_pp <- rlnorm(n = nmc,
               meanlog = int_Winf_jags_out$sims.list$mu_k,
               sdlog = int_Winf_jags_out$sims.list$sig_k)
plotdens(k_pp)
# normalcurve2(k_pp)
lognormalcurve2(k_pp)
k_mu <- median(log(k_pp))
k_tau <- sd(log(k_pp))^(-2)

Linf_mu <- NA
Linf_tau <- NA
for(j in which(morphometry$censor)) {
  logmu <- with(int_Winf_jags_out$sims.list,
                log(gam*(1-exp(-lam*(1+log(morphometry$SurfaceArea_h[j]))))))
  Linf_pp <- rlnorm(n = nmc,
                    meanlog = logmu,
                    sdlog = int_Winf_jags_out$sims.list$sig_LA)
  plotdens(Linf_pp)
  lognormalcurve2(Linf_pp)
  Linf_mu[j] <- median(log(Linf_pp))
  Linf_tau[j] <- sd(log(Linf_pp))^(-2)
}

# sig_L - also
# sig_L[j] <- pow((eta_L^2) + ((zeta_L^2)/nL[j]), 0.5)
sig_L_mu <- NA
sig_L_tau <- NA
for(j in which(morphometry$censor)) {
  sig_L_pp <- with(int_Winf_jags_out$sims.list,
            sqrt((eta_L^2) + ((zeta_L^2)/int_Winf_data$nL[j])))
  plotdens(sig_L_pp)
  lognormalcurve2(sig_L_pp)
  # normalcurve2(sig_L_pp)  ## not super thrilled with this either
  sig_L_mu[j] <- median(log(sig_L_pp))
  sig_L_tau[j] <- sd(log(sig_L_pp))^(-2)
}

# sig_W - also
sig_W_mu <- NA
sig_W_tau <- NA
for(j in which(morphometry$censor)) {
  sig_W_pp <- with(int_Winf_jags_out$sims.list,
                   sqrt((eta_W^2) + ((zeta_W^2)/int_Winf_data$nW[j])))
  if(int_Winf_data$nW[j]) {
    plotdens(sig_W_pp)
    lognormalcurve2(sig_W_pp)
    # normalcurve2(sig_W_pp)  ## not super thrilled with this either
  }
  sig_W_mu[j] <- median(log(sig_W_pp))
  sig_W_tau[j] <- sd(log(sig_W_pp))^(-2)
}

# b0 for each lake
# b1
b0_mu <- NA
b0_tau <- NA
for(j in which(morphometry$censor)) {
  mu_b0_pp <- with(int_Winf_jags_out$sims.list,
                   b0_int + b0_area*int_Winf_data$logareac[j])
  # plotdens(mu_b0_pp)
  b0_pp <- rnorm(n = nmc,
                 mean = mu_b0_pp,
                 sd = int_Winf_jags_out$sims.list$sig_b0)
  plotdens(b0_pp)
  normalcurve2(b0_pp)
  b0_mu[j] <- median(b0_pp)
  b0_tau[j] <- sd(b0_pp)^(-2)
}


b1_pp <- rnorm(n = nmc,
               mean = int_Winf_jags_out$sims.list$mu_b1,
               sd = int_Winf_jags_out$sims.list$sig_b1)
plotdens(b1_pp)
normalcurve2(b1_pp)
b1_mu <- median(b1_pp)
b1_tau <- sd(b1_pp)^(-2)


int_Winf_jags_newlakes <- tempfile()
cat('model {

  ## ---- Length ~ Age Component ---- ##
  for(i in whichdata_Lt) {   # fish-level data from lakes with paired lengths & ages
    L[i] ~ dlnorm(logmu_Lt[i], tau_Lt)
    Lpp[i] ~ dlnorm(logmu_Lt[i], tau_Lt)
    logmu_Lt[i] <- log(Linf[lake[i]]*(1-exp(-k[lake[i]]*(Age[i]-t0[lake[i]]))))
  }

  tau_Lt <- pow(sig_Lt, -2)
  # sig_Lt ~ dunif(0, 3)   #######
  # sig_Lt_prior ~ dunif(0, 3)  ######
  sig_Lt ~ dnorm(sig_Lt_mu, sig_Lt_tau)


  for(j in whichlakes_Lt) {   # lake level
    # t0[j] ~ dnorm(mu_t0, tau_t0)T(,1)
    # k[j] ~ dlnorm(mu_k, tau_k)
    t0[j] ~ dnorm(t0_mu, t0_tau)T(,1)
    k[j] ~ dlnorm(k_mu, k_tau)


    # making plottable LVB curve envelopes for each lake
    for(ifit in 1:nfit) {
      Lfit[ifit, j] <- Linf[j]*(1-exp(-k[j]*(Agefit[ifit]-t0[j])))
    }
  }

  # mu_t0 ~ dnorm(0, 1)   # was 1
  # mu_t0_prior ~ dnorm(0, 1)
  # # sig_t0 ~ dunif(0, 0.2)
  # sig_t0 ~ dunif(0, 10)
  # sig_t0_prior ~ dunif(0, 10)
  # tau_t0_prior <- pow(sig_t0_prior, -2)
  # # sig_t0 ~ dexp(10)
  # tau_t0 <- pow(sig_t0, -2)
  #
  # mu_k ~ dnorm(0, 0.1)
  # mu_k_prior ~ dnorm(0, 0.1)
  # sig_k ~ dunif(0, 3)
  # sig_k_prior ~ dunif(0, 3)
  # tau_k_prior <- pow(sig_k_prior, -2)
  # tau_k <- pow(sig_k, -2)



  ## ---- Linf ~ Area component ---- ##
  for(j in whichlakes_LA) {  # lake-level data where there are areas
    # Linf[j] ~ dlnorm(logmu_LA[j], tau_LA)
    # logmu_LA[j] <- log(gam * (1 - exp(-lam * (1 + log(Area[j])))))
    Linf[j] ~ dlnorm(Linf_mu[j], Linf_tau[j])
  }
  # for(j in whichlakes_LA_c) {  # lake-level priors where there are no areas
  #   Linf[j] ~ dnorm(600, 0.00001)T(1,)
  # }

  # tau_LA <- pow(sig_LA, -2)
  # # sig_LA ~ dunif(0, 300)
  # # sig_LA_prior ~ dunif(0, 300)
  # sig_LA ~ dunif(0, 3)
  # sig_LA_prior ~ dunif(0, 3)
  # # sig_LA ~ dexp(0.05)
  # # sig_LA_prior ~ dexp(0.05)

  # gam ~ dnorm(gam_lester, pow(sd_gam_lester, -2))T(0.01,)
  # gam_prior ~ dnorm(gam_lester, pow(sd_gam_lester, -2))T(0.01,)
  # lam ~ dnorm(lam_lester, pow(sd_lam_lester, -2))T(0.01,)
  # lam_prior ~ dnorm(lam_lester, pow(sd_lam_lester, -2))T(0.01,)



  ## ---- Length quantile ~ Linf component ---- ##
  for(j in whichlakes_L) {
    qL[j] ~ dnorm(Linf[j], tau_L[j])
    qLpp[j] ~ dnorm(Linf[j], tau_L[j])
    tau_L[j] <- pow(sig_L[j], -2)
    # sig_L[j] <- pow((eta_L^2) + ((zeta_L^2)/nL[j]), 0.5)
    sig_L[j] ~ dlnorm(sig_L_mu[j], sig_L_tau[j])
  }

  # eta_L ~ dunif(0, eta_L_cap)
  # eta_L_prior ~ dunif(0, eta_L_cap)
  # zeta_L ~ dunif(0, zeta_L_cap)
  # zeta_L_prior ~ dunif(0, zeta_L_cap)
  #
  # # eta_L ~ dexp(0.1)
  # # eta_L_prior ~ dexp(0.1)
  # # zeta_L ~ dexp(0.1)
  # # zeta_L_prior ~ dexp(0.1)



  ## ---- log Weight ~ log Length component ---- ##
  # loop over all WL data
  for(i in whichdata_WL) {
    logW[i] ~ dnorm(mu_WL[i], tau_WL)
    # logWpp[i] ~ dnorm(mu_WL[i], tau_WL)    # might take this out if model seems ok
    mu_WL[i] <- b0[lake[i]] + b1[lake[i]]*logLc[i]
  }

  # loop over lakes with lat & area data
  for(j in whichlakes_WL) {
    # b0[j] ~ dnorm(mu_b0[j], tau_b0)
    b0[j] ~ dnorm(b0_mu[j], b0_tau[j])
    b0_interp[j] <- b0[j] - b1[j]*meanlogLc
    # mu_b0[j] <- b0_int
    #              + b0_area*logareac[j]

    # b1[j] ~ dnorm(mu_b1[j], tau_b1)
    b1[j] ~ dnorm(b1_mu, b1_tau)
    # mu_b1[j] <- b1_int
    #             # + b1_lat*latc[j]
  }

  # # global priors
  # sig_b0 ~ dunif(0, 10)
  # sig_b0_prior ~ dunif(0, 10)
  # tau_b0 <- pow(sig_b0, -2)
  #
  # sig_b1 ~ dunif(0, 10)
  # sig_b1_prior ~ dunif(0, 10)
  # tau_b1 <- pow(sig_b1, -2)
  #
  # b0_int ~ dnorm(0, 0.001)
  # b0_area ~ dnorm(0, 0.001)
  # b0_int_prior ~ dnorm(0, 0.001)
  # b0_area_prior ~ dnorm(0, 0.001)
  #
  # b1_int ~ dnorm(0, 0.001)
  # b1_lat ~ dnorm(0, 0.001)
  # b1_int_prior ~ dnorm(0, 0.001)
  # b1_lat_prior ~ dnorm(0, 0.001)
  #
  tau_WL <- pow(sig_WL, -2)
  # sig_WL ~ dunif(0, 10)
  sig_WL ~ dnorm(sig_WL_mu, sig_WL_tau)
  # sig_WL_prior ~ dunif(0, 10)



  ## ---- Weight quantile ~ Winf component ---- ##
  for(j in whichlakes_W) {
    qW[j] ~ dnorm(Winf[j], tau_W[j])
    qWpp[j] ~ dnorm(Winf[j], tau_W[j])
    tau_W[j] <- pow(sig_W[j], -2)
    # sig_W[j] <- pow((eta_W^2) + ((zeta_W^2)/nW[j]), 0.5)
    sig_W[j] ~ dlnorm(sig_W_mu[j], sig_W_tau[j])
  }

  # eta_W ~ dunif(0, eta_W_cap)
  # eta_W_prior ~ dunif(0, eta_W_cap)
  # zeta_W ~ dunif(0, zeta_W_cap)
  # zeta_W_prior ~ dunif(0, zeta_W_cap)
  #
  # # eta_W ~ dexp(0.1)
  # # eta_W_prior ~ dexp(0.1)
  # # zeta_W ~ dexp(0.1)
  # # zeta_W_prior ~ dexp(0.1)



  ## ---- Winf from Linf! ---- ##
  for(j in whichlakes_WL) {  # alllakes
    Winf[j] <- exp(b0_interp[j] + b1[j]*log(Linf[j]))
  }


}', file=int_Winf_jags_newlakes)



# bundle data to pass into JAGS
q_input <- 0.95   # assumed quantile value for asymptotic size
int_Winf_data_newlakes <- list(

  # fish-level data
  Age = laketrout$Age,
  L = laketrout$ForkLength_mm,
  logW = log(laketrout$Weight_g/1000),
  logLc = log(laketrout$ForkLength_mm) - mean(log(laketrout$ForkLength_mm), na.rm=TRUE),
  lake = as.numeric(as.character(laketrout$LakeNum)),
  whichdata_WL = which(laketrout$use_fish &
                         # !laketrout$censor &
                         laketrout$censor &  #############
                         !is.na(laketrout$ForkLength_mm) &
                         !is.na(laketrout$SurfaceArea_h) &
                         !is.na(laketrout$Latitude_WGS84)),  # note: only 14 use_fish are missing surface area
  # ^^^ this may be vestigial
  # whichdata_WL = which(laketrout$use_fish),  # note: only 14 use_fish are missing surface area

  # lake-level data
  Area = laketrout_Winf$Area_ha,
  latc = morphometry$Latitude_WGS84 - mean(morphometry$Latitude_WGS84, na.rm=TRUE),
  logareac = log(morphometry$SurfaceArea_h) - mean(log(morphometry$SurfaceArea_h), na.rm=TRUE),

  whichlakes_L = which(morphometry$use_fish & laketrout_Winf$n_Length > 10 &
                         morphometry$censor),  ###########
  qL = tapply(laketrout$ForkLength_mm, laketrout$LakeNum, quantile, p=q_input, na.rm=TRUE),
  nL = laketrout_Winf$n_Length,
  whichlakes_W = which(morphometry$use_fish & laketrout_Winf$n_Weight > 10 &
                         # !(morphometry$LakeName %in% c("Skilak Lake", "Big Lake (Healy)"))
                         !(morphometry$LakeName %in% c("Big Lake (Healy)")) &
                         morphometry$censor),  ###########
  qW = tapply(laketrout$Weight_g/1000, laketrout$LakeNum, quantile, p=q_input, na.rm=TRUE),
  nW = laketrout_Winf$n_Weight,

  #############
  whichlakes_Lt = which(morphometry$censor & morphometry$use_fish & laketrout_Winf$n_Age > 0 & laketrout_Winf$n_Length > 0),
  whichlakes_LA = which(morphometry$censor & !is.na(laketrout_Winf$Area_ha)),# & laketrout_Winf$n_Length > 0),
  whichlakes_LA_c = which(morphometry$censor & !(!is.na(laketrout_Winf$Area_ha))),# & laketrout_Winf$n_Length > 0)),
  whichlakes_WL = which(morphometry$censor & !is.na(laketrout_Winf$Lat) & !is.na(laketrout_Winf$Area_ha)),
  whichlakes_WL_c = which(morphometry$censor & is.na(laketrout_Winf$Lat) | is.na(laketrout_Winf$Area_ha)),
  # alllakes = 1:nrow(laketrout_Winf),

  # inputs for fitted vals of L~Age
  Agefit = 1:max(laketrout$Age, na.rm=TRUE),
  nfit = max(laketrout$Age, na.rm=TRUE),

  # priors from Lester
  gam_lester = 957,
  # cv_gam_lester = 0.2*2.5,
  sd_gam_lester = 54,
  lam_lester = 0.14,
  # cv_lam_lester = 1*2.5,
  sd_lam_lester = 0.016,

  # global values
  eta_L_cap = 500,
  zeta_L_cap = 2000,
  eta_W_cap = 5,
  zeta_W_cap = 20,
  meanlogLc = mean(log(laketrout$ForkLength_mm), na.rm=TRUE),

  # priors from all-lakes model post pred
  sig_Lt_mu = sig_Lt_mu,
  sig_Lt_tau = sig_Lt_tau,
  gam_mu = gam_mu,
  gam_tau = gam_tau,
  lam_mu = lam_mu,
  lam_tau = lam_tau,
  sig_WL_mu = sig_WL_mu,
  sig_WL_tau = sig_WL_tau,
  t0_mu = t0_mu,
  t0_tau = t0_tau,
  k_mu = k_mu,
  k_tau = k_tau,
  Linf_mu = Linf_mu,
  Linf_tau = Linf_tau,
  sig_L_mu = sig_L_mu,
  sig_L_tau = sig_L_tau,
  sig_W_mu = sig_W_mu,
  sig_W_tau = sig_W_tau,
  b0_mu = b0_mu,
  b0_tau = b0_tau,
  b1_mu = b1_mu,
  b1_tau = b1_tau
)

# can i make this within the above?
int_Winf_data_newlakes$whichdata_Lt <- which((laketrout$LakeNum %in% int_Winf_data$whichlakes_Lt) &
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
niter <- 2*1000
# niter <- 20*1000
# niter <- 50*1000      # 50k in 9 minutes
# niter <- 100*1000  # 37 min
# niter <- 200*1000
# niter <- 500*1000  # 2.7 hrs
# niter <- 2000*1000  # 12 hrs on laptop
# niter <- 5000*1000  # 14 hrs on desktop

# ncores <- 3
ncores <- 8
# ncores <- min(10, parallel::detectCores()-1)


###### RUNNING THE MODEL #####
{
  tstart <- Sys.time()
  print(tstart)
  int_Winf_jags_out_newlakes <- jagsUI::jags(model.file=int_Winf_jags_newlakes, data=int_Winf_data_newlakes,
                                    parameters.to.save=parameters,
                                    n.chains=ncores, parallel=T, n.iter=niter,
                                    n.burnin=niter/2, n.thin=niter/2000)
  print(Sys.time() - tstart)
}




###### SAVING INTERIM RESULTS #####
save(int_Winf_jags_out,
     int_Winf_data,
     laketrout_Winf,
     lakenames,
     niter, ncores,
     file="interim_posts/int_Winf_modelrun.Rdata")
