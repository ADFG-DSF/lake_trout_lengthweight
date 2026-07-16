# library(jpeg)
# # Read the JPG file
# img <- readJPEG("raw_data/lester_linf_area.jpg", native = TRUE)
#
# pars <- par(no.readonly = TRUE)
#
# # Create an empty plot area
# plot(0:1, 0:1, par(mar=c(0,0,0,0)))#, type = "n", ann = FALSE, axes = FALSE)
#
# # Render the image inside the plot window
# rasterImage(img, 0, 0, 1, 1)
#
# # locator(n=2, type="p", pch="+")
#
# corners <- locator(n=2, type="p", pch="+")
# # should be x=log(c(5,500000)), y=c(1000,0)
#
# dots <- locator(n=129, type="p", pch="+")
#
#
# origx <- corners$x[1]
# origy <- corners$y[2]
# points(origx, origy)
#
# rangex <- abs(diff(corners$x))
# rangey <- abs(diff(corners$y))
#
# scalex <- abs(diff(log(c(5,500000))))/rangex
# scaley <- 1000/rangey
#
# datay <- (dots$y-origy)*scaley
# datax <- exp((dots$x-origx)*scalex + log(5))
#
# par(pars)
# plot(datax, datay, log="x")
# abline(v=c(50,50000, 500000))
# save(datax, datay, file="raw_data/dotsfromdots.Rdata")


load(file="raw_data/dotsfromdots.Rdata")



load(file="Rdata/laketrout_sampling_formodel.Rdata")
load(file="interim_posts/int_Winf_modelrun.Rdata")
# load(file="interim_posts/int_Winf_modelrun_26_05_30.Rdata")
library(magrittr)
library(jagshelper)
# plot(datax, datay, log="x")

par(mar=c(3, 3, 1, 1)+0.1)
par(family="serif")
par(mgp=c(2, 1, 0))

# y_cat2(int_Winf_jags_out$sims.list$Linf,
#        xval=morphometry$SurfaceArea_h,
#        highlight=morphometry$make_estimates,
#        x0=TRUE,
#        log="x",
#        xlab="Surface area (ha)",
#        ylab="Asymptotic length (mm FL)")

datalakes <- rep(FALSE, length(int_Winf_jags_out$q50$k))
datalakes[union(union(int_Winf_data$whichlakes_L,
      int_Winf_data$whichlakes_W),
      int_Winf_data$whichlakes_Lt)] <- TRUE

caterpillar(int_Winf_jags_out, p="Linf",
            x = morphometry$SurfaceArea_h,
            # col = plotcol,
            col = ifelse(!datalakes, adjustcolor("grey", alpha.f=.3), 1),
            log="x", main="", ylim_add=0,
            xlim=range(morphometry$SurfaceArea_h, datax, na.rm=TRUE),
            xlab="Surface area (ha)",
            # ylab="Asymptotic length (mm FL)",
            ylab=expression(paste(L[infinity], " (mm FL)")))
curve(int_Winf_jags_out$q50$gam * (1 - exp(-int_Winf_jags_out$q50$lam * (1 + log(x)))), add=TRUE, lty=2)
curve(int_Winf_data$gam_lester * (1 - exp(-int_Winf_data$lam_lester * (1 + log(x)))), add=TRUE, lty=3)

points(datax, datay, col=4)




# let's try re-fitting the Lester dots!!
# specify model, which is written to a temporary file
dots_jags <- tempfile()
cat('model {
  for(i in 1:n) {
    # y[i] ~ dlnorm(logmu[i], tau)
    # logmu[i] <- log(gam * (1 - exp(-lam * (1 + log(x[i])))))
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- gam * (1 - exp(-lam * (1 + log(x[i]))))
  }

  tau <- pow(sig, -2)
  # sig ~ dunif(0, 3)
  sig ~ dunif(0, 300)

  gam ~ dnorm(gam_lester, pow(cv_gam_lester*gam_lester, -2))T(0.01,)
  lam ~ dnorm(lam_lester, pow(cv_lam_lester*lam_lester, -2))T(0.01,)

}', file=dots_jags)

# bundle data to pass into JAGS
dots_data <- list(x=datax,
                  y=datay,
                  n=length(datax),
                  gam_lester=int_Winf_data$gam_lester,
                  lam_lester=int_Winf_data$lam_lester,
                  cv_gam_lester=5*int_Winf_data$cv_gam_lester,
                  cv_lam_lester=5*int_Winf_data$cv_lam_lester)

# JAGS controls
niter <- 100000
# ncores <- 3
ncores <- min(10, parallel::detectCores()-1)

{
  tstart <- Sys.time()
  print(tstart)
  dots_jags_out <- jagsUI::jags(model.file=dots_jags, data=dots_data,
                                parameters.to.save=c("sig","gam","lam"),
                                n.chains=ncores, parallel=T, n.iter=niter,
                                n.burnin=niter/2, n.thin=niter/2000)
  print(Sys.time() - tstart)
}

# nbyname(dots_jags_out)
plotRhats(dots_jags_out)
traceworstRhat(dots_jags_out, parmfrow = c(3, 3))

dots_jags_out$q50
dots_jags_out$sd



########## Trying full model run using the inferences above as priors for gamma & lambda ############

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
  tau_t0_prior <- pow(sig_t0_prior, -2)
  # sig_t0 ~ dexp(10)
  tau_t0 <- pow(sig_t0, -2)

  mu_k ~ dnorm(0, 0.1)
  mu_k_prior ~ dnorm(0, 0.1)
  sig_k ~ dunif(0, 3)
  sig_k_prior ~ dunif(0, 3)
  tau_k_prior <- pow(sig_k_prior, -2)
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

  # gam ~ dnorm(gam_lester, pow(cv_gam_lester*gam_lester, -2))T(0.01,)
  # gam_prior ~ dnorm(gam_lester, pow(cv_gam_lester*gam_lester, -2))T(0.01,)
  # lam ~ dnorm(lam_lester, pow(cv_lam_lester*lam_lester, -2))T(0.01,)
  # lam_prior ~ dnorm(lam_lester, pow(cv_lam_lester*lam_lester, -2))T(0.01,)
  gam ~ dnorm(gam_lester, pow(gam_lester_sd, -2))T(0.01,)
  gam_prior ~ dnorm(gam_lester, pow(gam_lester_sd, -2))T(0.01,)
  lam ~ dnorm(lam_lester, pow(lam_lester_sd, -2))T(0.01,)
  lam_prior ~ dnorm(lam_lester, pow(lam_lester_sd, -2))T(0.01,)



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
                 + b0_area*logareac[j]

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
                         !is.na(laketrout$Latitude_WGS84)),  # note: only 14 use_fish are missing surface area
  # ^^^ this may be vestigial
  # whichdata_WL = which(laketrout$use_fish),  # note: only 14 use_fish are missing surface area

  # lake-level data
  Area = laketrout_Winf$Area_ha,
  latc = morphometry$Latitude_WGS84 - mean(morphometry$Latitude_WGS84, na.rm=TRUE),
  logareac = log(morphometry$SurfaceArea_h) - mean(log(morphometry$SurfaceArea_h), na.rm=TRUE),

  whichlakes_L = which(morphometry$use_fish & laketrout_Winf$n_Length > 10 &
                         !(morphometry$LakeName %in% c("Skilak Lake", "Big Lake (Healy)"))),
  qL = tapply(laketrout$ForkLength_mm, laketrout$LakeNum, quantile, p=q_input, na.rm=TRUE),
  nL = laketrout_Winf$n_Length,
  whichlakes_W = which(morphometry$use_fish & laketrout_Winf$n_Weight > 10 &
                         !(morphometry$LakeName %in% c("Skilak Lake", "Big Lake (Healy)"))),
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
  # cv_gam_lester = 0.2*2.5,
  lam_lester = 0.14,
  # cv_lam_lester = 1*2.5,
  gam_lester_sd = 10*dots_jags_out$sd$gam,
  lam_lester_sd = 10*dots_jags_out$sd$lam,

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
niter <- 50*1000      # 50k in 9 minutes
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
  int_Winf_jags_out_dotprior <- jagsUI::jags(model.file=int_Winf_jags, data=int_Winf_data,
                                    parameters.to.save=parameters,
                                    n.chains=ncores, parallel=T, n.iter=niter,
                                    n.burnin=niter/2, n.thin=niter/2000)
  print(Sys.time() - tstart)
}

datalakes <- rep(FALSE, length(int_Winf_jags_out_dotprior$q50$k))
datalakes[union(union(int_Winf_data$whichlakes_L,
                      int_Winf_data$whichlakes_W),
                int_Winf_data$whichlakes_Lt)] <- TRUE


# caterpillar(int_Winf_jags_out, p="Linf",
#             x = morphometry$SurfaceArea_h,
#             # col = plotcol,
#             col = ifelse(!datalakes, adjustcolor("grey", alpha.f=.3), 1),
#             log="x", main="", ylim_add=0,
#             xlim=range(morphometry$SurfaceArea_h, datax, na.rm=TRUE),
#             xlab="Surface area (ha)",
#             # ylab="Asymptotic length (mm FL)",
#             ylab=expression(paste(L[infinity], " (mm FL)")))
# curve(int_Winf_jags_out$q50$gam * (1 - exp(-int_Winf_jags_out$q50$lam * (1 + log(x)))), add=TRUE, lty=2)
# curve(int_Winf_data$gam_lester * (1 - exp(-int_Winf_data$lam_lester * (1 + log(x)))), add=TRUE, lty=3)
#
# points(datax, datay, col=4)


caterpillar(int_Winf_jags_out_dotprior, p="Linf",
            x = morphometry$SurfaceArea_h,
            # col = plotcol,
            col = ifelse(!datalakes, adjustcolor("grey", alpha.f=.3), 1),
            log="x", main="", ylim_add=0,
            xlim=range(morphometry$SurfaceArea_h, datax, na.rm=TRUE),
            xlab="Surface area (ha)",
            # ylab="Asymptotic length (mm FL)",
            ylab=expression(paste(L[infinity], " (mm FL)")))
curve(int_Winf_jags_out_dotprior$q50$gam * (1 - exp(-int_Winf_jags_out_dotprior$q50$lam * (1 + log(x)))), add=TRUE, lty=2)
curve(int_Winf_data$gam_lester * (1 - exp(-int_Winf_data$lam_lester * (1 + log(x)))), add=TRUE, lty=3)

points(datax, datay, col=4)

legend("topleft", lwd=c(3,3,NA), col=c(1,"grey90",4), pch=c(NA,NA,1),
       legend=c("AK data lakes","AK no-data lakes","Lester lakes"))

lester_resids <- datay - int_Winf_data$gam_lester * (1 - exp(-int_Winf_data$lam_lester * (1 + log(datax))))
hist(lester_resids)
sd(lester_resids)

ak_resids <- int_Winf_jags_out_dotprior$q50$Linf[datalakes] -
  int_Winf_data$gam_lester * (1 - exp(-int_Winf_data$lam_lester *
                                        (1 + log(morphometry$SurfaceArea_h[datalakes]))))
hist(ak_resids)
sd(ak_resids, na.rm=TRUE)

par(mfrow=c(2,1))
hist(lester_resids, xlim=c(-300,300))
hist(ak_resids, xlim=c(-300,300))

par(mfrow=c(1,1))
resid_names <- morphometry$LakeName[datalakes]
plot(x=rank(ak_resids), y=ak_resids, pch=16, ylim=c(-200,400))
grid(nx=NA, ny=NULL)
text(x=rank(ak_resids), y=ak_resids, labels=resid_names, srt=90,
     pos=4)
