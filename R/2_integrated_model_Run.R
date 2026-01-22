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
                         !is.na(laketrout$Latitude_WGS84)),  # note: only 14 use_fish are missing surface area
  # ^^^ this may be vestigial
  # whichdata_WL = which(laketrout$use_fish),  # note: only 14 use_fish are missing surface area

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
niter <- 20*1000
# niter <- 50*1000      # 50k in 9 minutes
# niter <- 100*1000
# niter <- 200*1000
# niter <- 500*1000

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



### a few helper functions and color vectors for plotting
estcols <- ifelse(morphometry$make_estimates, 4, 3)
datacols <- ifelse(int_Winf_data$alllakes %in% int_Winf_data$whichlakes_W, 4,
                   ifelse(int_Winf_data$alllakes %in% int_Winf_data$whichlakes_Lt, 2,
                          ifelse(int_Winf_data$alllakes %in% int_Winf_data$whichlakes_L, 3,
                                 ifelse(int_Winf_data$alllakes %in% int_Winf_data$whichlakes_LA, 5, 1))))
datarank <- rank(ifelse(int_Winf_data$alllakes %in% int_Winf_data$whichlakes_W, 1,
                        ifelse(int_Winf_data$alllakes %in% int_Winf_data$whichlakes_Lt, 2,
                               ifelse(int_Winf_data$alllakes %in% int_Winf_data$whichlakes_L, 3,
                                      ifelse(int_Winf_data$alllakes %in% int_Winf_data$whichlakes_LA, 4, 5)))), ties.method = "first")
datalegend <- function(loc="topleft") {
  legend(loc, legend=c("Wq, then", "Ages, then", "Lq, then", "Area, then","nothing"),
         lwd=3, col=c(4,2,3,5,1), cex=.5)
}

### should try tabulating datarank stuff more fully
# datacols uses:
# - presence of Wq (whichlakes_W)
# - presence of Age relationship (whichlakes_Lt)
# - presence of Lq (whichlakes_L)
# - presence of Area (whichlakes_LA)
nn <- 1:nrow(morphometry)
with(int_Winf_data,
     table(nn %in% whichlakes_W,
           nn %in% whichlakes_Lt))
## ^^ not convinced this is useful



### ------------------- convergence diagnostics  ------------------ ###

# how many nodes (or dim) per each named parameter
nbyname(int_Winf_jags_out)

# how many named parameters
length(nbyname(int_Winf_jags_out))

# plotting the Rhat values of all parameters
par(mfrow=c(1,1))
plotRhats(int_Winf_jags_out)

# trace plots for the nodes with the worst Rhat value for each named parameter
traceworstRhat(int_Winf_jags_out, parmfrow=c(3,3))





### ------------------- exploring model appropriateness  ------------------ ###

# comparing distributions of priors to the associated parameter
comparepriors(int_Winf_jags_out, parmfrow=c(3,3))

# plotting the correlation between all parameters
par(mfrow=c(1,1))
plotcor_jags(int_Winf_jags_out, p=c("Linf", "mu_LA", "t0", "k", "b0", "b1", "sig","gam","lam","eta","zeta"))



# cross plots between parameters ------ FILL THIS IN BETTER
crossplot(int_Winf_jags_out, p=c("eta_L","zeta_L"), drawblob=TRUE)
crossplot(int_Winf_jags_out, p=c("gam","lam"), drawblob=TRUE)
crossplot(int_Winf_jags_out, p=c("b0","b1"), drawblob=FALSE)




### Posterior predictive for qL
par(mfrow=c(1,2))
qq_postpred(int_Winf_jags_out$sims.list$qLpp[,int_Winf_data$whichlakes_L],
            y=int_Winf_data$qL[int_Winf_data$whichlakes_L],
            main="Post predictive for qL")
ts_postpred(int_Winf_jags_out$sims.list$qLpp[,int_Winf_data$whichlakes_L],
            y=int_Winf_data$qL[int_Winf_data$whichlakes_L],
            x=int_Winf_data$nL[int_Winf_data$whichlakes_L], log="x",
            xlab="n Length for quantiles")
par(mfcol=c(3,3))
plot_postpred(int_Winf_jags_out$sims.list$qLpp[,int_Winf_data$whichlakes_L],
              y=int_Winf_data$qL[int_Winf_data$whichlakes_L],
              x=log(int_Winf_data$nL[int_Winf_data$whichlakes_L]))


### plot how sig_L varies by lake
par(mfrow=c(1,2))
caterpillar(int_Winf_jags_out, "sig_L", col=estcols)
envelope(int_Winf_jags_out$sims.list$sig_L[, int_Winf_data$whichlakes_L],
         x=int_Winf_data$nL[int_Winf_data$whichlakes_L], log="x",
         xlab="n Length for quantiles", ylab="sig_L")



### Posterior predictive for qW
par(mfrow=c(1,2))
qq_postpred(int_Winf_jags_out$sims.list$qWpp[,int_Winf_data$whichlakes_W],
            y=int_Winf_data$qW[int_Winf_data$whichlakes_W],
            main="Post predictive for qW")
ts_postpred(int_Winf_jags_out$sims.list$qWpp[,int_Winf_data$whichlakes_W],
            y=int_Winf_data$qW[int_Winf_data$whichlakes_W],
            x=int_Winf_data$nW[int_Winf_data$whichlakes_W], log="x",
            xlab="n Wength for quantiles")
par(mfcol=c(3,3))
plot_postpred(int_Winf_jags_out$sims.list$qWpp[,int_Winf_data$whichlakes_W],
              y=int_Winf_data$qW[int_Winf_data$whichlakes_W],
              x=log(int_Winf_data$nW[int_Winf_data$whichlakes_W]))


### plot how sig_W varies by lake
par(mfrow=c(1,2))
caterpillar(int_Winf_jags_out, "sig_W", col=estcols)
envelope(int_Winf_jags_out$sims.list$sig_W[, int_Winf_data$whichlakes_W],
         x=int_Winf_data$nW[int_Winf_data$whichlakes_W], log="x",
         xlab="n Weight for quantiles", ylab="sig_W")





### ------------------- deep-dive into results  ------------------ ###


### LVB parameters
par(mfrow=c(1,2))
caterpillar(int_Winf_jags_out, "t0",
            col=estcols)
caterpillar(int_Winf_jags_out, "k",
            col=estcols)





## length ~ age (plots for all lakes where there are lengths and ages)
par(mfrow=c(3,3))
for(j in int_Winf_data$whichlakes_Lt) {
  envelope(int_Winf_jags_out$sims.list$Lfit[,,j], main=lakenames[j],
           ylim=range(0,int_Winf_data$L[!is.na(int_Winf_data$Age)]),
           # col=ifelse(j %in% int_Winf_data$whichlakes_L, 2, 4),
           col=datacols[j])
  points(int_Winf_data$Age[int_Winf_data$lake==j], int_Winf_data$L[int_Winf_data$lake==j])
  # legend("topleft", legend=c("L ~ Age only", "qL also"),
  #        fill=adjustcolor(c(4,2), alpha.f=.3), border=c(4,2), cex=.6)
  datalegend()
}


caterpillar1 <- function (df, p = NULL, x = NA, row = NULL, column = NULL, median = TRUE,
          mean = FALSE, ci = c(0.5, 0.95), lwd = 1, col = 4, add = FALSE,
          xlab = "", ylab = "", main = NULL, ylim = NULL, xax = NA, do_xax=TRUE,
          transform = c("none", "exp", "expit"), medlwd = lwd, medwd = 1,
          ...)
{
  if (!inherits(df, "jagsUI") & !inherits(df, c("matrix", "data.frame",
                                                "numeric", "integer"))) {
    stop("Input must be a data.frame or output from jagsUI::jags() plus parameter name")
  }
  if (inherits(df, "jagsUI") & length(p) != 1)
    stop("Need single parameter name in p= argument")
  if (class(df)[1] == "jagsUI" & !is.null(p)) {
    simslist <- df$sims.list
    if (all(names(simslist) != p))
      stop("No parameters with matching names")
    theparm <- simslist[names(simslist) == p][[1]]
    if (is.null(row) & is.null(column)) {
      df <- theparm
    }
    else {
      theparm <- simslist[names(simslist) == p][[1]]
      if (!is.null(row))
        df <- theparm[, row, ]
      if (!is.null(column))
        df <- theparm[, , column]
    }
    if (is.null(main) & length(p) == 1) {
      if (is.null(row) & is.null(column)) {
        main <- p
      }
      else {
        if (!is.null(row))
          main <- paste0(p, "[", row, ",]")
        if (!is.null(column))
          main <- paste0(p, "[,", column, "]")
      }
    }
  }
  if (is.null(main))
    main <- ""
  df <- as.matrix(df)
  transform <- match.arg(transform)
  if (transform == "exp")
    df <- exp(df)
  if (transform == "expit")
    df <- expit(df)
  ci <- rev(sort(ci))
  loq <- apply(df, 2, quantile, p = (1 - ci)/2, na.rm = T)
  hiq <- apply(df, 2, quantile, p = 1 - (1 - ci)/2, na.rm = T)
  if (length(ci) == 1) {
    loq <- t(as.matrix(loq))
    hiq <- t(as.matrix(hiq))
  }
  med <- apply(df, 2, median, na.rm = T)
  if (all(is.na(x)))
    x <- 1:ncol(df)
  d <- ifelse(length(x) > 1, diff(range(x, na.rm = TRUE))/length(x),
              1)
  nn <- ncol(df)
  if (all(is.na(xax)))
    xax <- names(df)
  lwds <- (1 + 2 * (1:length(ci) - 1)) * lwd
  if (!add) {
    if (is.null(ylim))
      ylim <- range(loq, hiq, na.rm = T)
    xlims <- range(x - (0.5 * d * (1 + (length(x) == 1))),
                   x + (0.5 * d * (1 + (length(x) == 1))), na.rm=TRUE)
    plot(NA, type = "l", xlim = xlims, xlab = xlab, ylab = ylab,
         main = main, ylim = ylim, xaxt = "n", ... = ...)
    if(do_xax) axis(1, x, labels = xax, las = list(...)$las)
  }
  if (median) {
    segments(x0 = x - 0.2 * d * medwd, x1 = x + 0.2 * d *
               medwd, y0 = med, y1 = med, col = col, lwd = medlwd,
             lend = 1)
  }
  if (mean)
    points(x, colMeans(df, na.rm = T), pch = 16, col = col)
  for (i in 1:length(ci)) segments(x0 = x, x1 = x, y0 = loq[i,
  ], y1 = hiq[i, ], col = col, lwd = lwds[i], lend = 1)
}   ##### ADD THIS PATCH TO JAGSHELPER

caterpillar_plus <- function(df, p, x=NA, col=datacols, ...) {
  nn <- length(df$q50[[p]])
  if(all(is.na(x))) x <- 1:nn
  caterpillar1(df=df, p=p, col=col, x=x, do_xax=FALSE, ...=...)
  axis(side=1, axTicks(1))
  text(x=x[1:nn], y=df$q50[[p]],
       labels=lakenames[1:nn],
       cex=.6, col=adjustcolor(col[1:nn], alpha.f=.6, blue.f=.8, red.f=.8, green.f=.8), srt=90)
  datalegend()
}


##### Linf #####

# vs data availability
caterpillar_plus(p="Linf", x=datarank,
                 df=int_Winf_jags_out, col=datacols)

# vs qL value
caterpillar_plus(p="Linf", x=int_Winf_data$qL, xlab="qL",
                 df=int_Winf_jags_out, col=datacols)
abline(0, 1, lty=3)

# vs Area
plot(x=morphometry$SurfaceArea_h, y=int_Winf_jags_out$q50$Linf,
     pch=16, col=datacols, log="x", ylim=c(0, max(int_Winf_jags_out$q97.5$Linf)),
     xlab="Surface Area (HA)", ylab="Linf (mm)")
caterpillar_plus(p="Linf", x=morphometry$SurfaceArea_h, #x=int_Winf_data$logareac,
                 df=int_Winf_jags_out, col=datacols,
                 add=TRUE, median=FALSE)
curve(int_Winf_data$gam_lester * (1 - exp(-int_Winf_data$lam_lester * (1 + log(x)))), add=TRUE, lty=2)

# vs Lester Linf
caterpillar_plus(p="Linf", x=morphometry$L_inf_lester, xlab="Lester Linf",
                 df=int_Winf_jags_out, col=datacols)
abline(0, 1, lty=3)




##### Winf #####

# vs data availability
caterpillar_plus(p="Winf", x=datarank,
                 df=int_Winf_jags_out, col=datacols)

# vs qW value
caterpillar_plus(p="Winf", x=int_Winf_data$qW, xlab="qW",
                 df=int_Winf_jags_out, col=datacols)
abline(0, 1, lty=3)

# vs Area
plot(x=morphometry$SurfaceArea_h, y=int_Winf_jags_out$q50$Winf,
     pch=16, col=datacols, log="x", ylim=c(0, max(int_Winf_jags_out$q97.5$Winf, na.rm=TRUE)),
     xlab="Surface Area (HA)", ylab="Winf (kg)")
caterpillar_plus(p="Winf", x=morphometry$SurfaceArea_h, #x=int_Winf_data$logareac,
                 df=int_Winf_jags_out, col=datacols,
                 add=TRUE, median=FALSE)
curve(exp(-19.56 +
            3.2*log(int_Winf_data$gam_lester * (1 - exp(-int_Winf_data$lam_lester * (1 + log(x)))))),
      add=TRUE, lty=2)

# vs Lester Linf
caterpillar_plus(p="Winf", x=morphometry$W_inf_lester, xlab="Lester Winf",
                 df=int_Winf_jags_out, col=datacols)
abline(0, 1, lty=3)



## possibly try this with
## - k and t0
## - b0 and b1
# versus Area, Elevation, Latitude, Temp, Linf
par(mfrow=c(1,2))
caterpillar_plus(p="k", x=NA,
                 df=int_Winf_jags_out, col=datacols)
caterpillar_plus(p="t0", x=NA,
                 df=int_Winf_jags_out, col=datacols)

caterpillar_plus(p="k", x=log(morphometry$SurfaceArea_h),
                 df=int_Winf_jags_out, col=datacols)
caterpillar_plus(p="t0", x=log(morphometry$SurfaceArea_h),,
                 df=int_Winf_jags_out, col=datacols)

caterpillar_plus(p="k", x=morphometry$`Temp (C)`,
                 df=int_Winf_jags_out, col=datacols)
caterpillar_plus(p="t0", x=morphometry$`Temp (C)`,
                 df=int_Winf_jags_out, col=datacols)

caterpillar_plus(p="k", x=morphometry$Latitude_WGS84,
                 df=int_Winf_jags_out, col=datacols)
caterpillar_plus(p="t0", x=morphometry$Latitude_WGS84,
                 df=int_Winf_jags_out, col=datacols)

caterpillar_plus(p="b0", x=NA,
                 df=int_Winf_jags_out, col=datacols)
caterpillar_plus(p="b1", x=NA,
                 df=int_Winf_jags_out, col=datacols)

caterpillar_plus(p="b0", x=log(morphometry$SurfaceArea_h),
                 df=int_Winf_jags_out, col=datacols)
caterpillar_plus(p="b1", x=log(morphometry$SurfaceArea_h),,
                 df=int_Winf_jags_out, col=datacols)

caterpillar_plus(p="b0", x=morphometry$`Temp (C)`,
                 df=int_Winf_jags_out, col=datacols)
caterpillar_plus(p="b1", x=morphometry$`Temp (C)`,
                 df=int_Winf_jags_out, col=datacols)

caterpillar_plus(p="b0", x=morphometry$Latitude_WGS84,
                 df=int_Winf_jags_out, col=datacols)
caterpillar_plus(p="b1", x=morphometry$Latitude_WGS84,
                 df=int_Winf_jags_out, col=datacols)



hasLW <- which(laketrout_Winf$n_Length > 0 & laketrout_Winf$n_Weight > 0)

plot(int_Winf_jags_out$q50$b0[hasLW] ~ log(morphometry$SurfaceArea_h)[hasLW])
summary(lm(int_Winf_jags_out$q50$b0[hasLW] ~ log(morphometry$SurfaceArea_h)[hasLW]))
plot(int_Winf_jags_out$q50$b1[hasLW] ~ log(morphometry$SurfaceArea_h)[hasLW])
summary(lm(int_Winf_jags_out$q50$b1[hasLW] ~ log(morphometry$SurfaceArea_h)[hasLW]))

plot(int_Winf_jags_out$q50$b0[hasLW] ~ morphometry$`Temp (C)`[hasLW])
summary(lm(int_Winf_jags_out$q50$b0[hasLW] ~ morphometry$`Temp (C)`[hasLW]))
plot(int_Winf_jags_out$q50$b1[hasLW] ~ morphometry$`Temp (C)`[hasLW])
summary(lm(int_Winf_jags_out$q50$b1[hasLW] ~ morphometry$`Temp (C)`[hasLW]))

plot(int_Winf_jags_out$q50$b0[hasLW] ~ morphometry$Latitude_WGS84[hasLW])
summary(lm(int_Winf_jags_out$q50$b0[hasLW] ~ morphometry$Latitude_WGS84[hasLW]))
plot(int_Winf_jags_out$q50$b1[hasLW] ~ morphometry$Latitude_WGS84[hasLW])
summary(lm(int_Winf_jags_out$q50$b1[hasLW] ~ morphometry$Latitude_WGS84[hasLW]))

summary(lm(int_Winf_jags_out$q50$b0[hasLW] ~ log(morphometry$SurfaceArea_h)[hasLW] +
             morphometry$`Temp (C)`[hasLW] +
             morphometry$Latitude_WGS84[hasLW]))

summary(lm(int_Winf_jags_out$q50$b1[hasLW] ~ log(morphometry$SurfaceArea_h)[hasLW] +
             morphometry$`Temp (C)`[hasLW] +
             morphometry$Latitude_WGS84[hasLW]))

### ^^^^^ this stuff is handled much better elsewhere




## ------- figure out better what is happening here


### Linf by data availability
par(mfrow=c(2,2))
caterpillar(int_Winf_jags_out, "Linf", col=3+int_Winf_data$alllakes %in% int_Winf_data$whichlakes_Lt)
legend("topleft", legend=c("has AGES", "no AGES"), col=4:3, lwd=3, cex=0.5)
caterpillar(int_Winf_jags_out, "Linf", col=3+int_Winf_data$alllakes %in% int_Winf_data$whichlakes_LA)
legend("topleft", legend=c("has AREA", "no AREA"), col=4:3, lwd=3, cex=0.5)
caterpillar(int_Winf_jags_out, "Linf", col=3+int_Winf_data$alllakes %in% int_Winf_data$whichlakes_L)
legend("topleft", legend=c("has Lq", "no Lq"), col=4:3, lwd=3, cex=0.5)
caterpillar(int_Winf_jags_out, "Linf", col=3+int_Winf_data$alllakes %in% int_Winf_data$whichlakes_W)
legend("topleft", legend=c("has Wq", "no Wq"), col=4:3, lwd=3, cex=0.5)

# # ages and lq
# caterpillar(int_Winf_jags_out, "Linf", col=3+(int_Winf_data$alllakes %in%
#               intersect(int_Winf_data$whichlakes_L, int_Winf_data$whichlakes_Lt)))
# legend("topleft", legend=c("has Lq and AGES", "no"), col=4:3, lwd=3, cex=0.5)
#
# # ages and lq and area
# caterpillar(int_Winf_jags_out, "Linf", col=3+(int_Winf_data$alllakes %in%
#                                              intersect(int_Winf_data$whichlakes_L,
#                                                        intersect(int_Winf_data$whichlakes_Lt,
#                                                        int_Winf_data$whichlakes_LA))))
# legend("topleft", legend=c("has Lq and AGES and AREA", "no"), col=4:3, lwd=3, cex=0.5)

# ages, then lq, then area
par(mfrow=c(1,1))
theorder <- ifelse(int_Winf_data$alllakes %in% int_Winf_data$whichlakes_Lt, 1,
                   ifelse(int_Winf_data$alllakes %in% int_Winf_data$whichlakes_L, 2,
                          ifelse(int_Winf_data$alllakes %in% int_Winf_data$whichlakes_LA, 3, 4))) %>% rank(ties.method = "first")
caterpillar(int_Winf_jags_out, "Linf", x=theorder,
            col=ifelse(int_Winf_data$alllakes %in% int_Winf_data$whichlakes_Lt, 4,
                       ifelse(int_Winf_data$alllakes %in% int_Winf_data$whichlakes_L, 2,
                              ifelse(int_Winf_data$alllakes %in% int_Winf_data$whichlakes_LA, 3, 1))))
legend("topleft", legend=c("AGES, then","Lq, then","AREA, then","nothing"), col=c(4,2,3,1), lwd=3, cex=0.5)
points(x=theorder[int_Winf_data$whichlakes_W], y=int_Winf_jags_out$q50$Linf[int_Winf_data$whichlakes_W])#, col=int_Winf_data$alllakes %in% int_Winf_data$whichlakes_W)
points(x=theorder[morphometry$make_estimates], y=int_Winf_jags_out$q50$Linf[morphometry$make_estimates], pch=16)

points(x=theorder, y=morphometry$L_inf_lester, pch="x")

caterpillar(int_Winf_jags_out, "Linf", x=datarank, col=datacols)
datalegend()

## ---

plot(NA, xlim=range(int_Winf_data$qL, na.rm=T),
     ylim=range(int_Winf_jags_out$q2.5$Linf, int_Winf_jags_out$q97.5$Linf, na.rm=TRUE),
     main="Linf", xlab="qL", ylab="")
caterpillar(int_Winf_jags_out, "Linf", col=ifelse(int_Winf_data$alllakes %in% int_Winf_data$whichlakes_Lt, 4,
                                               ifelse(int_Winf_data$alllakes %in% int_Winf_data$whichlakes_L, 2,
                                                      ifelse(int_Winf_data$alllakes %in% int_Winf_data$whichlakes_LA, 3, 1))),
            x=int_Winf_data$qL,
            add=TRUE)
legend("topleft", legend=c("AGES, then","Lq, then","AREA, then","nothing"), col=c(4,2,3,1), lwd=3, cex=0.5)
abline(0,1, lty=2)
# points(x=int_Winf_data$qL[int_Winf_data$whichlakes_LA],
#        y=int_Winf_jags_out$q50$mu_LA[int_Winf_data$whichlakes_LA])


plot(NA, xlim=range(int_Winf_data$qL, na.rm=T),
     ylim=range(int_Winf_jags_out$q2.5$Linf, int_Winf_jags_out$q97.5$Linf, na.rm=TRUE),
     main="Linf", xlab="qL", ylab="")
caterpillar(int_Winf_jags_out, "Linf", x=int_Winf_data$qL, col=datacols, add=T)
abline(0,1, lty=2)
datalegend()

#---

plot(NA, xlim=range(int_Winf_data$qW, na.rm=T),
     ylim=range(int_Winf_jags_out$q2.5$Winf, int_Winf_jags_out$q97.5$Winf, na.rm=TRUE),
     main="Winf", xlab="qW", ylab="")
caterpillar(int_Winf_jags_out, "Winf", x=int_Winf_data$qW, col=datacols, add=T)
abline(0,1, lty=2)
datalegend()


# ----------

par(mfrow=c(1,1))
plot(NA,
     ylim=range(int_Winf_jags_out$q2.5$Linf, int_Winf_jags_out$q97.5$Linf, na.rm=TRUE),
     xlim=range(log(int_Winf_data$Area), -1, na.rm=TRUE),
     main="Linf", xlab="Area (ha)", ylab="", xaxt="n")
lbls <- c(1, 5, 10, 50, 100, 500, 1000, 5000)
axis(side=1, at=log(lbls), labels=lbls)
thecols <- ifelse(int_Winf_data$alllakes %in% int_Winf_data$whichlakes_Lt, 4,
                  ifelse(int_Winf_data$alllakes %in% int_Winf_data$whichlakes_L, 2,
                         ifelse(int_Winf_data$alllakes %in% int_Winf_data$whichlakes_LA, 3, 1)))
caterpillar(int_Winf_jags_out, "Linf", col=datacols,#col=thecols,
            x=ifelse(is.na(int_Winf_data$Area), -1, log(int_Winf_data$Area)),
            add=TRUE)
points(x=ifelse(is.na(int_Winf_data$Area), -1, log(int_Winf_data$Area)),
       y=int_Winf_jags_out$q50$Linf,
       col=datacols, pch=16)
curve(int_Winf_data$gam_lester * (1 - exp(-int_Winf_data$lam_lester * (1 + x))), add=TRUE, lty=2)
# legend("topleft", legend=c("AGES, then","Lq, then","AREA, then","nothing","Lester priors"),
#        col=c(4,2,3,1,1), lwd=c(rep(3,4),1), lty=c(rep(1,4),2), cex=0.5)
datalegend()

# ---

# this bit was commented out
points(x=(ifelse(is.na(int_Winf_data$Area), -1, log(int_Winf_data$Area)))[int_Winf_data$whichlakes_L],
       y=int_Winf_data$qL[int_Winf_data$whichlakes_L])
points(x=(ifelse(is.na(int_Winf_data$Area), -1, log(int_Winf_data$Area)))[morphometry$make_estimates],
       y=int_Winf_jags_out$q50$Linf[morphometry$make_estimates], pch=16)

points(x=ifelse(is.na(int_Winf_data$Area), -1, log(int_Winf_data$Area)), y=morphometry$L_inf_lester, pch="x")

# more of the same, maybe these thoughts can be combined?
caterpillar(int_Winf_jags_out$sims.list$Linf[,int_Winf_data$whichlakes_L],
            x=int_Winf_data$qL[int_Winf_data$whichlakes_L],
            xlab="qL", main="subset with qL")
abline(0,1,lty=3)
caterpillar(int_Winf_jags_out$sims.list$Linf[,int_Winf_data$whichlakes_Lt],
            x=int_Winf_data$qL[int_Winf_data$whichlakes_Lt],
            xlab="qL", main="subset with Ages")
abline(0,1,lty=3)
int <- intersect(int_Winf_data$whichlakes_Lt, int_Winf_data$whichlakes_L)
caterpillar(int_Winf_jags_out$sims.list$Linf[,int],
            x=int_Winf_data$qL[int],
            xlab="qL", main="subset with qL and Ages")
abline(0,1,lty=3)


## logW ~ logL for each lake
## b0 ~ lat & area
## b1 ~ lat & area
## mu_b0 ~ lat & area
## mu_b1 ~ lat & area



### ------------------- compiling and saving output  ------------------ ###


### filling in summary output to save to external file

ff <- function(x, n=nrow(laketrout_Winf)) { # f for fill
  c(x, rep(NA, n - length(x)))
}

laketrout_Winf$Winf_est <- int_Winf_jags_out$q50$Winf %>% ff
laketrout_Winf$Winf_se <- int_Winf_jags_out$sd$Winf %>% ff
laketrout_Winf$Winf_cilo <- int_Winf_jags_out$q2.5$Winf %>% ff
laketrout_Winf$Winf_cihi <- int_Winf_jags_out$q97.5$Winf %>% ff

laketrout_Winf$Linf_est <- int_Winf_jags_out$q50$Linf %>% ff
laketrout_Winf$Linf_se <- int_Winf_jags_out$sd$Linf %>% ff
laketrout_Winf$Linf_cilo <- int_Winf_jags_out$q2.5$Linf %>% ff
laketrout_Winf$Linf_cihi <- int_Winf_jags_out$q97.5$Linf %>% ff

laketrout_Winf$t0_est <- int_Winf_jags_out$q50$t0 %>% ff
laketrout_Winf$t0_se <- int_Winf_jags_out$sd$t0 %>% ff
laketrout_Winf$k_est <- int_Winf_jags_out$q50$k %>% ff
laketrout_Winf$k_se <- int_Winf_jags_out$sd$k %>% ff

laketrout_Winf$b0_est <- int_Winf_jags_out$q50$b0_interp %>% ff
laketrout_Winf$b0_se <- int_Winf_jags_out$sd$b0_interp %>% ff
laketrout_Winf$b1_est <- int_Winf_jags_out$q50$b1 %>% ff
laketrout_Winf$b1_se <- int_Winf_jags_out$sd$b1 %>% ff

save_results
if(save_results) {
  write.csv(laketrout_Winf, file="Winf_estimates.csv")
}

##### add convergence diagnostics etc here!!

## make sure data object is appropriate
## make sure model output is complete (plots etc?)
## - pull out plotting stuff into its own 3_ script
