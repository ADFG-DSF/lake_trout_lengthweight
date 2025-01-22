## NEED TO WRITE A HEADER



## loading packages
library(tidyverse)   # for data manipulation
library(jagsUI)      # for running JAGS
library(jagshelper)  # for plotting JAGS output & model diagnostics



## Whether to save results from long JAGS runs, etc.
save_results <- FALSE




## load data, filter bad observations BUT NOT observations with missing data

# read lake-level data
morphometry1 <- read_csv("flat_data/lake_morphometry3.csv", skip=1) %>%
  filter(!is.na(LakeName)) %>%
  filter(!(is.na(Latitude_WGS84) & is.na(`Elevation (m)`) & is.na(`Temp (C)`) & is.na(SurfaceArea_h))) %>%
  mutate(use_fish = `Include in "Alaskanizing" Modeling Exercise` %in% c("Yes","yes")) %>%
  mutate(make_estimates = `Potentially Include in Lake Trout Management Plan` != "No") %>%
  arrange(LakeName)

# # is lake name unique?  YES
# sum(!is.na(morphometry1$LakeName))
# length(unique(morphometry1$LakeName))




# read fish-level data
# match lake names to lake-level (morphometry) data as needed
# take out rows where !use_fish
laketrout_all <- read_csv("flat_data/length_weight3.csv", skip=1) %>%
  mutate(LakeName = ifelse(LakeName == "Donnelly Lake", "Donnelly Lake (Richardson Highway)", LakeName)) %>%
  mutate(LakeName = ifelse(LakeName == "Four Mile Lake", "Fourmile Lake (Taylor Highway)", LakeName)) %>%
  mutate(LakeName = ifelse(LakeName == "Lost Lake", "Lost Lake  (Chisholm Lake) (near Birch Lake)", LakeName)) %>%
  mutate(LakeName = ifelse(LakeName == "North Twin Lake", "North Twin Lake (Meadows Road)", LakeName)) %>%
  mutate(LakeName = ifelse(LakeName == "Paul's Pond", "Pauls Pond (Coal Mine Road)", LakeName)) %>%
  mutate(LakeName = ifelse(LakeName == "Rapids Lake", "Rapids Lake (Richardson Highway)", LakeName)) %>%
  mutate(LakeName = ifelse(LakeName == "Summit Lake (Paxson)", "Summit Lake (Richardson Highway near Paxson)", LakeName)) %>%
  left_join(morphometry1) %>%
  filter(make_estimates | use_fish)



# table(morphometry1$make_estimates, morphometry1$LakeName %in% laketrout_all$LakeName, morphometry1$use_fish)
morphometry <- morphometry1 %>%
  filter(make_estimates | use_fish)
morphometry$LakeNum <- 1:nrow(morphometry)
nrow(morphometry)  # 58


# filling in numerical lake identifier in fish-level data, consistent with lake-level data
laketrout_all$LakeNum <- NA
for(i in 1:nrow(morphometry)) {
  laketrout_all$LakeNum[laketrout_all$LakeName == morphometry$LakeName[i]] <- i
}
table(laketrout_all$LakeNum)

nrow(laketrout_all)  # 35691, was 35516
sapply(laketrout_all, function(x) sum(is.na(x)))

summary(laketrout_all$Year)  # 1960-2024
length(unique(laketrout_all$ProjectTitle))  # 118, was 148
length(unique(laketrout_all$LakeName))  # 41, was 84


# making sure that all lake names in laketrout_all are contained in morphometry
sort(unique(laketrout_all$LakeName))[!(sort(unique(laketrout_all$LakeName)) %in% morphometry$LakeName)]
morphometry$LakeName[!is.na(morphometry$LakeName)]

## Filtering informed by Weight ~ Length relationship
## - one problematic project
## - outlying residuals from a log(Weight) ~ log(Length) regression

laketrout1 <- laketrout_all %>%
  mutate(Weight_g = ifelse(Weight_g < 50, Weight_g*1000, Weight_g)) %>%
  filter(is.na(Weight_g) | Weight_g < 100000) %>%
  # filter(is.na(ForkLength_mm) | ForkLength_mm > 150) %>%
  # filter(is.na(Age) | Age < 50) %>%
  filter(!ProjectTitle %in% c("Mark-Recapture Event 1 - (September - 2002)",
                              "Mark-Recapture Event 1 - (September - 2003)",
                              "Mark-Recapture Event 1 - (September - 2004)",
                              "Mark-Recapture Event 2 - (May - 2003)"))
# filter(is.na(ForkLength_mm) | ForkLength_mm )
nrow(laketrout1) # 32714, was 32539

lm1 <- with(laketrout1, lm(log(Weight_g) ~ log(ForkLength_mm)))
resids1 <- log(laketrout1$Weight_g) - predict(lm1, newdata=laketrout1)

# how many residual SDs to include
sdcutoff <- 5
laketrout2 <- filter(laketrout1, is.na(resids1) | abs(resids1) < sdcutoff*sd(resids1, na.rm=TRUE))
with(laketrout1, plot(log(Weight_g) ~ log(ForkLength_mm),
                      pch=ifelse(is.na(resids1) | abs(resids1) < sdcutoff*sd(resids1, na.rm=TRUE), 1, 16)))

laketrout <- laketrout2 %>%
  filter(is.na(Age) | Age < 50) # %>%
  # filter(`Include in "Alaskanizing" Modeling Exercise` %in% c("yes", "Yes"))# %>%
  # mutate(LakeNum = as.numeric(as.factor(LakeName)))

nrow(laketrout) # 32697 # was 32516

# lakenames <- levels(as.factor(laketrout$LakeName))
length(levels(as.factor(laketrout$LakeName))) # 41 was 84

# sort(unique(laketrout$LakeName[laketrout$use_fish]))
# sort(unique(laketrout$LakeName[laketrout$make_estimates]))
# morphometry$LakeName[morphometry$`Include in "Alaskanizing" Modeling Exercise` %in% c("Yes","yes")]
# morphometry$LakeName[morphometry$`Potentially Include in Lake Trout Management Plan` != "No"]



### plotting weight ~ length
laketrout_all %>%
  # filter(ForkLength_mm > 100) %>%
  # filter(LakeName == "Paxson Lake") %>%
  ggplot(aes(x=ForkLength_mm, y=Weight_g, color=LakeName)) +
  # ggplot(aes(x=ForkLength_mm, y=Weight_g)) +
  # facet_wrap(~LakeName) +
  # facet_wrap(~ProjectTitle) +
  geom_point() +
  scale_y_log10() +
  scale_x_log10() +
  theme_bw() +
  theme(legend.position = 'none')

laketrout %>%
  ggplot(aes(x=ForkLength_mm, y=Weight_g, color=LakeName)) +
  # ggplot(aes(x=ForkLength_mm, y=Weight_g)) +
  geom_point() +
  scale_y_log10() +
  scale_x_log10(limits=c(100,1200)) +
  theme_bw() +
  theme(legend.position = 'none')



## Filtering informed by Length ~ Age relationship

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

# plot(laketrout_all$ForkLength_mm)
# plot(laketrout$ForkLength_mm)
# plot(laketrout_all$Weight_g)
# plot(laketrout$Weight_g)

sapply(laketrout, function(x) sum(is.na(x)))
# table(is.na(laketrout$Latitude_WGS84),
#       is.na(laketrout$SurfaceArea_h),
#       is.na(laketrout$Weight_g))

nrow(laketrout)

sum(!is.na(laketrout$ForkLength_mm) & !is.na(laketrout$Weight_g))  # 4607 paired weight~length obs
length(unique(laketrout$LakeNum[!is.na(laketrout$ForkLength_mm) & !is.na(laketrout$Weight_g)]))  # 26 lakes

sum(!is.na(laketrout$ForkLength_mm) & !is.na(laketrout$Age))  # 1517 paired length~age obs
length(unique(laketrout$LakeNum[!is.na(laketrout$ForkLength_mm) & !is.na(laketrout$Age)]))  # 19 lakes

# length(unique(laketrout$LakeNum[!is.na(laketrout$Latitude_WGS84) & !is.na(laketrout$SurfaceArea_h)]))

sapply(laketrout, function(x) sum(!is.na(x)))




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
                             Temperature = morphometry$`Temp (C)`)
lakenames <- morphometry$LakeName



cindy_jags <- tempfile()
cat('model {

  ## ---- Length ~ Age Component ---- ##
  for(i in whichdata_LAge) {   # fish-level data from lakes with paired lengths & ages
    # L[i] ~ dnorm(mu_LAge[i], tau_LAge)
    # Lpp[i] ~ dnorm(mu_LAge[i], tau_LAge)
    # mu_LAge[i] <- Linf[lake[i]]*(1-exp(-k[lake[i]]*(Age[i]-t0[lake[i]])))
    L[i] ~ dlnorm(logmu_LAge[i], tau_LAge)
    Lpp[i] ~ dlnorm(logmu_LAge[i], tau_LAge)
    logmu_LAge[i] <- log(Linf[lake[i]]*(1-exp(-k[lake[i]]*(Age[i]-t0[lake[i]]))))
  }

  tau_LAge <- pow(sig_LAge, -2)
  # sig_LAge ~ dunif(0, 100)
  # sig_LAge_prior ~ dunif(0, 100)
  sig_LAge ~ dunif(0, 3)
  sig_LAge_prior ~ dunif(0, 3)


  for(j in whichlakes_LAge) {   # lake level
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
  for(j in whichlakes_LArea) {  # lake-level data where there are areas
    # Linf[j] ~ dnorm(mu_LArea[j], tau_LArea)T(1,)
    # mu_LArea[j] <- b0_LArea * (1 - exp(-b1_LArea * (1 + log(Area[j]))))
    Linf[j] ~ dlnorm(logmu_LArea[j], tau_LArea)
    logmu_LArea[j] <- log(b0_LArea * (1 - exp(-b1_LArea * (1 + log(Area[j])))))
  }
  for(j in whichlakes_LArea_c) {  # lake-level priors where there are no areas
    Linf[j] ~ dnorm(600, 0.00001)T(1,)
    # Linf[j] ~ dnorm(mu_Linf_noarea, tau_Linf_noarea)
  }
  # mu_Linf_noarea ~ dnorm(600, 0.0001)
  # tau_Linf_noarea <- pow(sig_Linf_noarea, -2)
  # sig_Linf_noarea ~ dunif(0, 500)
  # sig_Linf_noarea_prior ~ dunif(0, 500)

  tau_LArea <- pow(sig_LArea, -2)
  # sig_LArea ~ dunif(0, 300)
  # sig_LArea_prior ~ dunif(0, 300)
  sig_LArea ~ dunif(0, 3)
  sig_LArea_prior ~ dunif(0, 3)
  # sig_LArea ~ dexp(0.05)
  # sig_LArea_prior ~ dexp(0.05)

  b0_LArea ~ dnorm(b0_LArea_lester, pow(cv_b0_LArea_lester*b0_LArea_lester, -2))T(0.01,)
  b0_LArea_prior ~ dnorm(b0_LArea_lester, pow(cv_b0_LArea_lester*b0_LArea_lester, -2))T(0.01,)
  b1_LArea ~ dnorm(b1_LArea_lester, pow(cv_b1_LArea_lester*b1_LArea_lester, -2))T(0.01,)
  b1_LArea_prior ~ dnorm(b1_LArea_lester, pow(cv_b1_LArea_lester*b1_LArea_lester, -2))T(0.01,)



  ## ---- Length quantile ~ Linf component ---- ##
  for(j in whichlakes_L) {
    qL[j] ~ dnorm(Linf[j], tau_Lj[j])
    qLpp[j] ~ dnorm(Linf[j], tau_Lj[j])
    tau_Lj[j] <- pow(sig_Lj[j], -2)
    sig_Lj[j] <- pow((sig_L0^2) + ((sig_L1^2)/nL[j]), 0.5)
  }

  sig_L0 ~ dunif(0, sig_L0_cap)
  sig_L0_prior ~ dunif(0, sig_L0_cap)
  sig_L1 ~ dunif(0, sig_L1_cap)
  sig_L1_prior ~ dunif(0, sig_L1_cap)

  # sig_L0 ~ dexp(0.1)
  # sig_L0_prior ~ dexp(0.1)
  # sig_L1 ~ dexp(0.1)
  # sig_L1_prior ~ dexp(0.1)



  ## ---- log Weight ~ log Length component ---- ##
  # loop over all WL data
  for(i in whichdata_WL) {
    logW[i] ~ dnorm(mu_logW[i], tau_logW)
    # logWpp[i] ~ dnorm(mu_logW[i], tau_logW)    # might take this out if model seems ok
    mu_logW[i] <- b0_logW[lake[i]] + b1_logW[lake[i]]*logLc[i]
  }

  # loop over lakes with lat & area data
  for(j in whichlakes_WL) {
    b0_logW[j] ~ dnorm(mu_b0_logW[j], tau_b0_logW)
    b0_logW_interp[j] <- b0_logW[j] - b1_logW[j]*meanlogLc
    mu_b0_logW[j] <- b0_logW_int
                + b0_logW_area*logareac[j]

    b1_logW[j] ~ dnorm(mu_b1_logW[j], tau_b1_logW)
    mu_b1_logW[j] <- b1_logW_int
                + b1_logW_lat*latc[j]
  }

  # loop over lakes without area & lat data (c for complement)
  # draw log(area) and lat from Normal distributions (mean & sd of actual data)
  # corresponding b0_logW, b0_logW_interp, and b1_logW will effectively be from the post pred
  for(jc in whichlakes_WL_c) {
    b0_logW[jc] ~ dnorm(mu_b0_logW[jc], tau_b0_logW)
    b0_logW_interp[jc] <- b0_logW[jc] - b1_logW[jc]*meanlogLc
    mu_b0_logW[jc] <- b0_logW_int
                + b0_logW_area*areasim[jc]

    b1_logW[jc] ~ dnorm(mu_b1_logW[jc], tau_b1_logW)
    mu_b1_logW[jc] <- b1_logW_int
                + b1_logW_lat*latsim[jc]
    areasim[jc] ~ dnorm(mean(logareac[whichlakes_WL]), pow(sd(logareac[whichlakes_WL]), -2))
    latsim[jc] ~ dnorm(mean(latc[whichlakes_WL]), pow(sd(latc[whichlakes_WL]), -2))
  }

  # global priors
  sig_b0_logW ~ dunif(0, 10)
  sig_b0_logW_prior ~ dunif(0, 10)
  tau_b0_logW <- pow(sig_b0_logW, -2)

  sig_b1_logW ~ dunif(0, 10)
  sig_b1_logW_prior ~ dunif(0, 10)
  tau_b1_logW <- pow(sig_b1_logW, -2)

  b0_logW_int ~ dnorm(0, 0.001)
  b0_logW_area ~ dnorm(0, 0.001)
  b0_logW_int_prior ~ dnorm(0, 0.001)
  b0_logW_area_prior ~ dnorm(0, 0.001)

  b1_logW_int ~ dnorm(0, 0.001)
  b1_logW_lat ~ dnorm(0, 0.001)
  b1_logW_int_prior ~ dnorm(0, 0.001)
  b1_logW_lat_prior ~ dnorm(0, 0.001)

  tau_logW <- pow(sig_logW, -2)
  sig_logW ~ dunif(0, 10)
  sig_logW_prior ~ dunif(0, 10)



  ## ---- Weight quantile ~ Winf component ---- ##
  for(j in whichlakes_W) {
    qW[j] ~ dnorm(Winf[j], tau_Wj[j])
    qWpp[j] ~ dnorm(Winf[j], tau_Wj[j])
    tau_Wj[j] <- pow(sig_Wj[j], -2)
    sig_Wj[j] <- pow((sig_W0^2) + ((sig_W1^2)/nW[j]), 0.5)
  }

  sig_W0 ~ dunif(0, sig_W0_cap)
  sig_W0_prior ~ dunif(0, sig_W0_cap)
  sig_W1 ~ dunif(0, sig_W1_cap)
  sig_W1_prior ~ dunif(0, sig_W1_cap)

  # sig_W0 ~ dexp(0.1)
  # sig_W0_prior ~ dexp(0.1)
  # sig_W1 ~ dexp(0.1)
  # sig_W1_prior ~ dexp(0.1)



  ## ---- Winf from Linf! ---- ##
  for(j in alllakes) {
    Winf[j] <- exp(b0_logW_interp[j] + b1_logW[j]*log(Linf[j]))
  }


}', file=cindy_jags)



# bundle data to pass into JAGS
q_input <- 0.95   # assumed quantile value for asymptotic size
cindy_data <- list(Age = laketrout$Age,
                   L = laketrout$ForkLength_mm,
                   lake = as.numeric(as.character(laketrout$LakeNum)),
                   Area = laketrout_Winf$Area_ha,
                   whichlakes_LAge = which(morphometry$use_fish &
                                             laketrout_Winf$n_Age > 0 & laketrout_Winf$n_Length > 0),
                   whichlakes_LArea = which(!is.na(laketrout_Winf$Area_ha)),# & laketrout_Winf$n_Length > 0),
                   whichlakes_LArea_c = which(!(!is.na(laketrout_Winf$Area_ha))),# & laketrout_Winf$n_Length > 0)),
                   alllakes = 1:nrow(laketrout_Winf),
                   Agefit = 1:max(laketrout$Age, na.rm=TRUE),
                   nfit = max(laketrout$Age, na.rm=TRUE),
                   b0_LArea_lester = 957,
                   cv_b0_LArea_lester = 0.2*2.5,
                   b1_LArea_lester = 0.14,
                   cv_b1_LArea_lester = 1*2.5,
                   qL = tapply(laketrout$ForkLength_mm, laketrout$LakeNum,
                               quantile, p=q_input, na.rm=TRUE),
                   nL = laketrout_Winf$n_Length,
                   whichlakes_L = which(morphometry$use_fish & laketrout_Winf$n_Length > 10),
                   sig_L0_cap = 500,
                   sig_L1_cap = 2000,

                   logW = log(laketrout$Weight_g/1000),
                   logLc = log(laketrout$ForkLength_mm) - mean(log(laketrout$ForkLength_mm), na.rm=TRUE),
                   meanlogLc = mean(log(laketrout$ForkLength_mm), na.rm=TRUE),
                   whichdata_WL = which(laketrout$use_fish &
                                          !is.na(laketrout$ForkLength_mm) &
                                       !is.na(laketrout$SurfaceArea_h) &
                                       !is.na(laketrout$Latitude_WGS84)),
                   whichlakes_WL = which(!is.na(laketrout_Winf$Lat) & !is.na(laketrout_Winf$Area_ha)),
                   whichlakes_WL_c = which(is.na(laketrout_Winf$Lat) | is.na(laketrout_Winf$Area_ha)),
                   latc = morphometry$Latitude_WGS84 - mean(morphometry$Latitude_WGS84, na.rm=TRUE),
                   logareac = log(morphometry$SurfaceArea_h) - mean(log(morphometry$SurfaceArea_h), na.rm=TRUE),
                   qW = tapply(laketrout$Weight_g/1000, laketrout$LakeNum,
                               quantile, p=q_input, na.rm=TRUE),
                   nW = laketrout_Winf$n_Weight,
                   whichlakes_W = which(morphometry$use_fish & laketrout_Winf$n_Weight > 10),
                   sig_W0_cap = 5,
                   sig_W1_cap = 20
                   )

cindy_data$whichdata_LAge <- which((laketrout$LakeNum %in% cindy_data$whichlakes_LAge) &
                                     (!is.na(laketrout$Age)) & (!is.na(laketrout$ForkLength_mm)))
# cindy_data$max_nL <- max(cindy_data$which_nL, na.rm=TRUE)
# cindy_data$max_nW <- max(cindy_data$which_nW, na.rm=TRUE)

# cindy_data$whichlakes_L <- cindy_data$whichlakes_L[cindy_data$whichlakes_L != 5]





# JAGS controls
# niter <- 2*1000
# niter <- 20*1000
# niter <- 50*1000      # 50k in 9 minutes
niter <- 100*1000
# ncores <- 3
ncores <- min(10, parallel::detectCores()-1)


{
  tstart <- Sys.time()
  print(tstart)
  parameters <- c("sig_LAge", "sig_LAge_prior",
                  "t0", "t0_prior",
                  "k", "k_prior",
                  "mu_t0","mu_t0_prior",
                  "sig_t0","sig_t0_prior",
                  "mu_k","mu_k_prior",
                  "sig_k","sig_k_prior",
                  "Lfit",
                  "Linf", "Winf",
                  "mu_LArea",
                  "sig_LArea", "sig_LArea_prior",
                  "b0_LArea", "b0_LArea_prior",
                  "b1_LArea", "b1_LArea_prior",
                  # "sig_L", "sig_L_prior",
                  "sig_L0", "sig_L0_prior",
                  "sig_L1", "sig_L1_prior",
                  "sig_Lj",
                  "Lpp","qLpp",
                  "mu_Linf_noarea", "sig_Linf_noarea", "sig_Linf_noarea_prior",
                  "b0_logW","b1_logW","sig_logW","sig_logW_prior",
                  "mu_b0_logW","mu_b1_logW",
                  "sig_b0_logW","sig_b0_logW_prior",
                  "sig_b1_logW","sig_b1_logW_prior",
                  "b0_logW_interp",
                  "b0_logW_int","b1_logW_int","b0_logW_int_prior","b1_logW_int_prior",
                  "b0_logW_lat","b0_logW_area","b0_logW_lat_prior","b0_logW_area_prior",
                  "b1_logW_lat","b1_logW_area","b1_logW_lat_prior","b1_logW_area_prior",
                  "sig_W0", "sig_W0_prior",
                  "sig_W1", "sig_W1_prior",
                  "sig_Wj",
                  "Wpp","qWpp")
  cindy_jags_out <- jagsUI::jags(model.file=cindy_jags, data=cindy_data,
                                 parameters.to.save=parameters,
                                 n.chains=ncores, parallel=T, n.iter=niter,
                                 n.burnin=niter/2, n.thin=niter/2000)
  print(Sys.time() - tstart)
}
cindy_jags_out$DIC

# winf_list <- list()
# i_list <- 1
i_list <- i_list+1
winf_list[[i_list]] <- cindy_jags_out$sims.list$Winf
par(mfrow=c(1,1))
comparecat(winf_list)

#### global model diagnostics and summary

# nbyname(cindy_jags_out)
par(mfrow=c(1,1))
plotRhats(cindy_jags_out)
traceworstRhat(cindy_jags_out, parmfrow=c(3,3))

comparepriors(cindy_jags_out, parmfrow=c(3,3))

par(mfrow=c(1,1))
plotcor_jags(cindy_jags_out, p=c("Linf", "mu_LArea", "t0", "k", "b0", "b1", "sig"))
# crossplot(cindy_jags_out, p=c("sig_L0","sig_L1"), drawblob=TRUE)




### Posterior predictive for qL
par(mfrow=c(1,2))
qq_postpred(cindy_jags_out$sims.list$qLpp[,cindy_data$whichlakes_L],
            y=cindy_data$qL[cindy_data$whichlakes_L],
            main="Post predictive for qL")
ts_postpred(cindy_jags_out$sims.list$qLpp[,cindy_data$whichlakes_L],
            y=cindy_data$qL[cindy_data$whichlakes_L],
            x=cindy_data$nL[cindy_data$whichlakes_L], log="x",
            xlab="n Length for quantiles")
par(mfcol=c(3,3))
plot_postpred(cindy_jags_out$sims.list$qLpp[,cindy_data$whichlakes_L],
              y=cindy_data$qL[cindy_data$whichlakes_L],
              x=log(cindy_data$nL[cindy_data$whichlakes_L]))


### plot how sig_Lj varies by lake
par(mfrow=c(1,2))
estcols <- ifelse(morphometry$make_estimates, 4, 3)
caterpillar(cindy_jags_out, "sig_Lj", col=estcols)
envelope(cindy_jags_out$sims.list$sig_Lj[, cindy_data$whichlakes_L],
         x=cindy_data$nL[cindy_data$whichlakes_L], log="x",
         xlab="n Length for quantiles")



### Posterior predictive for qW
par(mfrow=c(1,2))
qq_postpred(cindy_jags_out$sims.list$qWpp[,cindy_data$whichlakes_W],
            y=cindy_data$qW[cindy_data$whichlakes_W],
            main="Post predictive for qW")
ts_postpred(cindy_jags_out$sims.list$qWpp[,cindy_data$whichlakes_W],
            y=cindy_data$qW[cindy_data$whichlakes_W],
            x=cindy_data$nW[cindy_data$whichlakes_W], log="x",
            xlab="n Wength for quantiles")
par(mfcol=c(3,3))
plot_postpred(cindy_jags_out$sims.list$qWpp[,cindy_data$whichlakes_W],
              y=cindy_data$qW[cindy_data$whichlakes_W],
              x=log(cindy_data$nW[cindy_data$whichlakes_W]))


### plot how sig_Lj varies by lake
par(mfrow=c(1,2))
caterpillar(cindy_jags_out, "sig_Wj", col=estcols)
envelope(cindy_jags_out$sims.list$sig_Wj[, cindy_data$whichlakes_W],
         x=cindy_data$nW[cindy_data$whichlakes_W], log="x",
         xlab="n Length for quantiles")


### Posterior predictive for L
# par(mfrow=c(1,2))
# qq_postpred(cindy_jags_out, p="Lpp", y=cindy_data$L)
# ts_postpred(cindy_jags_out, p="Lpp", y=cindy_data$L)
# par(mfcol=c(3,3))
# plot_postpred(cindy_jags_out, p="Lpp", y=cindy_data$L, x=cindy_data$Age)



### LVB parameters
par(mfrow=c(1,2))
caterpillar(cindy_jags_out, "t0",
            col=estcols)
caterpillar(cindy_jags_out, "k",
            col=estcols)

## length ~ age (plots for all lakes where there are lengths and ages)
par(mfrow=c(3,3))
for(j in cindy_data$whichlakes_LAge) {
  envelope(cindy_jags_out$sims.list$Lfit[,,j], main=lakenames[j],
           ylim=range(cindy_data$L[!is.na(cindy_data$Age)]),
           col=ifelse(j %in% cindy_data$whichlakes_L, 2, 4))
  points(cindy_data$Age[cindy_data$lake==j], cindy_data$L[cindy_data$lake==j])
  legend("topleft", legend=c("L ~ Age only", "qL also"),
         fill=adjustcolor(c(4,2), alpha.f=.3), border=c(4,2), cex=.6)
}




### Linf by data availability
par(mfrow=c(2,2))
caterpillar(cindy_jags_out, "Linf", col=3+cindy_data$alllakes %in% cindy_data$whichlakes_LAge)
legend("topleft", legend=c("has AGES", "no AGES"), col=4:3, lwd=3, cex=0.5)
caterpillar(cindy_jags_out, "Linf", col=3+cindy_data$alllakes %in% cindy_data$whichlakes_LArea)
legend("topleft", legend=c("has AREA", "no AREA"), col=4:3, lwd=3, cex=0.5)
caterpillar(cindy_jags_out, "Linf", col=3+cindy_data$alllakes %in% cindy_data$whichlakes_L)
legend("topleft", legend=c("has Lq", "no Lq"), col=4:3, lwd=3, cex=0.5)
caterpillar(cindy_jags_out, "Linf", col=3+cindy_data$alllakes %in% cindy_data$whichlakes_W)
legend("topleft", legend=c("has Wq", "no Wq"), col=4:3, lwd=3, cex=0.5)

# # ages and lq
# caterpillar(cindy_jags_out, "Linf", col=3+(cindy_data$alllakes %in%
#               intersect(cindy_data$whichlakes_L, cindy_data$whichlakes_LAge)))
# legend("topleft", legend=c("has Lq and AGES", "no"), col=4:3, lwd=3, cex=0.5)
#
# # ages and lq and area
# caterpillar(cindy_jags_out, "Linf", col=3+(cindy_data$alllakes %in%
#                                              intersect(cindy_data$whichlakes_L,
#                                                        intersect(cindy_data$whichlakes_LAge,
#                                                        cindy_data$whichlakes_LArea))))
# legend("topleft", legend=c("has Lq and AGES and AREA", "no"), col=4:3, lwd=3, cex=0.5)

# ages, then lq, then area
par(mfrow=c(1,1))
theorder <- ifelse(cindy_data$alllakes %in% cindy_data$whichlakes_LAge, 1,
                   ifelse(cindy_data$alllakes %in% cindy_data$whichlakes_L, 2,
                          ifelse(cindy_data$alllakes %in% cindy_data$whichlakes_LArea, 3, 4))) %>% rank(ties.method = "first")
caterpillar(cindy_jags_out, "Linf", x=theorder,
            col=ifelse(cindy_data$alllakes %in% cindy_data$whichlakes_LAge, 4,
                                               ifelse(cindy_data$alllakes %in% cindy_data$whichlakes_L, 2,
                                                      ifelse(cindy_data$alllakes %in% cindy_data$whichlakes_LArea, 3, 1))))
legend("topleft", legend=c("AGES, then","Lq, then","AREA, then","nothing"), col=c(4,2,3,1), lwd=3, cex=0.5)
points(x=theorder[cindy_data$whichlakes_W], y=cindy_jags_out$q50$Linf[cindy_data$whichlakes_W])#, col=cindy_data$alllakes %in% cindy_data$whichlakes_W)
points(x=theorder[morphometry$make_estimates], y=cindy_jags_out$q50$Linf[morphometry$make_estimates], pch=16)

plot(NA, xlim=range(cindy_data$qL, na.rm=T),
     ylim=range(cindy_jags_out$q2.5$Linf, cindy_jags_out$q97.5$Linf, na.rm=TRUE),
     main="Linf", xlab="qL", ylab="")
caterpillar(cindy_jags_out, "Linf", col=ifelse(cindy_data$alllakes %in% cindy_data$whichlakes_LAge, 4,
                                               ifelse(cindy_data$alllakes %in% cindy_data$whichlakes_L, 2,
                                                      ifelse(cindy_data$alllakes %in% cindy_data$whichlakes_LArea, 3, 1))),
            x=cindy_data$qL,
            add=TRUE)
legend("topleft", legend=c("AGES, then","Lq, then","AREA, then","nothing"), col=c(4,2,3,1), lwd=3, cex=0.5)
abline(0,1, lty=2)
# points(x=cindy_data$qL[cindy_data$whichlakes_LArea],
#        y=cindy_jags_out$q50$mu_LArea[cindy_data$whichlakes_LArea])

par(mfrow=c(1,1))
plot(NA,
     ylim=range(cindy_jags_out$q2.5$Linf, cindy_jags_out$q97.5$Linf, na.rm=TRUE),
     xlim=range(log(cindy_data$Area), -1, na.rm=TRUE),
     main="Linf", xlab="Area (ha)", ylab="", xaxt="n")
lbls <- c(1, 5, 10, 50, 100, 500, 1000, 5000)
axis(side=1, at=log(lbls), labels=lbls)
thecols <- ifelse(cindy_data$alllakes %in% cindy_data$whichlakes_LAge, 4,
                  ifelse(cindy_data$alllakes %in% cindy_data$whichlakes_L, 2,
                         ifelse(cindy_data$alllakes %in% cindy_data$whichlakes_LArea, 3, 1)))
caterpillar(cindy_jags_out, "Linf", col=thecols,
            x=ifelse(is.na(cindy_data$Area), -1, log(cindy_data$Area)),
            add=TRUE)
points(x=ifelse(is.na(cindy_data$Area), -1, log(cindy_data$Area)),
       y=cindy_jags_out$q50$Linf,
       col=thecols, pch=16)
curve(cindy_data$b0_LArea_lester * (1 - exp(-cindy_data$b1_LArea_lester * (1 + x))), add=TRUE, lty=2)
legend("topleft", legend=c("AGES, then","Lq, then","AREA, then","nothing","Lester priors"),
       col=c(4,2,3,1,1), lwd=c(rep(3,4),1), lty=c(rep(1,4),2), cex=0.5)

points(x=(ifelse(is.na(cindy_data$Area), -1, log(cindy_data$Area)))[cindy_data$whichlakes_L],
       y=cindy_data$qL[cindy_data$whichlakes_L])
points(x=(ifelse(is.na(cindy_data$Area), -1, log(cindy_data$Area)))[morphometry$make_estimates],
       y=cindy_jags_out$q50$Linf[morphometry$make_estimates], pch=16)

# more of the same, maybe these thoughts can be combined?
# caterpillar(cindy_jags_out$sims.list$Linf[,cindy_data$whichlakes_L],
#             x=cindy_data$qL[cindy_data$whichlakes_L],
#             xlab="qL", main="subset with qL")
# abline(0,1,lty=3)
# caterpillar(cindy_jags_out$sims.list$Linf[,cindy_data$whichlakes_LAge],
#             x=cindy_data$qL[cindy_data$whichlakes_LAge],
#             xlab="qL", main="subset with Ages")
# abline(0,1,lty=3)
# int <- intersect(cindy_data$whichlakes_LAge, cindy_data$whichlakes_L)
# caterpillar(cindy_jags_out$sims.list$Linf[,int],
#             x=cindy_data$qL[int],
#             xlab="qL", main="subset with qL and Ages")
# abline(0,1,lty=3)


## logW ~ logL for each lake
## b0_logW ~ lat & area
## b1_logW ~ lat & area
## mu_b0_logW ~ lat & area
## mu_b1_logW ~ lat & area


## colors are differentiated by whether lat/area data and weight data are present
LatArea_present <- sort(unique(laketrout$LakeNum)) %in% cindy_data$whichlakes_WL
weight_present <- table(laketrout$LakeNum,!is.na(laketrout$Weight_g))[,2] > 0
cols <- ifelse(weight_present, 3, ifelse(LatArea_present, 4, 2))
par(mfrow=c(3,1))
caterpillar(cindy_jags_out, p="b1_logW", col=cols, ylim=c(2.3, 4.4))
legend("topright", legend=c("Length & Wt", "Lat & Area only", "no data"),
       lwd=3, col=c(3, 4, 2))
abline(h=3.2, lty=3)
caterpillar(cindy_jags_out, p="b0_logW", col=cols, ylim=c(-0.3, 0.7))
legend("topright", legend=c("Length & Wt", "Lat & Area only", "no data"),
       lwd=3, col=c(3, 4, 2))
caterpillar(cindy_jags_out, p="b0_logW_interp", col=cols, ylim=c(-25, -13))
legend("topright", legend=c("Length & Wt", "Lat & Area only", "no data"),
       lwd=3, col=c(3, 4, 2))
abline(h=-19.56, lty=3)

# plot regression bands for each lake, overlay data and lester line
xpredict <- seq(from = min(log(laketrout$ForkLength_mm[!is.na(laketrout$Weight_g)]), na.rm=TRUE),
                to = max(log(laketrout$ForkLength_mm[!is.na(laketrout$Weight_g)]), na.rm=TRUE),
                length.out=50)
par(mfrow=c(3,3))
for(ilake in 1:nrow(laketrout_Winf)) {   ### I don't think colors are right
  logweight_predict <- cindy_jags_out$sims.list$b0_logW_interp[,ilake] +
    outer(cindy_jags_out$sims.list$b1_logW[,ilake], xpredict)
  plot(NA, xlim=range(xpredict), ylim=range(log(laketrout$Weight_g/1000), na.rm=TRUE),
       xlab="log(Length)", ylab="log(Weight)", main=lakenames[ilake])
  if(all(!is.na(logweight_predict))) envelope(logweight_predict, x=xpredict, add=TRUE, col=cols[ilake])
  points(x = log(laketrout$ForkLength_mm[laketrout$LakeNum==ilake]),
         y = log(laketrout$Weight_g[laketrout$LakeNum==ilake]/1000))
  abline(a=-19.56, b=3.2, lty=3)
}

par(mfrow=c(2,2))
caterpillar(cindy_jags_out$sims.list$b0_logW[, weight_present & LatArea_present],
            x=cindy_data$latc[weight_present & LatArea_present],
            xlab="Lat c", main="b0")
caterpillar(cindy_jags_out$sims.list$b1_logW[, weight_present & LatArea_present],
            x=cindy_data$latc[weight_present & LatArea_present],
            xlab="Lat c", main="b1")
envelope(cindy_jags_out$sims.list$mu_b1_logW[, weight_present & LatArea_present],
         x=cindy_data$latc[weight_present & LatArea_present],
         add=TRUE)
caterpillar(cindy_jags_out$sims.list$b0_logW[, weight_present & LatArea_present],
            x=cindy_data$logareac[weight_present & LatArea_present],
            xlab="Log Area c", main="b0")
envelope(cindy_jags_out$sims.list$mu_b0_logW[, weight_present & LatArea_present],
         x=cindy_data$logareac[weight_present & LatArea_present],
         add=TRUE)
caterpillar(cindy_jags_out$sims.list$b1_logW[, weight_present & LatArea_present],
            x=cindy_data$logareac[weight_present & LatArea_present],
            xlab="Log Area c", main="b1")



par(mfrow=c(1,1))
caterpillar(cindy_jags_out, "Winf", col=estcols)

# wq, then ages, then lq, then area
par(mfrow=c(1,1))
theorder <- ifelse(cindy_data$alllakes %in% cindy_data$whichlakes_W, 0,
  ifelse(cindy_data$alllakes %in% cindy_data$whichlakes_LAge, 1,
                   ifelse(cindy_data$alllakes %in% cindy_data$whichlakes_L, 2,
                          ifelse(cindy_data$alllakes %in% cindy_data$whichlakes_LArea, 3, 4)))) %>% rank(ties.method = "first")
caterpillar(cindy_jags_out, "Winf", x=theorder,
            col=ifelse(cindy_data$alllakes %in% cindy_data$whichlakes_W, 5,
                       ifelse(cindy_data$alllakes %in% cindy_data$whichlakes_LAge, 4,
                              ifelse(cindy_data$alllakes %in% cindy_data$whichlakes_L, 2,
                              ifelse(cindy_data$alllakes %in% cindy_data$whichlakes_LArea, 3, 1)))))
legend("topleft", legend=c("Wq, then","AGES, then","Lq, then","AREA, then","nothing"), col=c(5,4,2,3,1), lwd=3, cex=0.5)
points(x=theorder[morphometry$make_estimates],
       y=cindy_jags_out$q50$Winf[morphometry$make_estimates],
       pch=16)

#
# plot(NA, xlim=range(cindy_data$qL, na.rm=T),
#      ylim=range(cindy_jags_out$q2.5$Winf, cindy_jags_out$q97.5$Winf, na.rm=TRUE),
#      main="Winf", xlab="qL", ylab="")
# caterpillar(cindy_jags_out, "Winf", col=ifelse(cindy_data$alllakes %in% cindy_data$whichlakes_LAge, 4,
#                                                ifelse(cindy_data$alllakes %in% cindy_data$whichlakes_L, 2,
#                                                       ifelse(cindy_data$alllakes %in% cindy_data$whichlakes_LArea, 3, 1))),
#             x=cindy_data$qL,
#             add=TRUE)
# legend("topleft", legend=c("AGES, then","Lq, then","AREA, then","nothing"), col=c(4,2,3,1), lwd=3, cex=0.5)
# abline(0,1, lty=2)
# # points(x=cindy_data$qL[cindy_data$whichlakes_LArea],
# #        y=cindy_jags_out$q50$mu_LArea[cindy_data$whichlakes_LArea])
#
# par(mfrow=c(1,1))
# plot(NA,
#      ylim=range(cindy_jags_out$q2.5$Winf, cindy_jags_out$q97.5$Winf, na.rm=TRUE),
#      xlim=range(log(cindy_data$Area), -1, na.rm=TRUE),
#      main="Linf", xlab="Area (ha)", ylab="", xaxt="n")
# lbls <- c(1, 5, 10, 50, 100, 500, 1000, 5000)
# axis(side=1, at=log(lbls), labels=lbls)
# thecols <- ifelse(cindy_data$alllakes %in% cindy_data$whichlakes_LAge, 4,
#                   ifelse(cindy_data$alllakes %in% cindy_data$whichlakes_L, 2,
#                          ifelse(cindy_data$alllakes %in% cindy_data$whichlakes_LArea, 3, 1)))
# caterpillar(cindy_jags_out, "Winf", col=thecols,
#             x=ifelse(is.na(cindy_data$Area), -1, log(cindy_data$Area)),
#             add=TRUE)
# points(x=ifelse(is.na(cindy_data$Area), -1, log(cindy_data$Area)),
#        y=cindy_jags_out$q50$Winf,
#        col=thecols, pch=16)
# curve(cindy_data$b0_LArea_lester * (1 - exp(-cindy_data$b1_LArea_lester * (1 + x))), add=TRUE, lty=2)
# legend("topleft", legend=c("AGES, then","Lq, then","AREA, then","nothing","Lester priors"),
#        col=c(4,2,3,1,1), lwd=c(rep(3,4),1), lty=c(rep(1,4),2), cex=0.5)
#
# points(x=(ifelse(is.na(cindy_data$Area), -1, log(cindy_data$Area)))[cindy_data$whichlakes_L],
#        y=cindy_data$qL[cindy_data$whichlakes_L])
#








# somehow compare sources of information for Winf
# compare to previous methods!!
load(file="Winf_all_comparison1.Rdata")
equivs <- c(100, 10, 1000, 10)
wts_vec <- 1/equivs
wts_nmat <- cbind(laketrout_Winf$n_Weight,
                  laketrout_Winf$n_Age,
                  laketrout_Winf$n_Length,
                  1*!is.na(laketrout_Winf$Area_ha))
wts <- wts_nmat*matrix(wts_vec, nrow=nrow(wts_nmat), ncol=ncol(wts_nmat), byrow=TRUE)
par(mfrow=c(3,4))
for(i in 1:length(Winf_method)) {
  cols <- c(rep(2,4),4,1)
  cols[Winf_method[i]] <- 3
  caterpillar(cbind(Winf_all_all[,i,], Winf_weighted[,i], cindy_jags_out$sims.list$Winf[,i]),
              col = cols,
              xax = c(paste(c(laketrout_Winf$n_Weight[i],
                              laketrout_Winf$n_Age[i],
                              laketrout_Winf$n_Length[i],
                              ""), round(wts[i,]/sum(wts[i,]), 2), sep=" / "),
                      "Wt","Int"),
              main = lakenames[i])
  abline(h=cindy_jags_out$q50$Winf[i], lty=2, lwd=1+(lakenames[i] %in% planlakes))
}



# still to do:
# - long runs with a vector of q_input, compare DICs
# - refine choice of priors & choice of error distributions

# to do when the new data comes...
# - investigate possible relationships:
#   . Elevation instead of Latitude in logW ~ logL
#   . t0 & k ~ elevation & area
#   . (Linf~area) ~ elevation





niter1 <- 1000*1000
ncores1 <- 8
q_inputs <- seq(from=0.9, to=0.99, by=0.01)   # length 10 x 50k with 5 cores took 1.4 hrs on laptop
                                              # length 10 x 500k with 8 cores took 20 hrs
length(q_inputs)
DICs <- mnL <- mnW <- NA*q_inputs
cindy_data1 <- cindy_data
for(i in seq_along(q_inputs)) {
  cindy_data1$qL = tapply(laketrout$ForkLength_mm, laketrout$LakeNum,
                          quantile, p=q_inputs[i], na.rm=TRUE)
  cindy_data1$qW = tapply(laketrout$Weight_g/1000, laketrout$LakeNum,
                          quantile, p=q_inputs[i], na.rm=TRUE)
  tstart <- Sys.time()
  print(tstart)
  cindy_jags_out1 <- jagsUI::jags(model.file=cindy_jags, data=cindy_data1,
                                 parameters.to.save=parameters,
                                 n.chains=ncores1, parallel=T, n.iter=niter1,
                                 n.burnin=niter1/2, n.thin=niter1/2000)
  print(Sys.time() - tstart)
  print(i)
  DICs[i] <- cindy_jags_out1$DIC
  mnL[i] <- mean(cindy_jags_out1$q50$qLpp > cindy_data1$qL[seq_along(cindy_jags_out1$q50$qLpp)], na.rm=TRUE)
  mnW[i] <- mean(cindy_jags_out1$q50$qWpp > cindy_data1$qW[seq_along(cindy_jags_out1$q50$qWpp)], na.rm=TRUE)
}
DICs - min(DICs)  # 22.07259 22.05110 10.91344  0.00000 19.77171 49.50245 68.75798 61.41945 34.86005 31.63211
mnL  # 0.6571429 0.6571429 0.6285714 0.6000000 0.5142857 0.4857143 0.4285714 0.4000000 0.3714286 0.2571429
mnW  # 0.7142857 0.7142857 0.7142857 0.6190476 0.6190476 0.6190476 0.5714286 0.5714286 0.3809524 0.3333333
par(mfrow=c(2,2))
plot(q_inputs, DICs - min(DICs))
plot(q_inputs, mnL, ylim=0:1)
abline(h=.5, lty=2)
plot(q_inputs, mnW, ylim=0:1)
abline(h=.5, lty=2)
