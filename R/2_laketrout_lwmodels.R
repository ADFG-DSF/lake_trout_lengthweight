### This script runs a set of candidate models describing the length-weight
### relationship for lake trout in several lakes.

### At the heart is a log-log regression, which is done using a Bayesian model.
### log(weight) = b0 + b1log(length) + residual error
### The intercept (b0) and slope (b1) parameters are allowed to vary between lakes.

### Lake-level parameters (slope and intercept) are allowed to vary according to
### the following:
### - common Normal distribution (i.e. b0[j] ~ N(mu_b0, sig_b0))
### - linear relationship with latitude, plus lake-level error
### - linear relationship with log area, plus lake-level error
### - linear relationship with latitude and log area, plus lake-level error

### A flexible model framework is defined below that uses data inputs to turn
### on/off model parameters and define the model type.  All models are then
### compared according to DIC score and resulting inferences.

### Matt Tyers, May 2024



# ## loading all data
source("R/1_laketrout_lwdata.R")



# JAGS controls
niter <- 50*1000#100*1000 #500*1000                # 50k in about an hour
ncores <- min(10, parallel::detectCores()-1)




# ## loading all data - consistent with 4_Winf methods!!
#
# # morphometry <- read_csv("flat_data/lake_morphometry.csv", skip=2)
# morphometry <- read_csv("flat_data/lake_morphometry2.csv", skip=2)
#
# # # is lake name unique?  YES
# # sum(!is.na(morphometry$LakeName))
# # length(unique(morphometry$LakeName))
#
# # laketrout_all <- read_csv("flat_data/length_weight.csv", skip=2) %>%
# #   left_join(morphometry)
# laketrout_all <- read_csv("flat_data/length_weight2.csv", skip=2) %>%
#   left_join(morphometry)
# nrow(laketrout_all)  # 35516
# sapply(laketrout_all, function(x) sum(is.na(x)))
#
# summary(laketrout_all$Year)  # 1960-2024
# length(unique(laketrout_all$ProjectTitle))  # 148
# length(unique(laketrout_all$LakeName))  # 84
#
#
# ## Filtering informed by Weight ~ Length relationship
# ## - one problematic project
# ## - outlying residuals from a log(Weight) ~ log(Length) regression
#
# laketrout1 <- laketrout_all %>%
#   mutate(Weight_g = ifelse(Weight_g < 50, Weight_g*1000, Weight_g)) %>%
#   filter(is.na(Weight_g) | Weight_g < 100000) %>%
#   # filter(is.na(ForkLength_mm) | ForkLength_mm > 150) %>%
#   # filter(is.na(Age) | Age < 50) %>%
#   filter(!ProjectTitle %in% c("Mark-Recapture Event 1 - (September - 2002)",
#                               "Mark-Recapture Event 1 - (September - 2003)",
#                               "Mark-Recapture Event 1 - (September - 2004)",
#                               "Mark-Recapture Event 2 - (May - 2003)"))
# # filter(is.na(ForkLength_mm) | ForkLength_mm )
# nrow(laketrout1) # 32539
#
# lm1 <- with(laketrout1, lm(log(Weight_g) ~ log(ForkLength_mm)))
# resids1 <- log(laketrout1$Weight_g) - predict(lm1, newdata=laketrout1)
# laketrout2 <- filter(laketrout1, is.na(resids1) | abs(resids1) < 4*sd(resids1, na.rm=TRUE))
#
# laketrout <- laketrout2 %>%
#   filter(is.na(Age) | Age < 50) %>%
#   mutate(LakeNum = as.numeric(as.factor(LakeName)))
#
# nrow(laketrout) # 32516
#
# ### plotting weight ~ length
# laketrout_all %>%
#   ggplot(aes(x=ForkLength_mm, y=Weight_g, color=LakeName)) +
#   # ggplot(aes(x=ForkLength_mm, y=Weight_g)) +
#   geom_point() +
#   scale_y_log10() +
#   scale_x_log10() +
#   theme_bw() +
#   theme(legend.position = 'none')
#
# laketrout %>%
#   ggplot(aes(x=ForkLength_mm, y=Weight_g, color=LakeName)) +
#   # ggplot(aes(x=ForkLength_mm, y=Weight_g)) +
#   geom_point() +
#   scale_y_log10() +
#   scale_x_log10(limits=c(100,1200)) +
#   theme_bw() +
#   theme(legend.position = 'none')
#
#
#
# ## Filtering informed by Length ~ Age relationship
#
# ## plotting length ~ age
# laketrout_all %>%
#   # ggplot(aes(y=ForkLength_mm, x=Age, color=LakeName)) +
#   ggplot(aes(y=ForkLength_mm, x=Age)) +
#   geom_point() +
#   theme_bw()
#
# laketrout %>%
#   # ggplot(aes(y=ForkLength_mm, x=Age, color=LakeName)) +
#   ggplot(aes(y=ForkLength_mm, x=Age)) +
#   geom_point() +
#   theme_bw()
#
# # plot(laketrout_all$ForkLength_mm)
# # plot(laketrout$ForkLength_mm)
# # plot(laketrout_all$Weight_g)
# # plot(laketrout$Weight_g)
#
# sapply(laketrout, function(x) sum(is.na(x)))
# table(is.na(laketrout$Latitude_WGS84),
#       is.na(laketrout$SurfaceArea_h),
#       is.na(laketrout$Weight_g))
#
# nrow(laketrout)
#
# sum(!is.na(laketrout$ForkLength_mm) & !is.na(laketrout$Weight_g))
# length(unique(laketrout$LakeNum[!is.na(laketrout$ForkLength_mm) & !is.na(laketrout$Weight_g)]))
#
# sum(!is.na(laketrout$ForkLength_mm) & !is.na(laketrout$Age))
# length(unique(laketrout$LakeNum[!is.na(laketrout$ForkLength_mm) & !is.na(laketrout$Age)]))
#
# length(unique(laketrout$LakeNum[!is.na(laketrout$Latitude_WGS84) & !is.na(laketrout$SurfaceArea_h)]))
#
# sapply(laketrout, function(x) sum(!is.na(x)))
#
#
#
# ## finally making laketrout_lw
# laketrout_lw <- filter(laketrout, )


# bundle data to pass into JAGS
# laketrout_lw <- filter(laketrout_lw, !is.na(SurfaceArea_h))

logarea <- log(with(laketrout_lw, tapply(SurfaceArea_h, LakeName, median, na.rm=T)))
lat <- with(laketrout_lw, tapply(Latitude_WGS84, LakeName, median))
lakenames <- names(lat)

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
## actually this is way easier!
controlmat <- expand.grid(0:1,0:1,0:1,0:1)

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



### Models were also constructed with common intercepts, slopes, and both
### but DIC scores were TERRIBLE!!  Thank you DIC for saving us from pseudoreplication.
### Not even going to consider these!!

# lt_jags_commonint <- tempfile()
# cat('model {
#   for(i in 1:n) {
#     y[i] ~ dnorm(mu[i], tau)
#     ypp[i] ~ dnorm(mu[i], tau)
#     mu[i] <- b0[lake[i]] + b1[lake[i]]*x[i]
#   }
#
#   for(j in 1:nlake) {
#     b0[j] <- b0_int
#     b0_interp[j] <- b0[j] - b1[j]*meanx
#     # b0[j] ~ dnorm(mu_b0[j], tau_b0)
#     # mu_b0[j] <- b0_int
#     #             + b0_lat*lat[j]   * intlat
#     #             + b0_area*area[j] * intarea
#
#     b1[j] ~ dnorm(mu_b1[j], tau_b1)
#     mu_b1[j] <- b1_int
#                 + b1_lat*lat[j]   * slopelat
#                 + b1_area*area[j] * slopearea
#   }
#
#   sig_b0 ~ dunif(0, 10)
#   tau_b0 <- pow(sig_b0, -2)
#
#   sig_b1 ~ dunif(0, 10)
#   tau_b1 <- pow(sig_b1, -2)
#
#   b0_int ~ dnorm(0, 0.001)
#   b0_lat ~ dnorm(0, 0.001)
#   b0_area ~ dnorm(0, 0.001)
#
#   b1_int ~ dnorm(0, 0.001)
#   b1_lat ~ dnorm(0, 0.001)
#   b1_area ~ dnorm(0, 0.001)
#
#   tau <- pow(sig, -2)
#   sig ~ dunif(0, 10)
#
# }', file=lt_jags_commonint)
#
#
# lt_jags_commonslope <- tempfile()
# cat('model {
#   for(i in 1:n) {
#     y[i] ~ dnorm(mu[i], tau)
#     ypp[i] ~ dnorm(mu[i], tau)
#     mu[i] <- b0[lake[i]] + b1[lake[i]]*x[i]
#   }
#
#   for(j in 1:nlake) {
#     b0[j] ~ dnorm(mu_b0[j], tau_b0)
#     b0_interp[j] <- b0[j] - b1[j]*meanx
#     mu_b0[j] <- b0_int
#                 + b0_lat*lat[j]   * intlat
#                 + b0_area*area[j] * intarea
#
#     b1[j] <- b1_int
#     # b1[j] ~ dnorm(mu_b1[j], tau_b1)
#     # mu_b1[j] <- b1_int
#     #             + b1_lat*lat[j]   * slopelat
#     #             + b1_area*area[j] * slopearea
#   }
#
#   sig_b0 ~ dunif(0, 10)
#   tau_b0 <- pow(sig_b0, -2)
#
#   sig_b1 ~ dunif(0, 10)
#   tau_b1 <- pow(sig_b1, -2)
#
#   b0_int ~ dnorm(0, 0.001)
#   b0_lat ~ dnorm(0, 0.001)
#   b0_area ~ dnorm(0, 0.001)
#
#   b1_int ~ dnorm(0, 0.001)
#   b1_lat ~ dnorm(0, 0.001)
#   b1_area ~ dnorm(0, 0.001)
#
#   tau <- pow(sig, -2)
#   sig ~ dunif(0, 10)
#
# }', file=lt_jags_commonslope)
#
#
# lt_jags_common <- tempfile()
# cat('model {
#   for(i in 1:n) {
#     y[i] ~ dnorm(mu[i], tau)
#     ypp[i] ~ dnorm(mu[i], tau)
#     mu[i] <- b0[lake[i]] + b1[lake[i]]*x[i]
#   }
#
#   for(j in 1:nlake) {
#     b0[j] <- b0_int
#     b0_interp[j] <- b0[j] - b1[j]*meanx
#     # b0[j] ~ dnorm(mu_b0[j], tau_b0)
#     # mu_b0[j] <- b0_int
#     #             + b0_lat*lat[j]   * intlat
#     #             + b0_area*area[j] * intarea
#
#     b1[j] <- b1_int
#     # b1[j] ~ dnorm(mu_b1[j], tau_b1)
#     # mu_b1[j] <- b1_int
#     #             + b1_lat*lat[j]   * slopelat
#     #             + b1_area*area[j] * slopearea
#   }
#
#   sig_b0 ~ dunif(0, 10)
#   tau_b0 <- pow(sig_b0, -2)
#
#   sig_b1 ~ dunif(0, 10)
#   tau_b1 <- pow(sig_b1, -2)
#
#   b0_int ~ dnorm(0, 0.001)
#   b0_lat ~ dnorm(0, 0.001)
#   b0_area ~ dnorm(0, 0.001)
#
#   b1_int ~ dnorm(0, 0.001)
#   b1_lat ~ dnorm(0, 0.001)
#   b1_area ~ dnorm(0, 0.001)
#
#   tau <- pow(sig, -2)
#   sig ~ dunif(0, 10)
#
# }', file=lt_jags_common)


## will compare
## - DIC
## - sig
## - b0_interp, b0, b1



par(mfrow=c(4,4))


### defining a list to put all model outputs in
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
                                                     "fit","pred","b0_interp",
                                                     "b0_int","b1_int",
                                                     "b0_lat","b0_area",
                                                     "b1_lat","b1_area"),#,"ypp"
                                n.chains=ncores, parallel=T, n.iter=niter,
                                n.burnin=niter/2, n.thin=niter/2000)
    print(Sys.time() - tstart)
    plotRhats(lt_jags_out)
    alloutputs[[imodel]] <- lt_jags_out
  }
}

theDICs <- sapply(alloutputs, function(x) x$DIC)
bestmodels <- which(theDICs %in% tapply(theDICs, rowSums(controlmat), min))
bestmodels

plot(bestmodels,theDICs[bestmodels])

## saving output for big model runs

########################################################
save(alloutputs, file="alloutputs.Rdata")




## post pred is currently commented out - I was satisfied with it
# for(i in 1:length(alloutputs)) {
#   qq_postpred(alloutputs[[i]], p="ypp", y=lt_data$y, main=i)
# }

## no longer needed
# trim_ypp <- function(aa) {
#   aa <- aa[names(aa) != "samples"]
#   aa$sims.list <- aa$sims.list[names(aa$sims.list) != "ypp"]
#   class(aa) <- "jagsUI"
#   return(aa)
# }
# alloutputs1 <- lapply(alloutputs, trim_ypp)
# save(alloutputs1, file="alloutputs1.Rdata")



### This is where we would have run the models with common intercepts & slopes
### Except we're not.

# ## appending models with common intercept
# controlmat <- rbind(controlmat,
#                     expand.grid(0, 0, 0:1, 0:1))
# for(imodel in (imodel+1):nrow(controlmat)) {
#   lt_data$intlat <- controlmat[imodel, 1]
#   lt_data$intarea <- controlmat[imodel, 2]
#   lt_data$slopelat <- controlmat[imodel, 3]
#   lt_data$slopearea <- controlmat[imodel, 4]
#
#   {
#     tstart <- Sys.time()
#     print(tstart)
#     lt_jags_out <- jagsUI::jags(model.file=lt_jags_commonint, data=lt_data,
#                                 parameters.to.save=c("b0","b1","sig",
#                                                      "mu_b0","sig_b0","mu_b1","sig_b1",
#                                                      "fit","pred","b0_interp",
#                                                      "b0_int","b1_int",
#                                                      "b0_lat","b0_area",
#                                                      "b1_lat","b1_area"),#,"ypp"
#                                 n.chains=ncores, parallel=T, n.iter=niter,
#                                 n.burnin=niter/2, n.thin=niter/2000)
#     print(Sys.time() - tstart)
#     plotRhats(lt_jags_out)
#     alloutputs[[imodel]] <- lt_jags_out
#   }
# }
#
# ## appending models with common slope
# controlmat <- rbind(controlmat,
#                     expand.grid(0:1, 0:1, 0, 0))
# for(imodel in (imodel+1):nrow(controlmat)) {
#   lt_data$intlat <- controlmat[imodel, 1]
#   lt_data$intarea <- controlmat[imodel, 2]
#   lt_data$slopelat <- controlmat[imodel, 3]
#   lt_data$slopearea <- controlmat[imodel, 4]
#
#   {
#     tstart <- Sys.time()
#     print(tstart)
#     lt_jags_out <- jagsUI::jags(model.file=lt_jags_commonslope, data=lt_data,
#                                 parameters.to.save=c("b0","b1","sig",
#                                                      "mu_b0","sig_b0","mu_b1","sig_b1",
#                                                      "fit","pred","b0_interp",
#                                                      "b0_int","b1_int",
#                                                      "b0_lat","b0_area",
#                                                      "b1_lat","b1_area"),#,"ypp"
#                                 n.chains=ncores, parallel=T, n.iter=niter,
#                                 n.burnin=niter/2, n.thin=niter/2000)
#     print(Sys.time() - tstart)
#     plotRhats(lt_jags_out)
#     alloutputs[[imodel]] <- lt_jags_out
#   }
# }
#
# ## appending model with common all
# controlmat <- rbind(controlmat,
#                     expand.grid(0, 0, 0, 0))
# for(imodel in (imodel+1):nrow(controlmat)) {
#   lt_data$intlat <- controlmat[imodel, 1]
#   lt_data$intarea <- controlmat[imodel, 2]
#   lt_data$slopelat <- controlmat[imodel, 3]
#   lt_data$slopearea <- controlmat[imodel, 4]
#
#   {
#     tstart <- Sys.time()
#     print(tstart)
#     lt_jags_out <- jagsUI::jags(model.file=lt_jags_common, data=lt_data,
#                                 parameters.to.save=c("b0","b1","sig",
#                                                      "mu_b0","sig_b0","mu_b1","sig_b1",
#                                                      "fit","pred","b0_interp",
#                                                      "b0_int","b1_int",
#                                                      "b0_lat","b0_area",
#                                                      "b1_lat","b1_area"),#,"ypp"
#                                 n.chains=ncores, parallel=T, n.iter=niter,
#                                 n.burnin=niter/2, n.thin=niter/2000)
#     print(Sys.time() - tstart)
#     plotRhats(lt_jags_out)
#     alloutputs[[imodel]] <- lt_jags_out
#   }
# }


## Extracting & summarizing DIC scores

theDIC <- sapply(alloutputs, function(x) x$DIC)[1:16]
theDIC
# [1] -3600.850 -3600.567 -3601.797 -3601.455
# [5] -3599.787 -3601.461 -3602.278 -3602.288
# [9] -3600.580 -3600.750 -3601.453 -3603.172
# [13] -3602.625 -3599.683 -3603.175 -3602.159

theDIC - min(theDIC)
# [1] 2.325315068 2.608856195 1.378083421
# [4] 1.720082566 3.388809277 1.713968200
# [7] 0.897552104 0.886937279 2.595716843
# [10] 2.424896407 1.722593753 0.003871507
# [13] 0.550545279 3.492543073 0.000000000
# [16] 1.016213227

modeldescription[which.min(theDIC)]

## vector of the number of additional parameters (used for description)
nparam <- rowSums(controlmat)
bestmodels <- which(theDIC %in% tapply(theDIC, rowSums(controlmat), min))


## plotting DIC for each model (a couple ways)
par(mfrow=c(1, 1))
plot(x=theDIC - min(theDIC), y=rank(theDIC),
     xlim=c(0, 10))
text(x=theDIC - min(theDIC), y=rank(theDIC),
     labels=paste(modeldescription, nparam), pos=4)

par(mfrow=c(1, 1))
plot(x=theDIC - min(theDIC), y=seq_along(theDIC),
     xlim=c(0, 10))
text(x=theDIC - min(theDIC), y=seq_along(theDIC),
     labels=paste(modeldescription, nparam), pos=4,
     col=rowSums(controlmat)+1)

# par(mfrow=c(1,1))
# comparecat(alloutputs, p="sig")
# thesigs <- sapply(alloutputs, function(x) jags_df(x, p="sig", exact=T))
# thesigs %>% as.data.frame %>% caterpillar

## making a table
model_df <- data.frame(delta_DIC = theDIC - min(theDIC),
                       Addl_Param = nparam)
rownames(model_df) <- modeldescription
highlight <- rep("", nrow(model_df))
highlight[model_df$delta_DIC %in% tapply(model_df$delta_DIC, model_df$Addl_Param, min)] <- "***"
model_df$highlight <- highlight

model_df <- model_df[order(model_df$delta_DIC),]


## comparison that doesn't do all that much
par(mfrow=c(2,2))
comparecat(alloutputs[controlmat[,1]==1], p="b0_lat")
abline(h=0, lty=2)
comparecat(alloutputs[controlmat[,2]==1], p="b0_area")
abline(h=0, lty=2)
comparecat(alloutputs[controlmat[,3]==1], p="b1_lat")
abline(h=0, lty=2)
comparecat(alloutputs[controlmat[,4]==1], p="b1_area")
abline(h=0, lty=2)


## plotting relationships between lake-level parameters and lake-level variables
## for minimal, preferred, and full models
par(mfrow=c(2,2))
caterpillar(alloutputs[[1]], p="b0", x=lat, xlab="latitude")
caterpillar(alloutputs[[1]], p="b1", x=lat, xlab="latitude")
caterpillar(alloutputs[[1]], p="b0", x=logarea, xlab="log area")
caterpillar(alloutputs[[1]], p="b1", x=logarea, xlab="log area")

caterpillar(alloutputs[[7]], p="b0", x=lat, xlab="latitude")
caterpillar(alloutputs[[7]], p="b1", x=lat, xlab="latitude")
envelope(alloutputs[[7]], p="mu_b1", x=lat, add=TRUE, dark=.2)
caterpillar(alloutputs[[7]], p="b0", x=logarea, xlab="log area")
envelope(alloutputs[[7]], p="mu_b0", x=logarea, add=TRUE, dark=.2)
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



caterpillar(alloutputs[[7]], p="b0", x=lat, xlab="latitude")
envelope(alloutputs[[7]], p="mu_b0", x=lat, add=TRUE, dark=.2)
caterpillar(alloutputs[[7]], p="b1", x=lat, xlab="latitude")
envelope(alloutputs[[7]], p="mu_b1", x=lat, add=TRUE, dark=.2)
caterpillar(alloutputs[[7]], p="b0", x=logarea, xlab="log area")
envelope(alloutputs[[7]], p="mu_b0", x=logarea, add=TRUE, dark=.2)
caterpillar(alloutputs[[7]], p="b1", x=logarea, xlab="log area")
envelope(alloutputs[[7]], p="mu_b1", x=logarea, add=TRUE, dark=.2)

caterpillar(alloutputs[[16]], p="b0", x=lat, xlab="latitude")
envelope(alloutputs[[16]], p="mu_b0", x=lat, add=TRUE, dark=.2)
caterpillar(alloutputs[[16]], p="b1", x=lat, xlab="latitude")
envelope(alloutputs[[16]], p="mu_b1", x=lat, add=TRUE, dark=.2)
caterpillar(alloutputs[[16]], p="b0", x=logarea, xlab="log area")
envelope(alloutputs[[16]], p="mu_b0", x=logarea, add=TRUE, dark=.2)
caterpillar(alloutputs[[16]], p="b1", x=logarea, xlab="log area")
envelope(alloutputs[[16]], p="mu_b1", x=logarea, add=TRUE, dark=.2)

data.frame(modeldescription,
           deltaDIC=theDIC - min(theDIC))


## This section will look a little weird
## I wanted to test the model 7 output reported in the markdown
## it's kind of a scratchpad of things I tried, and I think the markdown is ok
b0test <- -33.5 - (0.031*logarea) + (0.218*lat)
b1test <- 5.55 - (0.036*lat)

caterpillar(alloutputs[[7]], p="b0_interp", x=lat, xlab="latitude")
points(lat, b0test)
caterpillar(alloutputs[[7]], p="b1", x=lat, xlab="latitude")
points(lat, b1test)
envelope(alloutputs[[7]], p="mu_b1", x=lat, add=TRUE, dark=.2)
caterpillar(alloutputs[[7]], p="b0_interp", x=logarea, xlab="log area")
points(logarea, b0test)
envelope(alloutputs[[7]], p="mu_b0", x=logarea, add=TRUE, dark=.2)
caterpillar(alloutputs[[7]], p="b1", x=logarea, xlab="log area")
points(logarea, b1test)

plot(b0test, alloutputs[[7]]$q50$b0_interp)
plot(b1test, alloutputs[[7]]$q50$b1)


## post pred "new" lakes from post object
## compare lake inferences
g0 <- alloutputs[[7]]$sims.list$b0_int
g1 <- alloutputs[[7]]$sims.list$b1_int
t0A <- alloutputs[[7]]$sims.list$b0_area
t1L <- alloutputs[[7]]$sims.list$b1_lat
meanlogA <- mean(log(laketrout_lw$SurfaceArea_h))
meanlogL <- mean(log(laketrout_lw$ForkLength_mm))
meanlat <- mean(laketrout_lw$Latitude_WGS84)

p1 <- g0 - t0A*meanlogA - g1*meanlogL + t1L*meanlat*meanlogL
p2 <- t0A
p3 <- -t1L*meanlogL
p4 <- g1 - t1L*meanlat
p5 <- t1L

b0test <- p1 + outer(p2, logarea) + outer(p3, lat)
b1test <- p4 + outer(p5, lat)

caterpillar(b0test, x=logarea)
envelope(alloutputs[[7]], p="b0_interp", add=T, x=logarea)
caterpillar(b0test, x=lat)
envelope(alloutputs[[7]], p="b0_interp", add=T, x=lat)
caterpillar(b1test, x=logarea)
envelope(alloutputs[[7]], p="b1", add=T, x=logarea)
caterpillar(b1test, x=lat)
envelope(alloutputs[[7]], p="b1", add=T, x=lat)

## ok this is getting stupid, just plot ablines for all lakes
par(mfrow=c(3,4))
for(i in 1:length(logarea)) {
  plot(x=lt_data$x[lt_data$lake==i]+meanlogL, y=lt_data$y[lt_data$lake==i],
       xlim=range(lt_data$x + meanlogL), ylim=range(lt_data$y),
       col=adjustcolor(1, alpha.f=.4))
  abline(a=-33.5 - (0.031*logarea[i]) + (0.218*lat[i]), b=5.55 - (0.036*lat[i]),
         lty=3, lwd=2)
  abline(a=alloutputs[[7]]$q50$b0_interp[i], b=alloutputs[[7]]$q50$b1[i])

}
plot(NA,
     xlim=range(lt_data$x + meanlogL), ylim=range(lt_data$y))
for(i in 1:length(logarea)) {
  abline(a=-33.5 - (0.031*logarea[i]) + (0.218*lat[i]), b=5.55 - (0.036*lat[i]))
}
plot(NA,
     xlim=range(lt_data$x + meanlogL), ylim=range(lt_data$y))
for(i in 1:length(logarea)) {
  abline(a=alloutputs[[7]]$q50$b0_interp[i], b=alloutputs[[7]]$q50$b1[i])
}


## post predict a new lake (Crosswind!)
bestmodels_ind <- bestmodels # which(bestmodels)
newarea <- median(log(laketrout_all$SurfaceArea_h[laketrout_all$LakeName=="Crosswind Lake"])) - meanlogA
  #log(morphometry$SurfaceArea_h[morphometry$LakeName=="Crosswind Lake"])
newlat <- median(laketrout_all$Latitude_WGS84[laketrout_all$LakeName=="Crosswind Lake"]) - meanlat
b0_new <- b1_new <- matrix(nrow=length(alloutputs[[1]]$sims.list$b0_int),
                           ncol=length(bestmodels_ind))

for(imodel in seq_along(bestmodels_ind)) {
  mu_b0_pred <- alloutputs[[bestmodels_ind[imodel]]]$sims.list$b0_int +
    alloutputs[[bestmodels_ind[imodel]]]$sims.list$b0_lat * newlat *
       controlmat[bestmodels_ind[imodel], 1] +
    alloutputs[[bestmodels_ind[imodel]]]$sims.list$b0_area * newarea *
       controlmat[bestmodels_ind[imodel], 2]

  mu_b1_pred <- alloutputs[[bestmodels_ind[imodel]]]$sims.list$b1_int +
    alloutputs[[bestmodels_ind[imodel]]]$sims.list$b1_lat * newlat *
       controlmat[bestmodels_ind[imodel], 3] +
    alloutputs[[bestmodels_ind[imodel]]]$sims.list$b1_area * newarea *
       controlmat[bestmodels_ind[imodel], 4]

  b0_new[, imodel] <- rnorm(n=length(alloutputs[[1]]$sims.list$b0_int),
                            mean=mu_b0_pred,
                            sd=alloutputs[[bestmodels_ind[imodel]]]$sims.list$sig_b0)
  b1_new[, imodel] <- rnorm(n=length(alloutputs[[1]]$sims.list$b0_int),
                            mean=mu_b1_pred,
                            sd=alloutputs[[bestmodels_ind[imodel]]]$sims.list$sig_b1)
}
caterpillar(b0_new)
caterpillar(b1_new)
abline(h=3.2, lty=2)
caterpillar(b0_new - b1_new*lt_data$meanx)
abline(h=-19.56, lty=2)
caterpillar(b1_new)
abline(h=3.2, lty=2)

