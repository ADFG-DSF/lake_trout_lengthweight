library(tidyverse)


## Whether to write the results to an external file
## This file will be read by
## - the model script
## - the summary document

save_results <- FALSE
# save_results <- TRUE


## load data, filter bad observations BUT NOT observations with missing data

# read lake-level data
morphometry1 <- read_csv("flat_data/lake_morphometry_25_12_23.csv", skip=1) %>%
  filter(!is.na(LakeName)) %>%
  filter(!(is.na(Latitude_WGS84) & is.na(`Elevation (m)`) & is.na(`Temp (C)`) & is.na(SurfaceArea_h))) %>%
  mutate(use_fish = `Include in "Alaskanizing" Modeling Exercise` %in% c("Yes","yes")) %>%
  mutate(make_estimates = `Potentially Include in Lake Trout Management Plan` != "No") %>%
  arrange(LakeName)

# # is lake name unique?  YES
# sum(!is.na(morphometry1$LakeName))
# length(unique(morphometry1$LakeName))


## adding Lester equations for comparison (without model error)
morphometry1 <- morphometry1 %>%
  mutate(DR_lester = MaximumDepth_m/MeanDepth_m) %>%
  mutate(D_th_lester = 3.26*SurfaceArea_h^0.109*MeanDepth_m^0.213*exp(-0.0263*`Temp (C)`)) %>%
  mutate(pV_hy_lester = (1-D_th_lester/MaximumDepth_m)^DR_lester) %>%
  mutate(pV_eb_lester = exp(-4.63*pV_hy_lester)) %>%
  mutate(L_inf_lester = 957*(1-exp(-0.14*(1+log(SurfaceArea_h))))) %>%
  mutate(W_inf_lester = (L_inf_lester/451)^3.2) %>%
  mutate(S_lester = 1/(1+exp(2.47+0.386*`Temp (C)`-16.8*pV_hy_lester))) %>%
  mutate(B_msy_lester = 8.47*(MeanDepth_m*pV_eb_lester*S_lester)/W_inf_lester^1.33) %>%
  mutate(M_lester = 0.26*(exp(0.021*`Temp (C)`+0.0004*`Temp (C)`^2))/W_inf_lester^0.30) %>%
  mutate(msy_ha_lester = B_msy_lester*M_lester) %>%
  mutate(msy_lester = round(msy_ha_lester*SurfaceArea_h, 2))


# read fish-level data
# match lake names to lake-level (morphometry) data as needed
# take out rows where !use_fish
laketrout_raw <- read_csv("flat_data/length_weight_25_12_23.csv", skip=1)
laketrout1 <- laketrout_raw  %>%
  mutate(LakeName = ifelse(LakeName == "Donnelly Lake", "Donnelly Lake (Richardson Highway)", LakeName)) %>%
  mutate(LakeName = ifelse(LakeName == "Four Mile Lake", "Fourmile Lake (Taylor Highway)", LakeName)) %>%
  mutate(LakeName = ifelse(LakeName == "Lost Lake", "Lost Lake  (Chisholm Lake) (near Birch Lake)", LakeName)) %>%
  mutate(LakeName = ifelse(LakeName == "North Twin Lake", "North Twin Lake (Meadows Road)", LakeName)) %>%
  mutate(LakeName = ifelse(LakeName == "Paul's Pond", "Pauls Pond (Coal Mine Road)", LakeName)) %>%
  mutate(LakeName = ifelse(LakeName == "Rapids Lake", "Rapids Lake (Richardson Highway)", LakeName)) %>%
  mutate(LakeName = ifelse(LakeName == "Summit Lake (Paxson)", "Summit Lake (Richardson Highway near Paxson)", LakeName)) %>%
  mutate(LakeName = ifelse(LakeName == "Chelle Lake", "Chelle Lake (Chaleau Lake or Hot Dog Lake)", LakeName)) %>%
  filter(LakeName != "Other") %>%
  left_join(morphometry1) %>%
  filter(make_estimates | use_fish)

nrow(laketrout_raw)  # 38040
nrow(laketrout1)  # 36739
nrow(laketrout_raw) - nrow(laketrout1)  # 1301


# table(morphometry1$make_estimates, morphometry1$LakeName %in% laketrout1$LakeName, morphometry1$use_fish)
morphometry <- morphometry1 %>%
  filter(make_estimates | use_fish)
morphometry$LakeNum <- 1:nrow(morphometry)
nrow(morphometry)  # 62

nrow(morphometry1)  # 263
nrow(morphometry)   #  62

# filling in numerical lake identifier in fish-level data, consistent with lake-level data
laketrout1$LakeNum <- NA
for(i in 1:nrow(morphometry)) {
  laketrout1$LakeNum[laketrout1$LakeName == morphometry$LakeName[i]] <- i
}
table(laketrout1$LakeNum)
table(laketrout1$LakeName)
### is it possible to check that this worked and is consistent??
all(as.numeric(table(laketrout1$LakeNum)) == as.numeric(table(laketrout1$LakeName)))

sort(unique(paste(laketrout1$LakeNum, laketrout1$LakeName)))
sort(unique(paste(morphometry$LakeNum, morphometry$LakeName)))

all(sort(unique(paste(laketrout1$LakeNum, laketrout1$LakeName))) %in%
      sort(unique(paste(morphometry$LakeNum, morphometry$LakeName))))



sapply(laketrout1, function(x) sum(is.na(x)))

summary(laketrout1$Year)  # 1960-2025
length(unique(laketrout1$ProjectTitle))  # 125
length(unique(laketrout1$LakeName))  # 44


# making sure that all lake names in laketrout1 are contained in morphometry
sort(unique(laketrout1$LakeName))[!(sort(unique(laketrout1$LakeName)) %in% morphometry$LakeName)]
morphometry$LakeName[!is.na(morphometry$LakeName)]



## Filtering outliers in
## - Length
## - Weight
## - Age?
## - Weight ~ Length relationship
## - Length ~ Age relationship?

## - one problematic project
## - outlying residuals from a log(Weight) ~ log(Length) regression


laketrout1 %>%
  ggplot(aes(y=ForkLength_mm, x=seq_along(ForkLength_mm))) +
  facet_wrap(~LakeName) +
  geom_point() +
  geom_hline(yintercept=150) +
  scale_y_log10()

laketrout1 %>%
  filter(LakeName=="Sevenmile Lake (Denali Hwy)") %>%
  ggplot(aes(y=ForkLength_mm, x=seq_along(ForkLength_mm),
             colour=factor(Year))) +
  # facet_wrap(~LakeName) +
  geom_point() +
  geom_hline(yintercept=150) +
  scale_y_log10()

laketrout1 %>%
  ggplot(aes(y=Weight_g, x=ForkLength_mm)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10()

laketrout1 %>%
  ggplot(aes(y=Weight_g, x=ForkLength_mm)) +
  facet_wrap(~LakeName) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10()

laketrout1 %>%
  filter(LakeName=="Sevenmile Lake (Denali Hwy)") %>%
  ggplot(aes(y=Weight_g, x=ForkLength_mm,
             colour=factor(Year))) +
  # facet_wrap(~LakeName) +
  geom_point() +
  geom_vline(xintercept=150) +
  scale_y_log10()


# are there many weights without lengths? if not, W~L will be an adequate vis
with(subset(laketrout1, !is.na(Weight_g)), table(is.na(ForkLength_mm)))

# laketrout1 %>%
#   ggplot(aes(y=Weight_g, x=seq_along(Weight_g))) +
#   geom_point() +
#   scale_y_log10()


### FILTER LENGTHS >= 150mm
laketrout2 <- laketrout1 %>%
  filter(ForkLength_mm>=150)
nrow(laketrout1) # 36739
nrow(laketrout2) # 36075
nrow(laketrout1)- nrow(laketrout2) # 664

laketrout1 %>%
  mutate(keep = ifelse(ForkLength_mm>=150, "yes", "no")) %>%
  ggplot(aes(y=ForkLength_mm, x=seq_along(ForkLength_mm), col=keep)) +
  facet_wrap(~LakeName) +
  geom_point() +
  geom_hline(yintercept=150) +
  scale_y_log10()

laketrout2 %>%
  ggplot(aes(y=ForkLength_mm, x=seq_along(ForkLength_mm))) +
  facet_wrap(~LakeName) +
  geom_point() +
  scale_y_log10()


laketrout1 %>%
  mutate(keep = ifelse(ForkLength_mm>=150, "yes", "no")) %>%
  ggplot(aes(y=Weight_g, x=ForkLength_mm, col=keep)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10()

laketrout2 %>%
  ggplot(aes(y=Weight_g, x=ForkLength_mm)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10()

laketrout2 %>%
  ggplot(aes(y=Weight_g, x=ForkLength_mm, col=LakeName)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  theme(legend.position = "none")

laketrout1 %>%
  mutate(keep = ifelse(ForkLength_mm>=150, "yes", "no")) %>%
  ggplot(aes(y=Weight_g, x=ForkLength_mm, col=keep)) +
  facet_wrap(~LakeName) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10()

laketrout2 %>%
  ggplot(aes(y=Weight_g, x=ForkLength_mm)) +
  facet_wrap(~LakeName) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10()

laketrout2 %>%
  filter(LakeName %in% c("Fielding Lake", "Kagati Lake", "Paxson Lake", "Round Tangle Lake")) %>%
  ggplot(aes(y=Weight_g, x=ForkLength_mm, colour=ProjectTitle)) +
  facet_wrap(~LakeName) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10()+
  theme(legend.position = "none")

laketrout1 %>%
  filter(LakeName=="Sevenmile Lake (Denali Hwy)") %>%
  ggplot(aes(y=Weight_g, x=ForkLength_mm,
             colour=factor(Year))) +
  # facet_wrap(~LakeName) +
  geom_point() +
  geom_vline(xintercept=150) +
  scale_y_log10()


# are there many weights without lengths? if not, W~L will be an adequate vis
with(subset(laketrout1, !is.na(Weight_g)), table(is.na(ForkLength_mm)))
with(subset(laketrout2, !is.na(Weight_g)), table(is.na(ForkLength_mm)))


## FILTER: OUTLIERS FROM LOG(W) ~ LOG(L) REGRESSION
laketrout2_lm <- lm(log(laketrout2$Weight_g) ~ log(laketrout2$ForkLength_mm))
summary(laketrout2_lm)
coefs <- coefficients(laketrout2_lm)
rsd <- sd(laketrout2_lm$residuals)


laketrout2 %>%
  ggplot(aes(y=Weight_g, x=ForkLength_mm)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  geom_function(fun=\(x) exp(coefs[1]-2*rsd)*x^coefs[2]) +
  geom_function(fun=\(x) exp(coefs[1]+2*rsd)*x^coefs[2]) +
  geom_function(fun=\(x) exp(coefs[1]-3*rsd)*x^coefs[2]) +
  geom_function(fun=\(x) exp(coefs[1]+3*rsd)*x^coefs[2]) +
  geom_function(fun=\(x) exp(coefs[1]-4*rsd)*x^coefs[2]) +
  geom_function(fun=\(x) exp(coefs[1]+4*rsd)*x^coefs[2]) +
  geom_function(fun=\(x) exp(coefs[1]-5*rsd)*x^coefs[2]) +
  geom_function(fun=\(x) exp(coefs[1]+5*rsd)*x^coefs[2])

sapply(1:5,
       \(x) mean(abs(laketrout2_lm$residuals) > x*rsd))
# [1] 0.7820849 0.9659508 0.9858565 0.9944124 0.9966824
# [1] 0.217915139 0.034049240 0.014143531 0.005587568 0.003317618

### Will censor points with residual > 4 times the residual sd

resids <- log(laketrout2$Weight_g) - predict(laketrout2_lm, newdata = laketrout2)
sd(resids, na.rm=TRUE) - rsd
laketrout2 %>%
  mutate(keep=ifelse(abs(resids) <= 4*rsd, "yes", "no")) %>%
  ggplot(aes(y=Weight_g, x=ForkLength_mm, col=keep)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10()


laketrout3 <- laketrout2 %>%
  filter(is.na(resids) | abs(resids) <= 4*rsd)

nrow(laketrout2) # 36075
nrow(laketrout3) # 36043
nrow(laketrout2) - nrow(laketrout3) # 32


# are there any fish with ages but no lengths?
with(subset(laketrout3, !is.na(Age)), table(is.na(ForkLength_mm)))


laketrout3 %>%
  filter(is.na(Age) | Age < 50) %>%
  ggplot(aes(x=Age, y=ForkLength_mm)) +
  geom_point()
max(laketrout3$Age, na.rm=TRUE)


laketrout3 %>%
  filter(Age < 50) %>%
  ggplot(aes(x=Age, y=ForkLength_mm)) +
  facet_wrap(~LakeName) +
  # ggplot(aes(x=Age, y=ForkLength_mm, col=LakeName)) +
  geom_point() +
  theme(legend.position = "none")


## finalizing and saving the current filtration
laketrout <- laketrout3

save_results
if(save_results) {
  save(laketrout, morphometry, file="Rdata/laketrout_sampling_formodel.Rdata")
}

## quick summary of how many samples there are

## number of lakes with samples
length(unique(laketrout$LakeName))
# 44

## number of lakes with samples used by the model
sum(morphometry$use_fish)
# 40

## number of lakes we want inferences for
sum(morphometry$make_estimates)
# 45

## cross table
with(morphometry, table(use_fish, make_estimates))
#          make_estimates
# use_fish FALSE TRUE
# FALSE        0   22
# TRUE        17   23
## 22 lakes of management interest with no fish samples
## 17 lakes with fish samples but not of management interest
## 23 lakes with both

## how many weights
sum(!is.na(laketrout$Weight_g))
# 5695

## how many lengths
sum(!is.na(laketrout$ForkLength_mm))
# 36043

## how many ages
sum(!is.na(laketrout$Age))
# 1496

## how many fish with both length and weight
sum(!is.na(laketrout$ForkLength_mm & laketrout$Weight_g))
# 5695

## how many fish with both length and age
sum(!is.na(laketrout$ForkLength_mm & laketrout$Age))
# 1496

## how many lakes with length and weight samples
length(unique(laketrout$LakeName[!is.na(laketrout$ForkLength_mm & laketrout$Weight_g)]))
# 29

## how many lakes with length and age samples
length(unique(laketrout$LakeName[!is.na(laketrout$ForkLength_mm & laketrout$Age)]))
# 19


## tabulate availability of lake-level data and make_estimates and use_fish
with(morphometry, table(use_fish, is.na(msy_lester)))
# use_fish FALSE TRUE
# FALSE        7   15
# TRUE        21   19
## for lakes with fish samples, 21 have data to support Lester, 19 do not

with(morphometry, table(make_estimates, is.na(msy_lester)))
# make_estimates FALSE TRUE
# FALSE              4   13
# TRUE              24   21
## for lakes of management interest, 24 have data to support Lester, 21 do not

with(morphometry, table(use_fish, is.na(SurfaceArea_h)))
# use_fish FALSE TRUE
#    FALSE    22    0
#    TRUE     36    4
## for lakes with fish samples, 36 have lake-level data, 4 do not

with(morphometry, table(make_estimates, is.na(SurfaceArea_h)))
# make_estimates FALSE TRUE
#          FALSE    13    4
#          TRUE     45    0
## for lakes of management interest, all 45 have lake-level data
