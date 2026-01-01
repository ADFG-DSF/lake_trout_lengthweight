library(tidyverse)
library(jagsUI)
library(jagshelper)


## original laketrout_lw had 4075 rows

# morphometry <- read_csv("flat_data/lake_morphometry.csv", skip=2)
morphometry <- read_csv("flat_data/lake_morphometry2.csv", skip=2)

# # is lake name unique?  YES
# sum(!is.na(morphometry$LakeName))
# length(unique(morphometry$LakeName))

# laketrout_all <- read_csv("flat_data/length_weight.csv", skip=2) %>%
#   left_join(morphometry)
laketrout_all <- read_csv("flat_data/length_weight2.csv", skip=2) %>%
  left_join(morphometry)
nrow(laketrout_all)  # 35516
sapply(laketrout_all, function(x) sum(is.na(x)))

summary(laketrout_all$Year)  # 1960-2024
length(unique(laketrout_all$ProjectTitle))  # 148
length(unique(laketrout_all$LakeName))  # 84


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
nrow(laketrout1) # 32539

lm1 <- with(laketrout1, lm(log(Weight_g) ~ log(ForkLength_mm)))
resids1 <- log(laketrout1$Weight_g) - predict(lm1, newdata=laketrout1)
laketrout2 <- filter(laketrout1, is.na(resids1) | abs(resids1) < 4*sd(resids1, na.rm=TRUE))

laketrout <- laketrout2 %>%
  filter(is.na(Age) | Age < 50) %>%
  mutate(LakeNum = as.numeric(as.factor(LakeName)))

nrow(laketrout) # 32516

### plotting weight ~ length
laketrout_all %>%
  ggplot(aes(x=ForkLength_mm, y=Weight_g, color=LakeName)) +
  # ggplot(aes(x=ForkLength_mm, y=Weight_g)) +
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
table(is.na(laketrout$Latitude_WGS84),
      is.na(laketrout$SurfaceArea_h),
      is.na(laketrout$Weight_g))

nrow(laketrout)

sum(!is.na(laketrout$ForkLength_mm) & !is.na(laketrout$Weight_g))
length(unique(laketrout$LakeNum[!is.na(laketrout$ForkLength_mm) & !is.na(laketrout$Weight_g)]))

sum(!is.na(laketrout$ForkLength_mm) & !is.na(laketrout$Age))
length(unique(laketrout$LakeNum[!is.na(laketrout$ForkLength_mm) & !is.na(laketrout$Age)]))

length(unique(laketrout$LakeNum[!is.na(laketrout$Latitude_WGS84) & !is.na(laketrout$SurfaceArea_h)]))

sapply(laketrout, function(x) sum(!is.na(x)))
sapply(laketrout, function(x) sum(is.na(x)))



laketrout_lw <- laketrout %>%
  filter(!is.na(ForkLength_mm)) %>%
  filter(!is.na(Weight_g)) %>%
  filter(!is.na(Latitude_WGS84)) %>%
  filter(!is.na(SurfaceArea_h))
nrow(laketrout_lw)




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
#
# sum(!is.na(laketrout_all$ForkLength_mm))  # 34107
# sum(!is.na(laketrout_all$Weight_g))  # 3894
# sum(!is.na(laketrout_all$Weight_g) & !is.na(laketrout_all$ForkLength_mm))  # 3887
#
# boxplot(laketrout_all$ForkLength_mm ~ laketrout_all$LakeName, las=2, xlab="")
#
# # we have way more data for length than weight!  But a decent amount for weight too
#
#
# ## subsetting for length-weight relationship
# laketrout_lw_filter1 <- laketrout_all %>%
#   # mutate(Weight_g = ifelse(Weight_g > 100000, Weight_g/1000, Weight_g)) %>%
#   mutate(Weight_g = ifelse(Weight_g < 50, Weight_g*1000, Weight_g)) %>%
#   filter(Weight_g < 100000) %>%
#   filter(ForkLength_mm > 150)
#
#
#
# laketrout_lw_filter1 %>%
#   # filter(!ProjectTitle %in% c("Mark-Recapture Event 1 - (September - 2002)",
#   #                             "Mark-Recapture Event 1 - (September - 2003)",
#   #                             "Mark-Recapture Event 1 - (September - 2004)",
#   #                             "Mark-Recapture Event 2 - (May - 2003)")) %>%
#   ggplot(aes(x=ForkLength_mm, y=Weight_g, colour=LakeName)) +
#   geom_point(alpha=.3) +
#   # facet_wrap(facets=~LakeName) +
#   scale_y_log10() +
#   scale_x_log10()
#
# laketrout_lw_filter2 <- laketrout_lw_filter1 %>%
#   filter(!ProjectTitle %in% c("Mark-Recapture Event 1 - (September - 2002)",
#                               "Mark-Recapture Event 1 - (September - 2003)",
#                               "Mark-Recapture Event 1 - (September - 2004)",
#                               "Mark-Recapture Event 2 - (May - 2003)"))
#
# laketrout_lw_filter2 %>%
#   ggplot(aes(x=ForkLength_mm, y=Weight_g, colour=LakeName)) +
#   geom_point(alpha=.3) +
#   # facet_wrap(facets=~LakeName) +
#   scale_y_log10() +
#   scale_x_log10()
#
# # laketrout_lw_filter2 %>%
# #   ggplot(aes(x=ForkLength_mm, y=Weight_g, colour=LakeName)) +
# #   geom_point(alpha=.3) +
# #   facet_wrap(facets=~LakeName) +
# #   scale_y_log10() +
# #   scale_x_log10()
#
# lm_filter2 <- lm(log(Weight_g/1000)~log(ForkLength_mm), data=laketrout_lw_filter2)
# plot(lm_filter2$residuals/sd(lm_filter2$residuals))
#
# laketrout_lw <- laketrout_lw_filter2 %>%
#   filter(lm_filter2$residuals/sd(lm_filter2$residuals) > -4 &
#            lm_filter2$residuals/sd(lm_filter2$residuals) < 4)
#
# laketrout_lw %>%
#   ggplot(aes(x=ForkLength_mm, y=Weight_g)) +#, colour=LakeName
#   geom_point(alpha=.3) +
#   facet_wrap(facets=~LakeName) +
#   scale_y_log10() +
#   scale_x_log10()
#
#
# lm1 <- lm(log(Weight_g/1000) ~
#             log(ForkLength_mm),
#           data=laketrout_lw)
# summary(lm1)
# # lm1$residuals %>% plot(col=laketrout_lw$LakeID)
# # lm1$residuals %>% plot(col=as.factor(paste(laketrout_lw$LakeID,laketrout_lw$Year)))
# # plot(lm1$residuals[laketrout_lw$LakeName=="Paxson Lake"], pch=16,
# #      col=as.factor(laketrout_lw$ProjectTitle[laketrout_lw$LakeName=="Paxson Lake"]))
# AIC(lm1)
#
# lm2 <- lm(log(Weight_g/1000) ~
#             log(ForkLength_mm)*LakeName,
#           data=laketrout_lw)
# summary(lm2)
# AIC(lm2)
#
# lm3 <- lm(log(Weight_g/1000) ~
#             log(ForkLength_mm)*LakeName +
#             log(ForkLength_mm)*Year,
#           data=laketrout_lw)
# summary(lm3)
# AIC(lm3)
#
# # lm3 <- lm(log(Weight_g/1000)~log(ForkLength_mm)*LakeName + Year, data=laketrout_lw)
# # summary(lm3)
# # AIC(lm3)
#
# lm4 <- lm(log(Weight_g/1000) ~
#             log(ForkLength_mm) +LakeName +
#             log(ForkLength_mm)*Latitude_WGS84,
#           data=laketrout_lw)
# summary(lm4)
# AIC(lm4)
#
# lm5 <- lm(log(Weight_g/1000) ~
#             log(ForkLength_mm) +
#             log(ForkLength_mm)*Latitude_WGS84,
#           data=laketrout_lw)
# summary(lm5)
# AIC(lm5)
#
# anova(lm1, lm2, lm3)
#
# plot(lm2$residuals)
# plot(lm2$residuals, col=as.numeric(as.factor(laketrout_lw$LakeName)))
# plot(lm2$residuals, col=as.numeric(as.factor(laketrout_lw$ProjectTitle)))
# plot(lm2$residuals ~ laketrout_lw$Latitude_WGS84)
