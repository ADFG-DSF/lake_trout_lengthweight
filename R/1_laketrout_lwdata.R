library(tidyverse)
library(jagsUI)
library(jagshelper)

# morphometry <- read_csv("flat_data/lake_morphometry.csv", skip=2)
morphometry <- read_csv("flat_data/lake_morphometry2.csv", skip=2)

# # is lake name unique?  YES
# sum(!is.na(morphometry$LakeName))
# length(unique(morphometry$LakeName))

# laketrout_all <- read_csv("flat_data/length_weight.csv", skip=2) %>%
#   left_join(morphometry)
laketrout_all <- read_csv("flat_data/length_weight2.csv", skip=2) %>%
  left_join(morphometry)

sum(!is.na(laketrout_all$ForkLength_mm))  # 34107
sum(!is.na(laketrout_all$Weight_g))  # 3894
sum(!is.na(laketrout_all$Weight_g) & !is.na(laketrout_all$ForkLength_mm))  # 3887

boxplot(laketrout_all$ForkLength_mm ~ laketrout_all$LakeName, las=2, xlab="")

# we have way more data for length than weight!  But a decent amount for weight too


## subsetting for length-weight relationship
laketrout_lw_filter1 <- laketrout_all %>%
  # mutate(Weight_g = ifelse(Weight_g > 100000, Weight_g/1000, Weight_g)) %>%
  mutate(Weight_g = ifelse(Weight_g < 50, Weight_g*1000, Weight_g)) %>%
  filter(Weight_g < 100000) %>%
  filter(ForkLength_mm > 150)



laketrout_lw_filter1 %>%
  # filter(!ProjectTitle %in% c("Mark-Recapture Event 1 - (September - 2002)",
  #                             "Mark-Recapture Event 1 - (September - 2003)",
  #                             "Mark-Recapture Event 1 - (September - 2004)",
  #                             "Mark-Recapture Event 2 - (May - 2003)")) %>%
  ggplot(aes(x=ForkLength_mm, y=Weight_g, colour=LakeName)) +
  geom_point(alpha=.3) +
  # facet_wrap(facets=~LakeName) +
  scale_y_log10() +
  scale_x_log10()

laketrout_lw_filter2 <- laketrout_lw_filter1 %>%
  filter(!ProjectTitle %in% c("Mark-Recapture Event 1 - (September - 2002)",
                              "Mark-Recapture Event 1 - (September - 2003)",
                              "Mark-Recapture Event 1 - (September - 2004)",
                              "Mark-Recapture Event 2 - (May - 2003)"))

laketrout_lw_filter2 %>%
  ggplot(aes(x=ForkLength_mm, y=Weight_g, colour=LakeName)) +
  geom_point(alpha=.3) +
  # facet_wrap(facets=~LakeName) +
  scale_y_log10() +
  scale_x_log10()

# laketrout_lw_filter2 %>%
#   ggplot(aes(x=ForkLength_mm, y=Weight_g, colour=LakeName)) +
#   geom_point(alpha=.3) +
#   facet_wrap(facets=~LakeName) +
#   scale_y_log10() +
#   scale_x_log10()

lm_filter2 <- lm(log(Weight_g/1000)~log(ForkLength_mm), data=laketrout_lw_filter2)
plot(lm_filter2$residuals/sd(lm_filter2$residuals))

laketrout_lw <- laketrout_lw_filter2 %>%
  filter(lm_filter2$residuals/sd(lm_filter2$residuals) > -4 &
           lm_filter2$residuals/sd(lm_filter2$residuals) < 4)

laketrout_lw %>%
  ggplot(aes(x=ForkLength_mm, y=Weight_g)) +#, colour=LakeName
  geom_point(alpha=.3) +
  facet_wrap(facets=~LakeName) +
  scale_y_log10() +
  scale_x_log10()


lm1 <- lm(log(Weight_g/1000) ~
            log(ForkLength_mm),
          data=laketrout_lw)
summary(lm1)
# lm1$residuals %>% plot(col=laketrout_lw$LakeID)
# lm1$residuals %>% plot(col=as.factor(paste(laketrout_lw$LakeID,laketrout_lw$Year)))
# plot(lm1$residuals[laketrout_lw$LakeName=="Paxson Lake"], pch=16,
#      col=as.factor(laketrout_lw$ProjectTitle[laketrout_lw$LakeName=="Paxson Lake"]))
AIC(lm1)

lm2 <- lm(log(Weight_g/1000) ~
            log(ForkLength_mm)*LakeName,
          data=laketrout_lw)
summary(lm2)
AIC(lm2)

lm3 <- lm(log(Weight_g/1000) ~
            log(ForkLength_mm)*LakeName +
            log(ForkLength_mm)*Year,
          data=laketrout_lw)
summary(lm3)
AIC(lm3)

# lm3 <- lm(log(Weight_g/1000)~log(ForkLength_mm)*LakeName + Year, data=laketrout_lw)
# summary(lm3)
# AIC(lm3)

lm4 <- lm(log(Weight_g/1000) ~
            log(ForkLength_mm) +LakeName +
            log(ForkLength_mm)*Latitude_WGS84,
          data=laketrout_lw)
summary(lm4)
AIC(lm4)

lm5 <- lm(log(Weight_g/1000) ~
            log(ForkLength_mm) +
            log(ForkLength_mm)*Latitude_WGS84,
          data=laketrout_lw)
summary(lm5)
AIC(lm5)

anova(lm1, lm2, lm3)

plot(lm2$residuals)
plot(lm2$residuals, col=as.numeric(as.factor(laketrout_lw$LakeName)))
plot(lm2$residuals, col=as.numeric(as.factor(laketrout_lw$ProjectTitle)))
plot(lm2$residuals ~ laketrout_lw$Latitude_WGS84)
