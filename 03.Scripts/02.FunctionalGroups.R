
################################################################################
######################## Ashlynn Smith #########################################
################## Script to Visualize and Analyze #############################
#############  Functional Group Vegetation Data : EPA Project  #################

# Load Packages 
library(tidyverse)
library(glmmTMB)
library(ggpattern)
library(Rmisc)
library(ISLR)
library(report)
library(sjPlot)
library(flextable)
library(boot)
library(bbmle)
library(performance)
library(DHARMa)
library(agridat)
library(car)
library(aods3)
library(vegan)
library(plyr)
library(cowplot)
library(Matrix)
library(coin)
library(grid)
library(ggplot2)
library(ggpubr)
library(effects)
library(lmerTest)
library(MuMIn)
library(glmm.hp)
library(emmeans)
library(multcomp)


# Call in Data 
# This data was selected from just the Spring data each year 

Data <- read_csv("02.CleanData/CleanData_SpringOnly.csv",
                 col_types = cols(Year = col_factor(),
                                  Position = col_factor(),
                                  Treatment = col_factor(),
                                  Site = col_factor()))

############################## WOODY COVER #####################################
############ Visualize the Woody Data ##### Decide on Outliers NOW #############

ggplot(Data, aes(x=Treatment, y=Woody_Cover, fill = SeasonYear)) + 
  geom_boxplot(notch=FALSE)

# I don't think I am going to remove any outliers from the Woody Data 
# The reason is because this was data collected in the field, not by a machine, & 
# I think all of the data represents real data 

####################### Distribution of Woody Data #############################
 
ggplot(Data, aes(x=Woody_Cover)) +
  geom_histogram()

# The distribution of this data is not normal, so it makes me think there will be issues 
# with modeling the data. 

######################## Woody Model Testing Code ##############################

## Full Model 

TestMod1 <- glmmTMB(Woody_Cover ~ Treatment * Year + (1|Site:Position), data=Data)

TestMod2 <- glmmTMB(Woody_Cover ~ Treatment * Year + (1|Site:Position), family = Gamma(link='log'), data=Data)

# Create the log of Woody_Cover then run a new model 

Data <- Data %>%
  mutate(Woody_log = log(Woody_Cover))

max(Data$Woody_log)
min(Data$Woody_log)

TestMod3 <- glmmTMB(Woody_log ~ Treatment * Year + (1|Site:Position), data=Data)

################### Transforming percent data to a proportion ##################

# The Beta distribution requires values be between, but can't actually BE 0 or 1 
# So I need to transform the percents to proportions just over zero and just below 1

### Find Max Value and Create the Beta Distribution 

max(Data$Woody_Cover) #196.82 (over 100% because of how woody cover was calculated - this is okay)
                      # This high number indicates there was nearly 100% woody cover at both the ground 
                      # and in the canopy. Converting to a 0-1 scale required for the Beta distribution 
                      # really is going to eliminate these high values, but I am going to try it 

# Now create a beta distribution (all values between 0-1) based on the max value 

Data <- Data %>%
  mutate(Woody_B = (Woody_Cover/197.5)+0.001)

max(Data$Woody_B)
min(Data$Woody_B)

TestMod4 <- glmmTMB(Woody_B ~ Treatment * Year + (1|Site:Position), zi=~Woody_B, data=Data, family=beta_family())

# Now look at the residuals for TestMod1 
plot(simulateResiduals(TestMod1))
hist(simulateResiduals(TestMod1)) ## histogram should be flat
qqPlot(resid(TestMod1))  ## residuals should line up pretty closely to the blue line
hist(residuals(TestMod1))

# Now look at the residuals for TestMod2
plot(simulateResiduals(TestMod2))
hist(simulateResiduals(TestMod2)) ## histogram should be flat
qqPlot(resid(TestMod2))  ## residuals should line up pretty closely to the blue line
hist(residuals(TestMod2))

# Now look at the residuals for TestMod3
plot(simulateResiduals(TestMod3))
hist(simulateResiduals(TestMod3)) ## histogram should be flat
qqPlot(resid(TestMod3))  ## residuals should line up pretty closely to the blue line
hist(residuals(TestMod3))

# Now look at the residuals for TestMod4
plot(simulateResiduals(TestMod4))
hist(simulateResiduals(TestMod4)) ## histogram should be flat
qqPlot(resid(TestMod4))  ## residuals should line up pretty closely to the blue line
hist(residuals(TestMod4))

# Based on visual examination of the residual plots, the best model appears to be 
# TestMod3 where the response variable (Woody_Cover) was log transformed. 

###################### Getting Model Results ####################################
Anova(TestMod3, type='III')

# Pairwise Treatment x Year 
PairsWoodyMod_All <- emmeans(TestMod3, ~Treatment*Year, type='response') 
pairs(PairsWoodyMod_All)
CI_Letters <- cld(PairsWoodyMod_All, Letters=letters, sort=TRUE, decreasing=TRUE)

# Now I want to calculate woody cover means to see if they match up to the model \
# output & 'CI_Letters' 

WoodyCover <- Data %>%
  dplyr::group_by(Treatment, Year) %>%
  dplyr::summarise(MEAN = mean(Woody_Cover)) 

# CI Letters reports in the original data which is the log of the actual mean
# So now I need to transform the emmean values and confidence limits back to the 
# original values in order to graph and compare to the mean values the code above gives me. 

CI <- CI_Letters %>%
  mutate(CL_low = exp(lower.CL) - 1)

CI <- CI %>%
  mutate(CL_up = exp(upper.CL)-1)

CI <- CI %>%
  mutate(Mean = exp(emmean)-1)

# The mean values reported by emmeans column do not match the mean values 
# calculated when I group by Treatment then Year and calculate the average 

# Another Means Table 

WoodyTable <- summarySE(Data, measurevar="Woody_Cover", groupvars=c("Treatment", "Year")) 

WoodyTable <- WoodyTable %>%
  mutate(Treatment = fct_relevel(Treatment, "Control", "Reference", "Cleared","Scraped", "Diaspore", "Cleared + Burned",
                                 "Cleared + Scraped", "Cleared + Diaspore","Scraped + Diaspore", "Cleared + Scraped + Diaspore",
                                 "Cleared + Diaspore + Burned","Cleared + Scraped + Diaspore + Burned"))

# TN_Summary <- summarySE(TNData, measurevar="Total_N", groupvars=c("Treatment", "Year")) %>%
  # mutate(Time = case_when(Year == '2022' ~ '2-Weeks Post Treatment (2022)',
                          # Year == '2023' ~ '1-Year Post Treatment (2023)',
                          # Year == '2024' ~ '2-Years Post Treatment (2024)')) %>%
  # mutate(Time = fct_relevel(Time, "2-Weeks Post Treatment (2022)", "1-Year Post Treatment (2023)",
                            # "2-Years Post Treatment (2024)"))


WoodyFig1 <- ggplot(WoodyTable, aes(x=Treatment, y=Woody_Cover, group = Year))+ 
  geom_errorbar(aes(ymin=Woody_Cover-ci, ymax=Woody_Cover+ci), width=.2, 
                position=position_dodge(0.5)) +
  geom_point(aes(x=Treatment, y=Woody_Cover), size=3, 
             position=position_dodge(0.5)) +
  labs(x="Treatment", y = "Woody Cover (%)")+
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 14),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.title = element_text(size=20))

# Figure Note: I exported this figure and made adjustments (shading and line connection only)
# In PowerPoint 

############################## GRAMINOID COVER #################################
############ Visualize the Graminoid Data ##### Decide on Outliers NOW #############

ggplot(Data, aes(x=Treatment, y=Gram_Cover, fill = SeasonYear)) + 
  geom_boxplot(notch=FALSE)

# There is a very high value for the graminoid cover for one Cleared plot in 
# Spring 2023. I need to check the raw data sheet and see if it was an entry error 
# or could have possibly been a recording error in the field 

# For now, I will remove any values over 150% for graminoid cover to get rid of that 
# extremely high value. 

GramData <- subset(Data, Gram_Cover<150) 

# Now re graph it to make sure that worked 

ggplot(GramData, aes(x=Treatment, y=Gram_Cover, fill = SeasonYear)) + 
  geom_boxplot(notch=FALSE)

####################### Distribution of Graminoid Data #############################

ggplot(GramData, aes(x=Gram_Cover)) +
  geom_histogram()

# The distribution of this data is not normal, so it makes me think there will be issues 
# with modeling the data. The data is very zero inflated so I will test a 
# zero inflation model 

######################## Graminoid Model Testing Code ##############################

## Full Guass Model 

GramMod1 <- glmmTMB(Gram_Cover ~ Treatment * Year + (1|Site:Position), data=GramData)

GramMod2 <- glmmTMB(Gram_Cover ~ Treatment * Year + (1|Site:Position), family = Gamma(link='log'), data=GramData)

# Create the log of Woody_Cover then run a new model 

GramData <- GramData %>%
  mutate(Gram_log = log(Gram_Cover))

max(GramData$Gram_log)
min(GramData$Gram_log)

# log transformation produces '-Inf' values. This means the model will not work 

TestMod3 <- glmmTMB(Gram_log ~ Treatment * Year + (1|Site:Position), data=Data)

########### Full Model Testing Beta ############################################

#I want to test model selection with Beta distribution because this data is very zero inflated 

### Find Max Value and Create the Beta Distribution 

max(GramData$Gram_Cover) #125.1 (over 100% because of how graminoid cover was calculated - this is okay)

# Now create a beta distribution (all values between 0-1) based on the max value 

GramData <- GramData %>%
  mutate(Gram_B = (Gram_Cover/126.6)+0.001)

max(GramData$Gram_B)
min(GramData$Gram_B)

# Now I can run model with the beta distribution because all values are greater than zero but less than 1 
## Full Model with Beta distribution but different links specified 

GramMod3 <- glmmTMB(Gram_B ~ Treatment * Year + (1|Site:Position), zi=~Gram_B, data=GramData, family=beta_family())

# Now look at the residuals for GramMod1 
plot(simulateResiduals(GramMod1))
hist(simulateResiduals(GramMod1)) ## histogram should be flat
qqPlot(resid(GramMod1))  ## residuals should line up pretty closely to the blue line
hist(residuals(GramMod1))

# Now look at the residuals for GramMod3 
plot(simulateResiduals(GramMod3))
hist(simulateResiduals(GramMod3)) ## histogram should be flat
qqPlot(resid(GramMod3))  ## residuals should line up pretty closely to the blue line
hist(residuals(GramMod3))

#GramMod3b <- glmmTMB(Gram_B ~ Treatment * Year + (1|Site:Position), data=GramData, family=beta_family(link="log"))

GramMod3c <- glmmTMB(Gram_B ~ Treatment * Year + (1|Site:Position), data=GramData, family=beta_family())

GramMod3d <- glmmTMB(Gram_B ~ Treatment * Year + (1|Site:Position), data=GramData, family=beta_family(link="probit"))

#GramMod3e <- glmmTMB(Gram_B ~ Treatment * Year + (1|Site:Position), data=GramData, family=beta_family(link="inverse"))

GramMod3f <- glmmTMB(Gram_B ~ Treatment * Year + (1|Site:Position), data=GramData, family=beta_family(link="cloglog"))

#GramMod3g <- glmmTMB(Gram_B ~ Treatment * Year + (1|Site:Position), data=GramData, family=beta_family(link="identity"))

#GramMod3h <- glmmTMB(Gram_B ~ Treatment * Year + (1|Site:Position), data=GramData, family=beta_family(link="sqrt"))

AIC(GramMod3, GramMod3c, GramMod3d, GramMod3f)
BIC(GramMod3, GramMod3c, GramMod3d, GramMod3f)

# Lowest AIC and BIC values are produced for GramMod3d but the values are not really that different 

# Now look at the residuals for GramMod3d 
plot(simulateResiduals(GramMod3d))
hist(simulateResiduals(GramMod3d)) ## histogram should be flat
qqPlot(resid(GramMod3d))  ## residuals should line up pretty closely to the blue line
hist(residuals(GramMod3d))

### Just going to try and see whether a sqrt transformation works better 
GramData <- GramData %>%
  mutate(Gram_SQRT = sqrt(Gram_Cover))

max(GramData$Gram_SQRT)
min(GramData$Gram_SQRT)

GramMod4 <- glmmTMB(Gram_SQRT ~ Treatment * Year + (1|Site:Position), data=GramData)

# Now look at the residuals for GramMod4 
plot(simulateResiduals(GramMod4))
hist(simulateResiduals(GramMod4)) ## histogram should be flat
qqPlot(resid(GramMod4))  ## residuals should line up pretty closely to the blue line
hist(residuals(GramMod4))

# A sqrt transformation of the data greatly improved model fit. Model residuals are very close to 
# a normal distribution and qqplot of model residuals hugs the line tightly through the 
# first quantile 

###################### Getting Model Results ####################################
Anova(GramMod4, type='III')

# Pairwise Treatment x Year 
PairsGramMod_All <- emmeans(GramMod4, ~Treatment*Year, type='response') 
pairs(PairsGramMod_All)
CI_Letters <- cld(PairsGramMod_All, Letters=letters, sort=TRUE, decreasing=TRUE)

# Pairwise Treatment | Year 
PairsGramMod <- emmeans(GramMod4, ~Treatment|Year, type='response') 
pairs(PairsGramMod)
CI_Letters <- cld(PairsGramMod, Letters=letters, sort=TRUE, decreasing=TRUE)

# Graminoid Means Table 

GramTable <- summarySE(GramData, measurevar="Gram_Cover", groupvars=c("Treatment", "Year")) 

GramTable <- GramTable %>%
  mutate(Treatment = fct_relevel(Treatment, "Control", "Reference", "Cleared","Scraped", "Diaspore", "Cleared + Burned",
                                 "Cleared + Scraped", "Cleared + Diaspore","Scraped + Diaspore", "Cleared + Scraped + Diaspore",
                                 "Cleared + Diaspore + Burned","Cleared + Scraped + Diaspore + Burned"))

# TN_Summary <- summarySE(TNData, measurevar="Total_N", groupvars=c("Treatment", "Year")) %>%
# mutate(Time = case_when(Year == '2022' ~ '2-Weeks Post Treatment (2022)',
# Year == '2023' ~ '1-Year Post Treatment (2023)',
# Year == '2024' ~ '2-Years Post Treatment (2024)')) %>%
# mutate(Time = fct_relevel(Time, "2-Weeks Post Treatment (2022)", "1-Year Post Treatment (2023)",
# "2-Years Post Treatment (2024)"))

GramFig1 <- ggplot(GramTable, aes(x=Treatment, y=Gram_Cover, group = Year))+ 
  geom_errorbar(aes(ymin=Gram_Cover-ci, ymax=Gram_Cover+ci), width=.2, 
                position=position_dodge(0.5)) +
  geom_point(aes(x=Treatment, y=Gram_Cover), size=3, 
             position=position_dodge(0.5)) +
  labs(x="Treatment", y = "Graminoid Cover (%)")+
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 14),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.title = element_text(size=20))

############################## FORB COVER #################################
############ Visualize the FORB Data ##### Decide on Outliers NOW #############

ggplot(Data, aes(x=Treatment, y=Forb_Cover, fill = SeasonYear)) + 
  geom_boxplot(notch=FALSE)


######################## Model Testing Code ####################################

#I want to test a model selection with Beta distribution because that was the best option for 
#the woody & Graminoid data sets and this one is even more zero inflated! 

### Find Max Value and Create the Beta Distribution 

max(Data$Forb_Cover) #85.5 

# Now create a beta distribution (all values between 0-1) based on the max value 

Data <- Data %>%
  mutate(Forb_B = (Forb_Cover/86)+0.001)

max(Data$Forb_B)
min(Data$Forb_B)

# Now I can run model with the beta distribution because all values are greater than zero but less than 1 

# Now I can run model with the beta distribution because all values are greater than zero but less than 1 
## Full Model with Beta distribution but different links specified 

#ForbMod1 <- glmmTMB(Forb_B ~ Treatment * Year + (1|Site:Position), data=Data, family=beta_family(link="log"))

ForbMod2 <- glmmTMB(Forb_B ~ Treatment * Year + (1|Site:Position), data=Data, family=beta_family())

ForbMod2a <- glmmTMB(Forb_B ~ Treatment * Year + (1|Site:Position), zi=~Forb_B, data=Data, family=beta_family())

ForbMod3 <- glmmTMB(Forb_B ~ Treatment * Year + (1|Site:Position), data=Data, family=beta_family(link="probit"))

#ForbMod4 <- glmmTMB(Forb_B ~ Treatment * Year + (1|Site:Position), data=Data, family=beta_family(link="inverse"))

ForbMod5 <- glmmTMB(Forb_B ~ Treatment * Year + (1|Site:Position), data=Data, family=beta_family(link="cloglog"))

#ForbMod6 <- glmmTMB(Forb_B ~ Treatment * Year + (1|Site:Position), data=Data, family=beta_family(link="identity"))

#ForbMod7 <- glmmTMB(Forb_B ~ Treatment * Year + (1|Site:Position), data=Data, family=beta_family(link="sqrt"))

AIC(ForbMod2, ForbMod2a, ForbMod3, ForbMod5)
BIC(ForbMod2, ForbMod2a, ForbMod3, ForbMod5)

# Of the Beta models, ForbMod3 has the lowest AIC and BIC values 
# Now I will take a look at the model residuals 
# If not a good fit, then I will explor other model options 

# Now look at the residuals for ForbMod3 
plot(simulateResiduals(ForbMod3))
hist(simulateResiduals(ForbMod3)) ## histogram should be flat
qqPlot(resid(ForbMod3))  ## residuals should line up pretty closely to the blue line
hist(residuals(ForbMod3))

### Just going to try and see whether a sqrt transformation works better 
Data <- Data %>%
  mutate(Forb_SQRT = sqrt(Forb_Cover))

max(Data$Forb_SQRT)
min(Data$Forb_SQRT)

ForbMod8 <- glmmTMB(Forb_SQRT ~ Treatment * Year + (1|Site:Position), data=Data)

# Now look at the residuals for GramMod4 
plot(simulateResiduals(ForbMod8))
hist(simulateResiduals(ForbMod8)) ## histogram should be flat
qqPlot(resid(ForbMod8))  ## residuals should line up pretty closely to the blue line
hist(residuals(ForbMod8))

# The simulated model residuals still look pretty bad but this model produces 
# the best qqplot and histogram normal distribution compared to the other models tested 
# thus far. I am going to stick with this for now and play around more with zero inflation later 

# Forb Model Results 

Anova(ForbMod8, type='III')

# Pairwise Treatment x Year 
PairsForbMod_All <- emmeans(ForbMod8, ~Treatment*Year, type='response') 
pairs(PairsForbMod_All)
CI_Letters_Forb <- cld(PairsForbMod_All, Letters=letters, sort=TRUE, decreasing=TRUE)

# Forb Means Table 

ForbTable <- summarySE(Data, measurevar="Forb_Cover", groupvars=c("Treatment", "Year")) 

ForbTable <- ForbTable %>%
  mutate(Treatment = fct_relevel(Treatment, "Control", "Reference", "Cleared","Scraped", "Diaspore", "Cleared + Burned",
                                 "Cleared + Scraped", "Cleared + Diaspore","Scraped + Diaspore", "Cleared + Scraped + Diaspore",
                                 "Cleared + Diaspore + Burned","Cleared + Scraped + Diaspore + Burned"))

# TN_Summary <- summarySE(TNData, measurevar="Total_N", groupvars=c("Treatment", "Year")) %>%
# mutate(Time = case_when(Year == '2022' ~ '2-Weeks Post Treatment (2022)',
# Year == '2023' ~ '1-Year Post Treatment (2023)',
# Year == '2024' ~ '2-Years Post Treatment (2024)')) %>%
# mutate(Time = fct_relevel(Time, "2-Weeks Post Treatment (2022)", "1-Year Post Treatment (2023)",
# "2-Years Post Treatment (2024)"))

ForbFig1 <- ggplot(ForbTable, aes(x=Treatment, y=Forb_Cover, group = Year))+ 
  geom_errorbar(aes(ymin=Forb_Cover-ci, ymax=Forb_Cover+ci), width=.2, 
                position=position_dodge(0.5)) +
  geom_point(aes(x=Treatment, y=Forb_Cover), size=3, 
             position=position_dodge(0.5)) +
  labs(x="Treatment", y = "Forb Cover (%)")+
  theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 14),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.title = element_text(size=20))



