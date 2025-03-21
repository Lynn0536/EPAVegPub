
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


# Call in Means Data 
# This data was averaged across season within each year 

Data <- read_csv("02.CleanData/CleanData_SpringOnly.csv",
                 col_types = cols(Year = col_factor(),
                                  Position = col_factor(),
                                  Treatment = col_factor(),
                                  Site = col_factor()))


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

# The Beta distribution requires values be between -, but can't actually BE 0 or 1 
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
# 

###################### Getting Model Results ####################################
Anova(TestMod1, type='III')

Anova(TestMod2, type="III")

Anova(TestMod3, type='III')

# Pairwise Treatment x Year 
PairsWoodyMod_All <- emmeans(TestMod1, ~Treatment*Year, type='response') 
pairs(PairsWoodyMod_All)
CI_Letters <- cld(PairsWoodyMod_All, Letters=letters, sort=TRUE, decreasing=TRUE)

