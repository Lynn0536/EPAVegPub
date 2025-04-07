#################################################
###               Ashlynn Smith               ###
###        SNRE - University of Florida       ###
###           Diversity Analysis              ###
###                                           ###
###           Restoration Experiment          ### 
###            Funding: U.S. EPA              ###

# Load Packages 
library(Rmisc)
library(pals)
library(ggsci)
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
library(tidyverse)
library(vegan)
library(plyr)
library(cowplot)
library(glmmTMB)
library(Matrix)
library(lme4)
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

####################### Call In Data ###########################################

## This data was fully manipulated and means were calculated in R 
## I saved the file and am just calling in this dataset 

Data <- read_csv("02.CleanData/CleanData_AnnualMeans.csv",
                 col_types = cols(Year = col_factor(),
                                  Position = col_factor(),
                                  Treatment = col_factor(),
                                  Site = col_factor()))


# CReate a UniqueID 

Data <- Data %>%
  group_by(Treatment, Site, Position, Year) %>% 
  mutate(ID = cur_group_id(), .before = 'Woody_Cover') %>% 
  ungroup

#Now I need to remove unwanted columns. 

DiversityData <- dplyr::select(Data, - c(1:4,6:12))
DiversityData

# Creating metadata info 
Meta <- Data %>%
  dplyr::select(Treatment, Site, Position, Year, ID)

# These next 2 code chunks calculate Richness and Shannon Diversity 
# I am simultaneously adding the new values as columns onto the DiversityData 

DiversityData <- DiversityData %>%
  mutate(Richness=specnumber(DiversityData))

DiversityData <- DiversityData %>%
  mutate(Diversity=diversity(DiversityData))

# Now merge all data together based on ID 
DivData <- inner_join(DiversityData, Meta)  

############ Visualize the Richness Data ##### Decide on Outliers NOW ##########

ggplot(DivData, aes(x=Richness)) +
  geom_histogram()

ggplot(DivData, aes(x=Treatment, y=Richness, fill = Year)) + 
  geom_boxplot(notch=TRUE)

################# Full Model Testing - Richness ################################

# Because Richness is count data, it is a Poisson distribution 
# The richness values I have are actually not very zero inflated 

?glmmTMB::family_glmmTMB # Family options 

RichMod1 <- glmmTMB(Richness ~ Treatment * Year + (1|Site:Position), 
                    data=DivData, family=poisson)

# Now look at the residuals 
plot(simulateResiduals(RichMod1))
hist(simulateResiduals(RichMod1)) ## histogram should be flat
qqPlot(resid(RichMod1))  ## residuals should line up pretty closely to the blue line
hist(residuals(RichMod1))
plot(fitted(RichMod1), residuals(RichMod1))

############################## Model Results ###################################

Anova(RichMod1, type="III")

# Pairwise Treatment x Year 
PairsRichMod <- emmeans(RichMod1, ~Treatment*Year, type='response') 
pairs(PairsRichMod)
CI_Letters_Rich <- cld(PairsRichMod, Letters=letters, sort=TRUE, decreasing=TRUE)

############# Means Table + Visualization  ##################################### 

RichTable <- summarySE(DivData, measurevar="Richness", groupvars=c("Treatment", "Year")) 

RichTable <- RichTable %>%
  mutate(Treatment = fct_relevel(Treatment, "Control", "Reference", "Cleared","Scraped", "Diaspore", "Cleared + Burned",
                                 "Cleared + Scraped", "Cleared + Diaspore","Scraped + Diaspore", "Cleared + Scraped + Diaspore",
                                 "Cleared + Diaspore + Burned","Cleared + Scraped + Diaspore + Burned"))

# TN_Summary <- summarySE(TNData, measurevar="Total_N", groupvars=c("Treatment", "Year")) %>%
# mutate(Time = case_when(Year == '2022' ~ '2-Weeks Post Treatment (2022)',
# Year == '2023' ~ '1-Year Post Treatment (2023)',
# Year == '2024' ~ '2-Years Post Treatment (2024)')) %>%
# mutate(Time = fct_relevel(Time, "2-Weeks Post Treatment (2022)", "1-Year Post Treatment (2023)",
# "2-Years Post Treatment (2024)"))


RichFig <- ggplot(RichTable, aes(x=Treatment, y=Richness, group = Year))+ 
  geom_errorbar(aes(ymin=Richness-ci, ymax=Richness+ci), width=.2, 
                position=position_dodge(0.5)) +
  geom_point(aes(x=Treatment, y=Richness), size=3, 
             position=position_dodge(0.5)) +
  labs(x="Treatment", y = "Plant Species Richness")+
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


################################################################################
############################# Shannon Diversity Data ###########################
############ Visualize the Shannon Diversity Data ##### Decide on Outliers NOW #
ggplot(DivData, aes(x=Diversity)) +
  geom_histogram()

ggplot(DivData, aes(x=Treatment, y=Diversity, fill = Year)) + 
  geom_boxplot(notch=TRUE)

######################### Diversity Model Testing ##############################

DivMod1 <- glmmTMB(Diversity ~ Treatment * Year + (1|Site:Position), 
                   data=DivData, family=gaussian)

# Now look at the residuals 
plot(simulateResiduals(DivMod1))
hist(simulateResiduals(DivMod1)) ## histogram should be flat
qqPlot(resid(DivMod1))  ## residuals should line up pretty closely to the blue line
hist(residuals(DivMod1))
plot(fitted(DivMod1), residuals(DivMod1))

# Not Perfect but not Terrible
# Going to go with this model for now, but may be able to improve it later 

###################### Getting Diversity Model Results #########################

Anova(DivMod1, type="III")

PairsDivMod <- emmeans(DivMod1, ~Treatment*Year, type='response') 
#options(max.print=5000)
pairs(PairsDivMod)
CI_Letters_DivYr <- cld(PairsDivMod, Letters=letters, sort=TRUE, decreasing=TRUE)

############# Diversity Means Table + Visualization  ##################################### 

DivTable <- summarySE(DivData, measurevar="Diversity", groupvars=c("Treatment", "Year")) 

DivTable <- DivTable %>%
  mutate(Treatment = fct_relevel(Treatment, "Control", "Reference", "Cleared","Scraped", "Diaspore", "Cleared + Burned",
                                 "Cleared + Scraped", "Cleared + Diaspore","Scraped + Diaspore", "Cleared + Scraped + Diaspore",
                                 "Cleared + Diaspore + Burned","Cleared + Scraped + Diaspore + Burned"))

# TN_Summary <- summarySE(TNData, measurevar="Total_N", groupvars=c("Treatment", "Year")) %>%
# mutate(Time = case_when(Year == '2022' ~ '2-Weeks Post Treatment (2022)',
# Year == '2023' ~ '1-Year Post Treatment (2023)',
# Year == '2024' ~ '2-Years Post Treatment (2024)')) %>%
# mutate(Time = fct_relevel(Time, "2-Weeks Post Treatment (2022)", "1-Year Post Treatment (2023)",
# "2-Years Post Treatment (2024)"))


DivFig <- ggplot(DivTable, aes(x=Treatment, y=Diversity, group = Year))+ 
  geom_errorbar(aes(ymin=Diversity-ci, ymax=Diversity+ci), width=.2, 
                position=position_dodge(0.5)) +
  geom_point(aes(x=Treatment, y=Diversity), size=3, 
             position=position_dodge(0.5)) +
  labs(x="Treatment", y = "Shannon Diversity")+
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


