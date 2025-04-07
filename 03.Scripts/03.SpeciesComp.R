################################################################################
######################## Ashlynn Smith #########################################
################## Script to Visualize and Analyze #############################
############# Plant Species Composition Data : EPA Project  ####################

install.packages('ggvegan')


library(pals)
library(ggsci)
library(wesanderson)
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
#library(ggvegan)
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
library(Rmisc)
library(RColorBrewer)

################ 2024 ONLY  ######################################################
# The only assumption of perMANONVA is that all samples included are independent #
# Because of the repeated measures nature (sampling the same plots year after year) #
# I decided to analyze just the final year of sampling. For publication, it might # 
# add to the paper by analyzing each year separately, but for now, I'm only going to # 
# look at 2024. I have analyzed 2023 data, but included position as a factor. I will # 
# Not look at slope position as a factor here, only Treatment. 

########################### Data Prep ###########################################

Data <- read_csv("02.CleanData/CleanData_AnnualMeans.csv",
                 col_types = cols(Year = col_factor(),
                                  Position = col_factor(),
                                  Treatment = col_factor(),
                                  Site = col_factor()))

CompData <- Data %>%
  filter(Year != "2022") %>%
  filter(Year != "2023")

# Saved this file for creating a means table for Publication
write_csv(CompData, "02.CleanData/CompData_2024.csv")

# Create a UniqueID 
CompData <- CompData %>%
  group_by(Treatment, Site, Position, Year) %>% 
  mutate(ID = cur_group_id(), .before = 'Woody_Cover') %>% 
  ungroup

# Create a Sum column to find out which rows have no species in them 
CompData <- CompData %>%
  mutate(Sum = rowSums(.[13:183]))

# Delete those rows 
CompData <- CompData %>%
  filter(Sum != 0)

# For this data set, no rows are removed 

# Delete Sum column 
CompData <- CompData %>%
  dplyr::select(-c(184))

# Change column ARISTR to ARIBEY 
CompData <- CompData %>%
  dplyr::rename(ARIBEY=ARISTR) # Have to specify dplyr here because there are two 
                               # packages installed that have 'rename' as a function 


#Now I need to remove additional unwanted columns. 

DiversityData <- dplyr::select(CompData, - c(1:4,6:12))
DiversityData

# Creating metadata info 
MetaData <- CompData %>%
  dplyr::select(Treatment, Site, Position, Year, ID)

MetaData$ID <- as.character(MetaData$ID)

#This step is required to formally make the Data into the format needed for NMDS 
DistanceMatrix <- DiversityData %>%
  column_to_rownames("ID")

#DistanceMtx2 is the TRUE distance Matrix (Actually a Dissimilarity Matrix) 
set.seed(1986)
DistanceMtx2 <- vegdist(DistanceMatrix, method = 'bray', na.rm=TRUE)
DistanceMtx2

all.sites <- DistanceMtx2

trt <- MetaData$Treatment

# pos <- MetaData$Position

set.seed(1986)
all.mds <- metaMDS(all.sites)
all.mds

data.scores <- as.data.frame(scores(all.mds)) 
data.scores$ID <- rownames(data.scores)  
data.scores <- inner_join(data.scores, MetaData) 

adon.results<-adonis2(all.sites ~ trt, method="bray",perm=999)
print(adon.results)

# Attempt at running a pairwise comparison of perMANOVA results 

library(devtools)
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)

# When installing (04.06.2025), I did install with updates 

Pair_Results_TRT <- pairwise.adonis(all.sites, MetaData$Treatment)

P_TRT <- Pair_Results_TRT$p.value

Pair_Results_TRT <- Pair_Results_TRT %>%
  dplyr::mutate(BH_p = p.adjust(Pair_Results_TRT$p.value, method = "BH"))

# BH adjustment = Benjamini & Hochberg (1995) 
# The BH adjustment, along with every other adjustment method is less conservative than 
# Bonforonni 

########################## Don't think this code is Needed ######################
#dis <- vegdist(all.sites)
#mod <- betadisper(dis, trt)
#mod

#centroids<-data.frame(grps=rownames(mod$centroids),data.frame(mod$centroids))
#vectors<-data.frame(group=mod$group,data.frame(mod$vectors))

#seg.data<-cbind(vectors[,1:3],centroids[rep(1:nrow(centroids),as.data.frame(table(vectors$group))$Freq),2:3])
#names(seg.data)<-c("group","v.PCoA1","v.PCoA2","PCoA1","PCoA2")

################################################################################

#Fit vector to ordination  
Fit <- envfit(all.mds, DistanceMatrix)
Fit

autoplot(Fit)

arrows <- data.frame(Fit$vector$arrows,R=Fit$vectors$r,P=Fit$vectors$pvals)
arrows

arrows$Species <- rownames(arrows) 
arrows.p <- arrows[arrows$P<=0.001,] #to see arrows with Species, set to 0.001
arrows.p

############### Visualizing the Results ########################################

data.scores$Treatment <- factor(data.scores$Treatment, levels = c("Reference",
                                                                  "Cleared + Scraped + Diaspore + Burned",
                                                                  "Cleared + Scraped + Diaspore",
                                                                  "Cleared + Diaspore + Burned",
                                                                  "Scraped + Diaspore",
                                                                  "Cleared + Burned",
                                                                  "Cleared + Diaspore",
                                                                  "Cleared + Scraped",
                                                                  "Cleared","Scraped","Diaspore",
                                                                  "Control"))

PLOT_TRT <- ggplot(data=data.scores) + 
  stat_ellipse(aes(x=NMDS1,y=NMDS2,colour=Treatment),level = 0.35, size=1) +
  geom_point(aes(x=NMDS1,y=NMDS2,colour=Treatment),size=2) + 
  #stat_summary(aes(x=NMDS1,y=NMDS2,fun = "mean"), geom="point", shape=17, size=4, colour=Treatment) +
  ggrepel::geom_text_repel(data=arrows.p, aes(x = NMDS1/.8, y = NMDS2/.8, label = Species), 
                           size=4) +
  scale_color_manual(values=as.vector(cols25(12))) +
  theme_classic2()
























