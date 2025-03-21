

################################################################################
######################## Ashlynn Smith #########################################
################## Script to Clean and Visualize ###############################
################# Vegetation Data : EPA Project  ###############################

# Load Packages 
library(tidyverse)


# Call in Raw Data 
Data <- read_csv("01.RawData/RawVegData.csv",
                 col_types = cols(SeasonYear = col_factor(),
                                  Position = col_factor(),
                                  Treatment = col_factor(),
                                  SampleUnitID = col_factor(),
                                  Site = col_factor(),
                                  Event = col_factor()))


# Based on a preliminary look at the data, the most complete and consistent data 
# is the 'Spring' collected data. So I will remove the Fall and Summer collection 
# events from the dataset 

# Step 1: 
Data2 <- Data %>%
  filter(Season == "Spring")

# Check data replication 
# The column 'n' in the data table produced below represents the replication for that 
# year of particular treatment for each slope position (i.e. n = the number of sites) 

Replication <- Data2 %>%
  group_by(Year, Treatment, Position) %>% 
  tally()
 
# Treatments with not enough replication (n < or = 2) at a given slope position need to be removed. 
# Removed treatments include 
# Burned, 
# Scraped + Burned, 
# Diaspore + Burned, 
# Scraped + Diaspore + Burned

# Additionally, the Cleared + Scraped + Burned Treatment needed to be removed because it is not a true treatment 
# In the field, fire would not carry in plots that were Cleared + Scraped without 
# Diaspore present to carry the fire. This treatment would really be just additional replication of 
# Cleared + Scraped, so it needs to be taken out. 

############## Remove Treatments ###############################################

# The code below basically says "remove the rows where Treatment equals this value." 

Data3 <- Data2 %>%
  filter(Treatment != "Burned") %>%
  filter(Treatment != "Scraped + Burned") %>%
  filter(Treatment != "Diaspore + Burned") %>%
  filter(Treatment != "Cleared + Scraped + Burned") %>%
  filter(Treatment != "Scraped + Diaspore + Burned")

# Also, there are n=12 Reference Sampling Plots right now. 
# To get the data closer to = replication, I need to remove 2 transects (6 plots)
# From each reference prairie. I talked with CC about this, and we decided to keep 
# Transect #2 for each Reference Site because sampling fatigue would be most consistent 
# In the middle transect of Reference Sites 

################# Remove Transects 1 and 3 from Reference Prairies #############

Data4 <- Data3 %>%
  filter(SampleUnitID != "RP02XXU1") %>%
  filter(SampleUnitID != "RP02XXM1") %>%
  filter(SampleUnitID != "RP02XXB1") %>%
  filter(SampleUnitID != "RP02XXU3") %>%
  filter(SampleUnitID != "RP02XXM3") %>%
  filter(SampleUnitID != "RP02XXB3") %>%
  filter(SampleUnitID != "RP84U1") %>%
  filter(SampleUnitID != "RP84M1") %>%
  filter(SampleUnitID != "RP84B1") %>%
  filter(SampleUnitID != "RP84U3") %>%
  filter(SampleUnitID != "RP84M3") %>%
  filter(SampleUnitID != "RP84B3") %>%
  filter(SampleUnitID != "RP85U1") %>%
  filter(SampleUnitID != "RP85M1") %>%
  filter(SampleUnitID != "RP85B1") %>%
  filter(SampleUnitID != "RP85U3") %>%
  filter(SampleUnitID != "RP85M3") %>%
  filter(SampleUnitID != "RP85B3") %>%
  filter(SampleUnitID != "RP42U1") %>%
  filter(SampleUnitID != "RP42M1") %>%
  filter(SampleUnitID != "RP42B1") %>%
  filter(SampleUnitID != "RP42U3") %>%
  filter(SampleUnitID != "RP42M3") %>%
  filter(SampleUnitID != "RP42B3")

# Now check that this worked! 
# When sorted and summed this way, the smallest number (for TRUE replication) 
# should be 3, representing 3 different sites

TEST <- Data4 %>%
  group_by(Year, Treatment, Position) %>%
  tally() 

# This checks out after running the code above and looking at the data 

summary(Data4)

# Write Data4 final dataset for use in other code scripts later 

write_csv(Data4, "02.CleanData/CleanData_SpringOnly.csv")


# The code below is no longer needed because I removed Fall / Summer data rather than 
# averaging across all seasons for each year 


# Now average across all columns (except this does not work for Woody_Height) 
# Based on Treatment, Position, Site, and Year 
# What this does is averages across seasons and provides a single mean 
# value for all variables for each year  

#Data5 <- Data4 %>%
  #group_by(Treatment, Position, Site, Year) %>%
  #summarise_at(vars(Woody_Cover:WISFRU), mean, na.rm = TRUE)




