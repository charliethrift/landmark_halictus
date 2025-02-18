# Script: read landmarks and occurrences, create csv of landmark values
# Code: Wing Venation and Asymmetry in Halictus ligatus Bees
# Charles Thrift
# 18 February 2025

## load packages
library(tidyverse)
library(geomorph) #NOTE: download XQuartz if using a Mac
library(devtools)
library(RRPP) #used for MANOVA 
library(adegenet) # used for DAPC


## read TPS data
# The TPS data are in one file, which was created using tpsUtil32. Here, 
# we read these data and specify a few parameters for the _readland.tps_ function. 
# There are no curves to be read (readcurves=FALSE), and any negative values for 
# landmarks are negative in Cartesian space (negNA=FALSE). If negNA were true, 
# the function would read any negative landmarks as a missing landmark (this 
# is because while landmarking, if a character is obscured you can "skip" 
# landmarks within the landmarking software. We do not include any specimens 
# missing any landmarks in this analysis).
all_tps_data <- readland.tps("../data/h.ligatus_18feb25.TPS",
                             specID = c("ID"), readcurves = FALSE, 
                             warnmsg = TRUE,negNA = FALSE)

all_tps_data <- readland.tps("../data/ligatus_test.tps",
                             specID = c("imageID"), readcurves = FALSE, 
                             warnmsg = TRUE,negNA = FALSE)

## read bee occurrence csv data
occurrence_data <- read_csv("beedata26jun22.csv")


## Generalized Procrustes Analaysis (GPA)
# Run GPA on TPS data, to then generate coordinate data for each landmark 
# on each specimen. In this code, we start with tps data and end with a CSV 
# of coordinate data for each specimen. This can then be merged with the bee 
# data generated above.
all_tps_gpa <- gpagen(all_tps_data, PrinAxes = TRUE)
# write csv files with the landmark coordinate information and centroid size
## ????? maybe replace with a single file that includes both coords and Csize??????
write.csv(all_tps_gpa$coords, "../data/all_tps_test.csv", row.names = TRUE)
write.csv(all_tps_gpa$Csize, "all_Csize_tps_21mar23.csv", row.names = TRUE)




## Merge Landmark Data with Bee Occurence Data (FIRST, clean up files)
# Note, this all looks pretty messy, and could maybe be cleaned up 
# or done more simply with fewer lines of code. Let's look at this later
lmdata <- read.csv("../data/all_tps_test.csv") #read landmark data in
Csize <- read.csv("all_Csize_tps_21mar23.csv") #read Csize data in (note, not using it in this analysis)
Csize1 <- Csize 
Csize2 <- setNames(cbind(rownames(Csize1), Csize1, row.names = NULL), 
                   c("number_delete", "specimenID", "Csize")) #add row names
Csize3 <- Csize2[,c(2:3)] #delete the first column (which is just numbering)
df_transpose = t(lmdata) #transpose the landmark data
df_transpose1 <- df_transpose
df_transpose1 <- df_transpose1[-1,]
#Format of LM data is currently two rows per specimen
#with one row being X coordinate values and one row
#being Y coordinate values. Below, we split into two 
#data frames and then stitch them back together to get
#18 different variables for the 9 landmarks
#(9 X coordinates and 9 Y coordinates)
lm1 <- df_transpose1
lmX <- lm1
lmY <- lm1
lmX1 <- lmX
lmX1 <- data.frame(lmX1)
lmX2 <- setNames(cbind(rownames(lmX1), lmX1, row.names = NULL),
                 c("name", "LM1x", "LM2x", "LM3x",
                   "LM4x", "LM5x", "LM6x", "LM7x", "LM8x", "LM9x"))
lmXonly <- lmX2[str_detect(lmX2$name, "X.UCSB"), ]
#repeat with Y
lmY1 <- lmY
lmY1 <- data.frame(lmY1)
lmY2 <- setNames(cbind(rownames(lmY1), lmY1, row.names = NULL),
                 c("name", "LM1y", "LM2y", "LM3y",
                   "LM4y", "LM5y", "LM6y", "LM7y", "LM8y", "LM9y"))
lmYonly <- lmY2[str_detect(lmY2$name, "Y.UCSB"), ]
lmXonly1 <- lmXonly
lmYonly1 <- lmYonly
lmXonly2 <- lmXonly1 %>% 
  tidyr::separate(name,                      
                  c("X","UCSB", "barcode", "wing","species", 
                    "location", "wingSide"), extra='drop') %>%
  tidyr::unite('catalogNumber', c('UCSB','barcode')) 
#drop any wings that were Right instead of Left
lmXonly3 <- lmXonly2[lmXonly2$wingSide %in% c("ed", NA), ] #remove any "right" wings
lmXonly4 <- lmXonly3[lmXonly3$species %in% c("edited", "far",
                                             "lig", "tri"), ] #remove any "right" wings
lmXonly5 <- lmXonly4[,c(2,7:15)]
##now: lmXonly5 has x coordinate values for all 9 landmarks, and just the catalogNumber
###repeat for Y
lmYonly2 <- lmYonly1 %>% 
  tidyr::separate(name,                      
                  c("Y","UCSB", "barcode", "wing","species", 
                    "location", "wingSide"), extra='drop') %>%
  tidyr::unite('catalogNumber', c('UCSB','barcode')) 
#drop any wings that were Right instead of Left
lmYonly3 <- lmYonly2[lmYonly2$wingSide %in% c("ed", NA), ] #remove any "right" wings
lmYonly4 <- lmYonly3[lmYonly3$species %in% c("edited", "far",
                                             "lig", "tri"), ] #remove any "right" wings
lmYonly5 <- lmYonly4[,c(2,7:15)]
#Now: unite the Y and X coordinate dataframes into just one
lm_both <- merge(lmXonly5, lmYonly5, by=c("catalogNumber"))
####Final step: add in the Csize column
Csize4 <- Csize3
Csize5 <- Csize4 %>% 
  tidyr::separate(specimenID,                      
                  c("UCSB", "barcode", "wing","species", 
                    "location", "wingSide"), extra='drop') %>%
  tidyr::unite('catalogNumber', c('UCSB','barcode')) 
Csize6 <- Csize5[,c(1,6)]
lm_both_size <- merge(lm_both, Csize6, by=c("catalogNumber"))
####Done with landmark data. "lm" has each specimen and 18 variables for lm coordinates
####plus 1 variable for Csize




## Merge Landmark Data with Bee Occurence Data (SECOND, perform merge)
lm_data <- merge(lm_both_size, df1, by=c("catalogNumber"))



## Summarizing Specimens by Location
lm <- lm_data
lm$species[lm$scientificName == 
             "Halictus ligatus Say, 1837"]<- "H. ligatus"

lm$location <- NA #store NA values first, then populate below with the right group
# !! REPLACE THESE BOUNDING BOXES WITH THE RIGHT COORDS FOR GEORGIA AND CALIFORNIA
lm$location[lm$decimalLatitude <= 34.234 & 
              lm$decimalLatitude >= 33.86 & 
              lm$decimalLongitude >= -120.05 & 
              lm$decimalLongitude <= -119.45]<- "Santa Cruz Island"
lm$location[lm$decimalLatitude <= 36 & 
              lm$decimalLatitude >= 34.113887 & 
              lm$decimalLongitude >= -121 & 
              lm$decimalLongitude <= -116]<- "Mainland"
lm$location <- as.factor(lm$location)
summary(lm$location)
lm <- subset(lm, location != "NA") #if it didn't get a population assigned, drop it from pop analysis
## Check count of species and locations
table(lm$species, lm$location)


# Export this landmark data as a CSV for future analyses
write.csv(lm, "data/landmarks_DATE.csv")


