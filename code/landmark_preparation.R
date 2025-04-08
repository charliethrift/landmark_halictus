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

tps_data_18feb <- readland.tps("../data/CT_7apr_HligLM_18feb.tps",
                              specID = c("imageID"), readcurves = FALSE, 
                              warnmsg = TRUE,negNA = FALSE)
tps_data_27feb <- readland.tps("../data/CT_7apr_HligLM_27feb.tps",
                               specID = c("imageID"), readcurves = FALSE, 
                               warnmsg = TRUE,negNA = FALSE)
tps_data_4mar <- readland.tps("../data/CT_7apr_HligLM_4mar.tps",
                              specID = c("imageID"), readcurves = FALSE, 
                              warnmsg = TRUE,negNA = FALSE)
tps_data_17mar <- readland.tps("../data/CT_7apr_HligLM_17mar.tps",
                               specID = c("imageID"), readcurves = FALSE, 
                               warnmsg = TRUE,negNA = FALSE)

## read bee occurrence csv data
occurrence_data <- read_csv("../data/halictus_7apr25_occurrences.csv")
occurrence_data <- occurrence_data %>% 
  select("catalogNumber", "scientificName","year","sex","stateProvince",
         "decimalLatitude","decimalLongitude")


## Generalized Procrustes Analysis (GPA)
# Run GPA on TPS data, to then generate coordinate data for each landmark 
# on each specimen. In this code, we start with tps data and end with a CSV 
# of coordinate data for each specimen. This can then be merged with the bee 
# data generated above.
gpa_tps_data_18feb <- gpagen(tps_data_18feb, PrinAxes = TRUE)
gpa_tps_data_27feb <- gpagen(tps_data_27feb, PrinAxes = TRUE)
gpa_tps_data_4mar <- gpagen(tps_data_4mar, PrinAxes = TRUE)
gpa_tps_data_17mar <- gpagen(tps_data_17mar, PrinAxes = TRUE)

# write csv files with the landmark coordinate information
write.csv(gpa_tps_data_18feb$coords, "../data/gpa_18feb.csv", row.names = TRUE)
write.csv(gpa_tps_data_27feb$coords, "../data/gpa_27feb.csv", row.names = TRUE)
write.csv(gpa_tps_data_4mar$coords, "../data/gpa_4mar.csv", row.names = TRUE)
write.csv(gpa_tps_data_17mar$coords, "../data/gpa_17mar.csv", row.names = TRUE)

# clean up landmark data into a tidy dataframe
lmdata_18feb <- read.csv("../data/gpa_18feb.csv") # image name with path
lmdata_27feb <- read.csv("../data/gpa_27feb.csv") # image name only
lmdata_4mar <- read.csv("../data/gpa_4mar.csv") # image name with path
lmdata_17mar <- read.csv("../data/gpa_17mar.csv") # image name with path



## first, clean up 27feb, because this one doesn't contain image paths
df_transpose = t(lmdata_27feb) #transpose the landmark data
df_transpose1 <- df_transpose[-1,]
lm1 <- df_transpose1
lmX <- lm1
lmY <- lm1
lmX1 <- lmX
lmX1 <- data.frame(lmX1)
lmX2 <- setNames(cbind(rownames(lmX1), lmX1, row.names = NULL),
                 c("name", "LM1x", "LM2x", "LM3x",
                   "LM4x", "LM5x", "LM6x", "LM7x", "LM8x", "LM9x",
                   "LM10x", "LM11x", "LM12x", "LM13x", "LM14x", "LM15x",
                   "LM16x", "LM17x", "LM18x", "LM19x", "LM20x", "LM21x",
                   "LM22x", "LM23x", "LM24x", "LM25x"))
lmXonly <- lmX2[str_detect(lmX2$name, "X.UCSB"), ]
#repeat with Y
lmY1 <- lmY
lmY1 <- data.frame(lmY1)
lmY2 <- setNames(cbind(rownames(lmY1), lmY1, row.names = NULL),
                 c("name", "LM1y", "LM2y", "LM3y",
                   "LM4y", "LM5y", "LM6y", "LM7y", "LM8y", "LM9y",
                   "LM10y", "LM11y", "LM12y", "LM13y", "LM14y", "LM15y",
                   "LM16y", "LM17y", "LM18y", "LM19y", "LM20y", "LM21y",
                   "LM22y", "LM23y", "LM24y", "LM25y"))
lmYonly <- lmY2[str_detect(lmY2$name, "Y.UCSB"), ]
lmXonly1 <- lmXonly
lmYonly1 <- lmYonly
lmXonly2 <- lmXonly1 %>% 
  tidyr::separate(name,                      
                  c("X","UCSB", "barcode","wingSide"), extra='drop') %>%
  tidyr::unite('catalogNumber', c('UCSB','barcode')) 
lmXonly3 <- lmXonly2[,-1]
##now: lmXonly5 has x coordinate values for all 25 landmarks, and just the catalogNumber
###repeat for Y
lmYonly2 <- lmYonly1 %>% 
  tidyr::separate(name,                      
                  c("Y","UCSB", "barcode", "wingSide"), extra='drop') %>%
  tidyr::unite('catalogNumber', c('UCSB','barcode')) 
#drop any wings that were Right instead of Left
lmYonly3 <- lmYonly2[,-1]
#Now: unite the Y and X coordinate dataframes into just one
lmXonly3 <- lmXonly3 %>%  tidyr::unite('wingID',c('catalogNumber','wingSide'))
lmYonly3 <- lmYonly3 %>%  tidyr::unite('wingID',c('catalogNumber','wingSide'))

lm_both <- merge(lmXonly3, lmYonly3, by=c("wingID"))
lm_clean_27feb <- lm_both

## second, clean up 18feb, does have image paths
df_transpose1 = t(lmdata_18feb)
df_transpose1 <- df_transpose1[-1,]
lm1 <- df_transpose1
lmX <- lm1
lmY <- lm1
lmX1 <- lmX
lmX1 <- data.frame(lmX1)
lmX2 <- setNames(cbind(rownames(lmX1), lmX1, row.names = NULL),
                 c("name", "LM1x", "LM2x", "LM3x",
                   "LM4x", "LM5x", "LM6x", "LM7x", "LM8x", "LM9x",
                   "LM10x", "LM11x", "LM12x", "LM13x", "LM14x", "LM15x",
                   "LM16x", "LM17x", "LM18x", "LM19x", "LM20x", "LM21x",
                   "LM22x", "LM23x", "LM24x", "LM25x"))
lmXonly <- lmX2[str_detect(lmX2$name, "X.C"), ]
#repeat with Y
lmY1 <- lmY
lmY1 <- data.frame(lmY1)
lmY2 <- setNames(cbind(rownames(lmY1), lmY1, row.names = NULL),
                 c("name", "LM1y", "LM2y", "LM3y",
                   "LM4y", "LM5y", "LM6y", "LM7y", "LM8y", "LM9y",
                   "LM10y", "LM11y", "LM12y", "LM13y", "LM14y", "LM15y",
                   "LM16y", "LM17y", "LM18y", "LM19y", "LM20y", "LM21y",
                   "LM22y", "LM23y", "LM24y", "LM25y"))
lmYonly <- lmY2[str_detect(lmY2$name, "Y.C"), ]
lmXonly1 <- lmXonly
lmYonly1 <- lmYonly
lmXonly2 <- lmXonly1 %>% 
  tidyr::separate(name,                      
                  c("1",'2','3','4','5','6','7','8',
                    '9','10','11','12','13','14','15'), extra='drop') %>%
  tidyr::unite('catalogNumber', c('13','14')) 

lmXonly3 <- lmXonly2[,-c(1:12)]
lmXonly3 <- lmXonly3 %>%  tidyr::unite('wingID',c('catalogNumber','15'))
###repeat for Y
lmYonly2 <- lmYonly1 %>% 
  tidyr::separate(name,                      
                  c("1",'2','3','4','5','6','7','8',
                    '9','10','11','12','13','14','15'), extra='drop') %>%
  tidyr::unite('catalogNumber', c('13','14')) 
lmYonly3 <- lmYonly2[,-c(1:12)]
lmYonly3 <- lmYonly3 %>%  tidyr::unite('wingID',c('catalogNumber','15'))

#Now: unite the Y and X coordinate dataframes into just one
lm_both <- merge(lmXonly3, lmYonly3, by=c("wingID"))
lm_clean_18feb <- lm_both


## third, clean up 4mar, does have image paths
df_transpose1 = t(lmdata_4mar)
df_transpose1 <- df_transpose1[-1,]
lm1 <- df_transpose1
lmX <- lm1
lmY <- lm1
lmX1 <- lmX
lmX1 <- data.frame(lmX1)
lmX2 <- setNames(cbind(rownames(lmX1), lmX1, row.names = NULL),
                 c("name", "LM1x", "LM2x", "LM3x",
                   "LM4x", "LM5x", "LM6x", "LM7x", "LM8x", "LM9x",
                   "LM10x", "LM11x", "LM12x", "LM13x", "LM14x", "LM15x",
                   "LM16x", "LM17x", "LM18x", "LM19x", "LM20x", "LM21x",
                   "LM22x", "LM23x", "LM24x", "LM25x"))
lmXonly <- lmX2[str_detect(lmX2$name, "X.C"), ]
#repeat with Y
lmY1 <- lmY
lmY1 <- data.frame(lmY1)
lmY2 <- setNames(cbind(rownames(lmY1), lmY1, row.names = NULL),
                 c("name", "LM1y", "LM2y", "LM3y",
                   "LM4y", "LM5y", "LM6y", "LM7y", "LM8y", "LM9y",
                   "LM10y", "LM11y", "LM12y", "LM13y", "LM14y", "LM15y",
                   "LM16y", "LM17y", "LM18y", "LM19y", "LM20y", "LM21y",
                   "LM22y", "LM23y", "LM24y", "LM25y"))
lmYonly <- lmY2[str_detect(lmY2$name, "Y.C"), ]
lmXonly1 <- lmXonly
lmYonly1 <- lmYonly
lmXonly2 <- lmXonly1 %>% 
  tidyr::separate(name,                      
                  c("1",'2','3','4','5','6','7','8',
                    '9','10','11','12'), extra='drop') %>%
  tidyr::unite('catalogNumber', c('10','11')) 

lmXonly3 <- lmXonly2[,-c(1:9)]
lmXonly3 <- lmXonly3 %>%  tidyr::unite('wingID',c('catalogNumber','12'))
###repeat for Y
lmYonly2 <- lmYonly1 %>% 
  tidyr::separate(name,                      
                  c("1",'2','3','4','5','6','7','8',
                    '9','10','11','12'), extra='drop') %>%
  tidyr::unite('catalogNumber', c('10','11')) 
lmYonly3 <- lmYonly2[,-c(1:9)]
lmYonly3 <- lmYonly3 %>%  tidyr::unite('wingID',c('catalogNumber','12'))

#Now: unite the Y and X coordinate dataframes into just one
lm_both <- merge(lmXonly3, lmYonly3, by=c("wingID"))
lm_clean_4mar <- lm_both



#fourth, repeat for 17 march, with image path
df_transpose1 = t(lmdata_17mar)
df_transpose1 <- df_transpose1[-1,]
lm1 <- df_transpose1
lmX <- lm1
lmY <- lm1
lmX1 <- lmX
lmX1 <- data.frame(lmX1)
lmX2 <- setNames(cbind(rownames(lmX1), lmX1, row.names = NULL),
                 c("name", "LM1x", "LM2x", "LM3x",
                   "LM4x", "LM5x", "LM6x", "LM7x", "LM8x", "LM9x",
                   "LM10x", "LM11x", "LM12x", "LM13x", "LM14x", "LM15x",
                   "LM16x", "LM17x", "LM18x", "LM19x", "LM20x", "LM21x",
                   "LM22x", "LM23x", "LM24x", "LM25x"))
lmXonly <- lmX2[str_detect(lmX2$name, "X.C"), ]
#repeat with Y
lmY1 <- lmY
lmY1 <- data.frame(lmY1)
lmY2 <- setNames(cbind(rownames(lmY1), lmY1, row.names = NULL),
                 c("name", "LM1y", "LM2y", "LM3y",
                   "LM4y", "LM5y", "LM6y", "LM7y", "LM8y", "LM9y",
                   "LM10y", "LM11y", "LM12y", "LM13y", "LM14y", "LM15y",
                   "LM16y", "LM17y", "LM18y", "LM19y", "LM20y", "LM21y",
                   "LM22y", "LM23y", "LM24y", "LM25y"))
lmYonly <- lmY2[str_detect(lmY2$name, "Y.C"), ]
lmXonly1 <- lmXonly
lmYonly1 <- lmYonly
lmXonly2 <- lmXonly1 %>% 
  tidyr::separate(name,                      
                  c("1",'2','3','4','5','6','7','8',
                    '9','10','11','12'), extra='drop') %>%
  tidyr::unite('catalogNumber', c('10','11')) 

lmXonly3 <- lmXonly2[,-c(1:9)]
lmXonly3 <- lmXonly3 %>%  tidyr::unite('wingID',c('catalogNumber','12'))
###repeat for Y
lmYonly2 <- lmYonly1 %>% 
  tidyr::separate(name,                      
                  c("1",'2','3','4','5','6','7','8',
                    '9','10','11','12'), extra='drop') %>%
  tidyr::unite('catalogNumber', c('10','11')) 
lmYonly3 <- lmYonly2[,-c(1:9)]
lmYonly3 <- lmYonly3 %>%  tidyr::unite('wingID',c('catalogNumber','12'))

#Now: unite the Y and X coordinate dataframes into just one
lm_both <- merge(lmXonly3, lmYonly3, by=c("wingID"))
lm_clean_17mar <- lm_both




### now, combine the 4 dataframes of cleaned coordinates we have
lm_all <- rbind(lm_clean_18feb,lm_clean_27feb,lm_clean_4mar,lm_clean_17mar)

## Merge Landmark Data with Bee Occurrence Data (SECOND, perform merge)
lm_all_cat <- lm_all
lm_all_cat$catalogNumber_wing <- lm_all_cat$wingID
lm_all_cat <- lm_all_cat %>% 
  tidyr::separate(catalogNumber_wing,                      
                  c("UCSB","number","wing"))
lm_all_cat <- lm_all_cat %>% 
  select(-wing)
lm_all_cat <- lm_all_cat %>% tidyr::unite('catalogNumber',c('UCSB','number'))
lm_all_cat <- lm_all_cat %>% relocate(catalogNumber, .before = wingID)

write.csv(lm_all_cat, "../data/cleaned_lm_data_no_meta.csv")




occurrence_data <- occurrence_data %>% separate(catalogNumber, c("ucsb","number"))
occurrence_data <- occurrence_data %>% tidyr::unite("catalogNumber", c("ucsb","number"))


lm_data <- merge(lm_all_cat, occurrence_data, by=c("catalogNumber"))
lm_data <- lm_data %>% relocate(53:58, .before = wingID)




write.csv(lm_data, "../data/lm_data_7april2025.csv")






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


