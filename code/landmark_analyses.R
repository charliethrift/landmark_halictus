# Script: run analyses on landmark data
# Code: Wing Venation and Asymmetry in Halictus ligatus Bees
# Charles Thrift
# 18 February 2025

## load libraries
library(tidyverse)
library(geomorph) #NOTE: download XQuartz if using a Mac
library(devtools)
library(RRPP) #used for MANOVA 
library(adegenet) # used for DAPC

## read data
lm_data <- read_csv("../data/lm_data_7april2025.csv")
lm_data <- lm_data[,-1]


# MAIN ANALYSES BELOW

# First, population level analyses between the two populations. 
## Are these wings different?

## Multivariate Analysis of Variance (MANOVA) for 2 populations
pop <- lm_data

manova_pop <- pop 
manova_pop <- manova_pop[,c(5,9:58)] #keep only numeric, plus State
manova_pop$stateProvince <- as.factor(manova_pop$stateProvince)
manova_pop_data <- manova_pop %>% select(-stateProvince) #removing the population column
manova_pop_data <- as.matrix(manova_pop_data) 

fit <- lm.rrpp(manova_pop_data ~ stateProvince, SS.type = "I", 
               data = manova_pop, print.progress = FALSE) #run the linear model fit
fitm <- manova.update(fit, print.progress = FALSE, tol = 0) #run the manova update
summary(fitm, test = "Pillai") #summarize with a manova table


## DAPC (rather than a PCA, run DAPC to visualize diff between groups)
pop1 <- lm_data
pop1 <- pop1[,c(5,9:58)]
pop2 <- pop1 %>% select(-stateProvince) # Only preserves LM values

dapc_pop <- dapc(pop2, 
                 grp = pop1$stateProvince,
                 n.pca = 15,
                 n.da = 10) #run DAPC on 2 populations

population_col <- c("darkblue", "#56B4E9") #save colors for population

scatter(dapc_pop, 
        bg="white", 
        cstar=0, scree.pca=FALSE, posi.pca="topleft",
        legend=TRUE, posi.leg = "topleft",
        col = population_col) #plot the two populations




#### Cross Validation for two populations
#pop1 is the full population and LM values
grp2 <- pop1$stateProvince #assign grouping to be population variable
#pop2 is only LM values #only preserve LM values

xval_pop <- xvalDapc(pop2, grp2, n.pca.max = 300, n.da = NULL,
                     training.set = 0.9,
                     result = c("groupMean","overall"), 
                     center = TRUE, scale = FALSE,
                     n.pca = NULL, n.rep = 30, xval.plot = TRUE) #run cross validation

xval_pop$`Number of PCs Achieving Highest Mean Success`#which number of PCs retained has highest accuracy?
xval_pop$`Mean Successful Assignment by Number of PCs of PCA`#what is accuracy at various PCs retained?





####Below, try to analyze differences between left and right wings
wing <- lm_data
wing <- wing[,c()]

ggplot()+
  geom_point(data = pop1, aes(x=LM1x, y=LM1y, color = stateProvince))+
  geom_point(data = pop1, aes(x=LM2x, y=LM2y, color = stateProvince))+
  geom_point(data = pop1, aes(x=LM3x, y=LM3y, color = stateProvince))+
  geom_point(data = pop1, aes(x=LM4x, y=LM4y, color = stateProvince))+
  geom_point(data = pop1, aes(x=LM5x, y=LM5y, color = stateProvince))+
  geom_point(data = pop1, aes(x=LM6x, y=LM6y, color = stateProvince))+
  geom_point(data = pop1, aes(x=LM7x, y=LM7y, color = stateProvince))+
  geom_point(data = pop1, aes(x=LM8x, y=LM8y, color = stateProvince))+
  geom_point(data = pop1, aes(x=LM9x, y=LM9y, color = stateProvince))+
  geom_point(data = pop1, aes(x=LM10x, y=LM10y, color = stateProvince))+
  geom_point(data = pop1, aes(x=LM11x, y=LM11y, color = stateProvince))+
  geom_point(data = pop1, aes(x=LM12x, y=LM12y, color = stateProvince))+
  geom_point(data = pop1, aes(x=LM13x, y=LM13y, color = stateProvince))+
  geom_point(data = pop1, aes(x=LM14x, y=LM14y, color = stateProvince))+
  geom_point(data = pop1, aes(x=LM15x, y=LM15y, color = stateProvince))+
  geom_point(data = pop1, aes(x=LM16x, y=LM16y, color = stateProvince))+
  geom_point(data = pop1, aes(x=LM17x, y=LM17y, color = stateProvince))+
  geom_point(data = pop1, aes(x=LM18x, y=LM18y, color = stateProvince))+
  geom_point(data = pop1, aes(x=LM19x, y=LM19y, color = stateProvince))+
  geom_point(data = pop1, aes(x=LM20x, y=LM20y, color = stateProvince))+
  geom_point(data = pop1, aes(x=LM21x, y=LM21y, color = stateProvince))+
  geom_point(data = pop1, aes(x=LM22x, y=LM22y, color = stateProvince))+
  geom_point(data = pop1, aes(x=LM23x, y=LM23y, color = stateProvince))+
  geom_point(data = pop1, aes(x=LM24x, y=LM24y, color = stateProvince))+
  geom_point(data = pop1, aes(x=LM25x, y=LM25y, color = stateProvince))+
  theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        panel.grid.major = element_blank(), #remove major gridlines
        panel.grid.minor = element_blank(), #remove minor gridlines
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent'))

ggsave('myplot.png', p, bg='transparent')

  

