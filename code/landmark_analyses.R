# Script: run analyses on landmark data
# Code: Wing Venation and Asymmetry in Halictus ligatus Bees
# Charles Thrift
# 18 February 2025

## load libraries

## read data




# MAIN ANALYSES BELOW

# First, population level analyses between the two populations. 
## Are these wings different?

## Multivariate Analysis of Variance (MANOVA) for 2 populations
pop <- lm

manova_pop <- pop 
manova_pop <- manova_pop[,c(2:19,34)] #UPDATE FOR RIGHT COLUMNS -- keep only numeric
manova_pop$location <- as.factor(manova_pop$location)
manova_pop_data <- manova_pop[,c(1:18)] #removing the population column
manova_pop_data <- as.matrix(manova_pop_data) 

fit <- lm.rrpp(manova_pop_data ~ location, SS.type = "I", 
               data = manova_pop, print.progress = FALSE) #run the linear model fit
fitm <- manova.update(fit, print.progress = FALSE, tol = 0) #run the manova update
summary(fitm, test = "Pillai") #summarize with a manova table



