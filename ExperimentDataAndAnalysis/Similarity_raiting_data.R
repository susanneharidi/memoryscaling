
#########################################################################
# load packages
#########################################################################


library(tidyverse)
library(ggplot2)
theme_set(
  theme_bw() +
    theme(legend.position = "top")
)
library(plyr)
library(dplyr)
library(data.table)
library(gridExtra)
library(ggsignif)
library(lmerTest)
library(here)
library(jsonlite)
library(lme4)
library(readr)
library(ggthemes)
library(RWiener)

# SET WD:
# If you are not using the here package and R projects, this is the place to 
# set the working directory to the inside of the github project

############################################################
# Fuctions
############################################################
#standard error function
se <- function(x){sd(x)/sqrt(length(x))}

##############################################################################
# Load the data with the W2V Sim rating 
# created in the python script called AddW2VSimToSimilarityRatings.py
###############################################################################

current_directory <- dirname(rstudioapi::getSourceEditorContext()$path)
subdirectory = "/SimilarityRaitings"
setwd(paste(current_directory, subdirectory, sep = ""))
DataSim = read.csv("HumanAndW2VSims.csv")

#####################################################################
# calculate the correlation and plot the relationship
####################################################################

correlation <- cor(DataSim$meanSim, DataSim$W2VWordPairsSim, method = 'pearson')
correlation
#r = 0.8069761

ggplot(DataSim, aes(x = meanSim, y = W2VWordPairsSim, size = sd)) +
  geom_point(color = "blue", alpha = 0.5)+
  geom_smooth(method = lm)+
  theme_minimal()+
  labs(title = "Human vs W2V Similarity",
       x = "Human Rating",
       y = "W2V")

