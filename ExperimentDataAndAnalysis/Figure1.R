# Memory scaling Figure 1


# load packages ----------------------------------------------------------------

library(tidyverse)
library(ggplot2)
theme_set(
  theme_minimal() +
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
library(brms)
library(patchwork)

# SET WD:
# If you are not using the here package and R projects, this is the place to 
# set the working directory to the inside of the github project


# Functions --------------------------------------------------------------------

#standard error function
se <- function(x){sd(x)/sqrt(length(x))}


# Turning my txt and json data files into df -----------------------------------


current_directory <- dirname(rstudioapi::getSourceEditorContext()$path)
subdirectory = "/Experiment1"
setwd(paste(current_directory, subdirectory, sep = ""))

Data = read.csv("ExperimentDataExp1.csv")
DataSim = read.csv("ExperimentDataExp1W2VSim.csv")


# making some relevant sub-dataframes ------------------------------------------


# only get the response trials
dResponses <- subset(Data, rt < 20000 & validID == TRUE & Phase == "CuePresentation")

# same for sim dataframe
dResponsesSim <- subset(DataSim, Phase == "CuePresentation" & validID == "True" & rt < 20000)

dResponses$rt = dResponses$rt/1000
dResponsesSim$rt = dResponsesSim$rt/1000

# use only the data from the valid positions
validPos = c(2, 4, 5, 7, 9)
dResponsesValid <- subset(dResponses, is.element(queriedPosition, validPos))
dResponsesValidSim <- subset(dResponsesSim, is.element(queriedPosition, validPos))

# Plotting variables
textsize = 22
font = "sans"
default_color <- "#007BB8"
AccColor <- "#74c476"
RTColor <- "#56B4E9"
RegColor = "black"
RibbonAlpha = 0.2
RT_min_value = 0
RT_max_value = 4


# set size plots ---------------------------------------------------------------

# RTs 
meanRTSetSizeExp1 = ddply(dResponsesValid, ~ ListLength + correct,
                          summarise,
                          meanRT = mean(rt),
                          se = se(rt))

setSizeExp1_2 <- ggplot(subset(meanRTSetSizeExp1), aes(x = ListLength, y = meanRT, color = correct, fill = correct))+
  geom_line()+
  geom_ribbon(aes(ymin = meanRT - se, 
                  ymax = meanRT + se), alpha = RibbonAlpha, color = "white")+
  geom_smooth(method = lm, se = FALSE, linetype = "dashed")+
  labs(y = 'Mean RT in s +/- SE',
       x = "Set size",
       color = "Correct",
       fill = "Correct") +
  coord_cartesian(ylim = c(RT_min_value, 7))+
  theme(legend.position = "bottom")+
  scale_fill_manual(values = c("TRUE" = "forestgreen", "FALSE" = "red"))+
  scale_color_manual(values = c("TRUE" = "forestgreen", "FALSE" = "red"))+
  theme(text = element_text(size = textsize, family = "sans"))
setSizeExp1_2


# Accuracy

meanAccSetSizeExp1 = ddply(dResponsesValid, ~ ListLength,
                           summarise,
                           meanAcc = mean(correct_01),
                           se = se(correct_01))

setSizeExp1Acc <- ggplot(meanAccSetSizeExp1, aes(x = ListLength, y = meanAcc))+
  geom_line(color = AccColor)+
  geom_ribbon(aes(ymin = meanAcc - se, 
                  ymax = meanAcc + se), alpha = RibbonAlpha, fill = AccColor)+
  geom_smooth(data = dResponsesValid, 
              aes(x = ListLength, y = correct_01),
              color = AccColor, method = lm, se = FALSE, linetype = "dashed")+
  labs(y = 'Mean accuracy +/- SE',
       x = "Set size") +
  theme(legend.position = "bottom")+
  ylim(0,1)+
  theme(text = element_text(size = textsize, family = "sans"))
setSizeExp1Acc


# Similarity plot Exp1 --------------------------------------------------------

# RT

breaks = seq(0, 1, by = 0.1)
dResponsesValidSim$simBin <- cut(dResponsesValidSim$W2VWordPairsSim, 
                                 breaks = breaks, 
                                 labels = breaks[-1]-0.05,
                                 include.lowest = TRUE, 
                                 right = FALSE)

# aggregate data across sim bins and set size
AggragatedTimes <- ddply(dResponsesValidSim, 
                         ~simBin + correct, 
                         summarize, MeanRT = mean(rt), se = se(rt))
AggragatedTimes$simBin <- as.numeric(as.character(AggragatedTimes$simBin))


SimilarityRTExp1 <- ggplot(AggragatedTimes, aes(x = simBin, y = MeanRT, group = correct, color = correct)) +
  #geom_point() +
  geom_line() +
  geom_ribbon(aes(ymin = MeanRT - se, 
                  ymax = MeanRT + se, fill = correct), alpha = RibbonAlpha, color = "white")+
  # add x limit values
  scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
  labs(y = 'Mean RT in s +/- SE',
       x = "Cue target similarity",
       color = "Correct",
       fill = "Correct") +
  geom_smooth(method = lm, se = FALSE, linetype = "dashed")+
  coord_cartesian(ylim = c(RT_min_value, 7))+
  theme(legend.position = "bottom")+
  scale_fill_manual(values = c("TRUE" = "forestgreen", "FALSE" = "red"))+
  scale_color_manual(values = c("TRUE" = "forestgreen", "FALSE" = "red"))+
  theme(text = element_text(size = textsize, family = "sans"))
SimilarityRTExp1

# Accuracy

# aggregate data across sim bins and set size
AggragatedAcc <- ddply(dResponsesValidSim, 
                         ~simBin, 
                         summarize, MeanAcc = mean(correct_01), se = se(correct_01))
AggragatedAcc$simBin <- as.numeric(as.character(AggragatedAcc$simBin))

SimilarityAccExp1 <- ggplot(AggragatedAcc, aes(x = simBin, y = MeanAcc))+
  geom_line(color = AccColor) +
  geom_ribbon(aes(ymin = MeanAcc - se, 
                  ymax = MeanAcc + se), alpha = RibbonAlpha, color = "white", fill = AccColor)+
  # add x limit values
  scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
  labs(y = 'Mean accuracy +/- SE',
       x = "Cue target similarity",
       color = "Correct",
       fill = "Correct") +
  geom_smooth(method = lm, se = FALSE, linetype = "dashed", color = AccColor)+
  coord_cartesian(ylim = c(0, 1))+
  theme(legend.position = "bottom")+
  theme(text = element_text(size = textsize, family = "sans"))
SimilarityAccExp1

# Context Plot Exp1 -----------------------------------------------------------

# Define the color manual
color_manual <- c("1" = "#4E79A7", "2" = "#F28E2B", "3" = "#76B7B2", "4" = "#E15759")

# RT
meanRTSetSizeExp1Block = ddply(dResponsesValid, ~ ListLength + correct + Block,
                               summarise,
                               meanRT = mean(rt),
                               se = se(rt))

setSizeExp1Block <- ggplot(subset(meanRTSetSizeExp1Block, correct == TRUE), 
                           aes(x = ListLength, y = meanRT))+
  geom_line(aes(color = as.factor(Block)))+
  geom_ribbon(aes(ymin = meanRT - se, 
                  ymax = meanRT + se,
                  fill = as.factor(Block)), alpha = RibbonAlpha)+
  geom_smooth(data = subset(dResponsesValid, correct == TRUE), 
              aes(x = ListLength, y = rt, color = as.factor(Block)),
              method = lm,
              se = FALSE,
              linetype = "dashed")+
  labs(y = 'Mean RT in s +/- SE',
       x = "Set size",
       color = "Block",
       fill = "Block") +
  scale_color_manual(values = color_manual) +
  scale_fill_manual(values = color_manual) +
  coord_cartesian(ylim = c(RT_min_value, 6))+
  theme(text = element_text(size = textsize, family = "sans"))
setSizeExp1Block

# Acc 
meanAccSetSizeExp1Block = ddply(dResponsesValid, ~ ListLength + Block,
                                summarise,
                                meanAcc = mean(correct_01),
                                se = se(correct_01))

setSizeExp1AccBlock <- ggplot(meanAccSetSizeExp1Block, 
                              aes(x = ListLength, y = meanAcc))+
  geom_line(aes(color = as.factor(Block)))+
  geom_ribbon(aes(ymin = meanAcc - se, 
                  ymax = meanAcc + se,
                  fill = as.factor(Block)), alpha = RibbonAlpha)+
  geom_smooth(data = dResponsesValid, 
              aes(x = ListLength, y = correct_01, color = as.factor(Block)),
              method = lm,
              se = FALSE,
              linetype = "dashed")+
  labs(y = 'Mean accuracy in s +/- SE',
       x = "Set size",
       color = "Block",
       fill = "Block") +
  ylim(0,1)+
  scale_color_manual(values = color_manual) +
  scale_fill_manual(values = color_manual) +
  theme(text = element_text(size = textsize, family = "sans"))
setSizeExp1AccBlock


# Assembling the figure --------------------------------------------------------

# Remove axis labels and ticks from specific plots
setSizeExp1Block <- setSizeExp1Block + theme(axis.title.y = element_blank())
SimilarityRTExp1 <- SimilarityRTExp1 + theme(axis.title.y = element_blank())
SimilarityAccExp1 <- SimilarityAccExp1 + theme(axis.title.y = element_blank())
setSizeExp1AccBlock <- setSizeExp1AccBlock + theme(axis.title.y = element_blank(), legend.position = "None")
SimilarityRTExp1 <- SimilarityRTExp1 + theme(axis.title.y = element_blank(), legend.position = "None")

# Combine plots with shared legends and axes
Figure1 <- ((setSizeExp1_2 + SimilarityRTExp1 + setSizeExp1Block) / (setSizeExp1Acc + SimilarityAccExp1 + setSizeExp1AccBlock)) + 
  plot_annotation(tag_levels = "A")+ 
  plot_layout(guides = "collect")+
  plot_layout(heights = c(1, 1))  # Ensure equal heights for both rows
Figure1

savingPlace = "/Figures/Figure1.pdf"
ggsave(paste(current_directory, savingPlace, sep = ""), Figure1,  
       width = 30, 
       height = 10)
