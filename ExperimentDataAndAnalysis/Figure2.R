# Memory scaling Figure 1


# load packages ----------------------------------------------------------------

library(tidyverse)
library(ggplot2)
theme_set(
  theme_minimal() +
    theme(legend.position = "bottom")
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


# Getting Data from Experiment 3 ----------------------------------------------

current_directory <- dirname(rstudioapi::getSourceEditorContext()$path)
subdirectory = "/Experiment3"
setwd(paste(current_directory, subdirectory, sep = ""))

Data = read.csv("ExperimentDataExp3.csv")
dResponsesExp3 <- subset(Data, Phase == "CuePresentation" & rt < 20 & validity == TRUE)

# Getting Data from Experiment 2 ----------------------------------------------

current_directory <- dirname(rstudioapi::getSourceEditorContext()$path)
subdirectory = "/Experiment2"
setwd(paste(current_directory, subdirectory, sep = ""))

Data2 = read.csv("ExperimentDataExp2.csv")
dResponsesExp2 <- subset(Data2, Phase == "CuePresentation" & rt < 20 & validity == TRUE)


# Plotting variables -----------------------------------------------------------

textsize = 22
font = "sans"
default_color <- "#007BB8"
AccColor <- "#74c476"
RTColor <- "#56B4E9"
RegColor = "black"
RibbonAlpha = 0.2
RT_min_value = 0
RT_max_value = 7

# Context size plot Exp3 ------------------------------------------------------

# RT
exp3NotControl = subset(dResponsesExp3, ContextSize != "None" & validity_rt == TRUE)
exp3NotControl$ContextSize = as.numeric(exp3NotControl$ContextSize)
meanRTContextSizeExp3 <- ddply(exp3NotControl, 
                               ~correct + ContextSize, 
                               summarize, meanRT = mean(rt), se = se(rt))


meanRTContextSizeExp3Control <- ddply(subset(dResponsesExp3, ContextSize == "None" & validity_rt == TRUE), 
                                      ~correct + NoneListType, 
                                      summarise, meanRT = mean(rt), se = se(rt))

contextSizeExp3 <- ggplot(subset(meanRTContextSizeExp3, correct == TRUE ), aes(x = ContextSize, y = meanRT))+
  geom_line(color = RTColor)+
  geom_ribbon(aes(ymin = meanRT - se, 
                  ymax = meanRT + se), alpha = RibbonAlpha, fill = RTColor)+
  geom_hline(data = subset(meanRTContextSizeExp3Control, correct == TRUE), 
             aes(yintercept = meanRT, 
                 linetype = as.factor(NoneListType)),
             color = "black")+
  geom_smooth(data = subset(exp3NotControl, correct == TRUE), 
              aes(x = ContextSize, y = rt),
              color = RTColor, 
              method = lm,
              linetype = "dashed",
              se = FALSE)+
  labs(y = 'Mean RT in s +/- SE',
       x = "Context size",
       linetype = "Control condition") +
  scale_linetype_manual(values = c("solid", "dotted"),  # Adjust linetypes as needed
                        labels = c("High similarity pairs", "Low similarity list")) +
  theme(legend.position = "bottom")+
  coord_cartesian(ylim = c(RT_min_value, RT_max_value))+
  theme(text = element_text(size = textsize, family = "sans"))
contextSizeExp3

# Accuracy
exp3NotControl = subset(dResponsesExp3, ContextSize != "None" )
exp3NotControl$ContextSize = as.numeric(exp3NotControl$ContextSize)
meanRTContextSizeExp3 <- ddply(exp3NotControl, 
                               ~correct + ContextSize, 
                               summarize, meanRT = mean(rt), se = se(rt))


meanAccContextSizeExp3 <- ddply(exp3NotControl, 
                                ~ContextSize, 
                                summarize, meanAcc = mean(correct_01), se = se(correct_01))


meanAccContextSizeExp3Control <- ddply(subset(dResponsesExp3, ContextSize == "None"), 
                                       ~NoneListType, 
                                       summarise, meanAcc = mean(correct_01), se = se(correct_01))

contextSizeAccExp3 <- ggplot(meanAccContextSizeExp3, aes(x = ContextSize, y = meanAcc))+
  geom_line(color = AccColor)+
  geom_ribbon(aes(ymin = meanAcc - se, 
                  ymax = meanAcc + se), alpha = RibbonAlpha, fill = AccColor)+
  geom_hline(data = meanAccContextSizeExp3Control, 
             aes(yintercept = meanAcc, 
                 linetype = as.factor(NoneListType)),color = "black")+
  geom_smooth(data = exp3NotControl, 
              aes(x = ContextSize, y = correct_01),
              color = AccColor, 
              method = lm,
              linetype = "dashed",
              se = FALSE)+
  labs(y = 'Mean accuracy +/- SE',
       x = "Context size",
       linetype = "Control Condition") +
  scale_linetype_manual(values = c("solid", "dotted"),  # Adjust linetypes as needed
                        labels = c("High similarity pairs", "Low similarity list")) +
  coord_cartesian(ylim = c(0, 1))+
  theme(legend.position = "bottom")+
  theme(text = element_text(size = textsize, family = "sans"))
contextSizeAccExp3

# Context size plot Exp2 ------------------------------------------------------

# RT
exp2NotControl = subset(dResponsesExp2, ContextSize != "None")
exp2NotControl$ContextSize = as.numeric(exp2NotControl$ContextSize)
meanRTContextSizeExp2 <- ddply(exp2NotControl, 
                               ~correct + ContextSize + ListLength, 
                               summarize, meanRT = mean(rt), se = se(rt))


meanRTContextSizeExp2Control <- ddply(subset(dResponsesExp2, ContextSize == "None"), 
                                      ~correct + ListLength, 
                                      summarise, meanRT = mean(rt), se = se(rt))

contextSizeExp2 <- ggplot(subset(meanRTContextSizeExp2, correct == TRUE), 
                          aes(x = ContextSize, y = meanRT, color = as.factor(ListLength))) +
  geom_line() +
  geom_ribbon(aes(ymin = meanRT - se, 
                  ymax = meanRT + se, 
                  fill = as.factor(ListLength)), alpha = RibbonAlpha, color = NA) +
  geom_hline(data = subset(meanRTContextSizeExp2Control, correct == TRUE), 
             aes(yintercept = meanRT, color = as.factor(ListLength))) +
  geom_smooth(data = subset(exp2NotControl, correct == TRUE), 
              aes(x = ContextSize, y = rt, color = as.factor(ListLength)), 
              method = lm, 
              linetype = "dashed", 
              se = FALSE) +
  labs(y = 'Mean RT in s +/- SE',
       x = "Context size",
       linetype = "Control condition",
       color = "Set size",   # Add legend label for line color
       fill = "Set size") +  # Add legend label for fill color
  scale_color_manual(values = c("#E69F00", "#D55E00", "#56B4E9", "#009E73")) +  # Okabe-Ito palette
  scale_fill_manual(values = c("#E69F00",  "#D55E00", "#56B4E9", "#009E73")) +   # Okabe-Ito palette
  theme(legend.position = "bottom") +
  coord_cartesian(ylim = c(RT_min_value, RT_max_value)) +
  theme(text = element_text(size = textsize, family = "sans"))

contextSizeExp2



# Accuracy
exp2NotControl = subset(dResponsesExp2, ContextSize != "None" )
exp2NotControl$ContextSize = as.numeric(exp2NotControl$ContextSize)
meanRTContextSizeExp2 <- ddply(exp2NotControl, 
                               ~correct + ContextSize + ListLength, 
                               summarize, meanRT = mean(rt), se = se(rt))


meanAccContextSizeExp2 <- ddply(exp2NotControl, 
                                ~ContextSize + ListLength, 
                                summarize, meanAcc = mean(correct_01), se = se(correct_01))


meanAccContextSizeExp2Control <- ddply(subset(dResponsesExp2, ContextSize == "None"), 
                                       ~ContextSize + ListLength, 
                                       summarise, meanAcc = mean(correct_01), se = se(correct_01))

contextSizeAccExp2 <- ggplot(meanAccContextSizeExp2, 
                             aes(x = ContextSize, y = meanAcc, color = as.factor(ListLength))) +
  geom_line() +
  geom_ribbon(aes(ymin = meanAcc - se, 
                  ymax = meanAcc + se, 
                  fill = as.factor(ListLength)), alpha = RibbonAlpha, color = NA) +
  geom_hline(data = meanAccContextSizeExp2Control, 
             aes(yintercept = meanAcc, color = as.factor(ListLength))) +
  geom_smooth(data = exp2NotControl, 
              aes(x = ContextSize, y = correct_01, color = as.factor(ListLength)), 
              method = lm, 
              linetype = "dashed", 
              se = FALSE) +
  labs(y = 'Mean accuracy +/- SE',
       x = "Context size",
       color = "Set size",   # Add legend label for line color
       fill = "Set size") +  # Add legend label for fill color
  scale_color_manual(values = c("#E69F00", "#D55E00", "#56B4E9", "#009E73")) +  # Okabe-Ito palette
  scale_fill_manual(values = c("#E69F00",  "#D55E00", "#56B4E9", "#009E73")) +   # Okabe-Ito palette
  coord_cartesian(ylim = c(0, 1)) +
  theme(legend.position = "bottom") +
  theme(text = element_text(size = textsize, family = "sans"))

contextSizeAccExp2




# Assembling the figure --------------------------------------------------------

# Remove axis labels and legends from specific plots# Remove axis labels and legends from specific plots
contextSizeAccExp3 <- contextSizeAccExp3 + theme(legend.position = "none")
contextSizeAccExp2 <- contextSizeAccExp2 + theme(legend.position = "none")


# Group plots for Experiment 2 and Experiment 3 with titles
exp2_plots <- (contextSizeExp2 + contextSizeAccExp2) +
  plot_layout(guides = "collect") +
  plot_annotation(title = "Experiment 2") &
  theme(plot.title = element_text(hjust = 0.5, size = 14))

exp3_plots <- (contextSizeExp3 + contextSizeAccExp3) +
  plot_layout(guides = "collect") +
  plot_annotation(title = "Experiment 3") &
  theme(plot.title = element_text(hjust = 0.5, size = 14))

# Combine the two groups horizontally
Figure2 <- (exp2_plots | exp3_plots) +
  plot_annotation(tag_levels = "A")  # Add subplot labels (A, B, C, D)
Figure2


savingPlace = "/Figures/Figure2.pdf"
ggsave(paste(current_directory, savingPlace, sep = ""), Figure2,  
       width = 20, 
       height = 5)


