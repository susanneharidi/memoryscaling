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

# Comparison RT Exp1 vs 3 -----------------------------------------------------
# load data experiment1 -------------------------------------------------------

current_directory <- dirname(rstudioapi::getSourceEditorContext()$path)
subdirectory = "/Experiment1"
setwd(paste(current_directory, subdirectory, sep = ""))

Data1 = read.csv("ExperimentDataExp1.csv")
DataSim1 = read.csv("ExperimentDataExp1W2VSim.csv")
dResponsesExp1 <- subset(Data1, rt < 20000 & validID == TRUE & Phase == "CuePresentation")
dResponsesSimExp1 <- subset(DataSim1, Phase == "CuePresentation" & validID == "True" & rt < 20000)
dResponsesExp1$rt = dResponsesExp1$rt/1000
dResponsesSimExp1$rt = dResponsesSimExp1$rt/1000

# use only the data from the valid positions
validPos = c(2, 4, 5, 7, 9)
dResponsesValidExp1 <- subset(dResponsesExp1, is.element(queriedPosition, validPos))
dResponsesValidSimExp1 <- subset(dResponsesSimExp1, is.element(queriedPosition, validPos))

# RT

# Define the color manual
color_manual <- c("1" = "#007BB8", "3" = "#F28E2B")

RTSummaryExp1 <- ddply(dResponsesValidExp1, 
                       ~ListLength + correct, 
                       summarise, meanRT = mean(rt), se = se(rt))

RTSummaryExp3 <- ddply(subset(dResponsesExp3, ContextSize != "None" & validity_rt == TRUE), 
                       ~ContextSize + correct + ListLength, 
                       summarise, meanRT = mean(rt), se = se(rt))
RTSummaryExp3$ContextSize = as.numeric(RTSummaryExp3$ContextSize)
RTSummaryExp1$ContextSize = RTSummaryExp1$ListLength
RTSummaryExp1$Experiment = "1"
RTSummaryExp3$Experiment = "3"
meanRTContextSizeExp3Control <- ddply(subset(dResponsesExp3, ContextSize == "None" & validity_rt == TRUE), 
                                      ~correct + NoneListType, 
                                      summarise, meanRT = mean(rt), se = se(rt))
meanRTContextSizeExp3Control$Experiment = "3"

RTSummary = rbind(RTSummaryExp1, RTSummaryExp3)

setSizePlotCombi <- ggplot(data = RTSummary, aes(x = ContextSize, y = meanRT, color = Experiment))+
  geom_line()+
  geom_ribbon(data = RTSummary,
              aes(ymin = meanRT - se, 
                  ymax = meanRT + se,
                  fill = Experiment), alpha = 0.2, color = NA)+
  labs(color = "Experiment",
       linetype = "Control Condotion", 
       y = 'mean RT in s +/- SE',
       x = "Set size/Context size") +
  theme(text = element_text(size = textsize, family = "sans"))+
  geom_hline(data = meanRTContextSizeExp3Control, 
             aes(yintercept = meanRT, 
                 linetype = as.factor(NoneListType)),
             color = "black")+
  coord_cartesian(ylim = c(RT_min_value, RT_max_value))+
  scale_color_manual(values = color_manual) +
  scale_fill_manual(values = color_manual) +
  facet_grid(cols = vars(correct))
setSizePlotCombi


# Accuracy

AccSummaryExp1 <- ddply(dResponsesValidExp1, 
                        ~ListLength, 
                        summarise, meanAcc = mean(correct_01), se = se(correct_01))

AccSummaryExp3 <- ddply(subset(dResponsesExp3, ContextSize != "None"), 
                        ~ContextSize + ListLength, 
                        summarise, meanAcc = mean(correct_01), se = se(correct_01))
AccSummaryExp3$ContextSize = as.numeric(AccSummaryExp3$ContextSize)
AccSummaryExp1$ContextSize = AccSummaryExp1$ListLength
AccSummaryExp1$Experiment = "1"
AccSummaryExp3$Experiment = "3"

meanAccContextSizeExp3Control <- ddply(subset(dResponsesExp3, ContextSize == "None"), 
                                       ~NoneListType, 
                                       summarise, meanAcc = mean(correct_01), se = se(correct_01))
meanAccContextSizeExp3Control$Experiment = "3"

AccSummary = rbind(AccSummaryExp1, AccSummaryExp3)

setSizePlotCombiAcc <- ggplot(data = AccSummary, aes(x = ContextSize, y = meanAcc, color = Experiment))+
  geom_line()+
  geom_ribbon(data = AccSummary,
              aes(ymin = meanAcc - se, 
                  ymax = meanAcc + se, 
                  fill = Experiment), alpha = RibbonAlpha, color = NA)+
  labs(color = "Experiment",
       linetype = "Control Condotion", 
       y = 'mean Acc in s +/- SE',
       x = "Set size/Context size") +
  geom_hline(data = meanAccContextSizeExp3Control, 
             aes(yintercept = meanAcc, 
                 linetype = as.factor(NoneListType)),
             color = "black")+
  coord_cartesian(ylim = c(0, 1))+
  scale_color_manual(values = color_manual) +
  scale_fill_manual(values = color_manual) +
  theme(text = element_text(size = textsize, family = "sans"))
setSizePlotCombiAcc


# Assembling the figure --------------------------------------------------------

# Remove axis labels and ticks from specific plots

setSizePlotCombiAcc <- setSizePlotCombiAcc + theme(legend.position = "None")



# Combine plots with shared legends and axes
FigureCombi <- (setSizePlotCombi | setSizePlotCombiAcc) + 
  plot_annotation(tag_levels = "A")+ 
  plot_layout(guides = "collect") +
  plot_layout(heights = c(1, 1))  # Ensure equal heights for both rows
FigureCombi

savingPlace = "/Figures/ComparisonExperiment1And3.pdf"
ggsave(paste(current_directory, savingPlace, sep = ""), FigureCombi,  
       width = 20, 
       height = 5)


# Statistical analysis ---------------------------------------------------------

meanRTsExp1 = ddply(subset(dResponsesExp1, correct == "TRUE"), 
                          ~subject_id + ListLength, 
                          summarise, meanRT = mean(rt))
meanAccsExp1 = ddply(dResponsesExp1, 
                          ~subject_id + ListLength, 
                          summarise, meanAcc = mean(correct_01))

# Test for differences between the two experiments

meanRTsExp1_20 = subset(meanRTsExp1, ListLength == 20)
meanRTsExp3Context = ddply(subset(dResponsesExp3, correct == "TRUE" & ContextSize != "None"),
                                      ~subject_id, 
                                      summarise, meanRT = mean(rt))
BFRTExp1_vs_3_Context = ttestBF(meanRTsExp1_20$meanRT, 
                                meanRTsExp3Context$meanRT,
                     paired = FALSE)
BFRTExp1_vs_3_Context
# Bayes factor analysis
# --------------
#   [1] Alt., r=0.707 : 0.6560139 ±0.01%
# 
# Against denominator:
#   Null, mu1-mu2 = 0 
# ---
#   Bayes factor type: BFindepSample, JZS

# calculate mean difference
mean(meanRTsExp1_20$meanRT) - mean(meanRTsExp3Context$meanRT)
# -0.7627214


meanRTsExp3HighSimPairs = ddply(subset(dResponsesExp3, correct == "TRUE" & ContextSize == "None" & NoneListType == "highSimPairs"),
                                ~subject_id, 
                                summarise, meanRT = mean(rt))

BFRTExp1_vs_3_meanRTsExp3HighSimPairs = ttestBF(meanRTsExp1_20$meanRT, 
                                                meanRTsExp3HighSimPairs$meanRT,
                                paired = FALSE)
BFRTExp1_vs_3_meanRTsExp3HighSimPairs
# Bayes factor analysis
# --------------
#   [1] Alt., r=0.707 : 0.920687 ±0.01%
# 
# Against denominator:
#   Null, mu1-mu2 = 0 
# ---
#   Bayes factor type: BFindepSample, JZS

# calculate mean difference
mean(meanRTsExp1_20$meanRT) - mean(meanRTsExp3HighSimPairs$meanRT)
# 0.7949374



meanRTsExp3LowSimList = ddply(subset(dResponsesExp3, correct == "TRUE" & ContextSize == "None" & NoneListType == "lowSimList"),
                                ~subject_id, 
                                summarise, meanRT = mean(rt))

BFRTExp1_vs_3_meanRTsExp3LowSimList = ttestBF(meanRTsExp1_20$meanRT, 
                                                meanRTsExp3LowSimList$meanRT,
                                                paired = FALSE)
BFRTExp1_vs_3_meanRTsExp3LowSimList
# Bayes factor analysis
# --------------
#   [1] Alt., r=0.707 : 0.2482479 ±0.02%
# 
# Against denominator:
#   Null, mu1-mu2 = 0 
# ---
#   Bayes factor type: BFindepSample, JZS

# calculate mean difference
mean(meanRTsExp1_20$meanRT) - mean(meanRTsExp3LowSimList$meanRT)
# -0.1579138


# Do the same for accuracy

meanAccsExp1_20 = subset(meanAccsExp1, ListLength == 20)
meanAccsExp3Context = ddply(subset(dResponsesExp3, ContextSize != "None"),
                            ~subject_id, 
                            summarise, meanAcc = mean(correct_01))

BFAccExp1_vs_3_Context = ttestBF(meanAccsExp1_20$meanAcc, 
                                meanAccsExp3Context$meanAcc,
                                paired = FALSE)
BFAccExp1_vs_3_Context
# Bayes factor analysis
# --------------
#   [1] Alt., r=0.707 : 35578.06 ±0%
# 
# Against denominator:
#   Null, mu1-mu2 = 0 
# ---
#   Bayes factor type: BFindepSample, JZS

# calculate mean difference
mean(meanAccsExp1_20$meanAcc) - mean(meanAccsExp3Context$meanAcc)
#0.2974473

meanAccsExp3HighSimPairs = ddply(subset(dResponsesExp3, ContextSize == "None" & NoneListType == "highSimPairs"),
                                ~subject_id, 
                                summarise, meanAcc = mean(correct_01))

BFAccExp1_vs_3_meanAccsExp3HighSimPairs = ttestBF(meanAccsExp1_20$meanAcc, 
                                                meanAccsExp3HighSimPairs$meanAcc,
                                                paired = FALSE)
BFAccExp1_vs_3_meanAccsExp3HighSimPairs
# Bayes factor analysis
# --------------
#   [1] Alt., r=0.707 : 0.2590727 ±0.02%
# 
# Against denominator:
#   Null, mu1-mu2 = 0 
# ---
#   Bayes factor type: BFindepSample, JZS

# calculate mean difference
mean(meanAccsExp1_20$meanAcc) - mean(meanAccsExp3HighSimPairs$meanAcc)
# -0.01570524

meanAccsExp3LowSimList = ddply(subset(dResponsesExp3, ContextSize == "None" & NoneListType == "lowSimList"),
                                ~subject_id, 
                                summarise, meanAcc = mean(correct_01))

BFAccExp1_vs_3_meanAccsExp3LowSimList = ttestBF(meanAccsExp1_20$meanAcc, 
                                                meanAccsExp3LowSimList$meanAcc,
                                                paired = FALSE)
BFAccExp1_vs_3_meanAccsExp3LowSimList
# Bayes factor analysis
# --------------
#   [1] Alt., r=0.707 : 46.94 ±0%
# 
# Against denominator:
#   Null, mu1-mu2 = 0 
# ---
#   Bayes factor type: BFindepSample, JZS

# calculate mean difference
mean(meanAccsExp1_20$meanAcc) - mean(meanAccsExp3LowSimList$meanAcc)
# 0.213021

# Comparing only to small contexts

# RT

meanRTsExp3ContextSmall = subset(dResponsesExp3, correct == "TRUE" & ContextSize != "None")
meanRTsExp3ContextSmall$ContextSize = as.numeric(meanRTsExp3ContextSmall$ContextSize)
meanRTsExp3ContextSmall = ddply(subset(meanRTsExp3ContextSmall, ContextSize < 5),
                                 ~subject_id, 
                                 summarise, meanRT = mean(rt))

BFAccExp1_vs_3_meanRTsExp3ContextSmall = ttestBF(meanRTsExp1_20$meanRT, 
                                                  meanRTsExp3ContextSmall$meanRT,
                                                  paired = FALSE)
BFAccExp1_vs_3_meanRTsExp3ContextSmall
# Bayes factor analysis
# --------------
#   [1] Alt., r=0.707 : 0.2500235 ±0.02%
# 
# Against denominator:
#   Null, mu1-mu2 = 0 
# ---
#   Bayes factor type: BFindepSample, JZS

# calculate mean difference
mean(meanRTsExp1_20$meanRT) - mean(meanRTsExp3ContextSmall$meanRT)
# -0.251627



# Acc
meanAccsExp3ContextSmall = subset(dResponsesExp3, ContextSize != "None")
meanAccsExp3ContextSmall$ContextSize = as.numeric(meanAccsExp3ContextSmall$ContextSize)
meanAccsExp3ContextSmall = ddply(subset(meanAccsExp3ContextSmall, ContextSize < 5),
                            ~subject_id, 
                            summarise, meanAcc = mean(correct_01))

BFAccExp1_vs_3_meanAccsExp3ContextSmall = ttestBF(meanAccsExp1_20$meanAcc, 
                                                meanAccsExp3ContextSmall$meanAcc,
                                                paired = FALSE)
BFAccExp1_vs_3_meanAccsExp3ContextSmall
# Bayes factor analysis
# --------------
#   [1] Alt., r=0.707 : 10.0529 ±0%
# 
# Against denominator:
#   Null, mu1-mu2 = 0 
# ---
#   Bayes factor type: BFindepSample, JZS

# calculate mean difference
mean(meanAccsExp1_20$meanAcc) - mean(meanAccsExp3ContextSmall$meanAcc)
# 0.161234
