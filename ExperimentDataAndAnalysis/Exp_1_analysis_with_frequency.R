
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
library(brms)
library(BayesFactor)
library(readxl)

# SET WD:
# If you are not using the here package and R projects, this is the place to 
# set the working directory to the inside of the github project

########################
# Functions
#########################
#standard error function
se <- function(x){sd(x)/sqrt(length(x))}

seed = 2023

###############################################
# Read in the subtlex word freuqncy file
###############################################
current_directory <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(current_directory)
subtlex <- read_excel("SUBTLEX.xls")


####################################################################
# Turning my txt and json data files into df
###################################################################
current_directory <- dirname(rstudioapi::getSourceEditorContext()$path)
subdirectory = "/Experiment1"
setwd(paste(current_directory, subdirectory, sep = ""))

Data = read.csv("ExperimentDataExp1.csv")
DataSim = read.csv("ExperimentDataExp1W2VSim.csv")
Data$rt = Data$rt/1000
DataSim$rt = DataSim$rt/1000

###################################################
# prepare the dataframes
######################################################

# only get the response trials
dResponses <- subset(Data, rt < 20000 & validID == TRUE & Phase == "CuePresentation")
dResponsesSim <- subset(DataSim, Phase == "CuePresentation" & validID == "True" & rt < 20000)
validPos = c(2, 4, 5, 7, 9)
dResponses$validPos = is.element(dResponses$queriedPosition, validPos) | dResponses$ListLength < 10
dResponsesSim$validPos = is.element(dResponsesSim$queriedPosition, validPos) | dResponsesSim$ListLength < 10

# check if exclusion worked
ggplot(dResponses, aes(x = ListLength, y = queriedPosition, color = validPos))+
  geom_jitter(width = 0.4, height = 0)

dResponsesValid <- subset(dResponses, validPos == TRUE)
dResponsesValidSim <- subset(dResponsesSim, validPos == TRUE)

# meancentered variables
dResponsesValid$ListLengthMC <- scale(dResponsesValid$ListLength, scale = FALSE)
dResponsesValidSim$W2VWordPairsSimMC <- scale(dResponsesValidSim$W2VWordPairsSim, scale = FALSE)
dResponsesValidSim$ListLengthMC <- scale(dResponsesValidSim$ListLength, scale = FALSE)

# MC and sclaed variables
dResponsesValidSim$ListLengthMCScaled <- scale(dResponsesValidSim$ListLength, scale = TRUE)
dResponsesValidSim$W2VWordPairsSimMCScaled <- scale(dResponsesValidSim$W2VWordPairsSim, scale = TRUE)


###############################################
# add the word frequency to the data
###############################################
# Ensure both columns are lowercase for comparison
subtlex$Word <- tolower(subtlex$Word)
dResponsesValidSim$target <- tolower(dResponsesValidSim$target)

# Use match() to find indices
match_indices <- match(dResponsesValidSim$target, subtlex$Word)

# Assign the corresponding frequency values
dResponsesValidSim$WordFrequencyTarget <- subtlex$FREQcount[match_indices]

# Identify words that were not found
not_found <- is.na(dResponsesValidSim$WordFrequencyTarget)
if (any(not_found)) {
  print(paste("Words not found in SUBTLEX:", paste(dResponsesValidSim$target[not_found], collapse = ", ")))
}

# only show the target and thw word frequency column
dResponsesValidSimtemp <- dResponsesValidSim[, c("target", "WordFrequencyTarget")]
# MC and sclaed variables
dResponsesValidSim$WordFrequencyTargetScaled <- scale(dResponsesValidSim$WordFrequencyTarget, scale = TRUE)


##############################################
# calculate a brms model
##############################################

current_directory <- dirname(rstudioapi::getSourceEditorContext()$path)
subdirectory = "/brmModels"
setwd(paste(current_directory, subdirectory, sep = ""))

ggplot(dResponsesValidSim, aes(x = WordFrequencyTarget, y = rt))+
  geom_jitter(width = 0.4, height = 0)+
  coord_cartesian(ylim = c(0, 20))+
  xlim(0, 2000)+
  geom_smooth(method = "lm", se = FALSE)

ggplot(dResponsesValidSim, aes(x = WordFrequencyTarget, y = correct_01))+
  geom_jitter(width = 0.4, height = 0)+
  xlim(0, 2000)+
  geom_smooth(method = "lm", se = FALSE)


##############################################################################
# Hypothesis 2 RT ∼ set-size * similarity + (set-size + similarity | subject)
##############################################################################

RTModelSimilarityInteraction <- brm(rt ~ ListLengthMCScaled * W2VWordPairsSimMCScaled + WordFrequencyTargetScaled + 
                                      (ListLengthMCScaled + W2VWordPairsSimMCScaled + WordFrequencyTargetScaled|subject_id), 
                                    data = subset(dResponsesValidSim, correct == "TRUE"),
                                    iter = 10000,
                                    cores = 4,
                                    save_pars = save_pars(all = TRUE),
                                    seed = seed,
                                    control = list(max_treedepth = 15, adapt_delta = 0.99),
                                    family = lognormal(),
                                    file = "mRTExp1SimilarityWordFrequency")
summary(RTModelSimilarityInteraction)
#18.03.25
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept                                      0.65      0.05     0.56     0.74 1.00     1914     3750
# ListLengthMCScaled                             0.12      0.02     0.09     0.15 1.00     8325    13718
# W2VWordPairsSimMCScaled                       -0.06      0.01    -0.08    -0.04 1.00    32018    14237
# WordFrequencyTargetScaled                     -0.02      0.01    -0.04     0.01 1.00    22855    13342
# ListLengthMCScaled:W2VWordPairsSimMCScaled    -0.01      0.01    -0.03     0.01 1.00    32130    14692

# Check for convergence
mcmc_plot(RTModelSimilarityInteraction, type = "trace")
# Plot effects 
mcmc_plot(RTModelSimilarityInteraction, type = "intervals", prob = 0.95)
# only plot main populationlevel effects
mcmc_plot(RTModelSimilarityInteraction, type = "areas", prob = 0.95, variable = "^b_", regex = TRUE)+
  geom_vline(xintercept = 0, linetype = "longdash", color = "gray")+
  xlab("RT in s")+
  theme_minimal()+
  theme(text = element_text(size = 20, family = "sans"))+
  ggtitle("Experiment 1 similarity effect")
pp_check(RTModelSimilarityInteraction, type = "ecdf_overlay")
#visualize all effects
plot(conditional_effects(RTModelSimilarityInteraction), points = FALSE)


RTModelSimilarity <- brm(rt ~ ListLengthMCScaled + W2VWordPairsSimMCScaled + WordFrequencyTargetScaled + 
                           (ListLengthMCScaled + W2VWordPairsSimMCScaled + WordFrequencyTargetScaled|subject_id), 
                         data = subset(dResponsesValidSim, correct == "TRUE"),
                         iter = 10000,
                         cores = 4,
                         save_pars = save_pars(all = TRUE),
                         seed = seed,
                         control = list(max_treedepth = 15, adapt_delta = 0.99),
                         family = lognormal(),
                         file = "mRTExp1SimilarityNoInteractionFrequency")
summary(RTModelSimilarity)
# 18.03.25
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept                     0.65      0.05     0.56     0.74 1.00     2144     4335
# ListLengthMCScaled            0.12      0.02     0.09     0.15 1.00     9815    14449
# W2VWordPairsSimMCScaled      -0.06      0.01    -0.08    -0.04 1.00    31160    15873
# WordFrequencyTargetScaled    -0.02      0.01    -0.04     0.01 1.00    24350    13446

bayes_factor(RTModelSimilarityInteraction, RTModelSimilarity)
# 18.03.25
# BF = 0.05931

# Reverse: Model without interaction is clearly superior
bayes_factor(RTModelSimilarity, RTModelSimilarityInteraction)
# 18.03.25
# BF = 26.45887


RTModelSimilarityNull <- brm(rt ~ ListLengthMCScaled + WordFrequencyTargetScaled + (ListLengthMCScaled + W2VWordPairsSimMCScaled + WordFrequencyTargetScaled|subject_id), 
                             data = subset(dResponsesValidSim, correct == "TRUE"),
                             iter = 10000,
                             cores = 4,
                             save_pars = save_pars(all = TRUE),
                             seed = seed,
                             control = list(max_treedepth = 15, adapt_delta = 0.99),
                             family = lognormal(),
                             file = "mRTExp1SimilarityNullFrequency")
summary(RTModelSimilarityNull)
# 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept                     0.64      0.05     0.54     0.73 1.00     1695     3119
# ListLengthMCScaled            0.12      0.02     0.09     0.15 1.00     8514    12724
# WordFrequencyTargetScaled    -0.01      0.01    -0.04     0.02 1.00    26507    14614

# Reverse: Model without interaction is clearly superior
bayes_factor(RTModelSimilarity, RTModelSimilarityNull)
# 18.03.25
# BF = 2865.80542


mRTExp1FrequencyNullFrequency <- brm(rt ~ ListLengthMCScaled + W2VWordPairsSimMCScaled + 
                              (ListLengthMCScaled + W2VWordPairsSimMCScaled + WordFrequencyTargetScaled|subject_id), 
                             data = subset(dResponsesValidSim, correct == "TRUE"),
                             iter = 10000,
                             cores = 4,
                             save_pars = save_pars(all = TRUE),
                             seed = seed,
                             control = list(max_treedepth = 15, adapt_delta = 0.99),
                             family = lognormal(),
                             file = "mRTExp1FrequencyNullFrequency")
summary(mRTExp1FrequencyNullFrequency)

bayes_factor(RTModelSimilarity, mRTExp1FrequencyNullFrequency)
# BF = 0.07877



################################################################################
# Hypothesis 1 and 2 Accuracy ∼ set-size + similarity + (set-size + similarity| subject)
###############################################################################

#Model with only lag controlled RTs
AccuracyModelSetSizeSimLagControlled <- brm(correct_01 ~ ListLengthMCScaled + W2VWordPairsSimMCScaled + WordFrequencyTargetScaled +
                                              (ListLengthMCScaled + W2VWordPairsSimMCScaled + WordFrequencyTargetScaled|subject_id), 
                                            data = subset(dResponsesValidSim, ListLength > 9),
                                            iter = 10000,
                                            cores = 4,
                                            save_pars = save_pars(all = TRUE),
                                            seed = seed,
                                            control = list(max_treedepth = 15, adapt_delta = 0.99),
                                            family = bernoulli(),
                                            file = "mAccuracyExp1SetSizeSimLagControlledFrequency")

summary(AccuracyModelSetSizeSimLagControlled)
# results 18.03.25

# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept                     1.87      0.14     1.61     2.17 1.00    11145    13959
# ListLengthMCScaled           -0.23      0.09    -0.41    -0.04 1.00    32755    15418
# W2VWordPairsSimMCScaled       0.38      0.10     0.19     0.56 1.00    19943    15940
# WordFrequencyTargetScaled     0.24      0.19    -0.01     0.72 1.00     8271    10040


# Check for convergence
mcmc_plot(AccuracyModelSetSizeLagControlled, type = "trace")
# Plot effects 
mcmc_plot(AccuracyModelSetSizeLagControlled, type = "intervals", prob = 0.95)
# only plot main populationlevel effects
mcmc_plot(AccuracyModelSetSizeLagControlled, type = "areas", prob = 0.95, variable = "^b_", regex = TRUE)+
  geom_vline(xintercept = 0, linetype = "longdash", color = "gray")+
  xlab("RT in s")+
  theme_minimal()+
  theme(text = element_text(size = 20, family = "sans"))+
  scale_y_discrete(labels = c("Accuracy", "Set size"))+# 
  ggtitle("Experiment 1 set size effect")
pp_check(AccuracyModelSetSizeLagControlled, type = "ecdf_overlay")
#visualize all effects
plot(conditional_effects(AccuracyModelSetSizeLagControlled), points = FALSE)

# calulate Null
AccuracyModelSetSizeLagControlledNullSim <- brm(correct_01 ~ ListLengthMCScaled + WordFrequencyTargetScaled +
                                                  (ListLengthMCScaled + W2VWordPairsSimMCScaled + WordFrequencyTargetScaled|subject_id), 
                                                data = subset(dResponsesValidSim, ListLength > 9),
                                                iter = 10000,
                                                cores = 4,
                                                save_pars = save_pars(all = TRUE),
                                                seed = seed,
                                                control = list(max_treedepth = 15, adapt_delta = 0.99),
                                                family = bernoulli(),
                                                file = "mAccuracyExp1SetSizeLagControlledNullSimFrequency")
summary(AccuracyModelSetSizeLagControlledNullSim)

# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept                     2.01      0.14     1.75     2.31 1.00     8677    12410
# ListLengthMCScaled           -0.23      0.10    -0.42    -0.04 1.00    21521    14342
# WordFrequencyTargetScaled     0.23      0.20    -0.03     0.75 1.00     5965     7092

# calculate BF
bayes_factor(AccuracyModelSetSizeSimLagControlled, AccuracyModelSetSizeLagControlledNullSim)
# BF = 887.06073

# calulate Null
AccuracyModelSetSizeLagControlledNullSetSize <- brm(correct_01 ~ W2VWordPairsSimMCScaled + WordFrequencyTargetScaled + 
                                                      (ListLengthMCScaled + W2VWordPairsSimMCScaled + WordFrequencyTargetScaled|subject_id), 
                                                    data = subset(dResponsesValidSim, ListLength > 9),
                                                    iter = 10000,
                                                    cores = 4,
                                                    save_pars = save_pars(all = TRUE),
                                                    seed = seed,
                                                    control = list(max_treedepth = 15, adapt_delta = 0.99),
                                                    family = bernoulli(),
                                                    file = "mAccuracyExp1SetSizeLagControlledNullSetSizeFrequency")

# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept                     1.78      0.14     1.52     2.06 1.00     6950    10754
# W2VWordPairsSimMCScaled       0.38      0.10     0.19     0.57 1.00    12733    14149
# WordFrequencyTargetScaled     0.24      0.18    -0.01     0.71 1.00     6039     7048

# calculate BF
bayes_factor(AccuracyModelSetSizeSimLagControlled, AccuracyModelSetSizeLagControlledNullSetSize)



AccuracyModelFrequencyLagControlledNull <- brm(correct_01 ~ ListLengthMCScaled + W2VWordPairsSimMCScaled +
                                              (ListLengthMCScaled + W2VWordPairsSimMCScaled + WordFrequencyTargetScaled|subject_id), 
                                            data = subset(dResponsesValidSim, ListLength > 9),
                                            iter = 10000,
                                            cores = 4,
                                            save_pars = save_pars(all = TRUE),
                                            seed = seed,
                                            control = list(max_treedepth = 15, adapt_delta = 0.99),
                                            family = bernoulli(),
                                            file = "AccuracyModelFrequencyLagControlledNull")

summary(AccuracyModelFrequencyLagControlledNull)
# calculate BF
bayes_factor(AccuracyModelSetSizeSimLagControlled, AccuracyModelFrequencyLagControlledNull)
# BF = 1.29609