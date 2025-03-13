
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

# SET WD:
# If you are not using the here package and R projects, this is the place to 
# set the working directory to the inside of the github project

########################
# Functions
#########################
#standard error function
se <- function(x){sd(x)/sqrt(length(x))}

seed = 2023

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

dResponsesValid <- subset(dResponses, validPos == TRUE)
dResponsesValidSim <- subset(dResponsesSim, validPos == TRUE)

# meancentered variables
dResponsesValid$ListLengthMC <- scale(dResponsesValid$ListLength, scale = FALSE)
dResponsesValidSim$W2VWordPairsSimMC <- scale(dResponsesValidSim$W2VWordPairsSim, scale = FALSE)
dResponsesValidSim$ListLengthMC <- scale(dResponsesValidSim$ListLength, scale = FALSE)

# MC and sclaed variables
dResponsesValidSim$ListLengthMCScaled <- scale(dResponsesValidSim$ListLength, scale = TRUE)
dResponsesValidSim$W2VWordPairsSimMCScaled <- scale(dResponsesValidSim$W2VWordPairsSim, scale = TRUE)

##############################################
# calculate a brms model
##############################################

current_directory <- dirname(rstudioapi::getSourceEditorContext()$path)
subdirectory = "/brmModels"
setwd(paste(current_directory, subdirectory, sep = ""))

#####################################################
# Hypothesis 1 RT ∼ set-size + (set-size | subject)
#####################################################

#Model with only lag controlled RTs
RTModelSetSizeLagControlled <- brm(rt ~ ListLengthMCScaled + (ListLengthMCScaled|subject_id), 
                                          data = subset(dResponsesValidSim, correct == "TRUE" & ListLength > 9),
                                          iter = 10000,
                                          cores = 4,
                                          save_pars = save_pars(all = TRUE),
                                          seed = seed,
                                          control = list(max_treedepth = 15, adapt_delta = 0.99),
                                          family = lognormal(),
                                          file = "mRTExp1SetSizeLagControlled")
summary(RTModelSetSizeLagControlled)
# 26.09.24
# Regression Coefficients:
#                     Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept              0.66      0.05     0.57     0.76 1.00     2403     4719
# ListLengthMCScaled     0.09      0.02     0.05     0.13 1.00    17471    15032

# Check for convergence
mcmc_plot(RTModelSetSizeLagControlled, type = "trace")
# Plot effects 
mcmc_plot(RTModelSetSizeLagControlled, type = "intervals", prob = 0.95)
# only plot main populationlevel effects
mcmc_plot(RTModelSetSizeLagControlled, type = "areas", prob = 0.95, variable = "^b_", regex = TRUE)+
  geom_vline(xintercept = 0, linetype = "longdash", color = "gray")+
  xlab("RT in s")+
  theme_minimal()+
  theme(text = element_text(size = 20, family = "sans"))+
  scale_y_discrete(labels = c("Intercept", "Set size"))+# 
  ggtitle("Experiment 1 set size effect")
pp_check(RTModelSetSizeLagControlled, type = "ecdf_overlay")
#visualize all effects
plot(conditional_effects(RTModelSetSizeLagControlled), points = FALSE)


#Null
RTModelSetSizeLagControlledNull <- brm(rt ~ 1 + (ListLengthMCScaled|subject_id), 
                                   data = subset(dResponsesValidSim, correct == "TRUE" & ListLength > 9),
                                   iter = 10000,
                                   cores = 4,
                                   save_pars = save_pars(all = TRUE),
                                   seed = seed,
                                   control = list(max_treedepth = 15, adapt_delta = 0.99),
                                   family = lognormal(),
                                   file = "mRTExp1SetSizeLagControlledNull")
summary(RTModelSetSizeLagControlledNull)

bayes_factor(RTModelSetSizeLagControlled, RTModelSetSizeLagControlledNull)
# 26.09.24
# BF = 555.61425

#Model including RTs from smaller setsizes
RTModelSetSize <- brm(rt ~ ListLengthMCScaled + (ListLengthMCScaled|subject_id), 
                                          data = subset(dResponsesValidSim, correct == "TRUE"),
                                          iter = 10000,
                                          cores = 4,
                                          save_pars = save_pars(all = TRUE),
                                          seed = seed,
                                          control = list(max_treedepth = 15, adapt_delta = 0.99),
                                          family = lognormal(),
                                          file = "mRTExp1SetSize")
summary(RTModelSetSize)
# 26.09.24
# Regression Coefficients:
#                      Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept              0.66      0.05     0.57     0.76 1.00     2403     4719
# ListLengthMCScaled     0.09      0.02     0.05     0.13 1.00    17471    15032

# Check for convergence
mcmc_plot(RTModelSetSize, type = "trace")


# Null
RTModelSetSizeNull <- brm(rt ~  1 + (ListLengthMCScaled|subject_id), 
                      data = subset(dResponsesValidSim, correct == "TRUE"),
                      iter = 10000,
                      cores = 4,
                      save_pars = save_pars(all = TRUE),
                      seed = seed,
                      control = list(max_treedepth = 15, adapt_delta = 0.99),
                      family = lognormal(),
                      file = "mRTExp1SetSizeNull")
summary(RTModelSetSizeNull)

bayes_factor(RTModelSetSize, RTModelSetSizeNull)
# 26.09.24
# BF = 553.42379

# calculate mean increase in RT per unscaled set size increase
# get all unique list lengths sorted by size
setSizes = sort(unique(dResponsesValidSim$ListLength))
# Now calculate the mean RT increase for each set size in a loop
meanRTIncrease = c()
for (i in 1:(length(setSizes)-1)){
  meanRTIncrease = c(meanRTIncrease, mean(dResponsesValidSim$rt[dResponsesValidSim$ListLength == setSizes[i+1]]) - mean(dResponsesValidSim$rt[dResponsesValidSim$ListLength == setSizes[i]]))
}
# get the mean
mean(meanRTIncrease)/2

##############################################################################
# Hypothesis 2 RT ∼ set-size * similarity + (set-size + similarity | subject)
##############################################################################

RTModelSimilarityInteraction <- brm(rt ~ ListLengthMCScaled * W2VWordPairsSimMCScaled + (ListLengthMCScaled + W2VWordPairsSimMCScaled|subject_id), 
                                          data = subset(dResponsesValidSim, correct == "TRUE"),
                                          iter = 10000,
                                          cores = 4,
                                          save_pars = save_pars(all = TRUE),
                                          seed = seed,
                                          control = list(max_treedepth = 15, adapt_delta = 0.99),
                                          family = lognormal(),
                                          file = "mRTExp1SimilarityCorrection")
summary(RTModelSimilarityInteraction)
# 26.09.24
#                                               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept                                      0.65      0.05     0.56     0.74 1.00     1984     3897
# ListLengthMCScaled                             0.12      0.02     0.09     0.15 1.00     8790    14363
# W2VWordPairsSimMCScaled                       -0.06      0.01    -0.08    -0.04 1.00    32824    16048
# ListLengthMCScaled:W2VWordPairsSimMCScaled    -0.01      0.01    -0.03     0.01 1.00    28374    15280

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
  scale_y_discrete(labels = c("Intercept", "Set size", "Similarity", "Interaction"))+# 
  ggtitle("Experiment 1 similarity effect")
pp_check(RTModelSimilarityInteraction, type = "ecdf_overlay")
#visualize all effects
plot(conditional_effects(RTModelSimilarityInteraction), points = FALSE)


RTModelSimilarity <- brm(rt ~ ListLengthMCScaled + W2VWordPairsSimMCScaled + (ListLengthMCScaled + W2VWordPairsSimMCScaled|subject_id), 
                                          data = subset(dResponsesValidSim, correct == "TRUE"),
                                          iter = 10000,
                                          cores = 4,
                                          save_pars = save_pars(all = TRUE),
                                          seed = seed,
                                          control = list(max_treedepth = 15, adapt_delta = 0.99),
                                          family = lognormal(),
                                          file = "mRTExp1SimilarityNoInteractionCorrection")
summary(RTModelSimilarity)

bayes_factor(RTModelSimilarityInteraction, RTModelSimilarity)
# 26.09.24
# BF = 0.00058

# Reverse: Model without interaction is clearly superior
bayes_factor(RTModelSimilarity, RTModelSimilarityInteraction)
# 26.09.24
# BF = 1798.25769

RTModelSimilarityNull <- brm(rt ~ ListLengthMCScaled + (ListLengthMCScaled + W2VWordPairsSimMCScaled|subject_id), 
                         data = subset(dResponsesValidSim, correct == "TRUE"),
                         iter = 10000,
                         cores = 4,
                         save_pars = save_pars(all = TRUE),
                         seed = seed,
                         control = list(max_treedepth = 15, adapt_delta = 0.99),
                         family = lognormal(),
                         file = "mRTExp1SimilarityNull")
summary(RTModelSimilarityNull)

# Reverse: Model without interaction is clearly superior
bayes_factor(RTModelSimilarity, RTModelSimilarityNull)
# 26.09.24
# BF = 1798.25769

# caculate the difference between mean RTs for sims below 0.2 and sims above 0.8
mean(dResponsesValidSim$rt[dResponsesValidSim$W2VWordPairsSim < 0.2]) - mean(dResponsesValidSim$rt[dResponsesValidSim$W2VWordPairsSim > 0.8])


##############################################################################
# Hypothesis 3 RT ∼ accuracy + (accuracy | subject)
##############################################################################

RTModelAccuracy <- brm(rt ~ correct_01 + (correct_01|subject_id), 
                             data = subset(dResponses),
                             iter = 10000,
                             cores = 4,
                             save_pars = save_pars(all = TRUE),
                             seed = seed,
                             control = list(max_treedepth = 15, adapt_delta = 0.99),
                             family = lognormal(),
                             file = "mRTExp1Accuracy")
summary(RTModelAccuracy)
# 26.09.24
#             Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept      1.46      0.06     1.34     1.57 1.00     9287    13056
# correct_01    -0.76      0.04    -0.85    -0.68 1.00    15681    15257

# Check for convergence
mcmc_plot(RTModelAccuracy, type = "trace")
# Plot effects 
mcmc_plot(RTModelAccuracy, type = "intervals", prob = 0.95)
# only plot main populationlevel effects
mcmc_plot(RTModelAccuracy, type = "areas", prob = 0.95, variable = "^b_", regex = TRUE)+
  geom_vline(xintercept = 0, linetype = "longdash", color = "gray")+
  xlab("RT in s")+
  theme_minimal()+
  theme(text = element_text(size = 20, family = "sans"))+
  scale_y_discrete(labels = c("Accuracy", "Set size"))+# 
  ggtitle("Experiment 1 set size effect")
pp_check(RTModelAccuracy, type = "ecdf_overlay")
#visualize all effects
plot(conditional_effects(RTModelAccuracy), points = FALSE)


RTModelAccuracyNull <- brm(rt ~ 1 + (correct_01|subject_id), 
                       data = subset(dResponses),
                       iter = 10000,
                       cores = 4,
                       save_pars = save_pars(all = TRUE),
                       seed = seed,
                       control = list(max_treedepth = 15, adapt_delta = 0.99),
                       family = lognormal(),
                       file = "mRTExp1AccuracyNull")
summary(RTModelAccuracyNull)

bayes_factor(RTModelAccuracy, RTModelAccuracyNull)
# 26.09.24
# BF = 27423112360680199372501155840.00000

# calculate the mean RT difference between correct and incorrect trials
mean(dResponses$rt[dResponses$correct_01 == 1]) - mean(dResponses$rt[dResponses$correct_01 == 0])


################################################################################
# Hypothesis 1 and 2 Accuracy ∼ set-size + similarity + (set-size + similarity| subject)
###############################################################################

#Model with only lag controlled RTs
AccuracyModelSetSizeSimLagControlled <- brm(correct_01 ~ ListLengthMCScaled + W2VWordPairsSimMCScaled + (ListLengthMCScaled + W2VWordPairsSimMCScaled|subject_id), 
                                   data = subset(dResponsesValidSim, ListLength > 9),
                                   iter = 10000,
                                   cores = 4,
                                   save_pars = save_pars(all = TRUE),
                                   seed = seed,
                                   control = list(max_treedepth = 15, adapt_delta = 0.99),
                                   family = bernoulli(),
                                   file = "mAccuracyExp1SetSizeSimLagControlled")

summary(AccuracyModelSetSizeSimLagControlled)
# results 16.12.2024

# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept                   1.84      0.14     1.58     2.13 1.00     8946    11981
# ListLengthMCScaled         -0.23      0.10    -0.42    -0.04 1.00    29566    15779
# W2VWordPairsSimMCScaled     0.36      0.09     0.18     0.55 1.00    17561    15909



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
AccuracyModelSetSizeLagControlledNullSim <- brm(correct_01 ~ ListLengthMCScaled + (ListLengthMCScaled + W2VWordPairsSimMCScaled|subject_id), 
                                        data = subset(dResponsesValidSim, ListLength > 9),
                                        iter = 10000,
                                        cores = 4,
                                        save_pars = save_pars(all = TRUE),
                                        seed = seed,
                                        control = list(max_treedepth = 15, adapt_delta = 0.99),
                                        family = bernoulli(),
                                        file = "mAccuracyExp1SetSizeLagControlledNullSim")

# calculate BF
bayes_factor(AccuracyModelSetSizeSimLagControlled, AccuracyModelSetSizeLagControlledNullSim)

# calulate Null
AccuracyModelSetSizeLagControlledNullSetSize <- brm(correct_01 ~ W2VWordPairsSimMCScaled + (ListLengthMCScaled + W2VWordPairsSimMCScaled|subject_id), 
                                                data = subset(dResponsesValidSim, ListLength > 9),
                                                iter = 10000,
                                                cores = 4,
                                                save_pars = save_pars(all = TRUE),
                                                seed = seed,
                                                control = list(max_treedepth = 15, adapt_delta = 0.99),
                                                family = bernoulli(),
                                                file = "mAccuracyExp1SetSizeLagControlledNullSetSize")

# calculate BF
bayes_factor(AccuracyModelSetSizeSimLagControlled, AccuracyModelSetSizeLagControlledNullSetSize)

##############################################################################
# Exploration: is there a Block effect
##############################################################################

# Interaction effect
RTModelSetSizeLagControlledBlockInteraction <- brm(rt ~ ListLengthMCScaled * Block + (ListLengthMCScaled + Block|subject_id), 
                                                   data = subset(dResponsesValidSim, correct == "TRUE" & ListLength > 9),
                                                   iter = 10000,
                                                   cores = 4,
                                                   save_pars = save_pars(all = TRUE),
                                                   seed = seed,
                                                   control = list(max_treedepth = 15, adapt_delta = 0.99),
                                                   family = lognormal(),
                                                   file = "mRTExp1SetSizeLagControlledBlockInteraction")
summary(RTModelSetSizeLagControlledBlockInteraction)
# results 11.12.24
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept                    0.68      0.07     0.54     0.83 1.00     5847     9852
# ListLengthMCScaled           0.14      0.06     0.02     0.27 1.00     9523    13572
# Block                       -0.01      0.02    -0.05     0.03 1.00    11859    14398
# ListLengthMCScaled:Block    -0.02      0.02    -0.07     0.03 1.00     9447    12884

# calculate bayes factor for unteraction
bayes_factor(RTModelSetSizeLagControlledBlockInteraction, RTModelSetSizeLagControlledBlock)
#0.09318

#Model with only lag controlled RTs
RTModelSetSizeLagControlledBlock <- brm(rt ~ ListLengthMCScaled + Block + (ListLengthMCScaled + Block|subject_id), 
                                        data = subset(dResponsesValidSim, correct == "TRUE" & ListLength > 9),
                                        iter = 10000,
                                        cores = 4,
                                        save_pars = save_pars(all = TRUE),
                                        seed = seed,
                                        control = list(max_treedepth = 15, adapt_delta = 0.99),
                                        family = lognormal(),
                                        file = "mRTExp1SetSizeLagControlledBlock")
summary(RTModelSetSizeLagControlledBlock)
# 7.10.24

#Null
RTModelSetSizeLagControlledBlockNull <- brm(rt ~ ListLengthMCScaled + (ListLengthMCScaled + Block|subject_id), 
                                            data = subset(dResponsesValidSim, correct == "TRUE" & ListLength > 9),
                                            iter = 10000,
                                            cores = 4,
                                            save_pars = save_pars(all = TRUE),
                                            seed = seed,
                                            control = list(max_treedepth = 15, adapt_delta = 0.99),
                                            family = lognormal(),
                                            file = "mRTExp1SetSizeLagControlledBlockInteractionBull")
summary(RTModelSetSizeLagControlledBlockNull)

# calculate bayes factor for main effect of block
bayes_factor(RTModelSetSizeLagControlledBlock, RTModelSetSizeLagControlledBlockNull)
# 0.04646

#Test Block 1 directly against block 4
Block1 = subset(dResponsesValidSim, correct == "TRUE" & Block == 1)
Block4 = subset(dResponsesValidSim, correct == "TRUE" & Block == 4)

Block1MeanRTs = ddply(Block1, ~ subject_id, summarize, meanRT = mean(rt))
Block4MeanRTs = ddply(Block4, ~ subject_id, summarize, meanRT = mean(rt))

# combine the above dataframes, so that only subjects that have data in both blocks are compared
BlockMeanRTs = merge(Block1MeanRTs, Block4MeanRTs, by = "subject_id")


# Do baysian pairwise ttest
ttestBF = ttestBF(x = BlockMeanRTs$meanRT.x, y = BlockMeanRTs$meanRT.y, paired = TRUE)
# Bayes factor analysis
# --------------
#   [1] Alt., r=0.707 : 0.3132393 ±0.06%
# 
# Against denominator:
#   Null, mu = 0 
# ---
#   Bayes factor type: BFoneSample, JZS


# same for accy=uracy
Block1 = subset(dResponsesValidSim, Block == 1)
Block4 = subset(dResponsesValidSim, Block == 4)

Block1Mean = ddply(Block1, ~ subject_id, summarize, meanAcc = mean(correct_01))
Block4Mean = ddply(Block4, ~ subject_id, summarize, meanAcc = mean(correct_01))

# combine the above dataframes, so that only subjects that have data in both blocks are compared
BlockMeanAcc = merge(Block1Mean, Block4Mean, by = "subject_id")


# Do baysian pairwise ttest
ttestBF = ttestBF(x = BlockMeanAcc$meanAcc.x, y = BlockMeanAcc$meanAcc.y, paired = TRUE)
# Bayes factor analysis
# --------------
#   [1] Alt., r=0.707 : 0.2073064 ±0.09%
# 
# Against denominator:
#   Null, mu = 0 
# ---
#   Bayes factor type: BFoneSample, JZS



