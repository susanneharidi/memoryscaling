

# load packages --------------------------------------------------------------



library(tidyverse)
library(ggplot2)
theme_set(
  theme_minimal() +
    theme(legend.position = "top")
)

library(dplyr)
library(plyr)
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
library(sys)
library(brms)
library(BayesFactor)

# SET WD:
# If you are not using the here package and R projects, this is the place to 
# set the working directory to the inside of the github project


# Fuctions --------------------------------------------------------------------

# standard error function
se <- function(x){sd(x)/sqrt(length(x))}


# load all data from experiment2 ----------------------------------------------

current_directory <- dirname(rstudioapi::getSourceEditorContext()$path)
subdirectory = "/Experiment2"
setwd(paste(current_directory, subdirectory, sep = ""))

Data = read.csv("ExperimentDataExp2.csv")


# Get Response Data and Block exclusion ---------------------------------------

# only get the response trials
dResponses <- subset(Data, Phase == "CuePresentation" & rt < 20 & validity == TRUE)
dResponsesLetters <- subset(Data, Phase == "LetterRecallResponse" & validity == TRUE)

# Exclude Blocks with an accuracy below 0.3

# calculate average accuracy for each block
BlockAccuracy <- ddply(dResponses, ~subject_id+Block, 
                       summarize, block_accuracy = mean(correct_01))
BlockAccuracy <- BlockAccuracy %>%
  mutate(valid_block = block_accuracy >= 0.3)

# Merge the validity back into the original dataframe
dResponses <- merge(dResponses, BlockAccuracy, by = c("subject_id", "Block"))

dResponsesValidBlocks <- subset(dResponses, valid_block == TRUE)
# remove "None" Context size
dResponsesValidBlocksContext <- subset(dResponsesValidBlocks, ContextSize != "None")
# make context size numeric
dResponsesValidBlocksContext$ContextSize <- as.numeric(dResponsesValidBlocksContext$ContextSize)
# MC and sclaed variables
dResponsesValidBlocksContext$ContextSizeMCScaled <- scale(dResponsesValidBlocksContext$ContextSize, scale = TRUE)

# Calculate the brms models ---------------------------------------------------
current_directory <- dirname(rstudioapi::getSourceEditorContext()$path)
subdirectory = "/brmModels"
setwd(paste(current_directory, subdirectory, sep = ""))
seed = 2024

# Hypothesis 1 RT ∼ context-size + set-size + (context-size + set-size | subject) ---

RTModelContextSizeSetSize <- brm(rt ~ ContextSizeMCScaled + as.factor(ListLength) + (ContextSizeMCScaled|subject_id), 
                                 data = subset(dResponsesValidBlocksContext, correct == "TRUE"),
                                 iter = 10000,
                                 cores = 4,
                                 save_pars = save_pars(all = TRUE),
                                 seed = seed,
                                 control = list(max_treedepth = 15, adapt_delta = 0.99),
                                 family = lognormal(),
                                 file = "mRTExp2ContextSizeSetSizeScaled_2")
summary(RTModelContextSizeSetSize)
# plot model effects
# Redone Results
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept                 0.86      0.06     0.74     0.97 1.00     1921     4082
# ContextSizeMCScaled      -0.03      0.01    -0.05    -0.01 1.00    19454    15533
# as.factorListLength30     0.22      0.08     0.06     0.38 1.00     1645     2863

RTModelContextSizeNull <- brm(rt ~  as.factor(ListLength) + (ContextSizeMCScaled|subject_id), 
                                 data = subset(dResponsesValidBlocksContext, correct == "TRUE"),
                                 iter = 10000,
                                 cores = 4,
                                 save_pars = save_pars(all = TRUE),
                                 seed = seed,
                                 control = list(max_treedepth = 15, adapt_delta = 0.99),
                                 family = lognormal(),
                                 file = "mRTExp2ContextSizeNull")
summary(RTModelContextSizeNull)

#calculate bayes factor
BF <- bayes_factor(RTModelContextSizeSetSize, RTModelContextSizeNull)
BF
# results
# 0.68970

#
# Accuracy analysis -----------------------------------------------------------
#

AccModelContextSize <- brm(correct_01 ~ ContextSizeMCScaled + as.factor(ListLength) + (ContextSizeMCScaled|subject_id), 
                           data = dResponsesValidBlocksContext,
                           iter = 10000,
                           cores = 4,
                           save_pars = save_pars(all = TRUE),
                           seed = seed,
                           control = list(max_treedepth = 15, adapt_delta = 0.99),
                           file = "exp2AccContextSizeFactor",
                           family = "bernoulli")

summary(AccModelContextSize)
# # revised results
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS
# Intercept                 1.11      0.16     0.80     1.42 1.00     2977
# ContextSizeMCScaled      -0.08      0.05    -0.17     0.01 1.00    17018
# as.factorListLength30     0.08      0.23    -0.36     0.53 1.00     2534


AccModelContextSizeNull <- brm(correct_01 ~ as.factor(ListLength) + (ContextSizeMCScaled|subject_id), 
                               data = dResponsesValidBlocksContext,
                               iter = 10000,
                               cores = 4,
                               save_pars = save_pars(all = TRUE),
                               seed = seed,
                               control = list(max_treedepth = 15, adapt_delta = 0.99),
                               file = "exp2AccContextSizeFactorNull",
                               family = "bernoulli")

summary(AccModelContextSizeNull)

# calculate BF
BF <- bayes_factor(AccModelContextSize, AccModelContextSizeNull)
BF
# 0.46969


# Replication Exp 1 RT ∼ set-size + (set-size | subject) ----------------------


SummarisedRT = ddply(subset(dResponsesValidBlocks, correct == "TRUE"), ~ListLength + subject_id, summarize, meanRT = mean(rt))
# do pairwise comparisons of set size via ttestBF
ttestBF(subset(SummarisedRT, ListLength == 15)$meanRT, subset(SummarisedRT, ListLength == 30)$meanRT, paired = FALSE)
# Bayes factor analysis
# --------------
#   [1] Alt., r=0.707 : 4.567106 ±0.01%
# 
# Against denominator:
#   Null, mu1-mu2 = 0 
# ---
#   Bayes factor type: BFindepSample, JZS

# mean RT difference between set size 15 and 30
mean(subset(dResponsesValidBlocks, correct == "TRUE"& ListLength == 15)$rt) - mean(subset(dResponsesValidBlocks, correct == "TRUE"& ListLength == 30)$rt)
# -0.893854

# Replication Exp 1 RT ∼ accuracy + (accuracy | subject) ----------------------


SummarisedRTAcc = ddply(dResponsesValidBlocks, ~correct_01 + subject_id, summarize, meanRT = mean(rt))
# do pairwise comparisons of set size via ttestBF
ttestBF(subset(SummarisedRTAcc, correct_01 == 1)$meanRT, subset(SummarisedRTAcc, correct_01 == 0)$meanRT, paired = FALSE)
# Bayes factor analysis
# --------------
#   [1] Alt., r=0.707 : 2.370763e+13 ±0%
# 
# Against denominator:
#   Null, mu1-mu2 = 0 
# ---
#   Bayes factor type: BFindepSample, JZS

# mean RT difference between correct and incorrect responses
mean(subset(dResponsesValidBlocks, correct == "TRUE")$rt) - mean(subset(dResponsesValidBlocks, correct == "FALSE")$rt)
# -2.5160042



# Replicate Accuracy Liste Length Experiment 2 --------------------------------


SummarisedAcc = ddply(dResponsesValidBlocks, ~ListLength + subject_id, summarize, meanAcc = mean(correct_01))
# do pairwise comparisons of set size via ttestBF
ttestBF(subset(SummarisedAcc, ListLength == 15)$meanAcc, subset(SummarisedAcc, ListLength == 30)$meanAcc, paired = FALSE)
# Bayes factor analysis
# --------------
#   [1] Alt., r=0.707 : 0.2218371 ±0.03%
# 
# Against denominator:
#   Null, mu1-mu2 = 0 
# ---
#   Bayes factor type: BFindepSample, JZS

# mean RT difference between set size 15 and 30
mean(subset(dResponsesValidBlocks, ListLength == 15)$correct_01) - mean(subset(dResponsesValidBlocks, ListLength == 30)$correct_01)
# 
