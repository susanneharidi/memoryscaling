

# load packages --------------------------------------------------------------

library(tidyverse)
library(ggplot2)
theme_set(
  theme_minimal() +
    theme(legend.position = "top")
)

library(dplyr)
library(plyr)
library(lme4)
library(brms)

# SET WD:
# If you are not using the here package and R projects, this is the place to 
# set the working directory to the inside of the github project


# Fuctions --------------------------------------------------------------------

# standard error function
se <- function(x){sd(x)/sqrt(length(x))}

binary_se <- function(binary_array) {
  if (!is.vector(binary_array)) {
    stop("Input must be a vector")
  }
  
  if (!all(binary_array %in% c(0, 1))) {
    stop("Input vector must contain only 0s and 1s")
  }
  
  n <- length(binary_array)
  p <- mean(binary_array)
  
  se <- sqrt(p * (1 - p) / n)
  
  return(se)
}

# Variables ------------------------------------------------------------------

seed = 2024

# load all data from experiment2 ----------------------------------------------

current_directory <- dirname(rstudioapi::getSourceEditorContext()$path)
subdirectory = "/Experiment3"
setwd(paste(current_directory, subdirectory, sep = ""))

Data = read.csv("ExperimentDataExp3.csv")
                
# get some demographics -------------------------------------------------------

N = length(unique(Data$subject_id))
# should be 116
N_valid = length(unique(subset(Data, validity == TRUE)$subject_id))
# should be 112

# Get Response Data and Block exclusion ---------------------------------------

# only get the response trials that are valid
dResponses <- subset(Data, Phase == "CuePresentation" & rt < 20 & validity == TRUE)
dResponsesLetters <- subset(Data, Phase == "LetterRecallResponse" & validity == TRUE)

# Classify things into small and large context vs no context ------------------

dResponses <- dResponses %>%
  mutate(RoughContextSize = case_when(
    ContextSize == "None" & NoneListType == "highSimPairs" ~ "0_NoneHighSim",
    ContextSize == "None" & NoneListType == "lowSimList" ~ "0_NoneLowSim",
    as.numeric(ContextSize) < 5 ~ "1_small",
    as.numeric(ContextSize) > 15 ~ "3_large",
    TRUE ~ "2_medium"
  ))

# plot this as a boxplot
ggplot(subset(dResponses, correct == TRUE), 
       aes(x = RoughContextSize, y = rt, color = RoughContextSize))+
  geom_boxplot()+
  facet_grid(cols = vars(as.factor(ListLength)))


###############################################################################
# Main analyses -------
###############################################################################

###############################################################################
# RT -------
###############################################################################

# BRMS Models: RT Context Size ------------------------------------------------
current_directory <- dirname(rstudioapi::getSourceEditorContext()$path)
subdirectory = "/brmModels"
setwd(paste(current_directory, subdirectory, sep = ""))

correctResponses = subset(dResponses, correct == TRUE & 
                            ContextSize != "None" &
                            validity_rt == TRUE)
correctResponses$ContextSize = as.numeric(correctResponses$ContextSize)
correctResponses$ContextSizeMCScaled <- scale(correctResponses$ContextSize, scale = TRUE)

RTModelContextSize <- brm(rt ~ ContextSizeMCScaled + (ContextSizeMCScaled|subject_id), 
               data = correctResponses,
               iter = 10000,
               cores = 4,
               save_pars = save_pars(all = TRUE),
               seed = seed,
               control = list(max_treedepth = 15, adapt_delta = 0.99),
               file = "exp3RTContextSizelogScaled",
               family = lognormal()
               )

# Model summary
summary(RTModelContextSize)

# results log 27.06.2024
#                        Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept               1.15      0.05     1.04     1.26 1.00     1623     3427
# ContextSizeMCScaled     0.03      0.01     0.01     0.05 1.00    23261    15330

# Check for convergence
mcmc_plot(RTModelContextSize, type = "trace")
# Plot effects 
mcmc_plot(RTModelContextSize, type = "intervals", prob = 0.95)
# only plot main populationlevel effects
mcmc_plot(RTModelContextSize, type = "areas", prob = 0.95, variable = "^b_", regex = TRUE)+
  geom_vline(xintercept = 0, linetype = "longdash", color = "gray")+
  xlab("RT in s")+
  theme_minimal()+
  theme(text = element_text(size = 20, family = "sans"))+
  #scale_y_discrete(labels = c("Intercept", "Sequence Length", "Query Structure", "Sequence Structure", "Sort Task"))+# 
  ggtitle("Context Size RT analysis")
pp_check(RTModelContextSize, type = "ecdf_overlay")
pp_check(RTModelContextSize, type = "error_hist", ndraws = 11)
pp_check(RTModelContextSize, type = "scatter_avg", ndraws = 100)
#visualize all effects
plot(conditional_effects(RTModelContextSize), points = FALSE)
tab_model(RTModelContextSize, transform = NULL)


RTModelContextSizeRough <- brm(rt ~ RoughContextSize + (RoughContextSize|subject_id), 
               data = correctResponses,
               iter = 10000,
               cores = 4,
               save_pars = save_pars(all = TRUE),
               seed = seed,
               control = list(max_treedepth = 15, adapt_delta = 0.99),
               file = "exp3RTRoughContextSizeLog",
               family = lognormal())

summary(RTModelContextSizeRough)

# results Log  27.06.2024
#                           Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept                    1.05      0.06     0.94     1.17 1.00     2389     5270
# RoughContextSize2_medium     0.15      0.05     0.05     0.24 1.00    15915    15531
# RoughContextSize3_large      0.11      0.03     0.05     0.17 1.00    24969    15881

# Check for convergence
mcmc_plot(RTModelContextSizeRough, type = "trace")
# Plot effects 
mcmc_plot(RTModelContextSizeRough, type = "intervals", prob = 0.95)
# only plot main populationlevel effects
mcmc_plot(RTModelContextSizeRough, type = "areas", prob = 0.95, variable = "^b_", regex = TRUE)+
  geom_vline(xintercept = 0, linetype = "longdash", color = "gray")+
  xlab("RT in s")+
  theme_minimal()+
  theme(text = element_text(size = 20, family = "sans"))+
  #scale_y_discrete(labels = c("Intercept", "Sequence Length", "Query Structure", "Sequence Structure", "Sort Task"))+# 
  ggtitle("Rough Context Size RT Analysis")
pp_check(RTModelContextSizeRough, type = "ecdf_overlay")
pp_check(RTModelContextSizeRough, type = "error_hist", ndraws = 11)
pp_check(RTModelContextSizeRough, type = "scatter_avg", ndraws = 100)
#visualize all effects
plot(conditional_effects(RTModelContextSizeRough), points = FALSE)
tab_model(RTModelContextSizeRough, transform = NULL)

# Check if there is a diffrence between medium and large contexts

# Extract posterior samples
posterior <- as_draws_df(RTModelContextSizeRough)

# Compute the difference between the two levels
posterior$diff_medium_large <- posterior$b_RoughContextSize2_medium - 
  posterior$b_RoughContextSize3_large

# Summary of the difference
summary_diff <- summary(posterior$diff_medium_large)
print(summary_diff)

# Plot the posterior distribution of the difference
mcmc_areas(posterior, pars = "diff_medium_large", prob = 0.95) +
  ggtitle("Posterior Distribution of Difference (Medium - Large Context)")

# Calculate the probability that the difference is greater than zero
p_diff_positive <- mean(posterior$diff_medium_large > 0)
p_diff_negative <- mean(posterior$diff_medium_large < 0)

cat("Probability difference > 0:", p_diff_positive, "\n")
cat("Probability difference < 0:", p_diff_negative, "\n")

# Compute the Bayes Factor (BF) for the hypothesis that the difference > 0
BF <- p_diff_positive / p_diff_negative
cat("Bayes Factor (Difference > 0):", BF, "\n")

# Calculate the mean difference in RTs between the small and medium and the mdeium and large context sizes
meanRTsmall = mean(subset(correctResponses, RoughContextSize == "1_small")$rt)
meanRTmedium = mean(subset(correctResponses, RoughContextSize == "2_medium")$rt)
meanRTlarge = mean(subset(correctResponses, RoughContextSize == "3_large")$rt)

# Calculate the difference
diffSmallMedium = meanRTsmall - meanRTmedium
diffMediumLarge = meanRTmedium - meanRTlarge
cat("Mean RT difference small - medium:", diffSmallMedium, "\n")
cat("Mean RT difference medium - large:", diffMediumLarge, "\n")



# Comparison with Null model
RTModelContextSizeNull <- brm(rt ~ 1 + (ContextSizeMCScaled|subject_id), 
                          data = correctResponses,
                          iter = 10000,
                          cores = 4,
                          save_pars = save_pars(all = TRUE),
                          seed = seed,
                          control = list(max_treedepth = 15, adapt_delta = 0.99),
                          file = "exp3RTContextSizelogScaledNull",
                          family = lognormal()
)
# calculate Bayes Factors
bfRTContextSize = bayes_factor(RTModelContextSize, RTModelContextSizeNull)
bfRTContextSize
# 1.21595

RTModelContextSizeRoughNull <- brm(rt ~ 1 + (RoughContextSize|subject_id), 
                               data = correctResponses,
                               iter = 10000,
                               cores = 4,
                               save_pars = save_pars(all = TRUE),
                               seed = seed,
                               control = list(max_treedepth = 15, adapt_delta = 0.99),
                               file = "exp3RTRoughContextSizeLogNull",
                               family = lognormal())

summary(RTModelContextSizeRoughNull)

# caluclate Bayes Factors
bfRTContextSizeRough = bayes_factor(RTModelContextSizeRough, RTModelContextSizeRoughNull)
bfRTContextSizeRough

# BRMS RT: Control condition -------------------------------------------------

correctResponsesNone = subset(dResponses, correct == TRUE & 
                            ContextSize == "None" &
                            validity_rt == TRUE)

ModelRTControl <- brm(rt ~ NoneListType + (1|subject_id), 
                               data = correctResponsesNone,
                               iter = 10000,
                               cores = 4,
                               save_pars = save_pars(all = TRUE),
                               seed = seed,
                               control = list(max_treedepth = 15, adapt_delta = 0.99),
                               file = "exp3RTControlConditionsLog",
                               family = lognormal())

summary(ModelRTControl)

# results Log  27.06.2024
# Population-Level Effects: 
#                          Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept                  0.58      0.09     0.41     0.75 1.00     1670     2970
# NoneListTypelowSimList     0.38      0.12     0.14     0.62 1.00     1594     3190

# Check for convergence
mcmc_plot(ModelRTControl, type = "trace")
# Plot effects 
mcmc_plot(ModelRTControl, type = "intervals", prob = 0.95)
# only plot main populationlevel effects
mcmc_plot(ModelRTControl, type = "areas", prob = 0.95, variable = "^b_", regex = TRUE)+
  geom_vline(xintercept = 0, linetype = "longdash", color = "gray")+
  xlab("RT in s")+
  theme_minimal()+
  theme(text = element_text(size = 20, family = "sans"))+
  #scale_y_discrete(labels = c("Intercept", "Sequence Length", "Query Structure", "Sequence Structure", "Sort Task"))+# 
  ggtitle("Control Condition Size RT Analysis")
pp_check(ModelRTControl, type = "ecdf_overlay")
pp_check(ModelRTControl, type = "error_hist", ndraws = 11)
pp_check(ModelRTControl, type = "scatter_avg", ndraws = 100)
#visualize all effects
plot(conditional_effects(ModelRTControl), points = FALSE)
tab_model(ModelRTControl, transform = NULL)

# compare to null model
ModelRTControlNull <- brm(rt ~ 1 + (1|subject_id), 
                               data = correctResponsesNone,
                               iter = 10000,
                               cores = 4,
                               save_pars = save_pars(all = TRUE),
                               seed = seed,
                               control = list(max_treedepth = 15, adapt_delta = 0.99),
                               file = "exp3RTControlConditionsLogNull",
                               family = lognormal())

# calculate BF
bayes_factor(ModelRTControl, ModelRTControlNull)

# calculate mean diffrence
mean(subset(correctResponsesNone, NoneListType == "lowSimList")$rt) - 
  mean(subset(correctResponsesNone, NoneListType == "highSimPairs")$rt)

###############################################################################
# Accuracy -------
###############################################################################

# BRMS Models: Accuracy Context Size -----------------------------------------

Responses = subset(dResponses, ContextSize != "None")
Responses$ContextSize = as.numeric(Responses$ContextSize)
Responses$ContextSizeMCScaled <- scale(Responses$ContextSize, scale = TRUE)

AccModelContextSize <- brm(correct_01 ~ ContextSizeMCScaled + (ContextSizeMCScaled|subject_id), 
                          data = Responses,
                          iter = 10000,
                          cores = 4,
                          save_pars = save_pars(all = TRUE),
                          seed = seed,
                          control = list(max_treedepth = 15, adapt_delta = 0.99),
                          file = "exp3AccContextSize",
                          family = "bernoulli")

summary(AccModelContextSize)

# results 27.06.24
#                       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept               0.19      0.16   -0.13     0.51 1.01     1117     2483
# ContextSizeMCScaled    -0.27      0.04    -0.35    -0.19 1.00     4011     9120

# Check for convergence
mcmc_plot(AccModelContextSize, type = "trace")
# Plot effects 
mcmc_plot(AccModelContextSize, type = "intervals", prob = 0.95)
# only plot main populationlevel effects
mcmc_plot(AccModelContextSize, type = "areas", prob = 0.95, variable = "^b_", regex = TRUE)+
  geom_vline(xintercept = 0, linetype = "longdash", color = "gray")+
  xlab("RT in s")+
  theme_minimal()+
  theme(text = element_text(size = 20, family = "sans"))+
  #scale_y_discrete(labels = c("Intercept", "Sequence Length", "Query Structure", "Sequence Structure", "Sort Task"))+# 
  ggtitle("Context Size Accuracy Analysis")
pp_check(AccModelContextSize, type = "ecdf_overlay")
pp_check(AccModelContextSize, type = "error_hist", ndraws = 11)
pp_check(AccModelContextSize, type = "scatter_avg", ndraws = 100)
#visualize all effects
plot(conditional_effects(AccModelContextSize), points = FALSE)
tab_model(AccModelContextSize, transform = NULL)


AccModelRoughContextSize <- brm(correct_01 ~ RoughContextSize + (RoughContextSize|subject_id), 
                          data = Responses,
                          iter = 10000,
                          cores = 4,
                          save_pars = save_pars(all = TRUE),
                          seed = seed,
                          control = list(max_treedepth = 15, adapt_delta = 0.99),
                          file = "exp3AccRoughContextSize",
                          family = "bernoulli")
summary(AccModelRoughContextSize)

# results  27.11.2024
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept                    0.93      0.17     0.61     1.27 1.00     3773     7570
# RoughContextSize2_medium    -0.68      0.19    -1.05    -0.31 1.00    10810    13998
# RoughContextSize3_large     -0.85      0.12    -1.08    -0.60 1.00    11683    13132

# Check for convergence
mcmc_plot(AccModelRoughContextSize, type = "trace")
# Plot effects 
mcmc_plot(AccModelRoughContextSize, type = "intervals", prob = 0.95)
# only plot main populationlevel effects
mcmc_plot(AccModelRoughContextSize, type = "areas", prob = 0.95, variable = "^b_", regex = TRUE)+
  geom_vline(xintercept = 0, linetype = "longdash", color = "gray")+
  xlab("RT in s")+
  theme_minimal()+
  theme(text = element_text(size = 20, family = "sans"))+
  #scale_y_discrete(labels = c("Intercept", "Sequence Length", "Query Structure", "Sequence Structure", "Sort Task"))+# 
  ggtitle("Rough Context Size Accuracy Analysis")
pp_check(AccModelRoughContextSize, type = "ecdf_overlay")
pp_check(AccModelRoughContextSize, type = "error_hist", ndraws = 11)
pp_check(AccModelRoughContextSize, type = "scatter_avg", ndraws = 100)
#visualize all effects
plot(conditional_effects(AccModelRoughContextSize), points = FALSE)
tab_model(AccModelRoughContextSize, transform = NULL)


# Check if there is a diffrence between medium and large contexts

# Extract posterior samples
posterior <- as_draws_df(AccModelRoughContextSize)

# Compute the difference between the two levels
posterior$diff_medium_large <- posterior$b_RoughContextSize2_medium - 
  posterior$b_RoughContextSize3_large

# Summary of the difference
summary_diff <- summary(posterior$diff_medium_large)
print(summary_diff)

# Plot the posterior distribution of the difference
mcmc_areas(posterior, pars = "diff_medium_large", prob = 0.95) +
  ggtitle("Posterior Distribution of Difference (Medium - Large Context)")

# Calculate the probability that the difference is greater than zero
p_diff_positive <- mean(posterior$diff_medium_large > 0)
p_diff_negative <- mean(posterior$diff_medium_large < 0)

cat("Probability difference > 0:", p_diff_positive, "\n")
cat("Probability difference < 0:", p_diff_negative, "\n")

# Compute the Bayes Factor (BF) for the hypothesis that the difference > 0
BF <- p_diff_positive / p_diff_negative
cat("Bayes Factor (Difference > 0):", BF, "\n")

# Calculate the mean difference in RTs between the small and medium and the mdeium and large context sizes
meanAccsmall = mean(subset(Responses, RoughContextSize == "1_small")$correct_01)
meanAccmedium = mean(subset(Responses, RoughContextSize == "2_medium")$correct_01)
meanAcclarge = mean(subset(Responses, RoughContextSize == "3_large")$correct_01)

# Calculate the difference
diffSmallMedium = meanAccsmall - meanAccmedium
diffMediumLarge = meanAccmedium - meanAcclarge
cat("Mean accuracy difference small - medium:", diffSmallMedium, "\n")
cat("Mean accuracy difference medium - large:", diffMediumLarge, "\n")



AccModelContextSizeNull <- brm(correct_01 ~ 1 + (ContextSizeMCScaled|subject_id), 
                           data = Responses,
                           iter = 10000,
                           cores = 4,
                           save_pars = save_pars(all = TRUE),
                           seed = seed,
                           control = list(max_treedepth = 15, adapt_delta = 0.99),
                           file = "exp3AccContextSizeNull",
                           family = "bernoulli")

summary(AccModelContextSizeNull)

# Calculate Bayes Factors
bfAccContextSize = bayes_factor(AccModelContextSize, AccModelContextSizeNull)
bfAccContextSize

AccModelRoughContextSizeNull <- brm(correct_01 ~ 1 + (RoughContextSize|subject_id), 
                                data = Responses,
                                iter = 10000,
                                cores = 4,
                                save_pars = save_pars(all = TRUE),
                                seed = seed,
                                control = list(max_treedepth = 15, adapt_delta = 0.99),
                                file = "exp3AccRoughContextSizeNull",
                                family = "bernoulli")
summary(AccModelRoughContextSize)

bfAccContextSize2 = bayes_factor(AccModelRoughContextSize, AccModelRoughContextSizeNull)
bfAccContextSize2

# BRMS Models: Accuracy Control Condition--------------------------------------

ResponsesNone = subset(dResponses, ContextSize == "None")

ModelAccuracyControl <- brm(correct_01 ~ NoneListType + (1|subject_id), 
                          data = ResponsesNone,
                          iter = 10000,
                          cores = 4,
                          chains = 4,
                          save_pars = save_pars(all = TRUE),
                          seed = seed,
                          control = list(max_treedepth = 15, adapt_delta = 0.99),
                          file = "exp3AccControlCondition",
                          family = "bernoulli")
summary(ModelAccuracyControl)

# results 27.06.2024
#                         Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept                  2.04      0.22     1.63     2.47 1.00     5805     9709
# NoneListTypelowSimList    -1.47      0.29    -2.06    -0.90 1.00     4953     7944


# Check for convergence
mcmc_plot(ModelAccuracyControl, type = "trace")
# Plot effects 
mcmc_plot(ModelAccuracyControl, type = "intervals", prob = 0.95)
# only plot main populationlevel effects
mcmc_plot(ModelAccuracyControl, type = "areas", prob = 0.95, variable = "^b_", regex = TRUE)+
  geom_vline(xintercept = 0, linetype = "longdash", color = "gray")+
  xlab("RT in s")+
  theme_minimal()+
  theme(text = element_text(size = 20, family = "sans"))+
  #scale_y_discrete(labels = c("Intercept", "Sequence Length", "Query Structure", "Sequence Structure", "Sort Task"))+# 
  ggtitle("Control Condition Accuracy Analysis")
pp_check(ModelAccuracyControl, type = "ecdf_overlay")
pp_check(ModelAccuracyControl, type = "error_hist", ndraws = 11)
pp_check(ModelAccuracyControl, type = "scatter_avg", ndraws = 100)
#visualize all effects
plot(conditional_effects(ModelAccuracyControl), points = FALSE)
tab_model(ModelAccuracyControl, transform = NULL)

ModelAccuracyControlNull <- brm(correct_01 ~ 1 + (1|subject_id), 
                            data = ResponsesNone,
                            iter = 10000,
                            cores = 4,
                            chains = 4,
                            save_pars = save_pars(all = TRUE),
                            seed = seed,
                            control = list(max_treedepth = 15, adapt_delta = 0.99),
                            file = "exp3AccControlConditionNull",
                            family = "bernoulli")

bayes_factor(ModelAccuracyControl, ModelAccuracyControlNull)

# calculate the mean difference in accuracy 
mean(subset(ResponsesNone, NoneListType == "lowSimList")$correct_01) - 
  mean(subset(ResponsesNone, NoneListType == "highSimPairs")$correct_01)

###############################################################################
# Hypothesis 2 the benefits of context cues dependednt on context size -------
###############################################################################
# RT -------------------------------------------------------------------------

# Calculate mean RT for each subject and rough context size as well as the "none" 
# context size condition i.e. based on the RoughContextSize column
meanConditionRTs = ddply(subset(dResponses, correct == TRUE & 
                                  validity_rt == TRUE), 
                         c("subject_id", "RoughContextSize", "NoneListType"), summarise,
                         meanRT = mean(rt),
                         seRT = se(rt))

# calculate the difference between the mean RTs of the "none" condition and the small, medium and large context size conditions
# for each participan
# make empty dataframe with rows for subject id, control condition small, medium and large RT difference
RTDifferences = data.frame(subject_id = character(),
                           difference = numeric(),
                           type = character(),
                           NoneListType = character(),
                           ContextRT = numeric(),
                           NoneContextRT = numeric())
for (subject in unique(meanConditionRTs$subject_id)){
  noneListType = subset(meanConditionRTs, subject_id == subject)$NoneListType[1]
  if (noneListType == "lowSimList"){
    noneRT = subset(meanConditionRTs, subject_id == subject & RoughContextSize == "0_NoneLowSim")$meanRT
  }
  else{
    noneRT = subset(meanConditionRTs, subject_id == subject & RoughContextSize == "0_NoneHighSim")$meanRT
  }
  # check if RoughContextSize 1_small is present
  if ("1_small" %in% subset(meanConditionRTs, subject_id == subject)$RoughContextSize){
    smallRT = subset(meanConditionRTs, subject_id == subject & RoughContextSize == "1_small")$meanRT
    # calculate the difference
    diff = smallRT - noneRT
    # add row to RTDifferences
    RTDifferences = rbind(RTDifferences, data.frame(subject_id = subject,
                                                     difference = diff,
                                                     type = "small",
                                                     NoneListType = noneListType,
                                                     ContextRT = smallRT,
                                                     NoneContextRT = noneRT))
  }
  #check if RoughContextSize 2_medium is present
  if ("2_medium" %in% subset(meanConditionRTs, subject_id == subject)$RoughContextSize){
    mediumRT = subset(meanConditionRTs, subject_id == subject & RoughContextSize == "2_medium")$meanRT
    # calculate the difference
    diff = mediumRT - noneRT
    # add row to RTDifferences
    RTDifferences = rbind(RTDifferences, data.frame(subject_id = subject,
                                                     difference = diff,
                                                     type = "medium",
                                                     NoneListType = noneListType,
                                                     ContextRT = mediumRT,
                                                     NoneContextRT = noneRT))
  }
  #check if RoughContextSize 3_large is present
  if ("3_large" %in% subset(meanConditionRTs, subject_id == subject)$RoughContextSize){
    largeRT = subset(meanConditionRTs, subject_id == subject & RoughContextSize == "3_large")$meanRT
    # calculate the difference
    diff = largeRT - noneRT
    # add row to RTDifferences
    RTDifferences = rbind(RTDifferences, data.frame(subject_id = subject,
                                                     difference = diff,
                                                     type = "large",
                                                     NoneListType = noneListType,
                                                     ContextRT = largeRT,
                                                     NoneContextRT = noneRT))
  }

}

# Plot the Rt diffreneces by type and none list type as facet grid with ggplot
ggplot(RTDifferences, aes(x = type, y = difference, color = type))+
  geom_boxplot()+
  geom_jitter(width = 0.2)+
  facet_grid(cols = vars(NoneListType))

# Do the same for accuracy ---------------------------------------------------

# Calculate mean accuracy for each subject and rough context size as well as the "none"
# context size condition i.e. based on the RoughContextSize column
meanConditionAcc = ddply(dResponses, 
                         c("subject_id", "RoughContextSize", "NoneListType"), summarise,
                         meanAcc = mean(correct_01),
                         seAcc = binary_se(correct_01))

# calculate the difference between the mean accuracy of the "none" condition and the small, medium and large context size conditions
# for each participan
# make empty dataframe with rows for subject id, control condition small, medium and large RT difference
AccDifferences = data.frame(subject_id = character(),
                            difference = numeric(),
                            type = character(),
                            NoneListType = character(),
                            ContextAcc = numeric(),
                            NoneContextAcc = numeric())
for (subject in unique(meanConditionAcc$subject_id)){
  noneListType = subset(meanConditionAcc, subject_id == subject)$NoneListType[1]
  if (noneListType == "lowSimList"){
    noneAcc = subset(meanConditionAcc, subject_id == subject & RoughContextSize == "0_NoneLowSim")$meanAcc
  }
  else{
    noneAcc = subset(meanConditionAcc, subject_id == subject & RoughContextSize == "0_NoneHighSim")$meanAcc
  }
  # check if RoughContextSize 1_small is present
  if ("1_small" %in% subset(meanConditionAcc, subject_id == subject)$RoughContextSize){
    smallAcc = subset(meanConditionAcc, subject_id == subject & RoughContextSize == "1_small")$meanAcc
    # calculate the difference
    diff = smallAcc - noneAcc
    # add row to AccDifferences
    AccDifferences = rbind(AccDifferences, data.frame(subject_id = subject,
                                                       difference = diff,
                                                       type = "small",
                                                       NoneListType = noneListType,
                                                       ContextAcc = smallAcc,
                                                       NoneContextAcc = noneAcc))
  }
  #check if RoughContextSize 2_medium is present
  if ("2_medium" %in% subset(meanConditionAcc, subject_id == subject)$RoughContextSize){
    mediumAcc = subset(meanConditionAcc, subject_id == subject & RoughContextSize == "2_medium")$meanAcc
    # calculate the difference
    diff = mediumAcc - noneAcc
    # add row to AccDifferences
    AccDifferences = rbind(AccDifferences, data.frame(subject_id = subject,
                                                       difference = diff,
                                                       type = "medium",
                                                       NoneListType = noneListType,
                                                       ContextAcc = mediumAcc,
                                                       NoneContextAcc = noneAcc))
  }
  #check if RoughContextSize 3_large is present
  if ("3_large" %in% subset(meanConditionAcc, subject_id == subject)$RoughContextSize){
    largeAcc = subset(meanConditionAcc, subject_id == subject & RoughContextSize == "3_large")$meanAcc
    # calculate the difference
    diff = largeAcc - noneAcc
    # add row to AccDifferences
    AccDifferences = rbind(AccDifferences, data.frame(subject_id = subject,
                                                       difference = diff,
                                                       type = "large",
                                                       NoneListType = noneListType,
                                                       ContextAcc = largeAcc,
                                                       NoneContextAcc = noneAcc))
  }
}

# Plot the accuracy differences by type and none list type as facet grid with ggplot
ggplot(AccDifferences, aes(x = type, y = difference, color = type))+
  geom_boxplot()+
  geom_jitter(width = 0.2)+
  facet_grid(cols = vars(NoneListType))

###############################################################################
# Perform the t-tests and store the p-values ----------------------------------
###############################################################################
# Load the BayesFactor package
library(BayesFactor)

# Calculate Bayes Factors
BFs <- list(
  ttestBF(subset(RTDifferences, type == "small" & NoneListType == "highSimPairs")$ContextRT, 
          subset(RTDifferences, type == "small" & NoneListType == "highSimPairs")$NoneContextRT,
          paired = TRUE),
  ttestBF(subset(RTDifferences, type == "medium" & NoneListType == "highSimPairs")$ContextRT, 
          subset(RTDifferences, type == "medium" & NoneListType == "highSimPairs")$NoneContextRT,
          paired = TRUE),
  ttestBF(subset(RTDifferences, type == "large" & NoneListType == "highSimPairs")$ContextRT, 
          subset(RTDifferences, type == "large" & NoneListType == "highSimPairs")$NoneContextRT,
          paired = TRUE),
  ttestBF(subset(RTDifferences, type == "small" & NoneListType == "lowSimList")$ContextRT, 
          subset(RTDifferences, type == "small" & NoneListType == "lowSimList")$NoneContextRT,
          paired = TRUE),
  ttestBF(subset(RTDifferences, type == "medium" & NoneListType == "lowSimList")$ContextRT, 
          subset(RTDifferences, type == "medium" & NoneListType == "lowSimList")$NoneContextRT,
          paired = TRUE),
  ttestBF(subset(RTDifferences, type == "large" & NoneListType == "lowSimList")$ContextRT, 
          subset(RTDifferences, type == "large" & NoneListType == "lowSimList")$NoneContextRT,
          paired = TRUE),
  ttestBF(subset(AccDifferences, type == "small" & NoneListType == "highSimPairs")$ContextAcc, 
          subset(AccDifferences, type == "small" & NoneListType == "highSimPairs")$NoneContextAcc,
          paired = TRUE),
  ttestBF(subset(AccDifferences, type == "medium" & NoneListType == "highSimPairs")$ContextAcc, 
          subset(AccDifferences, type == "medium" & NoneListType == "highSimPairs")$NoneContextAcc,
          paired = TRUE),
  ttestBF(subset(AccDifferences, type == "large" & NoneListType == "highSimPairs")$ContextAcc, 
          subset(AccDifferences, type == "large" & NoneListType == "highSimPairs")$NoneContextAcc,
          paired = TRUE),
  ttestBF(subset(AccDifferences, type == "small" & NoneListType == "lowSimList")$ContextAcc, 
          subset(AccDifferences, type == "small" & NoneListType == "lowSimList")$NoneContextAcc,
          paired = TRUE),
  ttestBF(subset(AccDifferences, type == "medium" & NoneListType == "lowSimList")$ContextAcc, 
          subset(AccDifferences, type == "medium" & NoneListType == "lowSimList")$NoneContextAcc,
          paired = TRUE),
  ttestBF(subset(AccDifferences, type == "large" & NoneListType == "lowSimList")$ContextAcc, 
          subset(AccDifferences, type == "large" & NoneListType == "lowSimList")$NoneContextAcc,
          paired = TRUE)
)

# Extract numeric Bayes Factors
BF_values <- sapply(BFs, function(x) extractBF(x)$bf)


# Do pairwise t test on the differnces between the context size conditions and the "none" conditions
# just doing this as a sanity check
p_values <- c(
  t.test(subset(RTDifferences, type == "small" & NoneListType == "highSimPairs")$difference)$p.value,
  t.test(subset(RTDifferences, type == "medium" & NoneListType == "highSimPairs")$difference)$p.value,
  t.test(subset(RTDifferences, type == "large" & NoneListType == "highSimPairs")$difference)$p.value,
  t.test(subset(RTDifferences, type == "small" & NoneListType == "lowSimList")$difference)$p.value,
  t.test(subset(RTDifferences, type == "medium" & NoneListType == "lowSimList")$difference)$p.value,
  t.test(subset(RTDifferences, type == "large" & NoneListType == "lowSimList")$difference)$p.value,
  t.test(subset(AccDifferences, type == "small" & NoneListType == "highSimPairs")$difference)$p.value,
  t.test(subset(AccDifferences, type == "medium" & NoneListType == "highSimPairs")$difference)$p.value,
  t.test(subset(AccDifferences, type == "large" & NoneListType == "highSimPairs")$difference)$p.value,
  t.test(subset(AccDifferences, type == "small" & NoneListType == "lowSimList")$difference)$p.value,
  t.test(subset(AccDifferences, type == "medium" & NoneListType == "lowSimList")$difference)$p.value,
  t.test(subset(AccDifferences, type == "large" & NoneListType == "lowSimList")$difference)$p.value
)

mean_differences = c(
  mean(subset(RTDifferences, type == "small" & NoneListType == "highSimPairs")$difference),
  mean(subset(RTDifferences, type == "medium" & NoneListType == "highSimPairs")$difference),
  mean(subset(RTDifferences, type == "large" & NoneListType == "highSimPairs")$difference),
  mean(subset(RTDifferences, type == "small" & NoneListType == "lowSimList")$difference),
  mean(subset(RTDifferences, type == "medium" & NoneListType == "lowSimList")$difference),
  mean(subset(RTDifferences, type == "large" & NoneListType == "lowSimList")$difference),
  mean(subset(AccDifferences, type == "small" & NoneListType == "highSimPairs")$difference),
  mean(subset(AccDifferences, type == "medium" & NoneListType == "highSimPairs")$difference),
  mean(subset(AccDifferences, type == "large" & NoneListType == "highSimPairs")$difference),
  mean(subset(AccDifferences, type == "small" & NoneListType == "lowSimList")$difference),
  mean(subset(AccDifferences, type == "medium" & NoneListType == "lowSimList")$difference),
  mean(subset(AccDifferences, type == "large" & NoneListType == "lowSimList")$difference)
)


# Adjust the p-values using the desired method
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")  # You can also try "BH" for Benjamini-Hochberg

# Combine the results for easy viewing
results <- data.frame(
  Comparison = c(
    "Small & HighSimPairs RT",
    "Medium & HighSimPairs RT",
    "Large & HighSimPairs RT",
    "Small & LowSimList RT",
    "Medium & LowSimList RT",
    "Large & LowSimList RT",
    "Small & HighSimPairs Acc",
    "Medium & HighSimPairs Acc",
    "Large & HighSimPairs Acc",
    "Small & LowSimList Acc",
    "Medium & LowSimList Acc",
    "Large & LowSimList Acc"
    
  ),
  mean_difference = mean_differences,
  #P_Value = p_values,
  Adjusted_P_Value = adjusted_p_values,
  BayesFactor = BF_values
)

print(results)
###############################################################################


# Errortype analysis ---------------------------------------------------------
current_directory <- dirname(rstudioapi::getSourceEditorContext()$path)
subdirectory = "/Experiment3"
setwd(paste(current_directory, subdirectory, sep = ""))
tempData = read.csv("ExperimentDataExp3ErrorType.csv")
dResponses <- subset(tempData, Phase == "CuePresentation" & rt < 20 & validity == TRUE)
dResponsesErrorSummary = ddply(dResponses, ~ subject_id + NoneListType + RoughContextSize + Error_type, summarise, Number = length(Error_type))

# Calculate the proportion of each error type per subject and RoughContextSize
# Step 1: Summarize total responses per subject and RoughContextSize
dTotalResponses <- ddply(dResponses, ~ subject_id + RoughContextSize, summarise,
                         Total = length(Error_type))

# Step 2: Merge the total responses with the error summary
dResponsesErrorSummary <- merge(dResponsesErrorSummary, dTotalResponses, 
                                by = c("subject_id", "RoughContextSize"))

# Step 3: Calculate the proportion for each response type
dResponsesErrorSummary$Proportion <- dResponsesErrorSummary$Number / dResponsesErrorSummary$Total

# Now add intrusions if they dont already exist


# Step 1: Identify all unique combinations of subject_id, RoughContextSize
unique_subject_context <- unique(dResponsesErrorSummary[, c("subject_id", "RoughContextSize")])

# Add intrusion proportion (0) if they are missing
missing_responses <- data.frame(subject_id = character(),
                                 RoughContextSize = character(),
                                 Error_type = character(),
                                 Number = numeric(),
                                 Proportion = numeric(),
                                 Total = numeric(),
                                 NoneListType = character())
# Chekc if error type "Intrusion" is present for each unique combination
for (i in 1:nrow(unique_subject_context)) {
  subject_id_1 <- unique_subject_context$subject_id[i]
  RoughContextSize_1 <- unique_subject_context$RoughContextSize[i]
  relevantSubset <- subset(dResponsesErrorSummary, subject_id == subject_id_1 & RoughContextSize == RoughContextSize_1)
  if (!("Intrusion" %in% relevantSubset$Error_type)) {
    missing_responses <- rbind(missing_responses, data.frame(subject_id = subject_id_1,
                                                               RoughContextSize = RoughContextSize_1,
                                                               Error_type = "Intrusion",
                                                               Number = 0,
                                                               Proportion = 0,
                                                               Total = relevantSubset$Total[1],
                                                               NoneListType = relevantSubset$NoneListType[1]))
  }
  # now check for omissions
  if (!("Omission" %in% relevantSubset$Error_type)) {
    missing_responses <- rbind(missing_responses, data.frame(subject_id = subject_id_1,
                                                               RoughContextSize = RoughContextSize_1,
                                                               Error_type = "Omission",
                                                               Number = 0,
                                                               Proportion = 0,
                                                               Total = relevantSubset$Total[1],
                                                               NoneListType = relevantSubset$NoneListType[1]))
  } 
}

dResponsesErrorSummary <- combined_df <- rbind(dResponsesErrorSummary, missing_responses)

# do a pairwise baysian ttest for intrusions btween large contexts and none contexts
LargeContextIntrusion = subset(dResponsesErrorSummary, Error_type == "Intrusion" & RoughContextSize == "3_large" & NoneListType == "lowSimList")
LowSimContextIntrusion = subset(dResponsesErrorSummary, Error_type == "Intrusion" & RoughContextSize == "0_NoneLowSim" & NoneListType == "lowSimList")

intrusionBF = ttestBF(LargeContextIntrusion$Proportion, 
                      LowSimContextIntrusion$Proportion,
                      paired = TRUE)
intrusionBF

# Bayes factor analysis
# --------------
#   [1] Alt., r=0.707 : 122.7256 ±0%
# 
# Against denominator:
#   Null, mu = 0 
# ---
#   Bayes factor type: BFoneSample, JZS

# calculate average intrusion proportion
mean(subset(dResponsesErrorSummary, Error_type == "Intrusion" & RoughContextSize == "3_large" & NoneListType == "lowSimList")$Proportion)- 
mean(subset(dResponsesErrorSummary, Error_type == "Intrusion" & RoughContextSize == "0_NoneLowSim" & NoneListType == "lowSimList")$Proportion)

# now do the same for omissions
LargeContextOmissions = subset(dResponsesErrorSummary, Error_type == "Omission" & RoughContextSize == "3_large" & NoneListType == "lowSimList")
LowSimContextOmissions = subset(dResponsesErrorSummary, Error_type == "Omission" & RoughContextSize == "0_NoneLowSim" & NoneListType == "lowSimList")

omissionBF = ttestBF(LargeContextOmissions$Proportion, 
                      LowSimContextOmissions$Proportion,
                      paired = TRUE)
omissionBF
# Bayes factor analysis
# --------------
#   [1] Alt., r=0.707 : 0.3755447 ±0.04%
# 
# Against denominator:
#   Null, mu = 0 
# ---
#   Bayes factor type: BFoneSample, JZS

# calculate average omission proportion
mean(subset(dResponsesErrorSummary, Error_type == "Omission" & RoughContextSize == "3_large" & NoneListType == "lowSimList")$Proportion)-
mean(subset(dResponsesErrorSummary, Error_type == "Omission" & RoughContextSize == "0_NoneLowSim" & NoneListType == "lowSimList")$Proportion)



