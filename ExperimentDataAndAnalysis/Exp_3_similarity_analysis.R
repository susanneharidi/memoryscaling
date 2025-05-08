

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
library(broom)
library(patchwork)

# SET WD:
# If you are not using the here package and R projects, this is the place to 
# set the working directory to the inside of the github project


# Fuctions --------------------------------------------------------------------

# standard error function
se <- function(x){sd(x)/sqrt(length(x))}


# Variables ------------------------------------------------------------------

seed = 2024

textsize = 22
font = "sans"

# load all data from experiment2 ----------------------------------------------

current_directory <- dirname(rstudioapi::getSourceEditorContext()$path)
subdirectory = "/Experiment3"
setwd(paste(current_directory, subdirectory, sep = ""))

Data = read.csv("ExperimentDataExp3W2VSim.csv")
dResponses <- subset(Data, Phase == "CuePresentation" & rt < 20 & validity == "True")
# Plot semantic similarity against RTs ----------------------------------------
correctResponses = subset(dResponses, correct_01 == 1 & 
                            !is.na(ContextSize) &
                            validity_rt == "True")

ggplot(correctResponses, aes(x = W2VWordPairsSim, y = rt)) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Semantic similarity vs RTs",
       x = "Semantic similarity",
       y = "RTs") +
  theme_minimal() +
  coord_cartesian(ylim = c(0, 6)) +
  theme(legend.position = "top")



# Fit a hierarchical model: similarity predicting RT, varying by subject_id
brm_model <- brm(
  rt ~ W2VWordPairsSim * as.factor(ContextSize) + (1 + W2VWordPairsSim + as.factor(ContextSize)| subject_id),
  data = correctResponses,
  cores = 4, chains = 4, iter = 2000,
  control = list(adapt_delta = 0.95)
)

# Extract fixed effects into a data frame
fixed_effects <- as.data.frame(fixef(brm_model, summary = TRUE))

# Add term names as a column
fixed_effects$term <- rownames(fixed_effects)
rownames(fixed_effects) <- NULL

# Grab the baseline slope
baseline_slope <- fixed_effects %>%
  filter(term == "W2VWordPairsSim") %>%
  select(Estimate, Q2.5, Q97.5) %>%
  mutate(ContextSize = "main effect")  


# Adjust interaction terms by adding the baseline slope
context_slopes <- fixed_effects %>%
  filter(grepl("W2VWordPairsSim:as.factorContextSize", term)) %>%
  mutate(
    ContextSize = gsub("b_W2VWordPairsSim:as.factorContextSize", "", term),
    Estimate = Estimate + baseline_slope$Estimate,
    Q2.5 = Q2.5 + baseline_slope$Q2.5,
    Q97.5 = Q97.5 + baseline_slope$Q97.5
  ) %>%
  select(ContextSize, Estimate, Q2.5, Q97.5)

summary_slopes <- bind_rows(baseline_slope, context_slopes)


slopePlot = ggplot(context_slopes, aes(x = ContextSize, y = Estimate)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = 0.3) +
  labs(x = "Context size",
      y = "Estimated slope (with 95% CI) on RT"
  ) +
  # relable the x axis ticks
  geom_hline(yintercept = 0, linetype = "dotted", color = "black") +
  scale_x_discrete(labels = c("2", "4", "6", "9", "11", "14", "16", "18", "19", "20")) +
  geom_hline(yintercept = baseline_slope$Estimate, linetype = "dashed", color = "blue") +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = baseline_slope$Q2.5, ymax = baseline_slope$Q97.5),
            fill = "blue", alpha = 0.05, inherit.aes = FALSE) +
  annotate("text", x = 1, y = baseline_slope$Estimate, label = "Main effect of similarity", 
           vjust = -1, hjust = 0, color = "blue", size = 5.5) +
  theme_minimal() +
  theme(legend.position = "bottom")+
  theme(text = element_text(size = textsize, family = font))
slopePlot



# same for accuracy --------------------------------------------
relevantResponses = subset(dResponses, !is.na(ContextSize) &
                             validity_rt == "True")
brm_modelAccuracy <- brm(
  correct_01 ~ W2VWordPairsSim * as.factor(ContextSize) + (1 + W2VWordPairsSim + as.factor(ContextSize)| subject_id),
  data = relevantResponses,
  cores = 4, chains = 4, iter = 2000,
  control = list(adapt_delta = 0.95),
  family = "bernoulli"
)

# Extract fixed effects into a data frame
fixed_effectsAcc <- as.data.frame(fixef(brm_modelAccuracy, summary = TRUE))

# Add term names as a column
fixed_effectsAcc$term <- rownames(fixed_effectsAcc)
rownames(fixed_effectsAcc) <- NULL

# Grab the baseline slope
baseline_slopeAcc <- fixed_effectsAcc %>%
  filter(term == "W2VWordPairsSim") %>%
  select(Estimate, Q2.5, Q97.5) %>%
  mutate(ContextSize = "main effect")  


# Adjust interaction terms by adding the baseline slope
context_slopesAcc <- fixed_effectsAcc %>%
  filter(grepl("W2VWordPairsSim:as.factorContextSize", term)) %>%
  mutate(
    ContextSize = gsub("b_W2VWordPairsSim:as.factorContextSize", "", term),
    Estimate = Estimate + baseline_slopeAcc$Estimate,
    Q2.5 = Q2.5 + baseline_slopeAcc$Q2.5,
    Q97.5 = Q97.5 + baseline_slopeAcc$Q97.5
  ) %>%
  select(ContextSize, Estimate, Q2.5, Q97.5)


slopePlotAcc = ggplot(context_slopesAcc, aes(x = ContextSize, y = Estimate)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = 0.3) +
  labs(x = "Context size",
       y = "Estimated slope (with 95% CI) on accuracy"
  ) +
  # relable the x axis ticks
  geom_hline(yintercept = 0, linetype = "dotted", color = "black") +
  scale_x_discrete(labels = c("2", "4", "6", "9", "11", "14", "16", "18", "19", "20")) +
  geom_hline(yintercept = baseline_slopeAcc$Estimate, linetype = "dashed", color = "blue") +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = baseline_slopeAcc$Q2.5, ymax = baseline_slopeAcc$Q97.5),
            fill = "blue", alpha = 0.05, inherit.aes = FALSE) +
  annotate("text", x = 1, y = baseline_slopeAcc$Estimate, label = "Main effect of similarity", 
           vjust = -1, hjust = 0, color = "blue", size = 5.5) +
  theme_minimal() +
  theme(legend.position = "bottom")+
  theme(text = element_text(size = textsize, family = font))
slopePlotAcc
# plot histograms of similarities -------------------------------------------
# plot histogramm of the cue target similarities
histPlot = ggplot(relevantResponses, aes(x = W2VWordPairsSim, fill = correct)) +
  geom_histogram(binwidth = 0.05, alpha = 0.7, position = "dodge") +
  labs(x = "Cue target similarity",
       y = "Frequency") +
  theme_minimal() +
  coord_cartesian(xlim = c(0, 1)) +
  theme(legend.position = "bottom")+
  theme(text = element_text(size = textsize, family = font))
histPlot

# just for control lets check on the rest
correctResponsesControl = subset(dResponses, correct_01 == 1 & 
                                   is.na(ContextSize) &
                                   validity_rt == "True")
ggplot(correctResponsesControl, aes(x = W2VWordPairsSim, fill = NoneListType)) +
  geom_histogram(binwidth = 0.1, alpha = 0.7) +
  labs(x = "Cue target similarity",
       y = "Frequency") +
  theme_minimal() +
  coord_cartesian(xlim = c(0, 1)) +
  theme(legend.position = "top")


# make compined figure --------------------------------
combinedPlot = slopePlot + slopePlotAcc + histPlot +
  plot_layout(ncol = 3) +
  plot_annotation(tag_levels = "A")+ 
  plot_annotation(theme = theme(plot.title = element_text(hjust = 0.5)))
combinedPlot


savingPlace = "/Figures/Figure3_5.pdf"
ggsave(paste(current_directory, savingPlace, sep = ""), combinedPlot,  
       width = 30, 
       height = 10)
