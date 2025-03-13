
# Script with model predictions:

# Load libraries ---------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# set seed for reproducibility ------------------------------------------------

set.seed(123)

#SimSS predictions ------------------------------------------------------------

# create functions for the model


modelOneRun = function(sim_cue_target, p_correct, p_incorrect, p_omission, p_intrusion){
  Times = c()
  Response_types = c()
  
  for(i in 1:length(sim_cue_target)){
    # Sample 1 if I find the target or 0 if I dont
    time = 0
    response_type = "something went wrong"
    stop = 0
    found = 0
    while(stop == 0) {
      time = time + 1
      found = sample(c(1, 0), size = 1, prob = c(p_correct[i], p_incorrect[i]))
      if(found){
        stop = 1
        response_type = "Correct"
      }else{
        omission = sample(c(1, 0), size = 1, prob = c(p_omission, 1-p_omission))
        intrusion = sample(c(1, 0), size = 1, prob = c(p_intrusion, 1-p_intrusion))
        if(omission){
          stop = 1
          response_type = "Omission"
        }else if(intrusion){
          stop = 1
          response_type = "Intrusion"
        }
      }
    }
    Times = c(Times, time)
    Response_types = c(Response_types, response_type)
  }
  return(list(Times = Times, Response_types = Response_types))
}


modelSimulation = function(Participant_data, p_omission, p_intrusion, association_strength = 0, repetitions = 1, UseSimilarity = TRUE){
  sim_cue_target = Participant_data$SimCueTarget
  sim_correct = sim_cue_target + association_strength
  sim_distractors = Participant_data$SimSumDistractors
  if (UseSimilarity){
    p_correct = sim_correct/(sim_distractors+sim_correct)
  }else{
    p_correct = ((1/Participant_data$ListLength) + association_strength)/(1+association_strength)
  }
  p_incorrect = 1-p_correct
  ListLength = Participant_data$ListLength
  Condition = Participant_data$Condition
  Experiment = Participant_data$Experiment
  ContextSize = Participant_data$ContextSize
  CueTargetSimilarity = Participant_data$SimCueTarget
  
  # Initialize lists to store the results
  all_times <- list()
  all_response_types <- list()
  
  for (i in 1:repetitions) {
    out = modelOneRun(sim_cue_target, p_correct, p_incorrect, p_omission, p_intrusion)
    Times = out$Times
    Response_types = out$Response_types
    
    # Store the results
    all_times[[i]] <- Times
    all_response_types[[i]] <- Response_types
  }
  
  # Convert lists to matrices for easier calculation
  all_times_matrix <- do.call(rbind, all_times)
  all_response_types_matrix <- do.call(rbind, all_response_types)
  
  # Calculate the average time for each entry across repetitions
  average_times <- colMeans(all_times_matrix)
  
  # Randomly choose one response type for each entry
  random_response_types <- apply(all_response_types_matrix, 2, function(x) sample(x, 1))
  
  df = data.frame(ListLength = ListLength,
                  Condition = Condition,
                  Experiment = Experiment,
                  ContextSize = ContextSize,
                  ModelTimes = average_times,
                  ModelResponseType = random_response_types,
                  CueTargetSimilarity = CueTargetSimilarity)
  return(df)
}

# Define general parameters ---------------------------------------------------
reps = 1000
p_omission = 0.02 
p_intrusion = 0.02 
association_strength = 0.2 
UseSimilarity = TRUE


# Experiment 1 -----------------------------------------------------------------

# simulate some participant data that mimic the setup of experiment 1

# Set sizes between 0 and 32 in 1 increments
ListLength = rep(1:32, reps)

# prep dataframe
df = data.frame(ListLength = numeric(),
                SimCueTarget = numeric(),
                SimSumDistractors = numeric(),
                Condition = character(),
                Experiment = character(),
                ContextSize = numeric())
# now for each Set sizes lets simulate the cue target similarity and the sum of the distractors similarity
for (i in 1:length(ListLength)){
  # cue target similarity (sample from a uniform distribution between 0 and 1)
  sim_cue_target = runif(1, min = 0, max = 1)
  # distractors similarity
  if (ListLength[i] == 1){
    sim_sum_distractors = sim_cue_target
  }else{
    sim_sum_distractors = sum(runif(ListLength[i]-1, min = 0, max = 1)) + sim_cue_target
  }
  # condition
  condition = "None"
  # experiment
  experiment = "1"
  # context size
  context_size = 0
  
  # add to the dataframe
  df = rbind(df, data.frame(ListLength = ListLength[i],
                            SimCueTarget = sim_cue_target,
                            SimSumDistractors = sim_sum_distractors,
                            Condition = condition,
                            Experiment = experiment,
                            ContextSize = context_size))
}

# simulate the model
model_predictions = modelSimulation(df, p_omission, p_intrusion, association_strength, reps, UseSimilarity)

# plot predictions Experiment 1 ------------------------------------------------

textsize = 22
font = "sans"
default_color <- "#007BB8"
AccColor <- "#74c476"
RTColor <- "#56B4E9"
RegColor = "black"
RibbonAlpha = 0.2
RT_min_value = 0
RT_max_value = 4

# crete a coding for correct and incorrect correct is everything where the model response type is correct
model_predictions$Correct = ifelse(model_predictions$ModelResponseType == "Correct", TRUE, FALSE)
model_predictions$Correct_01 = ifelse(model_predictions$ModelResponseType == "Correct", 1, 0)

# calculate mean RT for each set size
meanRTs = ddply(model_predictions, ~ ListLength + Correct, summarise, MeanRT = mean(ModelTimes))
# 
Exp1LLRT = ggplot(meanRTs, aes(x = ListLength, y = MeanRT, group = Correct, color = Correct)) +
  #geom_point() +
  geom_line(se = FALSE, size = 2, method = lm) +
  # add mean RT
  # plot mean RTs
  labs(x = "Set size",
       y = "Model times") +
  theme_minimal()+
  scale_color_manual(values = c("TRUE" = "forestgreen", "FALSE" = "red"))+
  theme(legend.position = "bottom")+
  theme(text = element_text(size = textsize, family = "sans"))
Exp1LLRT

# make sim Bins

breaks = seq(0, 1, by = 0.1)
model_predictions$simBin <- cut(model_predictions$CueTargetSimilarity, 
                                breaks = breaks, 
                                labels = breaks[-1]-0.05,
                                include.lowest = TRUE, 
                                right = FALSE)

# aggregate data across sim bins and set size
AggragatedModelTimes <- ddply(model_predictions, 
                              ~simBin + Correct, 
                              summarize, MeanModelTime = mean(ModelTimes))
AggragatedModelTimes$simBin <- as.numeric(as.character(AggragatedModelTimes$simBin))

# Plot ModelTimes based on CueTargetSimilarity
Exp1SimRT = ggplot(AggragatedModelTimes, aes(x = simBin, y = MeanModelTime, group = Correct, color = Correct)) +
  #geom_point() +
  geom_line(size = 2) +
  # add x limit values
  scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
  labs(x = "Cue target similarity",
       y = "Model times") +
  theme_minimal()+
  theme(legend.position = "bottom")+
  scale_color_manual(values = c("TRUE" = "forestgreen", "FALSE" = "red"))+
  theme(text = element_text(size = textsize, family = "sans"))
Exp1SimRT


AggragatedModelTimes <- ddply(model_predictions, 
                              ~simBin + ListLength + Correct, 
                              summarize, MeanModelTime = mean(ModelTimes))
AggragatedModelTimes$simBin <- as.numeric(as.character(AggragatedModelTimes$simBin))

Exp1SimRTLL = ggplot(subset(AggragatedModelTimes, Correct == TRUE &
                              ListLength %in% c(1, 5, 10, 20, 32)),
                     aes(x = simBin, y = MeanModelTime, group = ListLength, color = as.factor(ListLength))) +
  #geom_point() +
  geom_line(size = 2) +
  labs(x = "Cue target similarity",
       y = "Model times",
       color = "Set size") +
  theme_minimal()+
  theme(legend.position = "bottom")+
  theme(text = element_text(size = textsize, family = "sans"))
Exp1SimRTLL


# Make figure that summarises the model predictions for experiment 1

Exp1ModelPredictions = Exp1LLRT + Exp1SimRT + Exp1SimRTLL +
  plot_annotation(title = "Model Predictions Experiment 1")
Exp1ModelPredictions

# Experiment 2 -----------------------------------------------------------------

# simulate some participant data that mimic the setup of experiment 2

# Set sizes between 0 and 32 in 1 increments
ContextSizes15 = rep(0:15, reps)
ContextSizes30 = rep(0:30, reps)

ContextBonus = 0.3


# prep dataframe
df2 = data.frame(ListLength = numeric(),
                 SimCueTarget = numeric(),
                 SimSumDistractors = numeric(),
                 Condition = character(),
                 Experiment = character(),
                 ContextSize = numeric())
# now for each Set sizes lets simulate the cue target similarity and the sum of the distractors similarity
for (i in 1:length(ContextSizes15)){
  # cue target similarity (sample from a uniform distribution between 0 and 1)
  sim_cue_target = runif(1, min = 0, max = 1)
  if (ContextSizes15[i] != 0){
    sim_cue_target = sim_cue_target + ContextBonus
  }
  # distractors similarity
  sim_sum_distractors = sum(runif(15-1, min = 0, max = 1)) + sim_cue_target + ContextBonus * ContextSizes15[i]
  # condition
  if (ContextSizes15[i] == 0){
    condition = "None"
  }else{
    condition = "Context"
  }
  # experiment
  experiment = "2"
  # context size
  context_size = ContextSizes15[i]
  # add to the dataframe
  df2 = rbind(df2, data.frame(ListLength = 15,
                              SimCueTarget = sim_cue_target,
                              SimSumDistractors = sim_sum_distractors,
                              Condition = condition,
                              Experiment = experiment,
                              ContextSize = context_size))
}

for (i in 1:length(ContextSizes30)){
  # cue target similarity (sample from a uniform distribution between 0 and 1)
  sim_cue_target = runif(1, min = 0, max = 1)
  if (ContextSizes30[i] != 0){
    sim_cue_target = sim_cue_target + ContextBonus
  }
  # distractors similarity
  sim_sum_distractors = sum(runif(30-1, min = 0, max = 1)) + sim_cue_target + ContextBonus * ContextSizes30[i]
  # condition
  if (ContextSizes30[i] == 0){
    condition = "None"
  }else{
    condition = "Context"
  }
  # experiment
  experiment = "2"
  # context size
  context_size = ContextSizes30[i]
  # add to the dataframe
  df2 = rbind(df2, data.frame(ListLength = 30,
                              SimCueTarget = sim_cue_target,
                              SimSumDistractors = sim_sum_distractors,
                              Condition = condition,
                              Experiment = experiment,
                              ContextSize = context_size))
}


# simulate the model
model_predictions2 = modelSimulation(df2, p_omission, p_intrusion, association_strength, reps, UseSimilarity)

# plot predictions Experiment 2 ------------------------------------------------

# create a coding for correct and incorrect correct is everything where the model response type is correct
model_predictions2$Correct = ifelse(model_predictions2$ModelResponseType == "Correct", TRUE, FALSE)
model_predictions2$Correct_01 = ifelse(model_predictions2$ModelResponseType == "Correct", 1, 0)
# male Set sizea factor
model_predictions2$ListLength = as.factor(model_predictions2$ListLength)

# calculate mean RT for the none context condition for each set size
meanRTsNone = ddply(subset(model_predictions2, Condition == "None"), ~ ListLength, summarise, MeanRT = mean(ModelTimes))

meanRTs2 = ddply(subset(model_predictions2, Condition != "None"), ~ ListLength + ContextSize + Correct + Condition, summarise, MeanRT = mean(ModelTimes))

# Plot ModelTimes based on ContextSize
Exp2ContextRT <- ggplot(
  subset(meanRTs2, Condition == "Context" & Correct == TRUE), 
  aes(x = ContextSize, y = MeanRT, group = ListLength, color = ListLength)) +
  # Add a smooth line
  geom_line(size = 2) +
  # Add labels
  labs(x = "Context size", y = "Model times", color = "Set size") +
  # Use a minimal theme
  theme_minimal() +
  # Add horizontal lines for meanRTsNone as dotted lines
  geom_hline(data = meanRTsNone, aes(yintercept = MeanRT, color = ListLength), linetype = "dotted", size = 2) +
  scale_color_manual(values = c("#E69F00", "#D55E00", "#56B4E9", "#009E73")) +  # Okabe-Ito palette
  scale_fill_manual(values = c("#E69F00",  "#D55E00", "#56B4E9", "#009E73")) +   # Okabe-Ito palette
  # Adjust legend position
  theme(legend.position = "bottom")+
  theme(text = element_text(size = textsize, family = "sans"))
Exp2ContextRT

# plot the same for accuracy
meanAccNone = ddply(subset(model_predictions2, Condition == "None"), ~ ListLength, summarise, MeanAcc = mean(Correct_01))

meanAcc2 = ddply(subset(model_predictions2, Condition == "Context"), ~ ListLength + ContextSize, summarise, MeanAcc = mean(Correct_01))

Exp2ContextAcc <- ggplot(
  meanAcc2, 
  aes(x = ContextSize, y = MeanAcc, group = ListLength, color = ListLength)) +
  # Add a smooth line
  geom_line(size = 2) +
  # Add labels
  labs(x = "Context size", y = "Accuracy", color = "Set size") +
  # Use a minimal theme
  theme_minimal() +
  geom_hline(data = meanAccNone, aes(yintercept = MeanAcc, color = ListLength), linetype = "dotted", size = 2) +
  scale_color_manual(values = c("#E69F00", "#D55E00", "#56B4E9", "#009E73")) +  # Okabe-Ito palette
  scale_fill_manual(values = c("#E69F00",  "#D55E00", "#56B4E9", "#009E73")) +   # Okabe-Ito palette
  # Adjust legend position
  theme(legend.position = "bottom")+
  theme(text = element_text(size = textsize, family = "sans"))
Exp2ContextAcc

# plot the mean RT for correct and incorrect responses
meanRTAcc2 = ddply(model_predictions2, ~ Correct, summarise, MeanRT = mean(ModelTimes))

# Plot mean RT for correct and incorrect responses as bar plot
Exp2RTAcc = ggplot(meanRTAcc2, aes(x = Correct, y = MeanRT, fill = Correct)) +
  geom_bar(stat = "identity", width = 0.5) +
  labs(x = "Correct", y = "Model times") +
  theme_minimal() +
  scale_fill_manual(values = c("TRUE" = "forestgreen", "FALSE" = "red"))+
  theme(legend.position = "none")+
  theme(text = element_text(size = textsize, family = "sans"))
Exp2RTAcc

# Make figure that summarises the model predictions for experiment 2

Exp2ModelPredictions = Exp2ContextRT + Exp2ContextAcc + Exp2RTAcc +
  plot_layout(guides = "collect")+
  plot_annotation(title = "Model Predictions Experiment 2")
Exp2ModelPredictions

# Experiment 3 -----------------------------------------------------------------

# simulate some participant data that mimic the setup of experiment 2

# Set sizes between 0 and 32 in 1 increments
ContextSizes20 = rep(1:20, reps)
ListLengthExp3 = 20


# prep dataframe
df3 = data.frame(ListLength = numeric(),
                 SimCueTarget = numeric(),
                 SimSumDistractors = numeric(),
                 Condition = character(),
                 Experiment = character(),
                 ContextSize = numeric())
# now for each Set sizes lets simulate the cue target similarity and the sum of the distractors similarity
for (i in 1:length(ContextSizes20)){
  # cue target similarity (sample from a uniform distribution between 0 and 1)
  sim_cue_target = runif(1, min = 0, max = 1)
  if (ContextSizes20[i] != 0){
    sim_cue_target = sim_cue_target + ContextBonus
  }
  # distractors similarity
  sim_sum_distractors = sum(c(runif(ListLengthExp3-ContextSizes20[i], min = 0, max = 0.3), runif(ContextSizes20[i]-1, min = 0.7, max = 1), min = 0, max = 1)) + sim_cue_target 
  # condition
  condition = "Context"
  # experiment
  experiment = "3"
  # context size
  context_size = ContextSizes20[i]
  # add to the dataframe
  df3 = rbind(df3, data.frame(ListLength = ListLengthExp3,
                              SimCueTarget = sim_cue_target,
                              SimSumDistractors = sim_sum_distractors,
                              Condition = condition,
                              Experiment = experiment,
                              ContextSize = context_size))
}

# Add the control conditions
# high pair similarity
for (i in 1:reps){
  # cue target similarity (sample from a uniform distribution between 0 and 1)
  sim_cue_target = runif(1, min = 0.7, max = 1)
  # distractors similarity
  sim_sum_distractors = sum(runif(ListLengthExp3-1, min = 0, max = 0.3)) + sim_cue_target
  # condition
  condition = "HighPairSimilarity"
  # experiment
  experiment = "3"
  # context size
  context_size = 0
  # add to the dataframe
  df3 = rbind(df3, data.frame(ListLength = ListLengthExp3,
                              SimCueTarget = sim_cue_target,
                              SimSumDistractors = sim_sum_distractors,
                              Condition = condition,
                              Experiment = experiment,
                              ContextSize = context_size))
}

for (i in 1:reps){
  # cue target similarity (sample from a uniform distribution between 0 and 1)
  sim_cue_target = runif(1, min = 0, max = 0.3)
  # distractors similarity
  sim_sum_distractors = sum(runif(ListLengthExp3-1, min = 0, max = 0.3)) + sim_cue_target
  # condition
  condition = "LowSimilarityList"
  # experiment
  experiment = "3"
  # context size
  context_size = 0
  # add to the dataframe
  df3 = rbind(df3, data.frame(ListLength = ListLengthExp3,
                              SimCueTarget = sim_cue_target,
                              SimSumDistractors = sim_sum_distractors,
                              Condition = condition,
                              Experiment = experiment,
                              ContextSize = context_size))
}

# simulate the model
model_predictions3 = modelSimulation(df3, p_omission, p_intrusion, association_strength, reps, UseSimilarity)

# plot predictions Experiment 3 ------------------------------------------------

# create a coding for correct and incorrect correct is everything where the model response type is correct
model_predictions3$Correct = ifelse(model_predictions3$ModelResponseType == "Correct", TRUE, FALSE)
model_predictions3$Correct_01 = ifelse(model_predictions3$ModelResponseType == "Correct", 1, 0)
# male Set sizea factor
model_predictions3$ListLength = as.factor(model_predictions3$ListLength)

# calculate mean RT for the none context condition for each set size
meanRTsNone = ddply(subset(model_predictions3, Condition != "Context"), ~ Condition, summarise, MeanRT = mean(ModelTimes))

meanRTs3 = ddply(subset(model_predictions3, Condition == "Context"), ~ Condition + Correct + ContextSize, summarise, MeanRT = mean(ModelTimes))

# Plot ModelTimes based on ContextSize
Exp3ContextRT <- ggplot(
  subset(meanRTs3, Condition == "Context" & Correct == TRUE), 
  aes(x = ContextSize, y = MeanRT)) +
  # Add a smooth line
  geom_line(size = 2, color = RTColor) +
  # Add labels
  labs(x = "Context size", y = "Model times") +
  # Use a minimal theme
  theme_minimal() +
  # Add horizontal lines for meanRTsNone as dotted lines
  geom_hline(data = meanRTsNone, aes(yintercept = MeanRT, linetype = Condition), size = 2, color = "black") +
  # Adjust legend position
  theme(legend.position = "bottom")+
  theme(text = element_text(size = textsize, family = "sans"))
Exp3ContextRT

# plot the same for accuracy

# calculate mean RT for the none context condition for each set size
meanAccNone = ddply(subset(model_predictions3, Condition != "Context"), ~ Condition, summarise, MeanAcc = mean(Correct_01))

# calculate mean RT for the none context condition for each set size
meanAcc3 = ddply(subset(model_predictions3, Condition == "Context"), ~ Condition + ContextSize, summarise, MeanAcc = mean(Correct_01))

# Plot ModelTimes based on ContextSize
Exp3ContextAcc <- ggplot(
  subset(meanAcc3, Condition == "Context"), 
  aes(x = ContextSize, y = MeanAcc)) +
  # Add a smooth line
  geom_line(size = 2, color = AccColor) +
  # Add labels
  labs(x = "Context size", y = "Model times") +
  # Use a minimal theme
  theme_minimal() +
  # Add horizontal lines for meanRTsNone as dotted lines
  geom_hline(data = meanAccNone, aes(yintercept = MeanAcc, linetype = Condition), color = "black", size = 2) +
  # Adjust legend position
  theme(legend.position = "bottom")+
  theme(text = element_text(size = textsize, family = "sans"))
Exp3ContextAcc

# plot the mean RT for correct and incorrect responses
meanRTAcc3 = ddply(model_predictions3, ~ Correct, summarise, MeanRT = mean(ModelTimes))

# Plot mean RT for correct and incorrect responses as bar plot
Exp3RTAcc = ggplot(meanRTAcc3, aes(x = Correct, y = MeanRT, fill = Correct)) +
  geom_bar(stat = "identity", width = 0.5) +
  labs(x = "Correct", y = "Model times") +
  theme_minimal() +
  scale_fill_manual(values = c("TRUE" = "forestgreen", "FALSE" = "red"))+
  theme(legend.position = "none")+
  theme(text = element_text(size = textsize, family = "sans"))
Exp3RTAcc


# make one summary figure with all of them -------------------------------------


# remove legend from Exp2RTAcc
Exp2RTAcc = Exp2RTAcc + theme(legend.position = "none")
Exp3RTAcc = Exp3RTAcc + theme(legend.position = "none")

# Combine plots for each experiment with shared legends and axes
Exp1Plots <- Exp1LLRT + Exp1SimRT + Exp1SimRTLL +
  plot_layout(guides = "collect")  & 
  theme(legend.position = "bottom")
Exp2Plots <- Exp2ContextRT + Exp2ContextAcc + Exp2RTAcc +
  plot_layout(guides = "collect")  & 
  theme(legend.position = "bottom")
Exp3Plots <- Exp3ContextRT + Exp3ContextAcc + Exp3RTAcc +
  plot_layout(guides = "collect")  & 
  theme(legend.position = "bottom")


# Arrange everything
ModelPredictions <- (
  Exp1Plots /
    Exp2Plots /
    Exp3Plots
) + plot_layout(heights = c(1, 1, 1)) +
  plot_annotation(tag_levels = "A")

ModelPredictions

current_directory <- dirname(rstudioapi::getSourceEditorContext()$path)
savingPlace = "/Figures/ModelPredictionFigure.pdf"
ggsave(paste(current_directory, savingPlace, sep = ""), ModelPredictions,  
       width = 20, 
       height = 10)

