
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


modelSimulation = function(Participant_data, p_omission, p_intrusion, association_strength, beta, repetitions = 1){
  sim_cue_target = Participant_data$SimCueTarget
  sim_distractors = Participant_data$SimsCueDistractors
  p_correct = c()
  for (i in 1:length(sim_cue_target)){
    sim_distractors_vec = as.numeric(strsplit(gsub("\\[|\\]", "", sim_distractors[i]), ",")[[1]])
    if (length(sim_distractors_vec) == 0){
      sim_distractors_vec = c(0)
    }
    p_correct = c(p_correct, min((sim_cue_target[i]**beta + association_strength)/
                                   (sim_cue_target[i]**beta + association_strength + sum(sim_distractors_vec**beta)), 1))
  }
  p_stopping = pmax((p_correct + (1-p_correct)*(min((p_omission + p_intrusion), 1))), 0.000000000001)
  p_incorrect = 1-p_correct
  ListLength = Participant_data$ListLength
  Condition = Participant_data$Condition
  Experiment = Participant_data$Experiment
  ContextSize = Participant_data$ContextSize
  CueTargetSimilarity = Participant_data$SimCueTarget
  
  
  # calculate the RTs and responsetzpes
  
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
  all_response_types_matrix <- do.call(rbind, all_response_types)
  all_model_rt_matrix <- do.call(rbind, all_times)
  
  # Calculate the average time for each entry across repetitions
  average_model_rt <- colMeans(all_model_rt_matrix)
  
  # Randomly choose one response type for each entry
  random_response_types <- apply(all_response_types_matrix, 2, function(x) sample(x, 1))
  
  
  
  df = data.frame(ListLength = ListLength,
                  Condition = Condition,
                  Experiment = Experiment,
                  ContextSize = ContextSize,
                  ModelTimes = average_model_rt,
                  ModelResponseType = random_response_types,
                  CueTargetSimilarity = CueTargetSimilarity)
  return(df)
}

# general varibale -----------------------------------------------------------
reps = 10

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
    sim_distractors = c()
  }else{
    sim_sum_distractors = sum(runif(ListLength[i]-1, min = 0, max = 1)) + sim_cue_target
    sim_distractors = runif(ListLength[i]-1, min = 0, max = 1)
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
                            ContextSize = context_size,
                            SimsCueDistractors = toString(sim_distractors)))
}

# Define general parameters ---------------------------------------------------

p_omissions = c(seq(0, 0.2, by = 0.02), seq(0.2, 0.5, by = 0.1))
p_intrusions = c(seq(0, 0.2, by = 0.02), seq(0.2, 0.5, by = 0.1))
association_strengths = seq(0, 6, by = 0.5) 
betas = seq(0, 1.5, by = 0.2)
UseSimilarity = TRUE

# get model predictions for robustness check -----------------------------------

# create en empty dataframe to store the model predictions
model_predictions = data.frame()

for (p_omission in p_omissions){
  print("p_omission")
  print(p_omission)
  for (p_intrusion in p_intrusions){
    for (association_strength in association_strengths){
      for (beta in betas){
        # simulate the model
        temp = modelSimulation(df, p_omission, p_intrusion, association_strength, beta, reps)
        
        # add parameter values to the dataframe
        temp$p_omission = p_omission
        temp$p_intrusion = p_intrusion
        temp$association_strength = association_strength
        temp$beta = beta
        
        # add to the model predictions
        model_predictions = rbind(model_predictions, temp)
        
      }
      
    }
  }
}

# save model predictions
#set diretory to current working directory
current_directory <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(current_directory)

#write.csv(model_predictions, "ModelPredictions_withbeta.csv")
# read the file in 
model_predictions = read.csv("ModelPredictions_withbeta.csv")



# plot predictions Experiment 1 ------------------------------------------------

textsize = 15
font = "sans"
default_color <- "#007BB8"
AccColor <- "#74c476"
RTColor <- "#56B4E9"
RegColor = "black"
RibbonAlpha = 0.2
RT_min_value = 0
RT_max_value = 4
trans = 1


# Define color mapping and readable names
effect_labels <- c(
  ListLength = "Set size",
  ContextSize = "Context size",
  ConditionLowSimilarityList = "Low similarity control condition",
  ConditionHighPairSimilarity = "High pair similarity control condition",
  CorrectTRUE = "Correct response",
  CueTargetSimilarity = "Cue-target similarity",
  `ListLength:CueTargetSimilarity` = "Set size × cue-target similarity"
)

effect_colors <- c(
  ListLength = "#66c2a5",
  ContextSize = "#fc8d62",
  ConditionLowSimilarityList = "#8da0cb",
  ConditionHighPairSimilarity = "#e78ac3",
  CorrectTRUE = "#a6d854",
  CueTargetSimilarity = "#ffd92f",
  `ListLength:CueTargetSimilarity` = "#e5c494"
)

# Shared color scale
shared_color_scale <- scale_color_manual(
  values = effect_colors,
  labels = effect_labels,
  name = "Effect Type"
)


# calculate linear regression mit ListLength and cue target similarity as main effects seperately for each value if p_omission, p_intrusion and association_strength
library(ggplot2)
library(dplyr)
library(purrr)
library(tidyr)
library(broom)

# Define variables to loop over
variables <- c("p_omission", "p_intrusion", "association_strength", "beta")
labels <- c("p(omission)", "p(intrusion)", "Association strength", "Beta")

# crete a coding for correct and incorrect correct is everything where the model response type is correct
model_predictions$Correct = ifelse(model_predictions$ModelResponseType == "Correct", TRUE, FALSE)
model_predictions$Correct_01 = ifelse(model_predictions$ModelResponseType == "Correct", 1, 0)

# Function to run regression and extract results
run_regression <- function(var) {
  temp = model_predictions %>%
    # only analyse beta values above 1.5
    filter(beta < 1.5) %>%
    mutate(ListLength = scale(ListLength, center = TRUE, scale = TRUE),
           CueTargetSimilarity = scale(CueTargetSimilarity, center = TRUE, scale = TRUE)) %>%
    group_by(across(all_of(var))) %>%
    do(tidy(lm(ModelTimes ~ ListLength * CueTargetSimilarity + Correct, data = .))) %>%
    filter(term %in% c("ListLength", "CueTargetSimilarity", "ListLength:CueTargetSimilarity", "CorrectTRUE")) %>%
    mutate(lower = estimate - 1.96 * std.error,
           upper = estimate + 1.96 * std.error,
           variable = var)  # Store x values for correct plotting
  temp[["x_value"]] = temp[[var]]
  return(temp)
}

# Run regressions for all variables
results_all <- map_dfr(variables, run_regression)

# Match variable names to readable labels
results_all <- results_all %>%
  mutate(variable = factor(variable, levels = variables, labels = labels))

# Create the combined plot
TimesExp1 = ggplot(results_all, aes(x = x_value, y = estimate, color = term)) +
  geom_point(size = 1) +
  geom_line() +  
  #geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.02, alpha = 0.5) + 
  labs(title = "Experiment 1 model times",
       x = "Parameter value", 
       y = "Effect size") +
  theme_minimal() +
  geom_hline(yintercept = 0, linetype = "dotted") +
  coord_cartesian(ylim = c(-0.5, 0.5)) +
  facet_wrap(~ variable, scales = "free_x", nrow = 1) +  # Row layout
  shared_color_scale +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))  # Place legend at the bottom
TimesExp1

# same for accuracy

run_logistic_regression <- function(var) {
  temp = model_predictions %>%
    filter(beta < 1.5) %>%
    mutate(ListLength = scale(ListLength, center = TRUE, scale = TRUE),
           CueTargetSimilarity = scale(CueTargetSimilarity, center = TRUE, scale = TRUE)) %>%
    group_by(across(all_of(var))) %>%
    do(tidy(glm(Correct ~ ListLength * CueTargetSimilarity, 
                data = ., 
                family = binomial(link = "logit")))) %>%
    filter(term %in% c("ListLength", "CueTargetSimilarity", "ListLength:CueTargetSimilarity")) %>%
    mutate(lower = estimate - 1.96 * std.error,
           upper = estimate + 1.96 * std.error,
           variable = var)  # Store x values for correct plotting
  temp[["x_value"]] = temp[[var]]
  return(temp)
}

# Run regressions for all variables
results_all_AC <- map_dfr(variables, run_logistic_regression)

# Match variable names to readable labels
results_all_AC <- results_all_AC %>%
  mutate(variable = factor(variable, levels = variables, labels = labels))

# Create the combined plot
AccExp1 = ggplot(results_all_AC, aes(x = x_value, y = estimate, color = term)) +
  geom_point(size = 1) +
  geom_line() +  
  #geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.02, alpha = 0.5) + 
  labs(title = "Experiment 1 acuracy",
       x = "Parameter value", 
       y = "Effect size") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  theme_minimal() +
  coord_cartesian(ylim = c(-0.8, 0.25)) +
  facet_wrap(~ variable, scales = "free_x", nrow = 1) +  # Row layout
  shared_color_scale +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))  # Place legend at the bottom
AccExp1

# Some other plots



# calculate mean RT for each set size
meanRTs = ddply(model_predictions, ~ ListLength + Correct + p_intrusion, summarise, MeanRT = mean(ModelTimes))
# 
Exp1LLRT_i = ggplot(meanRTs, aes(x = ListLength, y = MeanRT, color = p_intrusion, group = p_intrusion)) +
  #geom_point() +
  geom_line(size = 1, alpha = trans) +
  # add mean RT
  # plot mean RTs
  coord_cartesian(ylim = c(0, 5)) +
  # add mean RT
  # plot mean RTs
  labs(x = "Set size",
       y = "Model times",
       color = "p(intrusion)")+
  theme_minimal()+
  guides(color = guide_colorbar(ticks = FALSE, 
                                barwidth = 8, 
                                barheight = 1, 
                                label.position = "bottom")) +
  theme(legend.position = "bottom")+
  scale_color_viridis_c(option = "rocket") +
  theme(text = element_text(size = textsize, family = "sans"))+
  facet_wrap(~Correct)
Exp1LLRT_i



# calculate mean RT for each set size
meanRTs = ddply(model_predictions, ~ ListLength + Correct + p_omission, summarise, MeanRT = mean(ModelTimes))
# 
Exp1LLRT_o = ggplot(meanRTs, aes(x = ListLength, y = MeanRT, color = p_omission, group = p_omission)) +
  #geom_point() +
  geom_line(size = 1, alpha = trans) +
  coord_cartesian(ylim = c(0, 5)) +
  # add mean RT
  # plot mean RTs
  labs(x = "Set size",
       y = "Model times",
       color = "p(omission)") +
  guides(color = guide_colorbar(ticks = FALSE, 
                                barwidth = 8, 
                                barheight = 1, 
                                label.position = "bottom")) +
  theme_minimal()+
  scale_color_viridis_c(option = "viridis")  +
  theme(legend.position = "bottom")+
  theme(text = element_text(size = textsize, family = "sans"))+
  facet_wrap(~Correct)
Exp1LLRT_o

# calculate mean RT for each set size
meanRTs = ddply(model_predictions, ~ ListLength + Correct + association_strength, summarise, MeanRT = mean(ModelTimes))
# 
Exp1LLRT_as = ggplot(meanRTs, aes(x = ListLength, y = MeanRT, color = association_strength, group = association_strength)) +
  #geom_point() +
  geom_line(size = 1, alpha = trans) +
  # add mean RT
  # plot mean RTs
  coord_cartesian(ylim = c(0, 5)) +
  labs(x = "Set size",
       y = "Model times",
       color = "AS") +
  guides(color = guide_colorbar(ticks = FALSE, 
                                barwidth = 8, 
                                barheight = 1, 
                                label.position = "bottom")) +
  theme_minimal()+
  theme(legend.position = "bottom")+
  scale_color_viridis_c(option = "cividis")  +
  theme(text = element_text(size = textsize, family = "sans"))+
  facet_wrap(~Correct)
Exp1LLRT_as




# make sim Bins

breaks = seq(0, 1, by = 0.1)
model_predictions$simBin <- cut(model_predictions$CueTargetSimilarity, 
                                breaks = breaks, 
                                labels = breaks[-1]-0.05,
                                include.lowest = TRUE, 
                                right = FALSE)

# aggregate data across sim bins and set size
AggragatedModelTimes <- ddply(model_predictions, 
                              ~simBin + Correct + p_intrusion, 
                              summarize, MeanModelTime = mean(ModelTimes))
AggragatedModelTimes$simBin <- as.numeric(as.character(AggragatedModelTimes$simBin))

# Plot ModelTimes based on CueTargetSimilarity
Exp1SimRT = ggplot(AggragatedModelTimes, aes(x = simBin, y = MeanModelTime, color = p_intrusion, group = p_intrusion)) +
  #geom_point() +
  geom_line(size = 1, alpha = 0.5) +
  coord_cartesian(ylim = c(0, 5)) +
  # add x limit values
  scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
  labs(x = "Cue target similarity",
       y = "Model times",
       color = "p(intrusion)") +
  theme_minimal()+
  guides(color = guide_colorbar(ticks = FALSE, 
                                barwidth = 8, 
                                barheight = 1, 
                                label.position = "bottom")) +
  theme(legend.position = "bottom")+
  scale_color_viridis_c(option = "rocket") +
  theme(text = element_text(size = textsize, family = "sans"))+
  facet_wrap(~Correct)
Exp1SimRT

# aggregate data across sim bins and set size
AggragatedModelTimes <- ddply(model_predictions, 
                              ~simBin + Correct + p_omission, 
                              summarize, MeanModelTime = mean(ModelTimes))
AggragatedModelTimes$simBin <- as.numeric(as.character(AggragatedModelTimes$simBin))

# Plot ModelTimes based on CueTargetSimilarity
Exp1SimRT = ggplot(AggragatedModelTimes, aes(x = simBin, y = MeanModelTime, color = as.factor(p_omission))) +
  #geom_point() +
  geom_line(size = 1, alpha = 0.5) +
  coord_cartesian(ylim = c(0, 5)) +
  # add x limit values
  scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
  labs(x = "Cue target similarity",
       y = "Model times",
       color = "p(omission)") +
  theme_minimal()+
  theme(legend.position = "bottom")+
  theme(text = element_text(size = textsize, family = "sans"))+
  facet_wrap(~Correct)
Exp1SimRT

# aggregate data across sim bins and set size
AggragatedModelTimes <- ddply(model_predictions, 
                              ~simBin + Correct + association_strength, 
                              summarize, MeanModelTime = mean(ModelTimes))
AggragatedModelTimes$simBin <- as.numeric(as.character(AggragatedModelTimes$simBin))

# Plot ModelTimes based on CueTargetSimilarity
Exp1SimRT = ggplot(AggragatedModelTimes, aes(x = simBin, y = MeanModelTime, color = as.factor(association_strength))) +
  #geom_point() +
  geom_line(size = 1, alpha = 0.5) +
  # add x limit values
  scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
  coord_cartesian(ylim = c(2.5, 6)) +
  labs(x = "Cue target similarity",
       y = "Model times",
       color = "AS") +
  theme_minimal()+
  theme(legend.position = "bottom")+
  theme(text = element_text(size = textsize, family = "sans"))+
  facet_wrap(~Correct)
Exp1SimRT


# Experiment 2 -----------------------------------------------------------------

# simulate some participant data that mimic the setup of experiment 2

# Set sizes between 0 and 32 in 1 increments
ContextSizes15 = rep(0:15, reps)
ContextSizes30 = rep(0:30, reps)

ContextBonuses = seq(0, 0.5, by = 0.1)

# prep dataframe
df2 = data.frame(ListLength = numeric(),
                 SimCueTarget = numeric(),
                 SimSumDistractors = numeric(),
                 Condition = character(),
                 Experiment = character(),
                 ContextSize = numeric(),
                 ContextBonus = numeric())
# now for each Set sizes lets simulate the cue target similarity and the sum of the distractors similarity
for (ContextBonus in ContextBonuses){
  for (i in 1:length(ContextSizes15)){
    # cue target similarity (sample from a uniform distribution between 0 and 1)
    sim_cue_target = runif(1, min = 0, max = 1)
    if (ContextSizes15[i] != 0){
      sim_cue_target = sim_cue_target + ContextBonus
    }
    # distractors similarity
    if (ContextSizes15[i] < 1 ){
      sim_distractors = c(runif(15*2-1, min = 0, max = 1))
    }else{
      sim_distractors = c(runif(ContextSizes15[i]*2-1, min = 0, max = 1) + ContextBonus, 
                          runif(15*2-ContextSizes15[i]*2, min = 0, max = 1))
    }
    
    sim_sum_distractors = sum(sim_distractors) + sim_cue_target
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
                                ContextSize = context_size,
                                ContextBonus = ContextBonus,
                                SimsCueDistractors = toString(sim_distractors)))
  }
  
  for (i in 1:length(ContextSizes30)){
    # cue target similarity (sample from a uniform distribution between 0 and 1)
    sim_cue_target = runif(1, min = 0, max = 1)
    if (ContextSizes30[i] != 0){
      sim_cue_target = sim_cue_target + ContextBonus
    }
    # distractors similarity
    if (ContextSizes30[i] < 1 ){
      sim_distractors = c(runif(30*2-1, min = 0, max = 1))
    }else{
      sim_distractors = c(runif(ContextSizes30[i]*2-1, min = 0, max = 1) + ContextBonus, 
                          runif(30*2-ContextSizes30[i]*2, min = 0, max = 1))
    }
    sim_sum_distractors = sum(sim_distractors) + sim_cue_target
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
                                ContextSize = context_size,
                                ContextBonus = ContextBonus,
                                SimsCueDistractors = toString(sim_distractors)))
  }
}


# simulate the model
model_predictions2 = data.frame()

for (p_omission in p_omissions){
  print("p_omission")
  print(p_omission)
  for (p_intrusion in p_intrusions){
    for (association_strength in association_strengths){
      for (beta in betas){
        # simulate the model
        temp = modelSimulation(df2, p_omission, p_intrusion, association_strength, beta, 1)
        
        # add parameter values to the dataframe
        temp$p_omission = p_omission
        temp$p_intrusion = p_intrusion
        temp$association_strength = association_strength
        temp$ContextBonus = df2$ContextBonus
        temp$beta = beta
        
        # add to the model predictions
        model_predictions2 = rbind(model_predictions2, temp)
      }
    }
  }
}

# save model predictions
write.csv(model_predictions2, "ModelPredictions2_withbeta.csv")
# read the file in
model_predictions2 = read.csv("ModelPredictions2_withbeta.csv") 

# plot predictions Experiment 2 ------------------------------------------------

# create a coding for correct and incorrect correct is everything where the model response type is correct
model_predictions2$Correct = ifelse(model_predictions2$ModelResponseType == "Correct", TRUE, FALSE)
model_predictions2$Correct_01 = ifelse(model_predictions2$ModelResponseType == "Correct", 1, 0)

variables2 <- c("p_omission", "p_intrusion", "association_strength", "beta", "ContextBonus")
labels2 <- c("p(omission)", "p(intrusion)", "Association strength", "Beta", "Context bonus")

model_predictions2$ListLength = as.numeric(model_predictions2$ListLength)
# Function to run regression and extract results
run_regression2 <- function(var) {
  temp = model_predictions2 %>%
    mutate(ListLength = scale(ListLength, center = TRUE, scale = TRUE),
           ContextSize = scale(ContextSize, center = TRUE, scale = TRUE)) %>%
    group_by(across(all_of(var))) %>%
    do(tidy(lm(ModelTimes ~ ListLength + ContextSize + Correct, data = .))) %>%
    filter(term %in% c("ListLength", "ContextSize", "CorrectTRUE")) %>%
    mutate(lower = estimate - 1.96 * std.error,
           upper = estimate + 1.96 * std.error,
           variable = var)  # Store x values for correct plotting
  temp[["x_value"]] = temp[[var]]
  return(temp)
}

# Run regressions for all variables
results_all2 <- map_dfr(variables2, run_regression2)

# Match variable names to readable labels
results_all2 <- results_all2 %>%
  mutate(variable = factor(variable, levels = variables2, labels = labels2))

# Create the combined plot
TimesExp2 = ggplot(results_all2, aes(x = x_value, y = estimate, color = term)) +
  geom_point(size = 1) +
  geom_line() +  
  #geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.02, alpha = 0.5) + 
  labs(title = "Experiment 2 model times",
       x = "Parameter value", 
       y = "Effect size") +
  theme_minimal() +
  coord_cartesian(ylim = c(-0.2, 0.8)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  facet_wrap(~ variable, scales = "free_x", nrow = 1) +  # Row layout
  shared_color_scale +
  theme(legend.position = "top", 
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))  # Place legend at the bottom
TimesExp2

# same for accuracy

run_logistic_regression2 <- function(var) {
  temp = model_predictions2 %>%
    mutate(ListLength = scale(ListLength, center = TRUE, scale = TRUE),
           ContextSize = scale(ContextSize, center = TRUE, scale = TRUE)) %>%
    group_by(across(all_of(var))) %>%
    do(tidy(glm(Correct ~ ListLength + ContextSize, 
                data = ., 
                family = binomial(link = "logit")))) %>%
    filter(term %in% c("ListLength", "ContextSize")) %>%
    mutate(lower = estimate - 1.96 * std.error,
           upper = estimate + 1.96 * std.error,
           variable = var)  # Store x values for correct plotting
  temp[["x_value"]] = temp[[var]]
  return(temp)
  
}

# Run regressions for all variables
results_all_AC2 <- map_dfr(variables2, run_logistic_regression2)

# Match variable names to readable labels
results_all_AC2 <- results_all_AC2 %>%
  mutate(variable = factor(variable, levels = variables2, labels = labels2))

# Create the combined plot
AccExp2 = ggplot(results_all_AC2, aes(x = x_value, y = estimate, color = term)) +
  geom_point(size = 1) +
  geom_line() +  
  #geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.02, alpha = 0.5) + 
  labs(title = "Experiment 2 accuracy",
       x = "Parameter value", 
       y = "Effect size") +
  theme_minimal() +
  # add horizontal line at 0
  geom_hline(yintercept = 0, linetype = "dotted") +
  facet_wrap(~ variable, scales = "free_x", nrow = 1) +  # Separate by parameter
  shared_color_scale +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))  # Place legend at the bottom
AccExp2

# -----------------------------------------------------------------------------

# male Set sizea factor
model_predictions2$ListLength = as.factor(model_predictions2$ListLength)

# calculate mean RT for the none context condition for each set size
meanRTsNone = ddply(subset(model_predictions2, Condition == "None"), ~ ListLength + p_omission, summarise, MeanRT = mean(ModelTimes))

meanRTs2 = ddply(subset(model_predictions2, Condition != "None"), ~ ListLength + ContextSize + Correct + Condition + p_omission, summarise, MeanRT = mean(ModelTimes))

# Plot ModelTimes based on ContextSize
Exp2ContextRT <- ggplot(
  subset(meanRTs2, Condition == "Context" & Correct == TRUE), 
  aes(x = ContextSize, y = MeanRT, color = as.factor(p_omission))) +
  # Add a smooth line
  geom_line(size = 1) +
  # Add labels
  labs(x = "Context size", y = "Model times", color = "p(omission)") +
  # Use a minimal theme
  theme_minimal() +
  # Add horizontal lines for meanRTsNone as dotted lines
  geom_hline(data = meanRTsNone, aes(yintercept = MeanRT, color = as.factor(p_omission)), linetype = "dotted", size = 1) +
  # Adjust legend position
  theme(legend.position = "bottom")+
  theme(text = element_text(size = textsize, family = "sans"))+
  facet_wrap(~ListLength)
Exp2ContextRT

# plot the same for accuracy
meanAccNone = ddply(subset(model_predictions2, Condition == "None"), ~ ListLength + p_omission, summarise, MeanAcc = mean(Correct_01))

meanAcc2 = ddply(subset(model_predictions2, Condition == "Context"), ~ ListLength + ContextSize + p_omission, summarise, MeanAcc = mean(Correct_01))

Exp2ContextAcc <- ggplot(
  meanAcc2, 
  aes(x = ContextSize, y = MeanAcc, color = as.factor(p_omission))) +
  # Add a smooth line
  geom_line(size = 2) +
  # Add labels
  labs(x = "Context size", y = "Accuracy", color = "p(omission)") +
  # Use a minimal theme
  theme_minimal() +
  geom_hline(data = meanAccNone, aes(yintercept = MeanAcc, color = as.factor(p_omission)), linetype = "dotted", size = 2) +
  # Adjust legend position
  theme(legend.position = "bottom")+
  theme(text = element_text(size = textsize, family = "sans"))+
  facet_wrap(~ListLength)
Exp2ContextAcc


meanRTsNone = ddply(subset(model_predictions2, Condition == "None"), ~ ListLength + p_intrusion, summarise, MeanRT = mean(ModelTimes))

meanRTs2 = ddply(subset(model_predictions2, Condition != "None"), ~ ListLength + ContextSize + Correct + Condition + p_intrusion, summarise, MeanRT = mean(ModelTimes))

# Plot ModelTimes based on ContextSize
Exp2ContextRT <- ggplot(
  subset(meanRTs2, Condition == "Context" & Correct == TRUE), 
  aes(x = ContextSize, y = MeanRT, color = as.factor(p_intrusion))) +
  # Add a smooth line
  geom_line(size = 1) +
  # Add labels
  labs(x = "Context size", y = "Model times", color = "p(intrusion)") +
  # Use a minimal theme
  theme_minimal() +
  # Add horizontal lines for meanRTsNone as dotted lines
  geom_hline(data = meanRTsNone, aes(yintercept = MeanRT, color = as.factor(p_intrusion)), linetype = "dotted", size = 1) +
  # Adjust legend position
  theme(legend.position = "bottom")+
  theme(text = element_text(size = textsize, family = "sans"))+
  facet_wrap(~ListLength)
Exp2ContextRT

# plot the same for accuracy
meanAccNone = ddply(subset(model_predictions2, Condition == "None"), ~ ListLength + p_intrusion, summarise, MeanAcc = mean(Correct_01))

meanAcc2 = ddply(subset(model_predictions2, Condition == "Context"), ~ ListLength + ContextSize + p_intrusion, summarise, MeanAcc = mean(Correct_01))

Exp2ContextAcc <- ggplot(
  meanAcc2, 
  aes(x = ContextSize, y = MeanAcc, color = as.factor(p_intrusion))) +
  # Add a smooth line
  geom_line(size = 1) +
  # Add labels
  labs(x = "Context size", y = "Accuracy", color = "p(intrusion)") +
  # Use a minimal theme
  theme_minimal() +
  geom_hline(data = meanAccNone, aes(yintercept = MeanAcc, color = as.factor(p_intrusion)), linetype = "dotted", size = 2) +
  # Adjust legend position
  theme(legend.position = "bottom")+
  theme(text = element_text(size = textsize, family = "sans"))+
  facet_wrap(~ListLength)
Exp2ContextAcc


meanRTsNone = ddply(subset(model_predictions2, Condition == "None"), ~ ListLength + association_strength, summarise, MeanRT = mean(ModelTimes))

meanRTs2 = ddply(subset(model_predictions2, Condition != "None"), ~ ListLength + ContextSize + Correct + Condition + association_strength, summarise, MeanRT = mean(ModelTimes))

# Plot ModelTimes based on ContextSize
Exp2ContextRT <- ggplot(
  subset(meanRTs2, Condition == "Context" & Correct == TRUE), 
  aes(x = ContextSize, y = MeanRT, color = as.factor(association_strength))) +
  # Add a smooth line
  geom_line(size = 1) +
  # Add labels
  labs(x = "Context size", y = "Model times", color = "AS") +
  # Use a minimal theme
  theme_minimal() +
  # Add horizontal lines for meanRTsNone as dotted lines
  geom_hline(data = meanRTsNone, aes(yintercept = MeanRT, color = as.factor(association_strength)), linetype = "dotted", size = 1) +
  # Adjust legend position
  theme(legend.position = "bottom")+
  theme(text = element_text(size = textsize, family = "sans"))+
  facet_wrap(~ListLength)
Exp2ContextRT

# plot the same for accuracy
meanAccNone = ddply(subset(model_predictions2, Condition == "None"), ~ ListLength + association_strength, summarise, MeanAcc = mean(Correct_01))

meanAcc2 = ddply(subset(model_predictions2, Condition == "Context"), ~ ListLength + ContextSize + association_strength, summarise, MeanAcc = mean(Correct_01))

Exp2ContextAcc <- ggplot(
  meanAcc2, 
  aes(x = ContextSize, y = MeanAcc, color = as.factor(association_strength))) +
  # Add a smooth line
  geom_line(size = 1) +
  # Add labels
  labs(x = "Context size", y = "Accuracy", color = "AS") +
  # Use a minimal theme
  theme_minimal() +
  geom_hline(data = meanAccNone, aes(yintercept = MeanAcc, color = as.factor(association_strength)), linetype = "dotted", size = 2) +
  # Adjust legend position
  theme(legend.position = "bottom")+
  theme(text = element_text(size = textsize, family = "sans"))+
  facet_wrap(~ListLength)
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
  
  sim_distractors = c(runif(ListLengthExp3-ContextSizes20[i], min = 0, max = 0.3), runif(ContextSizes20[i]-1, min = 0.7, max = 1))
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
                              ContextSize = context_size,
                              SimsCueDistractors = toString(sim_distractors)))
}

# Add the control conditions
# high pair similarity
for (i in 1:reps){
  # cue target similarity (sample from a uniform distribution between 0 and 1)
  sim_cue_target = runif(1, min = 0.7, max = 1)
  # distractors similarity
  sim_sum_distractors = sum(runif(ListLengthExp3-1, min = 0, max = 0.3)) + sim_cue_target
  sim_distractors = runif(ListLengthExp3-1, min = 0, max = 0.3)
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
                              ContextSize = context_size,
                              SimsCueDistractors = toString(sim_distractors)))
}

for (i in 1:reps){
  # cue target similarity (sample from a uniform distribution between 0 and 1)
  sim_cue_target = runif(1, min = 0, max = 0.3)
  # distractors similarity
  sim_sum_distractors = sum(runif(ListLengthExp3-1, min = 0, max = 0.3)) + sim_cue_target
  sim_distractors = runif(ListLengthExp3-1, min = 0, max = 0.3)
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
                              ContextSize = context_size,
                              SimsCueDistractors = toString(sim_distractors)))
}

# simulate the model
model_predictions3 = data.frame()
for (p_omission in p_omissions){
  print("p_omission")
  print(p_omission)
  for (p_intrusion in p_intrusions){
    for (association_strength in association_strengths){
      for (beta in betas){
        # simulate the model
        temp = modelSimulation(df3, p_omission, p_intrusion, association_strength, beta, 1)
        
        # add parameter values to the dataframe
        temp$p_omission = p_omission
        temp$p_intrusion = p_intrusion
        temp$association_strength = association_strength
        temp$beta = beta
        
        # add to the model predictions
        model_predictions3 = rbind(model_predictions3, temp)
      }
    }
  }
}

# save
write.csv(model_predictions3, "ModelPredictions3_withbeta.csv")
# read in
model_predictions3 = read.csv("ModelPredictions3_withbeta.csv")

# plot predictions Experiment 3 ------------------------------------------------

# create a coding for correct and incorrect correct is everything where the model response type is correct
model_predictions3$Correct = ifelse(model_predictions3$ModelResponseType == "Correct", TRUE, FALSE)
model_predictions3$Correct_01 = ifelse(model_predictions3$ModelResponseType == "Correct", 1, 0)

model_predictions3$ListLength = as.numeric(model_predictions3$ListLength)
# Function to run regression and extract results
run_regression3 <- function(var) {
  temp = model_predictions3 %>%
    mutate(ContextSize = scale(ContextSize, center = TRUE, scale = TRUE)) %>%
    group_by(across(all_of(var))) %>%
    do(tidy(lm(ModelTimes ~ Condition + ContextSize + Correct, data = .))) %>%
    filter(term %in% c("ContextSize", "CorrectTRUE", "ConditionHighPairSimilarity", "ConditionLowSimilarityList")) %>%
    mutate(lower = estimate - 1.96 * std.error,
           upper = estimate + 1.96 * std.error,
           variable = var)  # Store x values for correct plotting
  temp[["x_value"]] = temp[[var]]
  return(temp)
}

# Run regressions for all variables
results_all3 <- map_dfr(variables, run_regression3)

# Match variable names to readable labels
results_all3 <- results_all3 %>%
  mutate(variable = factor(variable, levels = variables, labels = labels))

# Create the combined plot
TimesExp3 = ggplot(results_all3, aes(x = x_value, y = estimate, color = term)) +
  geom_point(size = 1) +
  geom_line() +  
  #geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.02, alpha = 0.5) + 
  labs(title = "Experiment 3 model times",
       x = "Parameter value", 
       y = "Effect size") +
  theme_minimal() +
  coord_cartesian(ylim = c(-0.2, 0.8)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  facet_wrap(~ variable, scales = "free_x", nrow = 1) +  # Separate by parameter
  shared_color_scale +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))  # Place legend at the bottom
TimesExp3

# same for accuracy

run_logistic_regression3 <- function(var) {
  temp = model_predictions3 %>%
    mutate(ContextSize = scale(ContextSize, center = TRUE, scale = TRUE)) %>%
    group_by(across(all_of(var))) %>%
    do(tidy(glm(Correct ~ ContextSize + Condition, 
                data = ., 
                family = binomial(link = "logit")))) %>%
    filter(term %in% c("ContextSize", "ConditionHighPairSimilarity", "ConditionLowSimilarityList")) %>%
    mutate(lower = estimate - 1.96 * std.error,
           upper = estimate + 1.96 * std.error,
           variable = var)  # Store x values for correct plotting
  temp[["x_value"]] = temp[[var]]
  return(temp)
}

# Run regressions for all variables
results_all_log3 <- map_dfr(variables, run_logistic_regression3)

# Match variable names to readable labels
results_all_log3 <- results_all_log3 %>%
  mutate(variable = factor(variable, levels = variables, labels = labels))


# Plot
AccExp3 = ggplot(results_all_log3, aes(x = x_value, y = estimate, color = term)) +
  geom_point(size = 1) +
  geom_line() +  
  #geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.02, alpha = 0.5) + 
  geom_hline(yintercept = 0, linetype = "dotted") +
  facet_wrap(~ variable, scales = "free_x", nrow = 1) +  # Row layout
  shared_color_scale +
  labs(title = "Experiment 3 accuracy",
       x = "Parameter value", 
       y = "Effect size") +
  theme_minimal() +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))   # Place legend at the bottom
AccExp3

# Combine the figure --------------------------------------

library(patchwork)
combined_plot <- (TimesExp1 / TimesExp2 / TimesExp3) +
  plot_annotation(tag_levels = 'A')+
  theme(legend.position = "top")
combined_plot

# save the plot
savingPlace = "/Figures/RobustnessCheckModelTimes.pdf"
ggsave(paste(current_directory, savingPlace, sep = ""), combined_plot,  
       width = 15, 
       height = 10)


combined_plot2 <- (AccExp1 / AccExp2 / AccExp3) +
  # Add ABC plot annotation
  plot_annotation(tag_levels = 'A')+
  theme(legend.position = "top")
combined_plot2

# save the plot
savingPlace = "/Figures/RobustnessCheckAC.pdf"
ggsave(paste(current_directory, savingPlace, sep = ""), combined_plot2,  
       width = 15, 
       height = 10)


#------------------------------------------------------------------------------


# male Set sizea factor
model_predictions3$ListLength = as.factor(model_predictions3$ListLength)

# calculate mean RT for the none context condition for each set size
meanRTsNone = ddply(subset(model_predictions3, Condition != "Context"), ~ Condition + p_omission, summarise, MeanRT = mean(ModelTimes))

meanRTs3 = ddply(subset(model_predictions3, Condition == "Context"), ~ Condition + Correct + ContextSize + p_omission, summarise, MeanRT = mean(ModelTimes))

# Plot ModelTimes based on ContextSize
Exp3ContextRT <- ggplot(subset(meanRTs3, Condition == "Context" & Correct == TRUE), aes(x = ContextSize, y = MeanRT, color = as.factor(p_omission))) +
  # Add a smooth line
  geom_line(size = 1) +
  # Add labels
  labs(x = "Context size", y = "Model times") +
  # Use a minimal theme
  theme_minimal() +
  # Add horizontal lines for meanRTsNone as dotted lines
  geom_hline(data = meanRTsNone, aes(yintercept = MeanRT, linetype = Condition, color = as.factor(p_omission)), size = 1) +
  # Adjust legend position
  theme(text = element_text(size = 10, family = "sans"))
Exp3ContextRT

# plot the same for accuracy

# calculate mean RT for the none context condition for each set size
meanAccNone = ddply(subset(model_predictions3, Condition != "Context"), ~ Condition + p_omission, summarise, MeanAcc = mean(Correct_01))

# calculate mean RT for the none context condition for each set size
meanAcc3 = ddply(subset(model_predictions3, Condition == "Context"), ~ Condition + ContextSize + p_omission, summarise, MeanAcc = mean(Correct_01))

# Plot ModelTimes based on ContextSize
Exp3ContextAcc <- ggplot(
  subset(meanAcc3, Condition == "Context"), 
  aes(x = ContextSize, y = MeanAcc, color = as.factor(p_omission))) +
  # Add a smooth line
  geom_line(size = 1) +
  # Add labels
  labs(x = "Context size", y = "Model times") +
  # Use a minimal theme
  theme_minimal() +
  # Add horizontal lines for meanRTsNone as dotted lines
  geom_hline(data = meanAccNone, aes(yintercept = MeanAcc, linetype = Condition, color = as.factor(p_omission)), size = 1) +
  # Adjust legend position
  theme(text = element_text(size = 10, family = "sans"))
Exp3ContextAcc


# calculate mean RT for the none context condition for each set size
meanRTsNone = ddply(subset(model_predictions3, Condition != "Context"), ~ Condition + association_strength, summarise, MeanRT = mean(ModelTimes))

meanRTs3 = ddply(subset(model_predictions3, Condition == "Context"), ~ Condition + Correct + ContextSize + association_strength, summarise, MeanRT = mean(ModelTimes))

# Plot ModelTimes based on ContextSize
Exp3ContextRT <- ggplot(subset(meanRTs3, Condition == "Context" & Correct == TRUE), aes(x = ContextSize, y = MeanRT, color = as.factor(association_strength))) +
  # Add a smooth line
  geom_line(size = 1) +
  # Add labels
  labs(x = "Context size", y = "Model times") +
  # Use a minimal theme
  theme_minimal() +
  # Add horizontal lines for meanRTsNone as dotted lines
  geom_hline(data = meanRTsNone, aes(yintercept = MeanRT, linetype = Condition, color = as.factor(association_strength)), size = 1) +
  # Adjust legend position
  theme(text = element_text(size = 10, family = "sans"))
Exp3ContextRT

# plot the same for accuracy

# calculate mean RT for the none context condition for each set size
meanAccNone = ddply(subset(model_predictions3, Condition != "Context"), ~ Condition + association_strength, summarise, MeanAcc = mean(Correct_01))

# calculate mean RT for the none context condition for each set size
meanAcc3 = ddply(subset(model_predictions3, Condition == "Context"), ~ Condition + ContextSize + association_strength, summarise, MeanAcc = mean(Correct_01))

# Plot ModelTimes based on ContextSize
Exp3ContextAcc <- ggplot(
  subset(meanAcc3, Condition == "Context"), 
  aes(x = ContextSize, y = MeanAcc, color = as.factor(association_strength))) +
  # Add a smooth line
  geom_line(size = 1) +
  # Add labels
  labs(x = "Context size", y = "Model times") +
  # Use a minimal theme
  theme_minimal() +
  # Add horizontal lines for meanRTsNone as dotted lines
  geom_hline(data = meanAccNone, aes(yintercept = MeanAcc, linetype = Condition, color = as.factor(association_strength)), size = 1) +
  # Adjust legend position
  theme(text = element_text(size = 10, family = "sans"))
Exp3ContextAcc


# calculate mean RT for the none context condition for each set size
meanRTsNone = ddply(subset(model_predictions3, Condition != "Context"), ~ Condition + p_intrusion, summarise, MeanRT = mean(ModelTimes))

meanRTs3 = ddply(subset(model_predictions3, Condition == "Context"), ~ Condition + Correct + ContextSize + p_intrusion, summarise, MeanRT = mean(ModelTimes))

# Plot ModelTimes based on ContextSize
Exp3ContextRT <- ggplot(subset(meanRTs3, Condition == "Context" & Correct == TRUE), aes(x = ContextSize, y = MeanRT, color = as.factor(p_intrusion))) +
  # Add a smooth line
  geom_line(size = 1) +
  # Add labels
  labs(x = "Context size", y = "Model times") +
  # Use a minimal theme
  theme_minimal() +
  # Add horizontal lines for meanRTsNone as dotted lines
  geom_hline(data = meanRTsNone, aes(yintercept = MeanRT, linetype = Condition, color = as.factor(p_intrusion)), size = 1) +
  # Adjust legend position
  theme(text = element_text(size = 10, family = "sans"))
Exp3ContextRT

# plot the same for accuracy

# calculate mean RT for the none context condition for each set size
meanAccNone = ddply(subset(model_predictions3, Condition != "Context"), ~ Condition + p_intrusion, summarise, MeanAcc = mean(Correct_01))

# calculate mean RT for the none context condition for each set size
meanAcc3 = ddply(subset(model_predictions3, Condition == "Context"), ~ Condition + ContextSize + p_intrusion, summarise, MeanAcc = mean(Correct_01))

# Plot ModelTimes based on ContextSize
Exp3ContextAcc <- ggplot(
  subset(meanAcc3, Condition == "Context"), 
  aes(x = ContextSize, y = MeanAcc, color = as.factor(p_intrusion))) +
  # Add a smooth line
  geom_line(size = 1) +
  # Add labels
  labs(x = "Context size", y = "Model times") +
  # Use a minimal theme
  theme_minimal() +
  # Add horizontal lines for meanRTsNone as dotted lines
  geom_hline(data = meanAccNone, aes(yintercept = MeanAcc, linetype = Condition, color = as.factor(p_intrusion)), size = 1) +
  # Adjust legend position
  theme(text = element_text(size = 10, family = "sans"))
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
