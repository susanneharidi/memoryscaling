

# A sequential model of memory retrieval in cued- recall modellled with a geometric distribution

# p(correct) = sim(cue|target)ˆb + a / (sim(cue|target)ˆb + sum_over_i(sim(cue|distractor_i)ˆb) + a)
# P(stopping) = p(correct) + p(incorrect) * (p(omission) + p(intrusion))
# P(continue) = 1-p(stopping)
# p(correct I stopping) = p(correct)/ p(stopping)

# Load packages ------------------------------------------------------
library(ggplot2)
theme_set(
  theme_minimal() +
    theme(legend.position = "top")
)

# Image Manipulation
library(stringr)
library(lme4) 
library(brms)
library(dplyr)
library(patchwork)
library(here)

# plotting parameter --------------------------------------------------

textsize = 22
font = "sans"

# load some necessary functions ------------------------------------
current_directory <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(current_directory)
source("FunctionsSimSSModelVariation.R")


# load the dataset 'DatasetPrepedForSimulationExp1.csv' --------------------

d1 = read.csv('DatasetPrepedForSimulationExp1.csv')
validData1 = subset(d1, validID == 1 & validTrial == 1 & RT < 20000)
setwd(current_directory)
d3 = read.csv('DatasetPrepedForSimulationExp3.csv')
d3$RT = d3$RT*1000
validData3 = subset(d3, validID == 1 & validTrial == 1 & RT < 20000 & SimCueTarget > 0)

valid_d_combined = rbind(validData1, validData3)
# plot distribution of RTs for each experiment
ggplot(valid_d_combined, aes(x = RT, fill = ListLength, color = ListLength, group = ListLength)) +
  #geom_histogram(bins = 100, alpha = 0.5) +
  # density plot
  geom_density(alpha = 0.5) +
  labs(title = "RT Distribution",
       x = "RT",
       y = "Count") +
  facet_wrap(~Experiment)


# define the starting parameters -------------------------------------------

# starting Params SimSS beta scaling
p_omission = -1     
p_intrusion = -1
associationstrength = ScaleParms(1, 0, 100)
beta = ScaleParms(0, -10, 10)
params = c(p_omission, p_intrusion, associationstrength, beta)

# starting parameters softmax
p_omission = -1
p_intrusion = -1
associationstrength = ScaleParms(1, 0, 100)
lamda = ScaleParms(1, 0, 100)
params_softmax = c(p_omission, p_intrusion, associationstrength, lamda)

# starting paramepters for previously learned lists
p_omission = -1
p_intrusion = -1
associationstrength = ScaleParms(10, 0, 100)
beta = ScaleParms(0, -10, 10)
context_devaluation = ScaleParms(1, 0, 10)
params_previously_learned = c(p_omission, p_intrusion, associationstrength, beta, context_devaluation)

#starting parameters correct acceptance
p_omission = -1
p_intrusion = -1
p_acceptance = 10
associationstrength = ScaleParms(10, 0, 100)
beta = ScaleParms(0, -10, 10)
params_correct_acceptance = c(p_omission, p_intrusion, associationstrength, beta, p_acceptance)


ModelVersions = c("SimSS no Sim",
                  "SimSS",
                  "SimSS no AS",
                  "SimSS no Sim/AS",
                  "Softmax", 
                  "OutOfContextDists", 
                  "CorrectAcceptance")
ModelFittingFunctions = list(Total_neg_LL, 
                             Total_neg_LL, 
                             Total_neg_LL, 
                             Total_neg_LL, 
                             Total_neg_LL_softmax, 
                             Total_neg_LL_OutOfListDistractors, 
                             Total_neg_LL_correct_acceptance)
ModelParametersN = c(length(params)-1, 
                     length(params), 
                     length(params)-1, 
                     length(params)-2, 
                     length(params_softmax), 
                     length(params_previously_learned), 
                     length(params_correct_acceptance))
AllParameters = list(params, 
                     params, 
                     params, 
                     params, 
                     params_softmax, 
                     params_previously_learned, 
                     params_correct_acceptance)

DoneModelVersions = c("SimSS no Sim",
                      "SimSS",
                      "SimSS no AS",
                      "SimSS no Sim/AS",
                      "Softmax", 
                      "OutOfContextDists")

#############################§
# Fit the different model versions -------------------------------------------
#############################§
stepsizes = c(100, 300, 500, 700, 1000)
intercepts = c(0, 200, 400, 600, 800)

# Only do the fitting if it hasnt already been done
FittingDone = FALSE

if(FittingDone){
  # Load the fitted parameters
  current_directory <- dirname(rstudioapi::getSourceEditorContext()$path)
  setwd(current_directory)
  fitted_parameters <- read.csv("FittedParametersSimSS_ModelVersions.csv")
  fitted_parameters <- subset(fitted_parameters, select = -c(X))
}else{
  fitted_parameters <- data.frame(Model = character(0),
                                  OmissionProb = numeric(0),
                                  IntrusionProb = numeric(0),
                                  AssociationStrength = numeric(0))
  for (stepsize in stepsizes){
    for (intercept in intercepts){
      # make an empty dataframe, to save the fitted parameters for each model in
      
      for (i in 1:length(ModelVersions)){
        # skip everything below if the model has already been fitted
        if (ModelVersions[i] %in% DoneModelVersions){
          next
        }
        ModelVersion = ModelVersions[i]
        print("started Fitting")
        print(ModelVersion)
        ModelFittingFunction = ModelFittingFunctions[[i]]
        n_params = ModelParametersN[i]
        Parameters = AllParameters[[i]]
        
        if (ModelVersion == "SimSS no Sim" || ModelVersion == "SimSS no Sim/AS"){
          UseSimilarity = FALSE
        }else{
          UseSimilarity = TRUE
        }
        if (ModelVersion == "SimSS no AS" || ModelVersion == "SimSS no Sim/AS"){
          UseAssociationStrength = FALSE
        }else{
          UseAssociationStrength = TRUE
        }
          
        
          
          
          # do Optimisation ------------------------------------------------------------
          
          out = optim(par = Parameters, 
                      fn = ModelFittingFunction, 
                      data = valid_d_combined,
                      UseSimilarity = UseSimilarity,
                      UseAssociationStrength = UseAssociationStrength,
                      stepsize = stepsize,
                      intercept = intercept)
          
          #retransform parameters
          print(out)
          fitted_p_omission = logistic_function(out$par[1])
          fitted_p_intrusion = logistic_function(out$par[2])
          if(UseAssociationStrength){
            fitted_association_strength = UnscaleParms(out$par[3], 0, 100)
          }else{
            fitted_association_strength = 0
          }
          if (ModelVersion != "Softmax" && UseSimilarity){
            fitted_beta = UnscaleParms(out$par[4], -10, 10)
            fitted_lambda = 0
          }else if (ModelVersion == "Softmax"){
            fitted_beta = 0
            fitted_lambda = UnscaleParms(out$par[4], 0, 100)
          }else{
            fitted_beta = 0
            fitted_lambda = 0
          }
          if (ModelVersion == "OutOfContextDists"){
            fitted_context_devaluation = UnscaleParms(out$par[5], 0, 10)
          }else{
            fitted_context_devaluation = 0
          }
          if (ModelVersion == "CorrectAcceptance"){
            fitted_p_acceptance = logistic_function(out$par[5])
          }else{
            fitted_p_acceptance = 1
          }
          print(ModelVersion)
          # print("omission prob")
          # print(fitted_p_omission)
          # print("intrusion pron")
          # print(fitted_p_intrusion)
          # print("association strength")
          # print(fitted_association_strength)
          # print("beta")
          # print(fitted_beta)
          # print("lambda")
          # print(fitted_lambda)
          # print("context devaluation")
          # print(fitted_context_devaluation)
          # print("p acceptance")
          # print(fitted_p_acceptance)
          print("stepsize")
          print(stepsize)
          print("intercept")
          print(intercept)
          
          # get the LogLikelihood
          neglogLik_fit <- out$value 
          
          # calculate BIC and AIC
          n_obs = nrow(valid_d_combined)
          BIC = 2*neglogLik_fit + n_params*log(n_obs)
          AIC = 2*neglogLik_fit + 2*n_params
          
          
          cat("Model:", ModelVersion, "BIC:", BIC, "AIC:", AIC, "NegLogLik:", neglogLik_fit,"\n")
          
          
          # save the parameters in the dataframe
          fitted_parameters <- rbind(fitted_parameters, data.frame(Model = ModelVersion, 
                                                                   OmissionProb = fitted_p_omission, 
                                                                   IntrusionProb = fitted_p_intrusion, 
                                                                   AssociationStrength = fitted_association_strength,
                                                                   Beta = fitted_beta,
                                                                   Lambda = fitted_lambda,
                                                                   ContextDevaluation = fitted_context_devaluation,
                                                                   PAcceptance = fitted_p_acceptance,
                                                                   StepSize = stepsize,
                                                                   Intercept = intercept,
                                                                   neglogLik_fit = neglogLik_fit,
                                                                   BIC = BIC,
                                                                   AIC = AIC))
          
          
          # Now this here is extra, to make the fitting easier. since the first model is a nested model in the following three, lets adjust the starting paramepters
          if (ModelVersion == "SimSS no Sim"){
            # if the model is SimSS, then we can use the parameters from the other models as starting parameters
            params_adjusted = c(out$par[1], out$par[2], ScaleParms(0.01, 0, 100), 0)
            AllParameters[[i+1]] = params_adjusted
            print("Adjusted parameters for")
            print(ModelVersions[i+1])
            AllParameters[[i+2]] = params_adjusted
            print("Adjusted parameters for")
            print(ModelVersions[i+2])
            AllParameters[[i+3]] = params_adjusted
            print("Adjusted parameters for")
            print(ModelVersions[i+3])
          }
      }
    }
  }
  
  # save the fitted parameters in a file csv file
  current_directory <- dirname(rstudioapi::getSourceEditorContext()$path)
  setwd(current_directory)
  write.csv(fitted_parameters, "FittedParametersSimSS_ModelVersions.csv")
}

ggplot(fitted_parameters, aes(x = StepSize, y = BIC, color = as.factor(Model), group = Model)) +
  geom_point(alpha = 0.5)+
  #coord_cartesian(ylim = c(209000, 212000)) +
  labs(title = "BIC values for the different models",
       x = "Model",
       y = "BIC") +
  facet_grid(~Intercept) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# lets make the above plot a heat map with each model getting its own heat map across stepsize and intercept
heatmapStepInterPlot = ggplot(fitted_parameters, aes(x = as.factor(StepSize), y = as.factor(Intercept), fill = BIC)) +
  geom_tile() +
  # use plasma color scale
  scale_fill_gradientn(colors = c("white", "yellow", "red")) +
  labs(title = "BIC values for the different models",
       x = "Step Size",
       y = "Intercept") +
  facet_wrap(~Model) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
heatmapStepInterPlot

# For each model lets find the step size and intercept, where the BIC is the lowest
bestModels = data.frame(Model = character(0), 
                          StepSize = numeric(0), 
                          Intercept = numeric(0), 
                          BIC = numeric(0))
for (model in ModelVersions){
  # get the model
  model_data = subset(fitted_parameters, Model == model)
  # get the minimum BIC
  min_BIC = min(model_data$BIC)
  # get the row with the minimum BIC
  min_BIC_row = model_data[which.min(model_data$BIC),]
  # print the row
  print(min_BIC_row)
  # add the row to the bestModels dataframe
  bestModels = rbind(bestModels, min_BIC_row)
}

ggplot(bestModels, aes(x = Model, y = BIC, color = as.factor(StepSize), group = StepSize)) +
  geom_point() +
  labs(title = "BIC values for the different models",
       x = "Model",
       y = "BIC") +
  facet_grid(~Intercept) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

bestModels$DeltaBIC = bestModels$BIC - subset(bestModels, Model == "SimSS")$BIC
BICPlot = ggplot(bestModels, aes(x = reorder(Model, DeltaBIC), y = DeltaBIC)) +
  geom_col(fill = "steelblue") +
  geom_text(aes(label = round(BIC, 1)), vjust = -0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(
    title = "BIC from base model",
    x = "Model",
    y = "BIC difference"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
BICPlot

bestModels$neglogLik_fitDelta = bestModels$neglogLik_fit - subset(bestModels, Model == "SimSS")$neglogLik_fit
ggplot(bestModels, aes(x = reorder(Model, neglogLik_fitDelta), y = neglogLik_fitDelta)) +
  geom_col(fill = "steelblue") +
  geom_text(aes(label = round(neglogLik_fit, 1)), vjust = -0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  # make x limit
  
  labs(
    title = "ΔNegLogLik from Best Model",
    x = "Model",
    y = "ΔNegLogLik"
  ) +
  theme_minimal()




#############################§
# Simulate the models -------------------------------------------------------
#############################§

Replications = 10
SimulationDone = FALSE
if(SimulationDone){
  # Load the fitted parameters
  current_directory <- dirname(rstudioapi::getSourceEditorContext()$path)
  setwd(current_directory)
  simulateddata <- read.csv("SimulatedData_geometric_modelVersions_best.csv")
  simulateddata <- subset(simulateddata, select = -c(X))
}else{
  First = TRUE
  for (j in 1:length(bestModels$Model)){
    modelVersion = bestModels$Model[j]
    stepsize_temp = bestModels$StepSize[j]
    print(modelVersion)
    
    # get the relevant parameters from the fitted_parameters dataframe
    relevant_fitted_parameters = bestModels[j,]
    
    # do the following a hundred times
    for (i in 1:Replications){
      
      if (First){
        simulateddata = modelSimulation(valid_d_combined,
                                        fitted_params = relevant_fitted_parameters, 
                                        stepsize = stepsize_temp,
                                        ModelVersion = modelVersion,
                                        repetitions = 1)
        
        
        #Transform the times into the approximately same space
        simulateddata$Model = modelVersion
        simulateddata$Run = i
        simulateddata$StepSize = stepsize_temp
        simulateddata$ModelRT = simulateddata$ModelRT + bestModels$Intercept[j]
        simulateddata$Intercept = bestModels$Intercept[j]
        simulateddata$Intercept = 0
        
        First = FALSE
        
      }else{
        temp = modelSimulation(valid_d_combined,
                               fitted_params = relevant_fitted_parameters, 
                               stepsize = stepsize_temp,
                               ModelVersion = modelVersion,
                               repetitions = 1)
        
        
        #Transform the times into the approximately same space
        temp$Model = modelVersion
        temp$Run = i
        temp$StepSize = stepsize_temp
        temp$ModelRT = temp$ModelRT + bestModels$Intercept[j]
        temp$Intercept = 0
        simulateddata = rbind(simulateddata, temp)
      }
      print(i)
    }
  }
  
  # save the simulated data in a csv file
  current_directory <- dirname(rstudioapi::getSourceEditorContext()$path)
  setwd(current_directory)
  write.csv(simulateddata, "SimulatedData_geometric_modelVersions_best.csv")
}

# add participant RT
temp = subset(simulateddata, Model == "SimSS")
temp$ModelRT = temp$ParticipantRT
temp$ModelResponseType = temp$ParticipantResponseType
temp$Model = "Participant"
simulateddata = rbind(simulateddata, temp)

# rename CorrectAcceptance into SimSS with correct acceptance
simulateddata$Model[simulateddata$Model == "CorrectAcceptance"] = "SimSS target rejection"
# rename OutOfContextDists into SimSS with out of context distractors
simulateddata$Model[simulateddata$Model == "OutOfContextDists"] = "SimSS out of list distractors"
# rename Softmax into SimSS with softmax
simulateddata$Model[simulateddata$Model == "Softmax"] = "SimSS softmax"

fitted_parameters$Model[fitted_parameters$Model == "CorrectAcceptance"] = "SimSS target rejection"
fitted_parameters$Model[fitted_parameters$Model == "OutOfContextDists"] = "SimSS out of list distractors"
fitted_parameters$Model[fitted_parameters$Model == "Softmax"] = "SimSS softmax"
# same for best models
bestModels$Model[bestModels$Model == "CorrectAcceptance"] = "SimSS target rejection"
bestModels$Model[bestModels$Model == "OutOfContextDists"] = "SimSS out of list distractors"
bestModels$Model[bestModels$Model == "Softmax"] = "SimSS softmax"

# redo plots with these names:
BICPlot = ggplot(bestModels, aes(x = reorder(Model, DeltaBIC), y = DeltaBIC)) +
  geom_col(fill = "steelblue") +
  geom_text(aes(label = round(BIC, 1)), vjust = -0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(
    x = "Model",
    y = "BIC difference to base model"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = textsize, family = font))
BICPlot

heatmapStepInterPlot = ggplot(fitted_parameters, aes(x = as.factor(StepSize), y = as.factor(Intercept), fill = BIC)) +
  geom_tile() +
  # use plasma color scale
  scale_fill_gradientn(colors = c("white", "steelblue", "midnightblue"),
                       breaks = range(fitted_parameters$BIC),  
                       labels = round(range(fitted_parameters$BIC))) +
  # remove tics for BIC values color scale
  labs(
       x = "Step size",
       y = "Intercept") +
  # with four plots per row
  facet_wrap(~Model, nrow = 2) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = textsize, family = font),
        legend.position = "right")
heatmapStepInterPlot

# save that plot
savingPlace = "/Figures/StepSizeInterceptHeatmap.pdf"
ggsave(paste(current_directory, savingPlace, sep = ""), heatmapStepInterPlot,  
       width = 30, 
       height = 10)




#############################§
# Plotting ---------------------
#############################§
# Define color mapping and readable names
model_colors <- c(
          "SimSS" = "#E69F00",
          "SimSS no AS" = "#56B4E9",
          "SimSS no Sim/AS" = "#009E73",
          "SimSS no Sim" = "#F0E442",
          "SimSS target rejection" = "#0072B2",
          "SimSS out of list distractors" = "#D55E00",
          "SimSS softmax" = "#CC79A7",
          "Participant" = "black")

# Shared color scale
shared_color_scale <- scale_color_manual(
  values = model_colors,
  name = "Model"
)

shared_fill_scale <- scale_fill_manual(
  values = model_colors,
  name = "Model"
)


# Plot distribution of RTs for each experiment -----------------------

# plot distribution of RTs for each experiment
ggplot(subset(simulateddata, Experiment == "Experiment 1"), aes(x = ModelRT, fill = ListLength, color = ListLength, group = ListLength)) +
  #geom_histogram(bins = 100, alpha = 0.5) +
  # density plot
  geom_density(alpha = 0.5) +
  # axis limits
  labs(title = "RT Distribution",
       x = "RT",
       y = "Desnity") +
  facet_wrap(~Model)

ggplot(subset(simulateddata, Experiment == "Experiment 3"), aes(x = ModelRT, fill = ContextSize, color = ContextSize, group = ContextSize)) +
  #geom_histogram(bins = 100, alpha = 0.5) +
  # density plot
  geom_density(alpha = 0.5) +
  labs(title = "RT Distribution",
       x = "RT",
       y = "Desnity") +
  facet_wrap(~Model)


# Plot RT Relationship ------------------------------------------------------
# summarise ModelRt over runs

ggplot(simulateddata, aes(x = ParticipantRT/1000, y = ModelRT/1000, color = Model)) +
  # make y limits
  geom_point(alpha = 0.01)+
  shared_color_scale +
  coord_cartesian(ylim = c(0, 20), xlim = c(0, 20)) +
  geom_smooth(se = FALSE, linetype = "dashed") +
  labs(x = "Participant RT in s",
       y = "Model RT in s") +
  theme(text = element_text(size = textsize, family = font))


# make trial variable which goes from the first trial in a run to the last trial
temp1 = ddply(simulateddata, .(Model, Experiment, ListLength, Condition, ContextSize, Trials, CueTargetSimilarity, ParticipantResponseType), 
             summarise, 
             ModelRT = mean(ModelRT), 
             ParticipantRT = mean(ParticipantRT))

ggplot(temp1, aes(x = ParticipantRT/1000, y = ModelRT/1000, color = Model)) +
  # make y limits
  geom_point(alpha = 0.1)+
  shared_color_scale+
  geom_smooth(se = FALSE, linetype = "dashed") +
  labs(title = "RT relationship",
       x = "Participant RT in s",
       y = "Model RT in s") +
  theme(text = element_text(size = textsize, family = font)) +
  facet_wrap(~Experiment)

breaks = seq(0, 1, by = 0.2)
breaks[length(breaks)] = breaks[length(breaks)] +0.01
temp1$simBin <- cut(temp1$CueTargetSimilarity, 
                            breaks = breaks, 
                            labels = c(0.1, 0.3, 0.5, 0.7, 0.9), 
                            include.lowest = TRUE, 
                            right = FALSE)

temp2 = ddply(temp1, .(Model, Experiment, ListLength, Condition, ContextSize, simBin), 
             summarise, 
             ModelRT = mean(ModelRT), 
             ParticipantRT = mean(ParticipantRT))

Experiment1aggregated = ggplot(subset(temp2, Experiment == "Experiment 1"), 
                               aes(x = ParticipantRT/1000, y = ModelRT/1000, color = Model, fill = Model)) +
  # make y limits
  geom_point(alpha = 0.1)+
  shared_color_scale +
  shared_fill_scale +
  geom_smooth(se = FALSE, linetype = "dashed", method = lm) +
  labs(x = "Participant RT in s",
       y = "Model RT in s") +
  theme(text = element_text(size = textsize, family = font),
        legend.position = "None") +
  facet_wrap(~Experiment)
Experiment1aggregated

temp3 = ddply(temp1, .(Model, Experiment, ListLength, Condition, ContextSize), 
              summarise, 
              ModelRT = mean(ModelRT), 
              ParticipantRT = mean(ParticipantRT))
Experiment3aggregated = ggplot(subset(temp3, Experiment == "Experiment 3"), 
                               aes(x = ParticipantRT/1000, y = ModelRT/1000, color = Model, fill = Model)) +
  # make y limits
  geom_point(alpha = 0.1)+
  shared_color_scale+
  shared_fill_scale +
  geom_smooth(se = FALSE, linetype = "dashed", method = lm) +
  labs(x = "Participant RT in s",
       y = "Model RT in s") +
  theme(text = element_text(size = textsize, family = font)) +
  facet_wrap(~Experiment)
Experiment3aggregated

joint_aggregated_plot = Experiment1aggregated /Experiment3aggregated +
  plot_layout(ncol = 1) +
  # collect guid
  plot_annotation(theme = theme(plot.title = element_text(size = textsize, family = font)))+
  plot_layout(guides = "collect")
joint_aggregated_plot


# Plot the response types ---------------------------------------------------
detach("package:plyr", unload=TRUE)
data_percent <- simulateddata %>%
  group_by(Experiment, Model, ModelResponseType) %>%
  summarise(Count = dplyr::n(), .groups = "drop") %>%
  group_by(Experiment, Model) %>%
  mutate(Percent = Count / sum(Count) * 100)


# Step 2: Plot
ResponsetypePlot = ggplot(data_percent, aes(x = ModelResponseType, y = Percent, fill = Model)) +
  geom_col(position = "dodge") +
  shared_color_scale +
  shared_fill_scale +
  scale_alpha_manual(values = c(0.5, 1)) +
  labs(x = "Response type",
       y = "Percentage") +
  theme(text = element_text(size = textsize, family = font),
        legend.position = "None") +
  facet_wrap(~Experiment)
ResponsetypePlot

simulateddata$Responsetypematching = simulateddata$ParticipantResponseType == simulateddata$ModelResponseType

# plot proportion of responsetype matching

data_percent <- simulateddata %>%
  group_by(Experiment, Model, Responsetypematching) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(Experiment, Model) %>%
  mutate(Percent = Count / sum(Count) * 100)

# Step 2: Plot
ggplot(data_percent, aes(x = Responsetypematching, y = Percent, fill = Model)) +
  geom_col(position = "dodge") +
  shared_color_scale +
  shared_fill_scale +
  scale_alpha_manual(values = c(0.5, 1)) +
  labs(x = "Matching responses",
       y = "Percentage") +
  theme(text = element_text(size = textsize, family = font)) +
  facet_wrap(~Experiment)


# Calculate correlations and plot them ---------------------------------------
library(plyr)
correlations <- ddply(simulateddata, .(Model, Experiment), summarise, cor = cor(ParticipantRT, ModelRT))
TrialCorrlationPlot = ggplot(subset(correlations, Model != "Participant"), aes(x = Model, y = cor, fill = Experiment)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Model",
       y = "Trial by trial correlation") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = textsize, family = font))
TrialCorrlationPlot

# Plot RTs against Similarity ------------------------------------------------

ggplot(simulateddata, aes(x = CueTargetSimilarity, y = ModelRT, color = Model)) +
  # make y limits
  #geom_point(alpha = 0.1) +
  geom_smooth(se = FALSE) +
  labs(x = "Cue target similarity",
       y = "Model RT") +
  shared_color_scale +
  shared_fill_scale +
  theme(text = element_text(size = textsize, family = font)) +
  facet_wrap(~Experiment)

ggplot(simulateddata, aes(x = ListLength, y = ModelRT, color = Model)) +
  #geom_point() +
  # make y limits
  geom_smooth(se = FALSE) +
  labs(title = "RT Relationship",
       x = "List Length",
       y = "Model RT") +
  facet_wrap(~Experiment)

ggplot(simulateddata, aes(x = ContextSize, y = ModelRT, color = Model)) +
  geom_smooth(se = FALSE) +
  labs(title = "RT Relationship",
       x = "Context Size",
       y = "Model RT") +
  facet_wrap(~Experiment)

# Plot set size effect ------------------------------------------------

ggplot(subset(simulateddata, Experiment == "Experiment 1"), aes(x = ListLength, y = ModelRT, color = Model)) +
  # make y limits
  #geom_point(alpha = 0.1) +
  geom_smooth(se = FALSE) +
  labs(x = "Set size",
       y = "Model RT") +
  shared_color_scale +
  shared_fill_scale +
  theme(text = element_text(size = textsize, family = font)) +
  facet_wrap(~Experiment)

# calculate mean accuracy across set sizes and models
dSetSize = ddply(simulateddata, 
                 .(Model, Experiment, ListLength), 
                 summarise, 
                 MeanAccuracy = mean(ModelResponseType == "Correct"))
# plot that
SetSizePlot = ggplot(subset(dSetSize, Experiment == "Experiment 1"), aes(x = ListLength, y = MeanAccuracy, color = Model)) +
  geom_smooth(se = FALSE) +
  labs(x = "Set size",
       y = "Mean accuracy") +
  shared_color_scale +
  shared_fill_scale +
  theme(text = element_text(size = textsize, family = font)) +
  facet_wrap(~Experiment)
SetSizePlot

# Plot RT for the diffrent reponsetypes --------------------------------------
#calculate mean Model Time for each model and responsetype
dResponsetypeRT = ddply(simulateddata, 
      .(Model, ModelResponseType, Experiment), 
      summarise, 
      MeanModelTime = mean(ModelRT), 
      MeanParticipantTime = mean(ParticipantRT)) 
# plot that 
participant_means = subset(dResponsetypeRT, Model == "Participant")

# Then, modify your plot
ResponsetypeRTPlot = ggplot(subset(dResponsetypeRT, Model != "Participant"), aes(x = Model, y = MeanModelTime/1000, fill = ModelResponseType)) +
  geom_hline(data = subset(participant_means, ModelResponseType == "Correct"),
             aes(yintercept = MeanModelTime/1000),
             linetype = "solid", color = "forestgreen", size = 1.3) +
  geom_hline(data = subset(participant_means, ModelResponseType == "Intrusion"),
             aes(yintercept = MeanModelTime/1000),
             linetype = "solid", color = "red", size = 1.3) +
  geom_hline(data = subset(participant_means, ModelResponseType == "Omission"),
             aes(yintercept = MeanModelTime/1000),
             linetype = "solid", color = "darkred", size = 1.3) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("forestgreen", "red", "darkred")) +
  labs(x = "Model",
       y = "Mean RT in s",
       fill = "Response") +
  # Here: Use geom_hline with data and aes
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = textsize, family = font)) +
  facet_wrap(~Experiment)
ResponsetypeRTPlot
# Plot similarity binned -----------------------------------------------------
breaks = seq(0, 1, by = 0.2)
breaks[length(breaks)] = breaks[length(breaks)] +0.01
simulateddata$simBin <- cut(simulateddata$CueTargetSimilarity, 
                            breaks = breaks, 
                            labels = c(0.1, 0.3, 0.5, 0.7, 0.9), 
                            include.lowest = TRUE, 
                            right = FALSE)
d1$simBin <- cut(d1$SimCueTarget, 
               breaks = breaks, 
               labels = c(0.1, 0.3, 0.5, 0.7, 0.9), 
               include.lowest = TRUE, 
               right = FALSE)
# aggregate data across sim bins and set size
aggregatedRTAndModelTimes <- ddply(simulateddata, 
                                   ~simBin + Model + Run + Experiment, 
                                   summarize, MeanModelTime = mean(ModelRT), MeanRT = mean(ParticipantRT))
aggregatedRTAndModelTimes$simBin <- as.numeric(as.character(aggregatedRTAndModelTimes$simBin))
ggplot(simulateddata, aes(x = CueTargetSimilarity, y = ParticipantRT)) +
  # make y limits
  #geom_point(alpha = 0.1) +
  geom_smooth(color = "black")+
  # add mean
  geom_line(data = aggregatedRTAndModelTimes, aes(x = simBin, y = MeanRT), color = "black") +
  geom_smooth(aes(x = CueTargetSimilarity, y = ModelRT, color = Model), se = FALSE) +
  labs(title = "RT Relationship",
       x = "Cue Target Similarity",
       y = "Model RT") +
  facet_wrap(~Experiment)

# Aggregate data across sim bins and set size -------------------------------
Exp1aggregatedRTAndModelTimes <- ddply(subset(simulateddata, Experiment == "Experiment 1"), 
                                       ~simBin + Model + Experiment + ListLength + Run,
                                       summarize, MeanModelRT = mean(ModelRT), MeanRT = mean(ParticipantRT))

Exp1aggregatedRTAndModelTimes$simBin <- as.numeric(as.character(aggregatedRTAndModelTimes$simBin))


# plot aggregated model times against participant times
ggplot(subset(Exp1aggregatedRTAndModelTimes, Model != "Participant"), aes(x = MeanRT, y = MeanModelRT, color = Model)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Aggregated RT Relationship. Experiment 1",
       x = "Participant RT",
       y = "Model RT")



Exp3aggregatedRTAndModelTimes <- ddply(subset(simulateddata, Experiment == "Experiment 3"), 
                                       ~ Model + Experiment + ContextSize + Condition  + Run,
                                       summarize, MeanModelRT = mean(ModelRT), MeanRT = mean(ParticipantRT))


# plot aggregated model times against participant times
ggplot(subset(Exp3aggregatedRTAndModelTimes, Model != "Participant"), aes(x = MeanRT, y = MeanModelRT, color = Model)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "aggregated RT Relationship Experiment 3",
       x = "Participant RT",
       y = "Model RT") +
  facet_wrap(~Experiment)


# calculate correlations
correlations1 <- ddply(Exp1aggregatedRTAndModelTimes, .(Model, Experiment), summarise, cor = cor(MeanRT, MeanModelRT))
correlations2 <- ddply(Exp3aggregatedRTAndModelTimes, .(Model, Experiment), summarise, cor = cor(MeanRT, MeanModelRT))

correlations = rbind(correlations1, correlations2)
AggrgatedCorrelationsPlot = ggplot(subset(correlations, Model != "Participant"), aes(x = Model, y = cor, fill = Experiment)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Model",
       y = "Correlation of aggregated data") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = textsize, family = font))
AggrgatedCorrelationsPlot

ggplot(correlations, aes(x = Experiment, y = cor, fill = Model)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Correlations between Participant RT and Model RT",
       x = "Model",
       y = "Correlation") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# remove x Axis ticks from TrialCorrlationPlot
TrialCorrlationPlot = TrialCorrlationPlot +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        text = element_text(size = textsize, family = font),
        legend.position = "None")

CorrelationFigure = TrialCorrlationPlot/AggrgatedCorrelationsPlot
CorrelationFigure

ResponsetypeFigure = (ResponsetypePlot/ResponsetypeRTPlot) 

combined_model_figure = BICPlot + CorrelationFigure+ joint_aggregated_plot + ResponsetypeFigure+
  plot_layout(ncol = 4) +
  plot_annotation(tag_levels = "A")
combined_model_figure
# save that plot
savingPlace = "/Figures/CombinedModelFigure.pdf"
ggsave(paste(current_directory, savingPlace, sep = ""), combined_model_figure,  
       width = 30, 
       height = 20)


