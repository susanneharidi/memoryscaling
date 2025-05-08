

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
library(plyr)
library(dplyr)
library(patchwork)
library(here)

# Define some necessary functions ------------------------------------

#standard error function
se <- function(x){sd(x)/sqrt(length(x))}


get_mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# Define the function for logistic transformation
logistic_function <- function(x) {
  return(1 / (1 + exp(-x)))
}

logistic_function_reversed <- function(y) {
  return(log(y / (1 - y)))
}


ScaleParms <- function(tpars,lb,ub){
  -log(((ub-lb)/(tpars-lb))-1)
}
UnscaleParms <- function(tpars,lb,ub){
  (ub-lb)/(1+exp(-tpars))+lb
}




Responsetype_log_likelihood <- function(p_correct, p_stopping, p_omission, p_intrusion, responsetypes, p_acceptance = 1) {
  LL_total = 0
  for (i in 1:length(responsetypes)){
    if(responsetypes[i] == "Correct"){
      LL_total = LL_total+log(max(p_correct[i]*p_acceptance/p_stopping[i],0))
    }else if (responsetypes[i] == "Omission"){
      LL_total = LL_total+log(max((1-p_correct[i])*p_omission/p_stopping[i],0))
    }else {
      LL_total = LL_total+log(max((1-p_correct[i])*p_intrusion/p_stopping[i],0))
    }
  }
  return(-LL_total)
}




Total_neg_LL <- function(params, data, stepsize, intercept = 0, listlength = 0, UseSimilarity = TRUE, UseAssociationStrength = TRUE) {
  responsetypes = data$ErrorType
  sim_cue_target = data$SimCueTarget
  sim_distractors = data$SimsCueDistractors
  rt = data$RT
  scaled_p_omission <- logistic_function(params[1])
  scaled_p_intrusion <- logistic_function(params[2])
  
  if(UseAssociationStrength){
    associationStrength = UnscaleParms(params[3], 0, 100)
  }else{
    associationStrength = 0
  }
  
  if(UseSimilarity){
    scaled_beta <- UnscaleParms(params[4], -10, 10)
  }else{
    scaled_beta <- 0
  }
  
  p_correct = c()
  for (i in 1:length(sim_cue_target)){
    sim_distractors_vec = as.numeric(strsplit(gsub("\\[|\\]", "", sim_distractors[i]), ",")[[1]])
    if (length(sim_distractors_vec) == 0){
      sim_distractors_vec = c(0)
    }
    p_correct = c(p_correct, min((sim_cue_target[i]**scaled_beta + associationStrength)/
                                   (sim_cue_target[i]**scaled_beta + associationStrength + sum(sim_distractors_vec**scaled_beta)), 1))
  }
  
  p_stopping = pmax((p_correct + (1-p_correct)*(min((scaled_p_omission + scaled_p_intrusion), 1))), 0.000000000001)
  if (max(p_stopping) > 1){
    print("something went wrong")
    print(scaled_p_omission)
    print(scaled_p_intrusion)
    print(max(p_correct))
    print(params[3])
  }
  
  LL_rt <- geometric_log_likelihood_corrected(p_stopping, rt, stepsize, intercept)
  
  LL_response_type <- Responsetype_log_likelihood(p_correct, p_stopping, scaled_p_omission, scaled_p_intrusion, responsetypes)
  
  TotalLL = LL_rt + LL_response_type
  
  return(TotalLL)
} 




softmax <- function(x, scale = 1) {
  exp_x = exp(x * scale)
  return(exp_x / sum(exp_x))
}

Total_neg_LL_softmax <- function(params, data, stepsize, intercept = 0, listlength = 0, UseSimilarity = TRUE, UseAssociationStrength = TRUE) {
  responsetypes = data$ErrorType
  sim_cue_target = data$SimCueTarget
  sim_distractors = data$SimsCueDistractors
  rt = data$RT
  scaled_p_omission <- logistic_function(params[1])
  scaled_p_intrusion <- logistic_function(params[2])
  
  if(UseAssociationStrength){
    associationStrength = UnscaleParms(params[3], 0, 100)
  }
  else{
    associationStrength = 0
  }
  
  if(UseSimilarity){
    scaled_lamda <- UnscaleParms(params[4], 0, 100)
  }
  else{
    scaled_lamda <- 0
  }
  
  p_correct = c()
  for (i in 1:length(sim_cue_target)){
    sim_distractors_vec = as.numeric(strsplit(gsub("\\[|\\]", "", sim_distractors[i]), ",")[[1]])
    if (length(sim_distractors_vec) == 0){
      sim_distractors_vec = c(0)
    }
    p_correct = c(p_correct, softmax(c(sim_cue_target[i] + associationStrength, sim_distractors_vec), scale = scaled_lamda)[1])
  }
  
  p_stopping = pmax((p_correct + (1-p_correct)*(min((scaled_p_omission + scaled_p_intrusion), 1))), 0.000000000001)
  if (max(p_stopping) > 1){
    print("something went wrong")
    print(scaled_p_omission)
    print(scaled_p_intrusion)
    print(max(p_correct))
    print(params[3])
  }
  
  LL_rt <- geometric_log_likelihood_corrected(p_stopping, rt, stepsize, intercept)
  
  LL_response_type <- Responsetype_log_likelihood(p_correct, p_stopping, scaled_p_omission, scaled_p_intrusion, responsetypes)
  
  TotalLL = LL_rt + LL_response_type
  
  return(TotalLL)
} 



Total_neg_LL_OutOfListDistractors <- function(params, data, stepsize, intercept = 0, listlength = 0, UseSimilarity = TRUE, UseAssociationStrength = TRUE) {
  responsetypes = data$ErrorType
  sim_cue_target = data$SimCueTarget
  sim_distractors = data$SimsCueDistractors
  rt = data$RT
  sim_out_of_list_distractors = data$SimsCueOutOfListDistractors
  scaled_p_omission <- logistic_function(params[1])
  scaled_p_intrusion <- logistic_function(params[2])
  
  if(UseAssociationStrength){
    associationStrength = UnscaleParms(params[3], 0, 100)
  }
  else{
    associationStrength = 0
  }
  
  if(UseSimilarity){
    scaled_beta <- UnscaleParms(params[4], -10, 10)
  }
  else{
    scaled_beta <- 0
  }
  
  scaled_context_devaluation <- UnscaleParms(params[5], 0, 10)
  
  p_correct = c()
  for (i in 1:length(sim_cue_target)){
    sim_distractors_vec = as.numeric(strsplit(gsub("\\[|\\]", "", sim_distractors[i]), ",")[[1]])
    sim_out_of_list_distractors_vec = as.numeric(strsplit(gsub("\\[|\\]", "", sim_out_of_list_distractors[i]), ",")[[1]])
    if (length(sim_distractors_vec) == 0){
      sim_distractors_vec = c(0)
    }
    if (length(sim_out_of_list_distractors_vec) == 0){
      sim_out_of_list_distractors_vec = c(0)
    }
    p_correct = c(p_correct, min((sim_cue_target[i]**scaled_beta + associationStrength)/
                                   (sim_cue_target[i]**scaled_beta + 
                                      associationStrength + 
                                      sum(sim_distractors_vec**scaled_beta) +
                                      sum(scaled_context_devaluation * (sim_out_of_list_distractors_vec**scaled_beta))), 1))
  }
  
  p_stopping = pmax((p_correct + (1-p_correct)*(min((scaled_p_omission + scaled_p_intrusion), 1))), 0.000000000001)
  if (max(p_stopping) > 1){
    print("something went wrong")
    print(scaled_p_omission)
    print(scaled_p_intrusion)
    print(max(p_correct))
    print(params[3])
  }
  
  LL_rt <- geometric_log_likelihood_corrected(p_stopping, rt, stepsize, intercept)
  
  LL_response_type <- Responsetype_log_likelihood(p_correct, p_stopping, scaled_p_omission, scaled_p_intrusion, responsetypes)
  
  TotalLL = LL_rt + LL_response_type
  
  return(TotalLL)
} 

softmax <- function(x, scale = 1) {
  exp_x = exp(x * scale)
  return(exp_x / sum(exp_x))
}

Total_neg_LL_correct_acceptance <- function(params, data, stepsize, intercept = 0, listlength = 0, UseSimilarity = TRUE, UseAssociationStrength = TRUE) {
  responsetypes = data$ErrorType
  sim_cue_target = data$SimCueTarget
  sim_distractors = data$SimsCueDistractors
  rt = data$RT
  scaled_p_omission <- logistic_function(params[1])
  scaled_p_intrusion <- logistic_function(params[2])
  scaled_p_acceptance <- logistic_function(params[5])
  
  if(UseAssociationStrength){
    associationStrength = UnscaleParms(params[3], 0, 100)
  }
  else{
    associationStrength = 0
  }
  
  if(UseSimilarity){
    scaled_beta <- UnscaleParms(params[4], -10, 10)
  }
  else{
    scaled_beta <- 0
  }
  
  p_correct = c()
  for (i in 1:length(sim_cue_target)){
    sim_distractors_vec = as.numeric(strsplit(gsub("\\[|\\]", "", sim_distractors[i]), ",")[[1]])
    if (length(sim_distractors_vec) == 0){
      sim_distractors_vec = c(0)
    }
    p_correct = c(p_correct, min((sim_cue_target[i]**scaled_beta + associationStrength)/
                                   (sim_cue_target[i]**scaled_beta + associationStrength + sum(sim_distractors_vec**scaled_beta)), 1))
  }
  
  p_stopping = pmax((p_correct*scaled_p_acceptance + (1-p_correct)*(min((scaled_p_omission + scaled_p_intrusion), 1))), 0.000000000001)
  if (max(p_stopping) > 1){
    print("something went wrong")
    print(scaled_p_omission)
    print(scaled_p_intrusion)
    print(max(p_correct))
    print(params[3])
  }
  
  LL_rt <- geometric_log_likelihood_corrected(p_stopping, rt, stepsize, intercept)
  
  LL_response_type <- Responsetype_log_likelihood(p_correct, p_stopping, scaled_p_omission, scaled_p_intrusion, responsetypes, 
                                                  p_acceptance = scaled_p_acceptance)
  
  TotalLL = LL_rt + LL_response_type
  
  return(TotalLL)
} 

# plot UnscaleParms(x, 0, 100) vs x, with x varzing between -inf and inf
x = seq(-10, 10, 0.1)
y = UnscaleParms(x, 0, 100)
plot(x, y, type = "l")

modelOneRun = function(p_correct, p_stopping, p_omission, p_intrusion){
  Times = mapply(rgeom, 1, p_stopping)
  Response_types = c()
  
  for(i in 1:length(p_correct)){
    final_p_correct = max(p_correct[i]/p_stopping[i],0)
    final_p_omission = max((1-p_correct[i])*p_omission/p_stopping[i],0)
    final_p_intrusion = max((1-p_correct[i])*p_intrusion/p_stopping[i],0)
    sum = final_p_correct + final_p_omission + final_p_intrusion
    response_type = sample(c("Correct", "Omission", "Intrusion"), size = 1, prob = c(final_p_correct/sum, final_p_omission/sum, final_p_intrusion/sum))
    Response_types = c(Response_types, response_type)
  }
  return(list(Times = Times, Response_types = Response_types))
}


modelSimulation = function(Participant_data, fitted_params, stepsize, ModelVersion = "Beta", repetitions = 1){
  p_omission = fitted_params[["OmissionProb"]]
  p_intrusion = fitted_params[["IntrusionProb"]]
  association_strength = fitted_params[["AssociationStrength"]]
  beta = fitted_params[["Beta"]]
  lamda = fitted_params[["Lambda"]]
  context_devaluation = fitted_params[["ContextDevaluation"]]
  p_acceptance = fitted_params[["PAcceptance"]]
  
  sim_cue_target = Participant_data$SimCueTarget
  sim_distractors = Participant_data$SimsCueDistractors
  sim_out_of_list_distractors = Participant_data$SimsCueOutOfListDistractors
  p_correct = c()
  for (i in 1:length(sim_cue_target)){
    sim_distractors_vec = as.numeric(strsplit(gsub("\\[|\\]", "", sim_distractors[i]), ",")[[1]])
    sim_out_of_list_distractors_vec = as.numeric(strsplit(gsub("\\[|\\]", "", sim_out_of_list_distractors[i]), ",")[[1]])
    if (length(sim_distractors_vec) == 0){
      sim_distractors_vec = c(0)
    }
    if (length(sim_out_of_list_distractors_vec) == 0){
      sim_out_of_list_distractors_vec = c(0)
    }
    
    
  # calculate p_correct based on the model version
    # check if the model version is one of the exceptions, alse use the base models way of calculating p_correct
    
    if (ModelVersion == "Softmax"){
      temp_p_correct = p_correct_softmax(sim_cue_target[i], association_strength, sim_distractors_vec, lamda)
    }else if (ModelVersion == "OutOfContextDists"){
      temp_p_correct = p_correct_OutOfListDistractors(sim_cue_target[i], 
                                                      association_strength, 
                                                      sim_distractors_vec, 
                                                      sim_out_of_list_distractors_vec, 
                                                      beta,
                                                      context_devaluation)
    }else{
      temp_p_correct = p_correct_beta(sim_cue_target[i], association_strength, sim_distractors_vec, beta)
    }
    
    
    p_correct = c(p_correct, min(temp_p_correct, 1))
  }
  p_incorrect = 1-p_correct
  if (ModelVersion =="CorrectAcceptance"){
    p_stopping = pmax((p_correct*p_acceptance + (1-p_correct)*(min((p_omission + p_intrusion), 1))), 0.000000000001)
  }else{
    p_stopping = pmax((p_correct + (1-p_correct)*(min((p_omission + p_intrusion), 1))), 0.000000000001)
  }
  
  ListLength = Participant_data$ListLength
  Condition = Participant_data$Condition
  Experiment = Participant_data$Experiment
  ContextSize = Participant_data$ContextSize
  ParticipantRT = Participant_data$RT
  ParticipantResponseType = Participant_data$ErrorType 
  CueTargetSimilarity = Participant_data$SimCueTarget
  
  
  # calculate the RTs and responsetzpes
  
  
  out = modelOneRun(p_correct, p_stopping, p_omission, p_intrusion)
  Times = out$Times * stepsize
  Response_types = out$Response_types
  Trials = seq(1, length(Times), 1)
  
  df = data.frame(ListLength = ListLength,
                  Condition = Condition,
                  Experiment = Experiment,
                  ContextSize = ContextSize,
                  ParticipantRT = ParticipantRT,
                  ParticipantResponseType = ParticipantResponseType,
                  ModelRT = Times,
                  ModelResponseType = Response_types,
                  CueTargetSimilarity = CueTargetSimilarity,
                  Trials = Trials)
  return(df)
}

p_correct_beta = function(sim_cue_target, association_strength, sim_distractors_vec, beta){
  p_correct = (sim_cue_target**beta + association_strength)/
                    (sim_cue_target**beta + association_strength + 
                       sum(sim_distractors_vec**beta))
  return(p_correct)
}


p_correct_softmax = function(sim_cue_target, association_strength, sim_distractors_vec, lamda){
  p_correct = softmax(c(sim_cue_target + association_strength, sim_distractors_vec), scale = lamda)[1]
  return(p_correct)
}


p_correct_OutOfListDistractors = function(sim_cue_target, association_strength, sim_distractors_vec, sim_out_of_list_distractors_vec, beta, context_devaluation = 1){
  p_correct = (sim_cue_target**beta + association_strength)/
                    (sim_cue_target**beta + association_strength + 
                       sum(sim_distractors_vec**beta) +
                       sum(context_devaluation * sim_out_of_list_distractors_vec**beta))
  return(p_correct)
}


modelSimulation_OutOfListDistractors = function(Participant_data, fitted_params, stepsize, repetitions = 1){
  p_omission = fitted_params["p_omission"]
  p_intrusion = fitted_params["p_intrusion"]
  association_strength = fitted_params["association_strength"]
  beta = fitted_params["beta"]
  context_devaluation = fitted_params["context_devaluation"]
  
  sim_cue_target = Participant_data$SimCueTarget
  sim_distractors = Participant_data$SimsCueDistractors
  sim
  p_correct = c()
  for (i in 1:length(sim_cue_target)){
    sim_distractors_vec = as.numeric(strsplit(gsub("\\[|\\]", "", sim_distractors[i]), ",")[[1]])
    if (length(sim_distractors_vec) == 0){
      sim_distractors_vec = c(0)
    }
    p_correct = c(p_correct, min((sim_cue_target[i]**beta + association_strength)/
                                   (sim_cue_target[i]**beta + association_strength + 
                                      sum(sim_distractors_vec**beta)), 1))
  }
  p_incorrect = 1-p_correct
  p_stopping = pmax((p_correct + (1-p_correct)*(min((p_omission + p_intrusion), 1))), 0.000000000001)
  
  ListLength = Participant_data$ListLength
  Condition = Participant_data$Condition
  Experiment = Participant_data$Experiment
  ContextSize = Participant_data$ContextSize
  ParticipantRT = Participant_data$RT
  ParticipantResponseType = Participant_data$ErrorType 
  CueTargetSimilarity = Participant_data$SimCueTarget
  
  
  # calculate the RTs and responsetzpes
  
  
  out = modelOneRun(p_correct, p_stopping, p_omission, p_intrusion)
  Times = out$Times * stepsize
  Response_types = out$Response_types
  
  df = data.frame(ListLength = ListLength,
                  Condition = Condition,
                  Experiment = Experiment,
                  ContextSize = ContextSize,
                  ParticipantRT = ParticipantRT,
                  ParticipantResponseType = ParticipantResponseType,
                  ModelRT = Times,
                  ModelResponseType = Response_types,
                  CueTargetSimilarity = CueTargetSimilarity)
  return(df)
}



# Function to calculate the likelihood of geometric distribution
geometric_log_likelihood_corrected <- function(p_stopping, rt, stepsize, intercept = 0) {
  
  # calculate the samplesetps as discrete values from the rts: rt/step and then round up and do minus 1 (because the geometric function starts at 0)
  
  # version no intercept
  #samplesteps = ceiling(rt/stepsize)-1
  # version with intercept
  samplesteps = pmax(1, ceiling((rt-intercept)/stepsize))-1
  
  # Now, using the same notation, the log-probability of r_t would be T_t is samplesteps
  # log p_G(T_t) - log(step) if T_t > 0
  # log p_G(T_t) - log(step+intercept) if T_t = 0
  # 
  add_ons = c()
  for (i in 1:length(rt)){
    if (samplesteps[i] == 0 ){
      add_ons = c(add_ons, log(stepsize + intercept))
    }else{
      add_ons = c(add_ons, log(stepsize))
    }
  }
  # Calculate the likelihood
  LL <- sum(dgeom(samplesteps, prob = p_stopping, log = TRUE) - add_ons)
  
  if(is.nan(LL)){
    print("something went wrong")
    for (i in 1:length(p_stopping)){
      if (is.nan(sum(dgeom(samplesteps, prob = p_stopping[i], log = TRUE)))){
        print(p_stopping[i])
      }
    }
  }
  
  return(-LL)
}

Total_neg_LL_intercept <- function(params, responsetypes, sim_cue_target, sim_distractors, rt, stepsize, intercept, listlength = 0, UseSimilarity = TRUE, UseAssociationStrength = TRUE) {
  scaled_p_omission <- logistic_function(params[1])
  scaled_p_intrusion <- logistic_function(params[2])
  
  if(UseAssociationStrength){
    associationStrength = UnscaleParms(params[3], 0, 100)
  }
  else{
    associationStrength = 0
  }
  
  if(UseSimilarity){
    scaled_beta <- UnscaleParms(params[4], -10, 10)
  }
  else{
    scaled_beta <- 0
  }
  
  p_correct = c()
  for (i in 1:length(sim_cue_target)){
    sim_distractors_vec = as.numeric(strsplit(gsub("\\[|\\]", "", sim_distractors[i]), ",")[[1]])
    if (length(sim_distractors_vec) == 0){
      sim_distractors_vec = c(0)
    }
    p_correct = c(p_correct, min((sim_cue_target[i]**scaled_beta + associationStrength)/
                                   (sim_cue_target[i]**scaled_beta + associationStrength + sum(sim_distractors_vec**scaled_beta)), 1))
  }
  
  p_stopping = pmax((p_correct + (1-p_correct)*(min((scaled_p_omission + scaled_p_intrusion), 1))), 0.000000000001)
  if (max(p_stopping) > 1){
    print("something went wrong")
    print(scaled_p_omission)
    print(scaled_p_intrusion)
    print(max(p_correct))
    print(params[3])
  }
  
  LL_rt <- geometric_log_likelihood_corrected(p_stopping, rt, stepsize, intercept)
  
  LL_response_type <- Responsetype_log_likelihood(p_correct, p_stopping, scaled_p_omission, scaled_p_intrusion, responsetypes)
  
  TotalLL = LL_rt + LL_response_type
  
  return(TotalLL) 
} 



Total_neg_LLOnlyRT <- function(params, data, stepsize, intercept = 0, listlength = 0, UseSimilarity = TRUE, UseAssociationStrength = TRUE) {
  responsetypes = data$ErrorType
  sim_cue_target = data$SimCueTarget
  sim_distractors = data$SimsCueDistractors
  rt = data$RT
  scaled_p_omission <- logistic_function(params[1])
  scaled_p_intrusion <- logistic_function(params[2])
  
  if(UseAssociationStrength){
    associationStrength = UnscaleParms(params[3], 0, 100)
  }else{
    associationStrength = 0
  }
  
  if(UseSimilarity){
    scaled_beta <- UnscaleParms(params[4], -10, 10)
  }else{
    scaled_beta <- 0
  }
  
  p_correct = c()
  for (i in 1:length(sim_cue_target)){
    sim_distractors_vec = as.numeric(strsplit(gsub("\\[|\\]", "", sim_distractors[i]), ",")[[1]])
    if (length(sim_distractors_vec) == 0){
      sim_distractors_vec = c(0)
    }
    p_correct = c(p_correct, min((sim_cue_target[i]**scaled_beta + associationStrength)/
                                   (sim_cue_target[i]**scaled_beta + associationStrength + sum(sim_distractors_vec**scaled_beta)), 1))
  }
  
  p_stopping = pmax((p_correct + (1-p_correct)*(min((scaled_p_omission + scaled_p_intrusion), 1))), 0.000000000001)
  if (max(p_stopping) > 1){
    print("something went wrong")
    print(scaled_p_omission)
    print(scaled_p_intrusion)
    print(max(p_correct))
    print(params[3])
  }
  
  LL_rt <- geometric_log_likelihood_corrected(p_stopping, rt, stepsize, intercept)
  
  return(LL_rt)
} 


Total_neg_LL_softmaxOnlyRT <- function(params, data, stepsize, intercept = 0, listlength = 0, UseSimilarity = TRUE, UseAssociationStrength = TRUE) {
  responsetypes = data$ErrorType
  sim_cue_target = data$SimCueTarget
  sim_distractors = data$SimsCueDistractors
  rt = data$RT
  scaled_p_omission <- logistic_function(params[1])
  scaled_p_intrusion <- logistic_function(params[2])
  
  if(UseAssociationStrength){
    associationStrength = UnscaleParms(params[3], 0, 100)
  }
  else{
    associationStrength = 0
  }
  
  if(UseSimilarity){
    scaled_lamda <- UnscaleParms(params[4], 0, 100)
  }
  else{
    scaled_lamda <- 0
  }
  
  p_correct = c()
  for (i in 1:length(sim_cue_target)){
    sim_distractors_vec = as.numeric(strsplit(gsub("\\[|\\]", "", sim_distractors[i]), ",")[[1]])
    if (length(sim_distractors_vec) == 0){
      sim_distractors_vec = c(0)
    }
    p_correct = c(p_correct, softmax(c(sim_cue_target[i] + associationStrength, sim_distractors_vec), scale = scaled_lamda)[1])
  }
  
  p_stopping = pmax((p_correct + (1-p_correct)*(min((scaled_p_omission + scaled_p_intrusion), 1))), 0.000000000001)
  if (max(p_stopping) > 1){
    print("something went wrong")
    print(scaled_p_omission)
    print(scaled_p_intrusion)
    print(max(p_correct))
    print(params[3])
  }
  
  LL_rt <- geometric_log_likelihood_corrected(p_stopping, rt, stepsize, intercept)
  
  return(LL_rt)
} 



Total_neg_LL_OutOfListDistractorsOnlyRT <- function(params, data, stepsize, intercept = 0, listlength = 0, UseSimilarity = TRUE, UseAssociationStrength = TRUE) {
  responsetypes = data$ErrorType
  sim_cue_target = data$SimCueTarget
  sim_distractors = data$SimsCueDistractors
  rt = data$RT
  sim_out_of_list_distractors = data$SimsCueOutOfListDistractors
  scaled_p_omission <- logistic_function(params[1])
  scaled_p_intrusion <- logistic_function(params[2])
  
  if(UseAssociationStrength){
    associationStrength = UnscaleParms(params[3], 0, 100)
  }
  else{
    associationStrength = 0
  }
  
  if(UseSimilarity){
    scaled_beta <- UnscaleParms(params[4], -10, 10)
  }
  else{
    scaled_beta <- 0
  }
  
  scaled_context_devaluation <- UnscaleParms(params[5], 0, 10)
  
  p_correct = c()
  for (i in 1:length(sim_cue_target)){
    sim_distractors_vec = as.numeric(strsplit(gsub("\\[|\\]", "", sim_distractors[i]), ",")[[1]])
    sim_out_of_list_distractors_vec = as.numeric(strsplit(gsub("\\[|\\]", "", sim_out_of_list_distractors[i]), ",")[[1]])
    if (length(sim_distractors_vec) == 0){
      sim_distractors_vec = c(0)
    }
    if (length(sim_out_of_list_distractors_vec) == 0){
      sim_out_of_list_distractors_vec = c(0)
    }
    p_correct = c(p_correct, min((sim_cue_target[i]**scaled_beta + associationStrength)/
                                   (sim_cue_target[i]**scaled_beta + 
                                      associationStrength + 
                                      sum(sim_distractors_vec**scaled_beta) +
                                      sum(scaled_context_devaluation * (sim_out_of_list_distractors_vec**scaled_beta))), 1))
  }
  
  p_stopping = pmax((p_correct + (1-p_correct)*(min((scaled_p_omission + scaled_p_intrusion), 1))), 0.000000000001)
  if (max(p_stopping) > 1){
    print("something went wrong")
    print(scaled_p_omission)
    print(scaled_p_intrusion)
    print(max(p_correct))
    print(params[3])
  }
  
  LL_rt <- geometric_log_likelihood_corrected(p_stopping, rt, stepsize, intercept)
  
  return(LL_rt)
} 



Total_neg_LL_correct_acceptanceOnlyRT <- function(params, data, stepsize, intercept = 0, listlength = 0, UseSimilarity = TRUE, UseAssociationStrength = TRUE) {
  responsetypes = data$ErrorType
  sim_cue_target = data$SimCueTarget
  sim_distractors = data$SimsCueDistractors
  rt = data$RT
  scaled_p_omission <- logistic_function(params[1])
  scaled_p_intrusion <- logistic_function(params[2])
  scaled_p_acceptance <- logistic_function(params[5])
  
  if(UseAssociationStrength){
    associationStrength = UnscaleParms(params[3], 0, 100)
  }
  else{
    associationStrength = 0
  }
  
  if(UseSimilarity){
    scaled_beta <- UnscaleParms(params[4], -10, 10)
  }
  else{
    scaled_beta <- 0
  }
  
  p_correct = c()
  for (i in 1:length(sim_cue_target)){
    sim_distractors_vec = as.numeric(strsplit(gsub("\\[|\\]", "", sim_distractors[i]), ",")[[1]])
    if (length(sim_distractors_vec) == 0){
      sim_distractors_vec = c(0)
    }
    p_correct = c(p_correct, min((sim_cue_target[i]**scaled_beta + associationStrength)/
                                   (sim_cue_target[i]**scaled_beta + associationStrength + sum(sim_distractors_vec**scaled_beta)), 1))
  }
  
  p_stopping = pmax((p_correct*scaled_p_acceptance + (1-p_correct)*(min((scaled_p_omission + scaled_p_intrusion), 1))), 0.000000000001)
  if (max(p_stopping) > 1){
    print("something went wrong")
    print(scaled_p_omission)
    print(scaled_p_intrusion)
    print(max(p_correct))
    print(params[3])
  }
  
  LL_rt <- geometric_log_likelihood_corrected(p_stopping, rt, stepsize, intercept)
  
  return(LL_rt)
} 