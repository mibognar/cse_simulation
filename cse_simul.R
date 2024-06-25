simulate_data = function(condition_parameters_data, participant_number, trial_number){
  library(tidyverse)
  library(EZ2)
  library(fdrtool)
  library(ez)
  library(rtdists)
  library(gamlss.dist)
  library(lme4)
  library(lmerTest)
  library(parallel)
### Create participants
source_data = data.frame(participant_id = c(1:participant_number),
                         rt_intercept = rnorm(participant_number, 0, 0.1),
                         rt_random_slope = rnorm(participant_number, 0, 0.01)) %>% 
  uncount(2) %>% 
  mutate(is_congruent = rep(0:1, each=1, length.out=n())) %>% 
  uncount(2) %>% 
  mutate(prev_congruent = rep(0:1, each=1, length.out=n())) %>% 
  left_join(condition_parameters_data, by=c("is_congruent","prev_congruent")) %>% 
  mutate(mean_rt = case_when(is_congruent == 1 ~ glob_rtm+rt_intercept-rt_random_slope,
                             is_congruent == 0 ~ glob_rtm+rt_intercept+rt_random_slope))

### Create diffusion model parameters
diffusion_data <- source_data %>%
  #filter(glob_rtc>.5) %>% 
  mutate(Result = pmap(list(
    ifelse(glob_rtc==1, 1-1/(N*2),glob_rtc),
    glob_rtv,
    mean_rt,
    s=1),
    Data2EZ),
    v = map_dbl(Result, "v"),
    a = map_dbl(Result, "a"),
    Ter = map_dbl(Result, "Ter")) %>% 
  dplyr::select(-Result)

### Generate trials based on model parameters
testdf = diffusion_data %>%
  mutate(Ter = ifelse(Ter<.1, .1, Ter)) %>% 
  mutate(Result=pmap(
    list(trial_number,
         a=a,
         v=v,
         t0=Ter),
    rdiffusion)) %>% 
  unnest(Result)

return(testdf)

}

test_simulation = function(condition_parameters_data, participant_number, trial_number, sd_filter, flag) {
  options(scipen = 999)
  options(dplyr.summarise.inform = FALSE)
  testdf = simulate_data(condition_parameters_data, participant_number, trial_number)
  
  ## Filtering data
  
  participant_summary = testdf %>% 
    mutate(correct = ifelse(response == "upper", 1, 0)) %>% 
    group_by(participant_id, is_congruent, prev_congruent) %>% 
    summarize(N = n(),
              participant_mean_rt = mean(rt),
              participant_var_rt = var(rt),
              participant_sd_rt = sd(rt),
              participant_correct_percent = mean(correct)) %>% 
    ungroup()
  
  diffusion_parameters = participant_summary %>% 
    mutate(Result = pmap(list(
      ifelse(participant_correct_percent==1, 1-1/(N*2),participant_correct_percent),
      participant_var_rt,
      participant_mean_rt,
      s=1),
      Data2EZ),
      v = map_dbl(Result, "v"),
      a = map_dbl(Result, "a"),
      Ter = map_dbl(Result, "Ter")) %>% 
    dplyr::select(-Result) %>% 
    mutate(is_congruent = as.factor(is_congruent),
           prev_congruent = as.factor(prev_congruent),
           participant_id = as.factor(participant_id))
  
  testdf = testdf %>% 
    left_join(participant_summary, by = c("participant_id", "is_congruent", "prev_congruent")) %>% 
    mutate(rt_zscore = ((rt-participant_mean_rt)/participant_sd_rt))
  
  testfilter = testdf %>% 
    filter(response == "upper",
           abs(rt_zscore)<sd_filter)
  
  #Fit model with glmer
  generalized_big_csemodel = glmer(rt ~ is_congruent*prev_congruent + (1+is_congruent*prev_congruent|participant_id), data = testfilter, family = inverse.gaussian(link = "log"))
  generalized_big_csemodel_summary = summary(generalized_big_csemodel)
  generalized_big_csemodel_p =generalized_big_csemodel_summary$coefficients[16]
  
  # Fit full linear model 
  full_csemodel = lmer(rt ~ is_congruent*prev_congruent + (1+is_congruent*prev_congruent | participant_id), data = testfilter, control = lmerControl(optimizer = "nlminbwrap"))
  full_csemodel_summary = summary(full_csemodel)
  full_csemodel_estimate = full_csemodel_summary$coefficients[4]*1000
  full_csemodel_p = full_csemodel_summary$coefficients[20]
  
  # Fit model with  intercept
  small_csemodel = lmer(rt ~ is_congruent*prev_congruent + (1|participant_id),data=testfilter, control = lmerControl(optimizer = "nlminbwrap"))
  small_csemodel_summary = summary(small_csemodel)
  small_csemodel_p = small_csemodel_summary$coefficients[20]
  
  #Fit ANOVA RT
  anova_summary = testfilter %>% 
    group_by(participant_id, is_congruent, prev_congruent) %>% 
    summarize(mean_rt = mean(rt, na.rm=T)) %>% 
    ungroup() %>% 
    group_by(participant_id) %>% 
    filter(n_distinct(is_congruent, prev_congruent)==4) %>% 
    mutate(is_congruent = as.factor(is_congruent),
           prev_congruent = as.factor(prev_congruent),
           participant_id = as.factor(participant_id))
  
  csenova <- ezANOVA(data = anova_summary,
                     dv = mean_rt,
                     wid = participant_id,
                     within = .(is_congruent, prev_congruent)
  )
  cse_table = testfilter %>% 
    group_by(prev_congruent, is_congruent) %>% 
    summarize(mean_rt = mean(rt, na.rm=T)) %>% 
    ungroup() %>% 
    pivot_wider(names_from = c("prev_congruent","is_congruent"), names_sep = "_", values_from = mean_rt)
  
  cse = (cse_table$`1_0`-cse_table$`1_1`)-(cse_table$`0_0`-cse_table$`0_1`)>0
  anova_p = csenova$ANOVA$p[3]
  
  #Fit ANOVA drift rate
  drift_anova = ezANOVA(data = diffusion_parameters,
                        dv = v,
                        wid = participant_id,
                        within = .(is_congruent, prev_congruent))
  
  drift_anova_p = drift_anova$ANOVA$p[3]
  
  return(c(full_csemodel_estimate,
           ifelse(cse,generalized_big_csemodel_p<.05, F),
           ifelse(cse,full_csemodel_p<.05,F),
           ifelse(cse,small_csemodel_p<.05,F),
           ifelse(cse,anova_p<.05,F),
           ifelse(cse, drift_anova_p<.05, F), flag))
  }

test_sequences = function(effect_table, runs, participants, trials){
  
  sd2list = vector("list", runs)
  sd3list = vector("list", runs)
  nonelist = vector("list", runs)
  for(j in 1:runs){
    sd2list[[j]] = test_simulation(effect_table, participants,trials,2.5, "2.5SD filter")
    print(j)
  }
  sd2df = as.data.frame(do.call(rbind, sd2list))
  
  for(j in 1:runs){
    sd3list[[j]] = test_simulation(effect_table, participants,trials,3, "3 SD filter")
    print(j)
  }
  sd3df = as.data.frame(do.call(rbind, sd3list))
  
  for(j in 1:runs){
    nonelist[[j]] = test_simulation(effect_table, participants,trials,1000, "No filter")
    print(j)
  }
  nonedf = as.data.frame(do.call(rbind, nonelist))
  
  alldf = rbind(sd2df,sd3df,nonedf)
  colnames(alldf) = c("estimate","gen_slope_p","slope_p", "intercept_p", "anova_p","drift_anova_p","filtering")
  
  testsummary = alldf %>%
    group_by(filtering) %>% 
    summarize(glmm_intercept_slope =mean(as.integer(as.logical(gen_slope_p)),na.rm=T),
              full_melr = mean(as.integer(as.logical(slope_p)), na.rm=T),
              intercept_melr = mean(as.integer(as.logical(intercept_p)), na.rm=T),
              anova_rt = mean(as.integer(as.logical(anova_p)), na.rm=T),
              anova_dr = mean(as.integer(as.logical(drift_anova_p)), na.rm = T))
  
  return(testsummary)
}

