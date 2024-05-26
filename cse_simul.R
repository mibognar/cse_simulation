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
  filter(Ter>.1) %>% 
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
    group_by(participant_id, is_congruent, prev_congruent) %>% 
    summarize(participant_mean_rt = mean(rt),
              participant_sd_rt = sd(rt))
  
  testdf = testdf %>% 
    left_join(participant_summary, by = c("participant_id", "is_congruent", "prev_congruent")) %>% 
    mutate(rt_zscore = (rt-participant_mean_rt/participant_sd_rt))
  
  testfilter = testdf %>% 
    filter(response == "upper",
           abs(rt_zscore)<sd_filter)
  
  # Fit model with slope and intercept
  big_csemodel = lmer(rt ~ is_congruent*prev_congruent + (1+is_congruent|participant_id),data=testfilter, control = lmerControl(optimizer = "Nelder_Mead"))
  big_csemodel_summary = summary(big_csemodel)
  big_csemodel_estimate = big_csemodel_summary$coefficients[4]*1000
  big_csemodel_p = big_csemodel_summary$coefficients[20]
  
  # Fit model with  intercept
  small_csemodel = lmer(rt ~ is_congruent*prev_congruent + (1|participant_id),data=testfilter, control = lmerControl(optimizer = "Nelder_Mead"))
  small_csemodel_summary = summary(small_csemodel)
  small_csemodel_p = small_csemodel_summary$coefficients[20]
  
  #Fit ANOVA
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
  
  anova_p = csenova$ANOVA$p[3]
  
  return(c(big_csemodel_estimate, big_csemodel_p<.05, small_csemodel_p<.05, anova_p<.05, flag))
}

test_sequences = function(runs, participants, trials){
  
  sd2list = vector("list", runs)
  sd3list = vector("list", runs)
  nonelist = vector("list", runs)
  for(j in 1:runs){
    sd2list[[j]] = test_simulation(flanker2020_summary, participants,trials,2, "2SD filter")
    print(j)
  }
  sd2df = as.data.frame(do.call(rbind, sd2list))
  
  for(j in 1:runs){
    sd3list[[j]] = test_simulation(flanker2020_summary, participants,trials,2.5, "2 and a half SD filter")
    print(j)
  }
  sd3df = as.data.frame(do.call(rbind, sd3list))
  
  for(j in 1:runs){
    nonelist[[j]] = test_simulation(flanker2020_summary, participants,trials,3, "3SD filter")
    print(j)
  }
  nonedf = as.data.frame(do.call(rbind, nonelist))
  
  alldf = rbind(sd2df,sd3df,nonedf)
  colnames(alldf) = c("estimate","slope_p", "intercept_p", "anova_p","filtering")
  
  testsummary = alldf %>%
    group_by(filtering) %>% 
    summarize(intercept_slope_melr = mean(as.integer(as.logical(slope_p)), na.rm=T),
              intercept_melr = mean(as.integer(as.logical(intercept_p)), na.rm=T),
              anova = mean(as.integer(as.logical(anova_p)), na.rm=T))
  
  return(testsummary)
}

