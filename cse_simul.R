### Read libraries
```{r}
library(tidyverse)
library(EZ2)
library(fdrtool)
library(ez)
library(rtdists)
library(gamlss.dist)
library(lme4)
library(lmerTest)
```

### Load empirical data summary

load("flanker2020.Rda")
flanker2020_summary = flanker2020_summary %>% 
  rename(glob_rtm = mean_rt,
         glob_rtv = var_rt,
         glob_rts = sd_rt,
         glob_rtc = correct_percent) %>% 
  dplyr::select(-"N")

### Create participants
source_data = data.frame(participant_id = c(1:100),
                         rt_intercept = rnorm(100, 0, 0.1),
                         rt_random_slope = rnorm(100, 0, 0.01)) %>% 
  uncount(2) %>% 
  mutate(is_congruent = rep(0:1, each=1, length.out=n())) %>% 
  uncount(2) %>% 
  mutate(prev_congruent = rep(0:1, each=1, length.out=n())) %>% 
  left_join(flanker2020_summary, by=c("is_congruent","prev_congruent")) %>% 
  mutate(mean_rt = case_when(is_congruent == 1 ~ glob_rtm+rt_intercept-rt_random_slope,
                             is_congruent == 0 ~ glob_rtm+rt_intercept+rt_random_slope))

### Create diffusion model parameters
diffusion_data <- source_data %>%
  filter(glob_rtc>.5) %>% 
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
    list(80,
         a=a,
         v=v,
         t0=Ter),
    rdiffusion)) %>% 
  unnest(Result)


