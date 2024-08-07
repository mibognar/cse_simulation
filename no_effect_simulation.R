library(tidyverse)
library(EZ2)
library(fdrtool)
library(ez)
library(rtdists)
library(gamlss.dist)
library(lme4)
library(lmerTest)
library(parallel)
source("cse_simul.R")
load("no_effect.Rda")

args = commandArgs(trailingOnly = T)
#effect_table = args[1]
runs = args[1]
participants = args[2]
trials_per_condition = args[3]
sd_filter = args[4]

filename = paste0("data/simulated/no_effect_",participants,"_",trials_per_condition,"_",runs,"runs", sd_filter,"filter", ".csv")
print(filename)

tested_sequence <- test_sequences(no_effect, runs, participants, trials_per_condition, sd_filter)
result_table <- do_test_summary(tested_sequence)

write_csv(result_table,filename)
