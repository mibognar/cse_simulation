library(dplyr)
library(readr)
library(EZ2)
library(fdrtool)
library(ez)
library(rtdists)
library(gamlss.dist)
library(lme4)
library(lmerTest)
library(parallel)
source("cse_simul.R")
load("small_effect.Rda")

args = commandArgs(trailingOnly = T)
#effect_table = args[1]
runs = args[1]
participants = args[2]
trials_per_condition = args[3]
sd_filter = args[4]

filename = paste0("data/simulated/small_effect_",participants,"_",trials_per_condition,"_",runs,"runs", sd_filter, "filter", ".csv")
raw_filename = paste0("data/simulated/small_effect_",participants,"_",trials_per_condition,"_",runs,"runs", sd_filter, "filter_RAW", ".tsv")
print(filename)

tested_sequence <- test_sequences(small_effect, runs, participants, trials_per_condition, sd_filter)
result_table <- do_test_summary(tested_sequence)

write_tsv(tested_sequence, raw_filename)
write_csv(result_table,filename)
