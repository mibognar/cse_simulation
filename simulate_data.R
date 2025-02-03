library(dplyr)
library(tidyr)
library(purrr)
library(furrr)
library(future)
library(rtdists)
library(EZ2)

# Loading data parameters

load("large_effect.Rda")
load("small_effect.Rda")
load("no_effect.Rda")

set.seed(42)  # For reproducibility

simulate_data_optimized <- function(condition_parameters_data, participant_number, trial_number) {

  # Generate random slopes more efficiently
  random_slopes <- MASS::mvrnorm(
    n = participant_number,
    mu = c(0, 0),
    Sigma = matrix(c(0.02^2, 0.8 * 0.02 * 0.005, 0.8 * 0.02 * 0.005, 0.005^2), nrow = 2)
  ) %>%
    as.data.frame() %>%
    setNames(c("congruency_random_slope", "interaction_random_slope"))

  # Create base participant data
  participants <- tibble(
    participant_id = 1:participant_number,
    rt_intercept = rnorm(participant_number, 0, 0.2)
  ) %>%
    bind_cols(random_slopes)

  # Expand participants with conditions
  expanded_data <- participants %>%
    crossing(is_congruent = 0:1, prev_congruent = 0:1) %>%
    left_join(condition_parameters_data, by = c("is_congruent", "prev_congruent"))

  # Compute random slopes and mean RT
  processed_data <- expanded_data %>%
    mutate(
      congruency_random_slope = ifelse(is_congruent == 1, congruency_random_slope, -congruency_random_slope),
      interaction_random_slope = case_when(
        is_congruent == 1 & prev_congruent == 1 ~ interaction_random_slope,
        is_congruent == 1 & prev_congruent == 0 ~ -interaction_random_slope,
        is_congruent == 0 & prev_congruent == 1 ~ -interaction_random_slope,
        is_congruent == 0 & prev_congruent == 0 ~ interaction_random_slope,
        TRUE ~ 0
      ),
      mean_rt = glob_rtm + rt_intercept + congruency_random_slope + interaction_random_slope
    )

# Generate diffusion model parameters using EZ2
  diffusion_data <- processed_data %>%
    mutate(
      pc = ifelse(glob_rtc == 1, 1 - 1 / (participant_number * 2), glob_rtc),
      Result = pmap(list(pc, glob_rtv, mean_rt, s = 1), Data2EZ),
      v = map_dbl(Result, "v"),
      a = map_dbl(Result, "a"),
      Ter = map_dbl(Result, "Ter"),
      Ter = pmax(Ter, 0.1)
    ) %>%
    select(-Result, -pc)

  # Generate trials using vectorized operations
  generate_trials <- function(data) {
    n <- nrow(data)
    tibble(
      participant_id = rep(data$participant_id, each = trial_number),
      is_congruent = rep(data$is_congruent, each = trial_number),
      prev_congruent = rep(data$prev_congruent, each = trial_number),
      rt = rdiffusion(
        n = n * trial_number,
        a = rep(data$a, each = trial_number),
        v = rep(data$v, each = trial_number),
        t0 = rep(data$Ter, each = trial_number)
      )
    )
  }

  # Apply trial generation and combine results
  trials <- diffusion_data %>%
    group_by(participant_id, is_congruent, prev_congruent) %>%
    group_modify(~ generate_trials(.x)) %>%
    ungroup()

  return(trials)
}


# Function to run the simulation and measure time

# run_simulation <- function() {
#   start_time <- Sys.time()
#   result <- simulate_data_optimized(large_effect, 10, 100)
#   end_time <- Sys.time()
#   return(as.numeric(end_time - start_time, units = "secs"))
# }

# Run the simulation 10 times and calculate average runtime

# n_runs <- 10
# runtimes <- numeric(n_runs)

# for (i in 1:n_runs) {
#   cat("Run", i, "of", n_runs, "\n")
#   runtimes[i] <- run_simulation()
# }

# average_runtime <- mean(runtimes)
# cat("\nAverage runtime over", n_runs, "runs:", round(average_runtime, 2), "seconds\n")

# # Optional: Display all runtimes
# print(data.frame(Run = 1:n_runs, Runtime = round(runtimes, 2)))

# original_data <- readr::read_csv("~/Downloads/df1.csv")
# new_data <- simulate_data_optimized(large_effect, 10, 100)

# debrief <- data.frame(
#   original_mean_rt = mean(original_data$rt),
#   new_mean_rt = mean(new_data$rt$rt),
#   original_sd_rt = sd(original_data$rt),
#   new_sd_rt = sd(new_data$rt$rt)
# )

# print(debrief)


### Save out 3x5x1000 Rds objects (effect-size * participant_number * 1000)

condition_parameters <- list(
  no_effect = no_effect,
  small_effect = small_effect,
  large_effect = large_effect
)

participant_numbers <- c(25, 50, 100, 200, 400)
num_runs <- 1000

num_cores <- parallel::detectCores() - 1
plan(multisession, workers = num_cores)

run_and_save <- function(effect, n, run) {
  effect_data <- condition_parameters[[effect]]

  sim_data <- simulate_data_optimized(
    condition_parameters_data = effect_data,
    participant_number = n,
    trial_number = 100
  )
  
  output_dir <- file.path("data/simulated", effect)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  file_name <- sprintf("%s_%03d.rds", n, run)
  saveRDS(sim_data, file.path(output_dir, file_name))

}


simulation_params <- expand_grid(
  effect = names(condition_parameters),
  n = participant_numbers,
  run = 1:num_runs
)

results <- future_pmap(
  simulation_params,
  ~ tryCatch(
    run_and_save(..1, ..2, ..3),
    error = function(e) {
      message("Error in ", ..1, "-", ..2, "-", ..3)
      return(NULL)
    }
  ),
  .progress = TRUE,
  .options = furrr_options(seed = TRUE)
)

# Print summary

cat("Total simulations run:", length(results), "\n")
cat("First few file names: \n")
print(head(results))

# Serial alternative:

# for (effect in names(condition_parameters)) {
#   for (n in participant_numbers) {
#     for (run in 1:num_runs) {
#       sim_data <- simulate_data_optimized(
#         condition_parameters[[effect]],
#         participant_number = n,
#         trial_number = 100
#       )
#       
#       # Save the file
#       saveRDS(
#         sim_data,
#         file = paste0("data/simulated/", effect, "/", n, "_participants_run", run, ".Rds")
#       )
#     }
#   }
# }





