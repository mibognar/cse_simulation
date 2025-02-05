#!/usr/bin/env Rscript

# sim_data.R
# authors: Miklos Bognar & Marton A. Varga
# affiliations: ELTE Eotvos Lorand University
# -------------------------------------------------

library(dplyr)
library(tidyr)
library(purrr)
library(furrr)
library(future)
library(rtdists)
library(EZ2)
library(qs)
library(matrixStats)

load("no_effect.Rda")
load("small_effect.Rda")
load("large_effect.Rda")

# Configuration ---------------------------------------------------------------
condition_parameters <- list(
  no_effect = no_effect,
  small_effect = small_effect,
  large_effect = large_effect
)

participant_numbers <- c(25, 50, 100, 200, 400)
num_runs <- 1000
trial_number <- 100  # Fixed trial count
num_cores <- parallel::detectCores() - 1
registry_file <- "job_registry.qs"

# Job Registry Management -----------------------------------------------------
initialize_registry <- function() {
  expand_grid(
    effect = names(condition_parameters),
    n = participant_numbers,
    run = 1:num_runs
  ) %>%
    mutate(
      status = "pending",
      file_path = NA_character_,
      attempts = 0L,
      last_error = NA_character_,
      timestamp = Sys.time()
    )
}

load_or_create_registry <- function() {
  if(file.exists(registry_file)) {
    qs::qread(registry_file)
  } else {
    registry <- initialize_registry()
    qs::qsave(registry, registry_file)
    registry
  }
}

# Optimized Simulation Function -----------------------------------------------
simulate_data_optimized <- function(condition_parameters_data, participant_number, trial_number) {

  # Generate random slopes efficiently
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


  expanded_data <- participants %>% 
    expand_grid(
      is_congruent = 0:1, 
      prev_congruent = 0:1
    ) %>% 
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
      pc = ifelse(glob_rtc == 1, 1 - 1/(participant_number*2), glob_rtc),
      ez_params = future_pmap(list(pc = pc, vrt = glob_rtv, mrt = mean_rt), function(pc, vrt, mrt) {
        tryCatch({
          result <- EZ2::Data2EZ(Pc = pc,  VRT = vrt, MRT = mrt, s = 1)
          list(v = result$v, a = result$a, Ter = pmax(result$Ter, 0.1))
        }, error = function(e) {
            message("Error in EZ2 params for pc=", pc, " vrt=", vrt, " mrt=", mrt)
            list(v = NA_real_, a = NA_real_, Ter = NA_real_)
          })
        },
        .options = furrr_options(
          globals = c("participant_number", "EZ2", "Data2EZ"),
          packages = "EZ2",
          seed = TRUE
        ),
        .progress = TRUE
      )
    ) %>% 
    unnest_wider(ez_params)


  generate_trials_vectorized <- function(a, v, Ter, n_trials) {
    n_conditions <- length(a)
    tibble(
      participant_id = rep(diffusion_data$participant_id, each = n_trials),
      is_congruent = rep(diffusion_data$is_congruent, each = n_trials),
      prev_congruent = rep(diffusion_data$prev_congruent, each = n_trials),
      rt = rdiffusion(
        n = n_conditions * n_trials,
        a = rep(a, each = n_trials),
        v = rep(v, each = n_trials),
        t0 = rep(Ter, each = n_trials)
      )
    )
  }

  trials <- generate_trials_vectorized(
    diffusion_data$a, 
    diffusion_data$v, 
    diffusion_data$Ter, 
    trial_number
  )

  return(trials)
}

# Job Execution Framework -----------------------------------------------------
run_job <- function(effect, n, run) {
  tryCatch({
    # Check existing status
    current_status <- job_registry %>%
      filter(effect == !!effect, n == !!n, run == !!run) %>%
      pull(status)
    
    if(current_status == "completed") return(TRUE)
    
    # Generate data
    sim_data <- simulate_data_optimized(
      condition_parameters[[effect]],
      participant_number = n,
      trial_number = trial_number
    )
    
    # Save output
    output_dir <- file.path("data/simulated", effect)
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    file_name <- sprintf("%s_%04d.qs", n, run)
    file_path <- file.path(output_dir, file_name)
    qs::qsave(sim_data, file_path, preset = "fast")
    
    # Update registry
    job_registry <<- job_registry %>%
      mutate(
        status = ifelse(
          effect == !!effect & n == !!n & run == !!run,
          "completed",
          status
        ),
        file_path = ifelse(
          effect == !!effect & n == !!n & run == !!run,
          file_path,
          file_path
        ),
        attempts = attempts + 1L,
        timestamp = Sys.time()
      )
    
    TRUE
  }, error = function(e) {
    # Record error details
    job_registry <<- job_registry %>%
      mutate(
        status = ifelse(
          effect == !!effect & n == !!n & run == !!run,
          "failed",
          status
        ),
        last_error = ifelse(
          effect == !!effect & n == !!n & run == !!run,
          conditionMessage(e),
          last_error
        ),
        attempts = attempts + 1L,
        timestamp = Sys.time()
      )
    FALSE
  })
}

# Parallel Execution Controller ------------------------------------------------
execute_simulations <- function() {
  plan(multisession, workers = num_cores)
  
  # Create job batches for better error handling
  job_batches <- job_registry %>%
    filter(status %in% c("pending", "failed")) %>%
    mutate(batch = (row_number() - 1) %/% 100) %>%
    group_split(batch)
  
  for(batch in job_batches) {
    results <- future_pmap(
      batch %>% select(effect, n, run),
      ~ run_job(..1, ..2, ..3),
      .progress = TRUE,
      .options = furrr_options(seed = TRUE)
    )
    
    # Save registry after each batch
    qs::qsave(job_registry, registry_file)
  }
}

# Job Recovery and Inspection --------------------------------------------------
retry_failed_jobs <- function(max_attempts = 3) {
  job_registry <<- job_registry %>%
    mutate(
      status = ifelse(
        status == "failed" & attempts < max_attempts,
        "pending",
        status
      )
    )
  execute_simulations()
}

show_status <- function() {
  job_registry %>%
    group_by(status) %>%
    summarise(
      n = n(),
      last_run = max(timestamp),
      .groups = "drop"
    )
}

# Main Execution Flow ---------------------------------------------------------
# Initialize or load existing registry
job_registry <- load_or_create_registry()

# First run
execute_simulations()

# Check status
show_status()

# Retry failed jobs (if needed)
retry_failed_jobs()

# Final status check
show_status()

# SERIAL BENCHMARK, COMMENTED OUT ----------------------------------
# library(microbenchmark)


# # Modified timing function with warmup runs and garbage collection control
# time_simulation <- function(effect_types = names(condition_parameters),
#                            participants = 155, #that is the mean
#                            trials = 100,
#                            iterations = 10) {
#   
#   # Warmup JIT compiler
#   invisible(simulate_data_optimized(condition_parameters[[1]], 10, 10))
#   
#   # Timing benchmark
#   results <- map_dfr(effect_types, function(effect) {
#     gc()  # Force garbage collection before timing
#     
#     mb <- microbenchmark(
#       {
#         simulate_data_optimized(
#           condition_parameters[[effect]],
#           participant_number = participants,
#           trial_number = trials
#         )
#       },
#       times = iterations,
#       unit = "s"
#     )
#     
#     tibble(
#       effect = effect,
#       mean_time = mean(mb$time / 1e9),
#       min_time = min(mb$time / 1e9),
#       max_time = max(mb$time / 1e9),
#       sd_time = sd(mb$time / 1e9)
#     )
#   })
#   
#   # Estimate total duration
#   total_estimates <- results %>%
#     mutate(
#       estimated_total = mean_time * 5000  # 5k runs altogether for each effect
#     ) %>%
#     group_by(effect) %>%
#     summarise(
#       avg_time_per_run = mean(mean_time),
#       total_estimated_hours = sum(estimated_total) / 3600,
#       .groups = "drop"
#     )
#   
#   list(
#     timing_results = results,
#     total_estimates = total_estimates
#   )
# }

# # Run benchmark with 10 iterations for statistical reliability
# set.seed(42)
# benchmark_results <- time_simulation(iterations = 10)

# # Print formatted results
# cat("Individual Run Timing (seconds):\n")
# print(benchmark_results$timing_results)

# cat("\nTotal Projected Duration:\n")
# print(benchmark_results$total_estimates)
