#!/mnt/st04pool/users/usumusu/local/bin/Rscript

#SBATCH --job-name=model.R
#SBATCH --output=out_model.log
#SBATCH --error=error_model.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=40G
#SBATCH --partition=hpc2019

# model.R
# authors: Miklos Bognar & Marton A. Varga
# affiliations: ELTE Eotvos Lorand University
# -------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(lme4)
  library(lmerTest)
  library(EZ2)
  library(ez)
  library(future)
  library(future.batchtools)
  library(furrr)
  library(qs)
})


# Configure SLURM cluster

plan(list(
  tweak(
    batchtools_slurm,
    template = "batchtools.slurm.tmpl",
    resources = list(
      memory = 5000,
      ncpus = 10,
      partition = "hpc2019",
      work_dir = getwd()
    )
  ),
  multisession
))


load_precomputed_data <- function(param_set) {
  file_path <- file.path("data/simulated", param_set$effect_size, paste0(param_set$id, ".qs"))
  qs::qread(file_path, strict = TRUE)
}

parse_filter_params <- function(filter_type) {
  if (filter_type == "no_filter") {
    return(list(type = "no_filter"))
  } else if (startsWith(filter_type, "sd_")) {
    threshold <- as.numeric(sub("sd_", "", filter_type))
    return(list(type = "sd", threshold = threshold))
  } else if (startsWith(filter_type, "mad_")) {
    threshold <- as.numeric(sub("mad_", "", filter_type))
    return(list(type = "mad", threshold = threshold))
  } else if (startsWith(filter_type, "time_")) {
    upper_limit <- as.numeric(sub("time_", "", filter_type))
    return(list(type = "time", lower = 0.2, upper = upper_limit))
  } else {
    stop("Unknown filter type: ", filter_type)
  }
}

ensure_complete_data <- function(data, participant_col, condition_cols) {
  # Ensure columns are factors
  data <- data %>%
    mutate(across(all_of(c(participant_col, condition_cols)), as.factor))

  # Count observations per participant per condition combination
  condition_counts <- data %>%
    group_by(across(all_of(c(participant_col, condition_cols)))) %>%
    summarise(n = n(), .groups = "drop")

  # Generate all possible combinations of participants and conditions
  all_combinations <- expand.grid(
    lapply(data[c(participant_col, condition_cols)], levels)
  ) %>%
    as_tibble()

  colnames(all_combinations) <- c(participant_col, condition_cols)

  # Identify missing combinations
  missing_combinations <- anti_join(
    all_combinations, condition_counts,
    by = c(participant_col, condition_cols)
  )

  if (nrow(missing_combinations) > 0) {
    incomplete_participants <- unique(missing_combinations[[participant_col]])
    message("Removing incomplete participants: ", paste(incomplete_participants, collapse = ", "))

    # Remove incomplete participants
    data_complete <- data %>%
      filter(!(!!sym(participant_col) %in% incomplete_participants)) %>%
      mutate(across(all_of(condition_cols), as.numeric))

  } else {
    message("Data is already complete across all conditions.")
    data_complete <- data %>%
      mutate(across(all_of(condition_cols), as.numeric))
  }

  return(data_complete)
}


run_model <- function(formula, data, family = NULL) {
  ctrl <- if (is.null(family)) {
    lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))
  } else {
    glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))
  }

  result <- tryCatch({
    model <- if (is.null(family)) {
      lmer(formula, data, control = ctrl, REML = FALSE)
    } else {
      glmer(formula, data, family = family, control = ctrl)
    }

    list(model = model, error = FALSE)

  }, error = function(e) list(model = NULL, error = TRUE, message = conditionMessage(e)))

  return(result)
}

# calculate_cse <- function(data) {
#   data %>%
#     group_by(prev_congruent, is_congruent) %>%
#     summarize(mean_rt = mean(rt, na.rm = TRUE)) %>%
#     pivot_wider(
#       names_from = c(prev_congruent, is_congruent),
#       values_from = mean_rt
#     ) %>%
#     mutate(cse = (`1_0` - `1_1`) - (`0_0` - `0_1`)) %>%
#     pull(cse)
# }

# Checkpoint System -----------------------------------------------------------
initialize_checkpoint <- function() {
  checkpoint_file <- "model_checkpoint.qs"
  if (file.exists(checkpoint_file)) {
    checkpoint <- qs::qread(checkpoint_file)
  } else {
    checkpoint <- list(
      completed_jobs = character(0),
      failed_jobs = list(),
      start_time = Sys.time(),
      last_updated = Sys.time()
    )
  }
  return(checkpoint)
}

update_checkpoint <- function(checkpoint, job_id, status) {
  checkpoint$last_updated <- Sys.time()

  if (status == "completed") {
    checkpoint$completed_jobs <- union(checkpoint$completed_jobs, job_id)
    checkpoint$failed_jobs[job_id] <- NULL
  } else if (status == "failed") {
    checkpoint$failed_jobs[[job_id]] <- list(
      timestamp = Sys.time(),
      attempts = length(checkpoint$failed_jobs[[job_id]]$attempts) + 1
    )
  }

  qs::qsave(checkpoint, "model_checkpoint.qs")
  invisible(checkpoint)
}

# Model fitting functions -----------------------------------------------------------

fit_full_glmer <- function(test_data) {
  run_model(
    diffusion_rt ~ is_congruent + is_congruent:prev_congruent + (1 + is_congruent | participant_id),
    test_data,
    inverse.gaussian(link = "log")
  )
}

fit_simple_glmer <- function(test_data) {
  run_model(
    diffusion_rt ~ is_congruent + is_congruent:prev_congruent + (1 | participant_id),
    test_data,
    inverse.gaussian(link = "log")
  )
}

fit_full_null_glmer <- function(test_data) {
  run_model(
    diffusion_rt ~ is_congruent + (1 + is_congruent | participant_id),
    test_data,
    inverse.gaussian(link = "log")
  )
}

fit_simple_null_glmer <- function(test_data) {
  run_model(
    diffusion_rt ~ is_congruent + (1 | participant_id),
    test_data,
    inverse.gaussian(link = "log")
  )
}

fit_full_lmer <- function(test_data) {
  run_model(
    diffusion_rt ~ is_congruent + is_congruent:prev_congruent + (1 + is_congruent | participant_id),
    test_data
  )
}

fit_simple_lmer <- function(test_data) {
  run_model(
    diffusion_rt ~ is_congruent + is_congruent:prev_congruent + (1 | participant_id),
    test_data
  )
}

fit_anova <- function(test_data) {
  #tryCatch({
    anova_model <- ezANOVA(
      data = test_data,
      dv = .(diffusion_rt),
      within = .(is_congruent, prev_congruent),
      wid = .(participant_id),
      detailed = TRUE
    )

    list(model = anova_model, error = FALSE)
  # },
  #   error = function(e) list(model = NULL, error = TRUE, message = conditionMessage(e))
  # )
  # anova_model <- aov(
  #   diffusion_rt ~ is_congruent * prev_congruent + Error(participant_id / (is_congruent * prev_congruent)),
  #   data = test_data
  # )
}

fit_full_null_lmer <- function(test_data) {
  run_model(
    diffusion_rt ~ is_congruent + (1 + is_congruent | participant_id),
    test_data
  )
}

fit_simple_null_lmer <- function(test_data) {
  run_model(
    diffusion_rt ~ is_congruent + (1 | participant_id),
    test_data
  )
}


# Simulation workflow -----------------------------------------------------
process_parameter_set <- function(param_set, checkpoint) {

  if (param_set$job_id %in% checkpoint$completed_jobs) {
    return(invisible())
  }

  raw_data <- load_precomputed_data(param_set)
  filter_params <- parse_filter_params(param_set$filter_type)


  filtered_data <- raw_data %>%
    mutate(
      correct = as.integer(diffusion_response == "upper")
    ) %>%
    group_by(participant_id, is_congruent, prev_congruent) %>%
    summarise(
      N = n(),
      participant_mean_rt = mean(diffusion_rt),
      participant_var_rt = var(diffusion_rt),
      participant_sd_rt = sd(diffusion_rt),
      participant_median_rt = median(diffusion_rt),
      participant_mad_rt = mad(diffusion_rt),
      .groups = "drop"
    ) %>%
    as_tibble()

  test_data <- raw_data %>%
    left_join(filtered_data, by = c("participant_id", "is_congruent", "prev_congruent")) %>%
    mutate(
      rt_zscore = (diffusion_rt - participant_mean_rt) / participant_sd_rt
    ) %>%
    filter(diffusion_response == "upper")

  if (filter_params$type == "sd") {
    test_data <- test_data %>% filter(abs(rt_zscore) < filter_params$threshold)
  } else if (filter_params$type == "mad") {
    test_data <- test_data %>%
      mutate(
        lower_bound = participant_median_rt - filter_params$threshold * participant_mad_rt,
        upper_bound = participant_median_rt + filter_params$threshold * participant_mad_rt
      ) %>%
      filter(diffusion_rt >= lower_bound & diffusion_rt <= upper_bound)
  } else if (filter_params$type == "time") {
    test_data <- test_data %>% filter(diffusion_rt >= filter_params$lower & diffusion_rt <= filter_params$upper)
  } else if (filter_params$type == "no_filter") {
    # Do nothing explicitly
  }

  test_data <- ensure_complete_data(
    data = test_data,
    participant_col = "participant_id",
    condition_cols = c("is_congruent", "prev_congruent")
  )


  rm(raw_data, filtered_data)
  gc()

  model_functions <- list(
    simple_glmer = fit_simple_glmer,
    full_glmer = fit_full_glmer,
    simple_null_glmer = fit_simple_null_glmer,
    full_null_glmer = fit_full_null_glmer,
    simple_lmer = fit_simple_lmer,
    full_lmer = fit_full_lmer,
    simple_null_lmer = fit_simple_null_lmer,
    full_null_lmer = fit_full_null_lmer,
    anova = fit_anova
  )

  # Fit models
  model_results <- future_map(
    model_functions,
    ~ .x(test_data),
    .options = furrr_options(seed = TRUE)
  )

  results <- list(
    params = param_set,
    models = model_results
  )


  output_dir <- file.path("data/results", param_set$effect_size, param_set$filter_type)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  file_path <- file.path(output_dir, paste0(param_set$id, ".qs"))

  qs::qsave(
    results,
    file_path,
    preset = "fast"
  )

  # Update checkpoint
  update_checkpoint(checkpoint, param_set$job_id, "completed")

  rm(results)
  gc()

}

# Parameter setup ---------------------------------------------------------
parameter_grid <- expand.grid(
  effect_size = c("no_effect", "small_effect", "large_effect"),
  filter_type = c(
    "no_filter",
    "sd_2.0", "sd_2.5", "sd_3.0",
    "mad_2.0", "mad_2.5", "mad_3.0",
    "time_1", "time_1.25", "time_1.5"
  ),
  participants = c(25, 50, 100, 200, 400),
  df_id = 1:1000,
  stringsAsFactors = FALSE
) %>%
  as_tibble() %>%
  mutate(
    id = sprintf("%s_%04d", participants, df_id),
    job_id = sprintf("%s_%s_%04d", participants, filter_type, df_id)
  )

run_jobs <- function() {
  checkpoint <- initialize_checkpoint()

  # Filter unprocessed jobs
  pending_jobs <- parameter_grid %>%
    filter(!job_id %in% checkpoint$completed_jobs)

  # Check if there are any pending jobs
  if (nrow(pending_jobs) == 0) {
    message("No pending jobs to process.")
    return()
  }

  # Process in optimized chunks using future_pmap
  future_map(
    .x = seq_len(nrow(pending_jobs)),
    .f = function(i) {
      param_row <- pending_jobs[i, ]
      process_parameter_set(param_row, checkpoint)
      qs::qsave(checkpoint, "model_checkpoint.qs")

      rm(param_row)
      gc()
    },
    .options = furrr_options(
      seed = TRUE,
    )
  )

}

# Main Execution --------------------------------------------------------------
run_jobs()
Sys.sleep(5)
message("Modeling completed successfully!")
