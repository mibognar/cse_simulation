#!/mnt/st04pool/users/usumusu/local/bin/Rscript

#SBATCH --job-name=model.R
#SBATCH --output=out_model.log
#SBATCH --error=error_model.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --partition=hpc2019

# main.R
# authors: Miklos Bognar & Marton A. Varga
# affiliations: ELTE Eotvos Lorand University
# -------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(lme4)
  library(EZ2)
  library(future)
  library(future.batchtools)
  library(furrr)
  library(qs)
  library(validate)
  library(futile.logger)
})

Sys.setenv(TZ = "UTC")
flog.appender(appender.file("model_internal.log"))

# Configure SLURM cluster


plan(list(
  tweak(batchtools_slurm,
        template = "batchtools.slurm.tmpl",
        resources = list(
          memory = 16000,
          ncpus = 16,
          partition = "hpc2019",
          work_dir = getwd(),
          chunks.as.array.jobs = TRUE
        )),
  tweak(multisession, workers = 8)
))

# number of cores
# num_cores <- parallel::detectCores() - 1
# plan(multisession, workers = num_cores)

# Helper functions --------------------------------------------------------
load_precomputed_data <- function(param_set) {
  # Input validation
  check <- check_that(param_set,
    is.character(param_set$effect_size),
    is.numeric(param_set$sd_filter),
    is.numeric(param_set$participants),
    is.character(param_set$id)
  )
  if (any(failing(check))) stop("Invalid parameter set structure")

  file_path <- file.path("data/simulated", param_set$effect_size, paste0(param_set$id, ".qs"))


  qs::qread(file_path, strict = TRUE)
}

run_model <- function(formula, data, family = NULL) {
  ctrl <- if (is.null(family)) {
    lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))
  } else {
    glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))
  }

  result <- tryCatch({
    model <- if (is.null(family)) {
      lmer(formula, data, control = ctrl)
    } else {
      glmer(formula, data, family = family, control = ctrl)
    }

    model@frame <- model@frame[0, ]
    model@resp <- new("glmResp", ...)

    list(model = model, error = FALSE)
  }, error = function(e) list(model = NULL, error = TRUE, message = conditionMessage(e)))

  return(result)
}

calculate_cse <- function(data) {
  data %>%
    group_by(prev_congruent, is_congruent) %>%
    summarize(mean_rt = mean(rt, na.rm = TRUE)) %>%
    pivot_wider(
      names_from = c(prev_congruent, is_congruent),
      values_from = mean_rt
    ) %>%
    mutate(cse = (`1_0` - `1_1`) - (`0_0` - `0_1`)) %>%
    pull(cse)
}

# Checkpoint System -----------------------------------------------------------
initialize_checkpoint <- function() {
  checkpoint_file <- "simulation_checkpoint.qs"
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

  qs::qsave(checkpoint, "simulation_checkpoint.qs")
  invisible(checkpoint)
}

# Model fitting functions -----------------------------------------------------------

fit_glmer <- function(test_data) {
  run_model(
    rt ~ is_congruent*prev_congruent + (1 + is_congruent*prev_congruent | participant_id),
    test_data,
    inverse.gaussian(link = "log")
  )
}

fit_full_lmer <- function(test_data) {
  run_model(
    rt ~ is_congruent*prev_congruent + (1 + is_congruent*prev_congruent | participant_id),
    test_data
  )
}

fit_simple_lmer <- function(test_data) {
  run_model(
    rt ~ is_congruent*prev_congruent + (1 | participant_id),
    test_data
  )
}

fit_anova <- function(test_data) {
  ezANOVA(
    data = test_data,
    dv = participant_mean_rt,
    wid = participant_id,
    within = .(is_congruent, prev_congruent)
  )
}


# Simulation workflow -----------------------------------------------------
process_parameter_set <- function(param_set, checkpoint) {
  flog.info("Starting job processing for parameter set %s", param_set$id)
  options(scipen = 999)
  options(dplyr.summarise.inform = FALSE)

  if (param_set$id %in% checkpoint$completed_jobs) {
    return(invisible())
  }

  tryCatch({
    raw_data <- load_precomputed_data(param_set)
    flog.debug("Raw data loaded for %s", param_set$id)

    filtered_data <- raw_data %>%
      dtplyr::lazy_dt() %>%
      mutate(
        correct = as.integer(response == "upper")
      ) %>%
      group_by(participant_id, is_congruent, prev_congruent) %>%
      summarise(
        N = n(),
        participant_mean_rt = mean(rt),
        participant_var_rt = var(rt),
        participant_sd_rt = sd(rt),
        participant_correct_percent = mean(correct),
        .groups = "drop"
      ) %>%
      as_tibble()

    test_data <- raw_data %>%
      inner_join(filtered_data, by = c("participant_id", "is_congruent", "prev_congruent")) %>%
      mutate(
        rt_zscore = (rt - participant_mean_rt) / participant_sd_rt,
        across(c(is_congruent, prev_congruent, participant_id), as.factor)
      ) %>%
      filter(response == "upper", abs(rt_zscore) < param_set$sd_filter)

    # Fit models
    model_results <- future_map(
      list(
        glmer = fit_glmer,
        full_lmer = fit_full_lmer,
        simple_lmer = fit_simple_lmer,
        anova = fit_anova
      ),
      ~ future(.x(test_data)),
      .options = furrr_options(seed = TRUE)
    )

    results <- list(
      params = param_set,
      models = model_results
    )

    flog.info("Models fitted for %s", param_set$id)

    # Save results incrementally
    result_file <- tempfile(pattern = "results_", tmpdir = getwd(), fileext = ".qs")
    qs::qsave(
      list(params = param_set, models = results),
      result_file,
      preset = "fast"
    )

    # Atomic move to final location
    final_path <- file.path("data/results", param_set$effect_size, paste0(param_set$id, ".qs"))
    file.rename(result_file, final_path)

    # Update checkpoint
    update_checkpoint(checkpoint, param_set$id, "completed")

  }, error = function(e) {
    update_checkpoint(checkpoint, param_set$id, "failed")
    stop("Error processing ", param_set$id, ": ", e$message)
  })
}

# Parameter setup ---------------------------------------------------------
parameter_grid <- expand.grid(
  effect_size = c("no_effect", "small_effect", "large_effect"),
  sd_filter = c(2.5, 3.0, Inf),
  participants = c(25, 50, 100, 200, 400),
  df_id = 1:1000,
  stringsAsFactors = FALSE
) %>%
  mutate(
    id = paste0(participants, "_", df_id)
  )


# Submit jobs -------------------------------------------------------------
run_jobs <- function() {

  checkpoint <- initialize_checkpoint()

  # Filter unprocessed jobs
  pending_jobs <- parameter_grid %>%
    filter(!id %in% checkpoint$completed_jobs)

  # Process in optimized chunks
  results <- pending_jobs %>%
    future_map(
      ~ tryCatch(
        process_parameter_set(.x, checkpoint),
        error = function(e) {
          message("Critical error: ", e$message)
          flog.error("Error processing %s: %s", pending_jobs$id, e$message)
        }
      ),
      .options = furrr_options(
        seed = TRUE,
        scheduling = 8,  # Process 8 jobs per worker
        chunk_size = 100, # Optimized for SLURM array jobs
        globals = c("checkpoint", "pending_jobs", "process_parameter_set") # Reduce memory overhead
      )
    )

  # Final checkpoint update
  qs::qsave(checkpoint, "simulation_checkpoint.qs")
}

# Main Execution --------------------------------------------------------------
if (!interactive()) {
  run_jobs()
  message("Simulation completed successfully!")
}
