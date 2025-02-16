suppressPackageStartupMessages({
  library(tidyverse)
  library(lme4)
  library(EZ2)
  library(ez)
  library(future)
  library(future.batchtools)
  library(furrr)
  library(qs)
})


plan(multisession, workers = 2)

load_precomputed_data <- function(param_set) {
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


    list(model = model, error = FALSE)
  }, error = function(e) list(model = NULL, error = TRUE, message = conditionMessage(e)))

  return(result)
}

fit_glmer <- function(test_data) {
  run_model(
    rt ~ is_congruent * prev_congruent + (1 + is_congruent * prev_congruent | participant_id),
    test_data,
    inverse.gaussian(link = "log")
  )
}

fit_full_lmer <- function(test_data) {
  run_model(
    rt ~ is_congruent * prev_congruent + (1 + is_congruent * prev_congruent | participant_id),
    test_data
  )
}

fit_simple_lmer <- function(test_data) {
  run_model(
    rt ~ is_congruent * prev_congruent + (1 | participant_id),
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


param_set <- expand.grid(
  effect_size = c("no_effect", "small_effect", "large_effect"),
  sd_filter = c(2.5, 3.0, Inf),
  participants = c(25, 50, 100, 200, 400),
  df_id = 1:1000,
  stringsAsFactors = FALSE
) %>%
  as_tibble() %>%
  mutate(
    id = sprintf("%s_%04d", participants, df_id)
  )



process_parameter_set <- function(param_set, checkpoint) {

  if (param_set$id %in% checkpoint$completed_jobs) {
    return(invisible())
  }

  tryCatch({
    raw_data <- load_precomputed_data(param_set) %>%
      unnest(rt)


    filtered_data <- raw_data %>%
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
      left_join(filtered_data, by = c("participant_id", "is_congruent", "prev_congruent")) %>%
      mutate(
        rt_zscore = (rt - participant_mean_rt) / participant_sd_rt,
        across(c(is_congruent, prev_congruent, participant_id), as.factor)
      ) %>%
      filter(response == "upper", abs(rt_zscore) < param_set$sd_filter)

    model_functions <- list(
      glmer = fit_glmer,
      full_lmer = fit_full_lmer,
      simple_lmer = fit_simple_lmer,
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

    return(results)

  }, error = function(e) {
    stop("Error processing ", param_set$id, ": ", e$message)
  })
}

checkpoint <- NULL

for (i in seq_len(nrow(param_set))) {
  parami_set <- param_set[i, ]
  message(sprintf("Processing job %d: %s,", i, param_set$id))
  results <- process_parameter_set(parami_set, checkpoint)
  print(results)
}
