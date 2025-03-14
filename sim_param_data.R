#!/mnt/st04pool/users/usumusu/local/bin/Rscript

#SBATCH --job-name=sim_param_data.R
#SBATCH --output=out_sim2.log
#SBATCH --error=error_sim2.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=50G
#SBATCH --partition=hpc2019

# sim_param_data.R
# authors: Miklos Bognar & Marton A. Varga
# affiliations: ELTE Eotvos Lorand University
# -------------------------------------------------

# CSE Simulation Pipeline ----------------------------------------------------

packages <- c(
  "tibble", "dplyr", "furrr", "future",
  "EZ2", "rtdists", "purrr", "tidyr",
  "readr", "lme4", "data.table", "future.apply",
  "future.batchtools", "fs"
)

loaded_pkgs <- lapply(packages, library, character.only = TRUE)

initialize_directories <- function(base_path) {
  fs::dir_create(base_path, recurse = TRUE, mode = "0775")
}

validate_parameters <- function(params) {
  required <- c("fixed", "vcov", "residual", "accuracy_model")
  if (!all(required %in% names(params))) {
    stop(sprintf("Invalid parameters for params '%s'", params))
  }
}

contrast_data <- function(empirical_data) {
  as.data.table(empirical_data)[
    , .(participant_id, rt, is_congruent, prev_congruent, correct)
  ][
     ,`:=`(
      is_congruent = as.integer(ifelse(is_congruent == 1, 1, -1)),
      prev_congruent = as.integer(ifelse(prev_congruent == 1, 1, -1))
    )
  ]
}

estimate_null_interaction <- function(empirical_data, accuracy_model) {
  cse_model <- lmer(
    rt ~ is_congruent + is_congruent:prev_congruent + (1 + is_congruent | participant_id),
    data = empirical_data,
    control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))
  )

  empirical_data$fitted_values <- fitted(cse_model)
  empirical_data$residuals <- residuals(cse_model)

  empirical_data$null_interaction <- sample(empirical_data$is_congruent * empirical_data$prev_congruent)

  interaction_effect <- fixef(cse_model)["is_congruent:prev_congruent"] * 
    (empirical_data$is_congruent * empirical_data$prev_congruent)

  empirical_data$rt_null <- empirical_data$fitted_values - interaction_effect +
    (fixef(cse_model)["is_congruent:prev_congruent"] * empirical_data$null_interaction) + 
    empirical_data$residuals

  masked_model <- lmer(
    rt_null ~ is_congruent + is_congruent:prev_congruent +
      (1 + is_congruent | participant_id),
    data = empirical_data,
    control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))
  )


  rand_eff <- VarCorr(masked_model)[["participant_id"]]
  rand_eff_cov <- as.matrix(rand_eff)

  return(list(
    fixed = fixef(masked_model),
    residual = sigma(masked_model),
    vcov = rand_eff_cov,
    accuracy_model = accuracy_model
  ))
}

estimate_parameters <- function(empirical_data) {

  cse_model <- lmer(
    rt ~ is_congruent + is_congruent:prev_congruent + (1 + is_congruent | participant_id),
    data = empirical_data,
    control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))
  )

  rand_eff <- VarCorr(cse_model)[["participant_id"]]
  rand_eff_cov <- as.matrix(rand_eff)

  accuracy_model <- glmer(
    correct ~ is_congruent + is_congruent:prev_congruent + (1 + is_congruent | participant_id),
    family = binomial,
    control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)),
    data = empirical_data
  )

  list(
    fixed = fixef(cse_model),
    residual = sigma(cse_model),
    vcov = rand_eff_cov,
    accuracy_model = accuracy_model
  )
}

simulate_cse_responses <- function(n_participants, n_trials, params) {

  validate_parameters(params)

  simulate_participant <- function(p) {

    participant_re <- MASS::mvrnorm(
      n = 1,
      mu = rep(0, nrow(params$vcov)),
      Sigma = params$vcov
    )

    rand_eff_acc <- VarCorr(params$accuracy_model)[["participant_id"]]
    rand_eff_cov_acc <- as.matrix(rand_eff_acc)

    participant_re_acc <- MASS::mvrnorm(
      n = 1,
      mu = rep(0, nrow(rand_eff_cov_acc)),
      Sigma = rand_eff_cov_acc
    )

    trials <- tibble(
      participant_id = as.character(p),
      trial = 1:n_trials,
      is_congruent = sample(c(-1, 1), n_trials, replace = TRUE),
      prev_congruent = lag(is_congruent, default = 1)
    ) %>%
      mutate(
        fixed_effect = as.numeric(
          params$fixed["(Intercept)"] +
          params$fixed["is_congruent"] * is_congruent +
          params$fixed["is_congruent:prev_congruent"] * is_congruent * prev_congruent
        ),
        random_effect = as.numeric(
          participant_re["(Intercept)"] +
          participant_re["is_congruent"] * is_congruent
        ),

        rt = fixed_effect + random_effect + rnorm(n_trials, 0, params$residual)
      )

    x_fixed <- model.matrix(~ is_congruent + is_congruent:prev_congruent, data = trials)

    x_random <- model.matrix(~ is_congruent, data = trials)

    fixed_logit <- x_fixed %*% fixef(params$accuracy_model)

    random_logit <- x_random %*% participant_re_acc

    logit <- fixed_logit + random_logit

    trials$correct <- rbinom(n_trials, 1, plogis(logit))

    return(trials)
  }

  results <- future_map_dfr(
    1:n_participants,
    simulate_participant,
    .options = furrr_options(seed = TRUE)
  )

  diffusion_data <- results %>%
    group_by(participant_id, is_congruent, prev_congruent) %>%
    summarise(
      pc = mean(correct, na.rm = TRUE),
      vrt = var(rt[correct == 1] / 1000, na.rm = TRUE),
      mrt = mean(rt[correct == 1] / 1000, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      pc = case_when(
        pc == 1 ~ 1 - 1 / (n_trials + 1),
        pc == 0 ~ 1 / (n_trials + 1),
        TRUE ~ pc
      ),
      vrt = ifelse(is.na(vrt), var(results$rt), vrt),
      mrt = ifelse(is.na(mrt), mean(results$rt), mrt)
    ) %>%
    mutate(
      ez_params = future_pmap(
        list(pc, vrt, mrt),
        function(p, v, m) {
          tryCatch({
            res <- EZ2::Data2EZ(Pc = p, VRT = v, MRT = m, s = 1)
            list(
              v = res$v,
              a = res$a,
              Ter = pmax(res$Ter, 0.01)
            )
          }, error = function(e) {
            message("EZ2 Error: ", e$message, " [Participant: ", pick()$participant_id, "]")
            list(v = NA_real_, a = NA_real_, Ter = NA_real_)
          })
        },
        .options = furrr_options(globals = "EZ2", packages = "EZ2", seed = TRUE)
      )
    ) %>%
    unnest_wider(ez_params)


  final_trials <- results %>%
    left_join(diffusion_data, by = c("participant_id", "is_congruent", "prev_congruent")) %>%
    filter(!is.na(v), !is.na(a), !is.na(Ter)) %>% # remove any problematic rows
    mutate(
      diffusion = future_pmap(
        list(a, v, Ter),
        ~ rdiffusion(n = 1, a = ..1, v = ..2, t0 = ..3),
        .options = furrr_options(seed = TRUE)
      )
    ) %>%
    unnest_wider(diffusion, names_sep = "_") %>%
    select(participant_id, trial, is_congruent, prev_congruent, diffusion_rt, diffusion_response)

  return(final_trials)
}

save_results <- function(data, effect_name, n_participants, run_number) {
  output_dir <- file.path("data/simulated", effect_name)
  initialize_directories(output_dir)

  file_path <- file.path(output_dir, sprintf("%s_%04d.qs", n_participants, run_number))
  qs::qsave(data, file_path, preset = "fast")

  return(file_path)
}

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

load_or_create_registry <- function(registry_file) {
  if (file.exists(registry_file)) {
    qs::qread(registry_file)
  } else {
    registry <<- initialize_registry()
    qs::qsave(registry, registry_file)
    registry
  }
}


update_registry <- function(effect, n, run, status, file_path = NA_character_, error_msg = NA_character_) {
  job_registry <<- job_registry %>%
    mutate(
      status = if_else(effect == !!effect & n == !!n & run == !!run, status, status),
      file_path = if_else(effect == !!effect & n == !!n & run == !!run & !is.na(file_path), file_path, file_path),
      last_error = if_else(effect == !!effect & n == !!n & run == !!run & !is.na(error_msg), error_msg, last_error),
      attempts = if_else(effect == !!effect & n == !!n & run == !!run, attempts + 1L, attempts),
      timestamp = if_else(effect == !!effect & n == !!n & run == !!run, Sys.time(), timestamp)
    )

  qs::qsave(job_registry, registry_file)
}

execute_simulations <- function(condition_parameters) {
  plan(multisession, workers = availableCores() - 2)
  # Configure SLURM cluster
  plan(list(
     tweak(
       batchtools_slurm,
       template = "batchtools.slurm.tmpl",
       resources = list(
         memory = 6000,
         ncpus = 10,
         partition = "hpc2019",
         work_dir = getwd()
       )
     ),
     multisession
  ))

  pending_jobs <-job_registry %>%
    filter(status == "pending") %>%
    select(effect, n, run)

  if (nrow(pending_jobs) == 0) {
    message("No pending jobs found.")
    return(NULL)
  }

  future_pwalk(
    .l = pending_jobs,
    .f = function(effect, n, run, ...) {
      tryCatch({
        message("Started simulation for effect: ", effect, ", participants: ", n, ", run: ", run)
        sim_data <- simulate_cse_responses(n, trial_number, condition_parameters[[effect]])
        path <- save_results(sim_data, effect, n, run)
        message("Saved results to:", path)
        update_registry(effect, n, run, "completed", file_path = path)
        message("Updated registry for effect: ", effect)
      }, error = function(e) {
        update_registry(effect, n, run, "failed", error_msg = e$message)
        message("Error encountered: ", e$message)
      })
    },
    .options = furrr_options(
      seed = TRUE
    )
  )
}


system.time({
  small_effect <- read_csv("./data/empirical/flanker_processed.csv") |>
    contrast_data()

  large_effect <- read_csv("./data/empirical/primeprobe_processed.csv") |>
    contrast_data()

  condition_parameters <- list(
    small_effect = estimate_parameters(small_effect),
    large_effect = estimate_parameters(large_effect)
  )

  condition_parameters$no_effect <- estimate_null_interaction(
    small_effect,
    condition_parameters$small_effect$accuracy_model
  )

  participant_numbers <- c(25, 50, 100, 200, 400)
  num_runs <- 1000
  trial_number <- 100  # Fixed trial count
  registry_file <<- "job_registry.qs"

  job_registry <<- load_or_create_registry(registry_file)

  execute_simulations(condition_parameters)

})
