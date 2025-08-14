#!/mnt/st04pool/users/usumusu/local/bin/Rscript

#SBATCH --job-name=raw.R
#SBATCH --output=out_raw.log
#SBATCH --error=error_raw.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=10G
#SBATCH --partition=hpc2019

# raw.R
# extract raw cse values for every condition
# authors: Miklos Bognar & Marton A. Varga
# affiliations: ELTE Eotvos Lorand University
# -------------------------------------------------

library(tidyverse)
library(future)
library(future.batchtools)
library(qs)
library(furrr)
library(fs)

plan(list(
  tweak(
    batchtools_slurm,
    template = "batchtools.slurm.tmpl",
    resources = list(
      memory = 500,
      ncpus = 1,
      partition = "hpc2019",
      work_dir = getwd()
    )
  ),
  multisession
))

ensure_complete_data <- function(data, participant_col, condition_cols) {
  # Ensure columns are factors
  data <- data %>%
    dplyr::mutate(dplyr::across(tidyselect::all_of(c(participant_col, condition_cols)), as.factor))

  # Count observations per participant per condition combination
  condition_counts <- data %>%
    dplyr::group_by(dplyr::across(tidyselect::all_of(c(participant_col, condition_cols)))) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop")

  # Generate all possible combinations of participants and conditions
  all_combinations <- expand.grid(
    lapply(data[c(participant_col, condition_cols)], levels)
  ) %>%
    tibble::as_tibble()

  colnames(all_combinations) <- c(participant_col, condition_cols)

  # Identify missing combinations
  missing_combinations <- dplyr::anti_join(
    all_combinations, condition_counts,
    by = c(participant_col, condition_cols)
  )

  if (nrow(missing_combinations) > 0) {
    incomplete_participants <- unique(missing_combinations[[participant_col]])
    message("Removing incomplete participants: ", paste(incomplete_participants, collapse = ", "))

    # Remove incomplete participants
    data_complete <- data %>%
      dplyr::filter(!(!!rlang::sym(participant_col) %in% incomplete_participants)) %>%
      dplyr::mutate(dplyr::across(tidyselect::all_of(condition_cols), as.numeric)) # Switching back to numeric

  } else {
    message("Data is already complete across all conditions.")
    data_complete <- data %>%
      dplyr::mutate(dplyr::across(tidyselect::all_of(condition_cols), as.numeric)) # Switching back to numeric
  }

  return(data_complete)
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

process_job <- function(effect_size, filter_type, participants, df_id, job_id, ...) {
  tryCatch({
    # Construct file path from parameters
    file_name <- sprintf("%04d", df_id)
    data_path <- path("data/simulated", effect_size,
                      paste0(participants, "_", file_name), ext = "qs")

    # Load and process data
    raw_data <- qread(data_path)

    # Apply filtering pipeline
    filter_params <- parse_filter_params(filter_type)

    filtered_data <- raw_data %>%
      mutate(correct = as.integer(diffusion_response == "upper")) %>%
      group_by(participant_id, is_congruent, prev_congruent) %>%
      summarise(
        participant_mean_rt = mean(diffusion_rt),
        participant_sd_rt = sd(diffusion_rt),
        participant_median_rt = median(diffusion_rt),
        participant_mad_rt = mad(diffusion_rt),
        .groups = "drop"
      )

    test_data <- raw_data %>%
      left_join(filtered_data, by = c("participant_id", "is_congruent", "prev_congruent")) %>%
      mutate(rt_zscore = (diffusion_rt - participant_mean_rt)/participant_sd_rt) %>%
      filter(diffusion_response == "upper", trial != 1)

    # Apply specific filter
    test_data <- switch(filter_params$type,
      "sd" = filter(test_data, abs(rt_zscore) < filter_params$threshold),
      "mad" = filter(test_data,
        between(diffusion_rt,
          participant_median_rt - filter_params$threshold * participant_mad_rt,
          participant_median_rt + filter_params$threshold * participant_mad_rt)),
      "time" = filter(test_data,
        between(diffusion_rt, filter_params$lower, filter_params$upper)),
      "no_filter" = test_data
    )

    # Ensure complete data and calculate CSE
    test_data <- ensure_complete_data(
      test_data %>% mutate(diffusion_rt = diffusion_rt * 1000),
      "participant_id",
      c("is_congruent", "prev_congruent")
    )

     test_data <- test_data %>%
      mutate(
        prev_type = if_else(prev_congruent == 1, "c", "i"),
        curr_type = if_else(is_congruent == 1, "C", "I")
      ) %>%
      group_by(prev_type, curr_type) %>%
      summarize(mean_rt = mean(diffusion_rt, na.rm = TRUE), .groups = "drop") %>%
      pivot_wider(
        names_from = c(prev_type, curr_type),
        values_from = mean_rt,
        names_sep = ""  # Creates columns like cI, cC, iI, iC
      )

    cse <- test_data %>%
      mutate(
        cse = (cI - cC) - (iI - iC),
        .keep = "none"
      ) %>%
      pull(cse)

    # Return formatted results
    tibble(
      job_id = job_id,
      effect_size = effect_size,
      filter_type = filter_type,
      participants = participants,
      df_id = df_id,
      cse = cse
    )
  }, error = function(e) {
    message("Error processing job ", job_id, ": ", e$message)
    tibble(
      job_id = job_id,
      effect_size = effect_size,
      filter_type = filter_type,
      participants = participants,
      df_id = df_id,
      cse = NA_real_
    )
  })
}

parameter_grid <- expand.grid(
  effect_size = c("small_no_effect", "large_no_effect", "small_effect", "large_effect"),
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
  tibble::as_tibble() %>%
  dplyr::mutate(
    id = sprintf("%s_%04d", participants, df_id),
    job_id = sprintf("%s_%s_%04d", participants, filter_type, df_id)
  )

# Process all combinations
results <- parameter_grid %>%
  select(-id) %>%
  future_pmap_dfr(process_job, .progress = TRUE)

# Save output
write_csv(results, "all_cse_results.csv")
