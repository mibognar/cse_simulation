library(tidyverse)
library(qs)

# Function to load all result files
load_all_results <- function(results_dir) {
  list.files(results_dir, pattern = "\\.qs$", full.names = TRUE, recursive = TRUE) %>%
    map(qs::qread) %>%
    bind_rows()
}

# Load all results
results <- load_all_results("data/results")

# Analyze model fitness
analyze_fitness <- function(results) {
  results %>%
    mutate(
      glmer_converged = map_lgl(models.glmer, ~!.x$error),
      full_lmer_converged = map_lgl(models.full_lmer, ~!.x$error),
      simple_lmer_converged = map_lgl(models.simple_lmer, ~!.x$error)
    ) %>%
    group_by(params.effect_size, params.sd_filter, params.participants) %>%
    summarize(
      glmer_convergence_rate = mean(glmer_converged),
      full_lmer_convergence_rate = mean(full_lmer_converged),
      simple_lmer_convergence_rate = mean(simple_lmer_converged),
      .groups = "drop"
    )
}

# Analyze model estimates
analyze_estimates <- function(results) {
  results %>%
    filter(map_lgl(models.full_lmer, ~!.x$error)) %>%
    mutate(
      interaction_estimate = map_dbl(models.full_lmer, 
                                     ~fixef(.x$model)["is_congruent1:prev_congruent1"])
    ) %>%
    group_by(params.effect_size, params.sd_filter, params.participants) %>%
    summarize(
      mean_interaction = mean(interaction_estimate),
      sd_interaction = sd(interaction_estimate),
      .groups = "drop"
    )
}

# Analyze false positive rates
analyze_false_positives <- function(results, alpha = 0.05) {
  results %>%
    filter(params.effect_size == "no_effect") %>%
    mutate(
      glmer_p = map_dbl(models.glmer, ~summary(.x$model)$coefficients["is_congruent1:prev_congruent1", "Pr(>|z|)"]),
      full_lmer_p = map_dbl(models.full_lmer, ~summary(.x$model)$coefficients["is_congruent1:prev_congruent1", "Pr(>|t|)"]),
      simple_lmer_p = map_dbl(models.simple_lmer, ~summary(.x$model)$coefficients["is_congruent1:prev_congruent1", "Pr(>|t|)"]),
      anova_p = map_dbl(models.anova, ~.x$ANOVA$`Pr(>F)`[3])
    ) %>%
    group_by(params.sd_filter, params.participants) %>%
    summarize(
      glmer_fp_rate = mean(glmer_p < alpha, na.rm = TRUE),
      full_lmer_fp_rate = mean(full_lmer_p < alpha, na.rm = TRUE),
      simple_lmer_fp_rate = mean(simple_lmer_p < alpha, na.rm = TRUE),
      anova_fp_rate = mean(anova_p < alpha, na.rm = TRUE),
      .groups = "drop"
    )
}

# Run analyses
fitness_results <- analyze_fitness(results)
estimate_results <- analyze_estimates(results)
false_positive_results <- analyze_false_positives(results)

# Print results
print(fitness_results)
print(estimate_results)
print(false_positive_results)

# Optionally, save results
write_csv(fitness_results, "fitness_results.csv")
write_csv(estimate_results, "estimate_results.csv")
write_csv(false_positive_results, "false_positive_results.csv")
