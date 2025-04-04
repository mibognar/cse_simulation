#!/mnt/st04pool/users/usumusu/local/bin/Rscript

#SBATCH --job-name=analysis.R
#SBATCH --output=out_anal.log
#SBATCH --error=error_anal.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=15G
#SBATCH --partition=hpc2019

# analysis.R
# authors: Miklos Bognar & Marton A. Varga
# affiliations: ELTE Eotvos Lorand University
# -------------------------------------------------------

library(tidyverse)
library(broom.mixed)
library(broom)
library(lme4)
library(lmerTest)
library(performance)
library(future.batchtools)
library(furrr)

list_all_files <- function(path) {
  list.files(path = path, full.names = TRUE, recursive = TRUE)
}

files <- list_all_files("data/results")

process_file <- function(file_path) {
  f <- qs::qread(file_path)

  result <- process_models(f)
  rm(f)
  gc()

  return(result)
}

plan(list(tweak(
  batchtools_slurm,
  template = "batchtools.slurm.tmpl",
  resources = list(
    memory = 8000,
    ncpus = 2,
    partition = "hpc2019",
    work_dir = getwd()
  )
),
multisession)
)

process_models <- function(f) {

  convergence_check <- function(model) {
    if (!inherits(model, "merMod")) return(list(status = "N/A", messages = NA_character_))

    optinfo <- model@optinfo
    messages <- optinfo$conv$lme4$messages %||% character(0)

    # Gradient analysis
    grad_metrics <- tryCatch({
      relgrad <- with(optinfo$derivs, solve(Hessian, gradient))
      list(
        max_gradient = max(abs(optinfo$derivs$gradient)),
        max_relgrad = max(abs(relgrad))
      )
    }, error = function(e) list(max_gradient = NA_real_, max_relgrad = NA_real_))

    verbal_status <- case_when(
      !is.null(optinfo$conv$lme4$messages) &&
        grad_metrics$max_relgrad < 0.001 ~
        "Convergence warning with acceptable gradient",

      !is.null(optinfo$conv$lme4$messages) &&
        grad_metrics$max_relgrad >= 0.001 ~
        "High gradient convergence warning",

      is.null(optinfo$conv$lme4$messages) &&
        grad_metrics$max_relgrad > 0.001 ~
        "Hidden convergence issues (high gradient)",

      optinfo$conv$opt != 0 ~ "Optimizer failure",

      TRUE ~ "No issue"
    )


    list(
      status = case_when(
        !is.null(optinfo$conv$lme4$messages) ~ "Convergence Warnings",
        optinfo$conv$opt != 0 ~ "Optimizer Failure",
        TRUE ~ "Converged"
      ),
      messages = paste(messages, collapse = "; "),
      max_gradient = grad_metrics$max_gradient,
      max_relgrad = grad_metrics$max_relgrad,
      verbal = verbal_status
    )
  }

  params_tibble <- as_tibble(f$params)

  imap_dfr(f$models, ~ {
    model_name <- .y
    current <- .x

    # Base template with error handling
    result <- params_tibble %>%
      mutate(
        model_id = model_name,
        error_message = if (!is.null(current$error)) current$error else NA_character_
      )

    model_object <- current$model

    null_model_id <- case_when(
      model_name == "full_glmer" ~ "full_null_glmer",
      model_name == "simple_glmer" ~ "simple_null_glmer",
      model_name == "full_lmer" ~ "full_null_lmer",
      model_name == "simple_lmer" ~ "simple_null_lmer",
      TRUE ~ NA
    )

    null_model <- f$models[[null_model_id]]$model

    null_model <- if (!is.na(null_model_id) && !is.null(f$models[[null_model_id]])) {
      f$models[[null_model_id]]$model
    } else {
      NULL
    }

    if (inherits(model_object, "merMod")) {
      if (current$error == TRUE) {
        result <- result %>%
          mutate(
            AIC = NA_real_, BIC = NA_real_,
            logLik = NA_real_, deviance = NA_real_,
            term = NA_character_, convergence_status = "Error",
            verbal = "Model fitting failed"
          )
      } else {
        convergence <- convergence_check(model_object)
        fitness <- tibble(
          AIC = AIC(model_object),
          BIC = BIC(model_object),
          logLik = as.numeric(logLik(model_object)),
          deviance = deviance(model_object, REML = FALSE),
          convergence_status = convergence$status,
          convergence_messages = convergence$messages,
          max_gradient = convergence$max_gradient,
          max_relgrad = convergence$max_relgrad,
          optimizers_converged = convergence$optimizers_converged,
          max_fe_sd = convergence$max_fe_sd,
          verbal = convergence$verbal
        )

        tidy_output <- broom.mixed::tidy(model_object)

        if (!is.null(null_model)) {
          null_r2 <- performance::r2(null_model)$R2_marginal
          alternative_r2 <- performance::r2(model_object)$R2_marginal
          cohens_f2 <- (alternative_r2 - null_r2) / (1 - alternative_r2)

          aic_diff <- AIC(null_model) - AIC(model_object)
          bic_diff <- BIC(null_model) - BIC(model_object)
          lrt_p <- anova(model_object, null_model)[["Pr(>Chisq)"]][2]

          interaction_term <- tidy_output %>%
            filter(term == "is_congruent:prev_congruent")

          interaction_p <- interaction_term$p.value
          interaction_estimate <- interaction_term$estimate

          direction_consistent <- interaction_estimate < 0

          # Evidence criteria evaluation
          evidence <- all(
            (
              (aic_diff >= 2) ||
                (bic_diff >= 6)
            ),
            lrt_p < 0.05,
            interaction_p < 0.05,
            direction_consistent
          )

          soft_evidence <- all(
            interaction_p < 0.05,
            direction_consistent
          )

        } else {
          cohens_f2 <- NA_real_
          evidence <- FALSE
          soft_evidence <- FALSE
        }


        result <- result %>%
          bind_cols(fitness) %>%
          bind_cols(tidy_output) %>%
          mutate(cohens_f2 = cohens_f2, evidence = evidence, soft_evidence = soft_evidence) %>%
          filter(term == "is_congruent:prev_congruent" | error_message == TRUE)
      }

    } else if ("ANOVA" %in% names(model_object)) {
      anova_table <- as_tibble(model_object$ANOVA)

      anova_terms <- anova_table %>%
        filter(Effect == "is_congruent:prev_congruent") %>%
        rename(term = Effect, p.value = p, statistic = `F`, eta_sq = ges)

      result <- result %>% bind_cols(
        tibble(
          AIC = NA_real_,
          BIC = NA_real_,
          logLik = NA_real_,
          deviance = NA_real_,
          convergence_status = NA_character_,
          convergence_messages = NA_character_,
          max_gradient = NA_real_,
          max_relgrad = NA_real_,
          optimizers_converged = NA_character_,
          max_fe_sd = NA_real_,
          verbal = NA_character_
        )
      ) %>%
        bind_cols(anova_terms) %>%
        mutate(cohens_f2 = eta_sq / (1 - eta_sq), evidence = p.value < 0.05 & term < 0)

    } else {
      result <- result %>%
        mutate(
          convergence_status = "Unsupported model type",
          verbal = "Convergence diagnostics not available for this model type"
        )
    }

    return(result)
  })
}

results <- future_map_dfr(files, process_file)
print(results, n = 100)
write_csv(results, "combined_results.csv")
