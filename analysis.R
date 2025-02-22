#!/mnt/st04pool/users/usumusu/local/bin/Rscript

#SBATCH --job-name=analysis.R
#SBATCH --output=out_anal.log
#SBATCH --error=error_anal.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
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

f <- qs::qread("data/results/small_effect/25_0001.qs")

process_models <- function(f) {
  # Enhanced convergence checker with verbal explanations
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

    # Multi-optimizer validation
    allfit_metrics <- tryCatch({
      if (length(messages) > 0) {
        af <- allFit(model, verbose = FALSE)
        converged <- map_lgl(af, ~is.null(.x@optinfo$conv$lme4$messages))
        if (any(converged)) {
          fe <- map_df(af[converged], ~as_tibble(t(fixef(.x))))
          list(
            optimizers_converged = sum(converged),
            max_fe_sd = max(map_dbl(fe, sd))
          )
        } else {
          list(optimizers_converged = 0L, max_fe_sd = NA_real_)
        }
      } else {
        list(optimizers_converged = NA_integer_, max_fe_sd = NA_real_)
      }
    }, error = function(e) list(optimizers_converged = NA_integer_, max_fe_sd = NA_real_))

    # Verbal status generator
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

      TRUE ~ "Proper convergence"
    )

    p_values <- tryCatch({
      summary(model)$coefficients[,5]
    }, error = function(e) rep(NA_real_, length(fixef(model))))

    list(
      status = case_when(
        !is.null(optinfo$conv$lme4$messages) ~ "Convergence Warnings",
        optinfo$conv$opt != 0 ~ "Optimizer Failure",
        TRUE ~ "Converged"
      ),
      messages = paste(messages, collapse = "; "),
      max_gradient = grad_metrics$max_gradient,
      max_relgrad = grad_metrics$max_relgrad,
      optimizers_converged = allfit_metrics$optimizers_converged,
      max_fe_sd = allfit_metrics$max_fe_sd,
      verbal = verbal_status,
      p_values = p_values
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
        error_message = if (!is.null(current$error)) as.character(current$error) else NA_character_
      )

    model_object <- current$model

    if (inherits(model_object, "merMod")) {
      if (current$error == TRUE) {
        return(
          result %>%
            mutate(
               AIC = NA_real_, BIC = NA_real_,
               logLik = NA_real_, deviance = NA_real_,
               term = NA_character_, convergence_status = "Error",
               verbal = "Model fitting failed"
            )
        )
      }

      # Enhanced metrics collection
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

      result %>%
        bind_cols(fitness) %>%
        bind_cols(broom.mixed::tidy(model_object)) %>%
        mutate(
          p.value = convergence$p_values[term]
        )

    } else if ("ANOVA" %in% names(current)) {
      anova_terms <- as_tibble(current$ANOVA) %>%
        rename(term = Effect, p.value = p, statistic = F)

      result %>%
        bind_cols(
          tibble(
            term = anova_terms$term,
            statistic = anova_terms$statistic,
            p.value = anova_terms$p.value,
            convergence_status = "N/A",
            verbal = "Convergence checks not applicable for ANOVA"
          )
        )
    } else {
      result %>%
        mutate(
          convergence_status = "Unsupported model type",
          verbal = "Convergence diagnostics not available for this model type"
        )
    }
  })
}

result <- process_models(f)
print(result, n = 100)
write_csv(result, "result.csv")
