library(broom.mixed)
library(tidyverse)
library(qs)
library(purrr)

small_emp <- qs::qread("data/small_effect_empirical_params.qs")
small_null <- qs::qread("data/small_effect_null_params.qs")

large_emp <- qs::qread("data/large_effect_empirical_params.qs")
large_null <- qs::qread("data/large_effect_null_params.qs")

param_list <- list(
  small_emp = small_emp,
  small_null = small_null,
  large_emp = large_emp,
  large_null = large_null
)

params_to_tibble <- function(params) {
  # Fixed effects
  fixed_tbl <- tibble(
    parameter_type = "fixed_effect",
    name = names(params$fixed),
    value = as.numeric(params$fixed)
  )

  # Residual SD
  residual_tbl <- tibble(
    parameter_type = "residual_sd",
    name = "sigma",
    value = params$residual
  )

  # Random effects covariance matrix
  vcov_tbl <- as.data.frame(as.table(params$vcov)) %>%
    filter(Var1 != "" & Var2 != "") %>%
    transmute(
      parameter_type = "random_effect_covariance",
      name = paste(Var1, Var2, sep = "_"),
      value = Freq
    ) %>%
    as_tibble()

  # Combine all into one tibble
  bind_rows(fixed_tbl, residual_tbl, vcov_tbl)
}


result <- imap_dfr(param_list,
  ~ {
    df <- params_to_tibble(.x)
    df$id <- .y
    return(df)
  }
)

write_csv(result, "empirical_model_parameters.csv")
