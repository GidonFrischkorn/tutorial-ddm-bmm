#' Extract fixed effects from a fitted bmm/brms model as a tidy tibble
#'
#' Shared helper for the comparison scripts (15_cross_model_comparison.R,
#' 21_appendix_parameter_recovery.R). For each requested parameter, pulls
#' all matching rows from fixef(), applies the link-inverse transform
#' (e.g., exp for log-linked bound/ndt, plogis for logit-linked zr),
#' and returns one row per coefficient with estimate and 95% CrI.

extract_fixef <- function(fit, model_name, params, transforms = list()) {
  fe <- brms::fixef(fit)
  purrr::map_dfr(params, function(par_name) {
    rows <- grep(paste0("^", par_name, "_"), rownames(fe), value = TRUE)
    purrr::map_dfr(rows, function(r) {
      est <- fe[r, ]
      # %||%: use the parameter's transform if supplied, identity otherwise
      tfm <- transforms[[par_name]] %||% identity
      tibble::tibble(
        model = model_name,
        parameter = par_name,
        level = sub(paste0("^", par_name, "_"), "", r),
        estimate = tfm(est["Estimate"]),
        lower = tfm(est["Q2.5"]),
        upper = tfm(est["Q97.5"])
      )
    })
  })
}

# Strip the predictor prefix from coefficient names (e.g., "conditionspeed"
# -> "speed") so the levels can be used directly as axis labels
clean_levels <- function(df) {
  df |>
    dplyr::mutate(
      level = stringr::str_remove(level, "^(stim_cat|condition|freq)"),
      level = stringr::str_remove(level, "Intercept")
    )
}
