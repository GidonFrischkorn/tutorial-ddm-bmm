# Population-level estimates of a bmm fit as a tibble: one row per coefficient,
# estimate and 95% CrI, back-transformed with the function given per parameter

extract_fixef <- function(fit, model_name, params, transforms = list()) {
  fe <- brms::fixef(fit)
  purrr::map_dfr(params, function(par_name) {
    rows <- grep(paste0("^", par_name, "_"), rownames(fe), value = TRUE)
    purrr::map_dfr(rows, function(r) {
      est <- fe[r, ]
      # identity if no transform is given for this parameter
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

# Strip predictor prefixes from coefficient names ("conditionspeed" -> "speed")
clean_levels <- function(df) {
  df |>
    dplyr::mutate(
      level = stringr::str_remove(level, "^(stim_cat|condition|freq)"),
      level = stringr::str_remove(level, "Intercept")
    )
}
