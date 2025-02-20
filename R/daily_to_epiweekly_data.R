daily_to_epiweekly_data <- function(dfall, forecast_date) {
  df_forecasts <- dfall |>
    mutate(
      epiweek = epiweek(date),
      year = year(date),
      reference_date = ymd(forecast_date) +
        (7 - wday(ymd(forecast_date), week_start = 7)),
      target_end_date = ymd(date) + (7 - wday(date, week_start = 7))
    ) |>
    arrange(date)

  df_weekly <- df_forecasts |>
    group_by(
      reference_date,
      target_end_date,
      location,
      draw
    ) |>
    summarize(
      n_days_data = n(),
      count = sum(count, na.rm = TRUE),
      obs_data = sum(obs_data)
    ) |>
    ungroup() |>
    dplyr::mutate(
      horizon = floor(as.integer(target_end_date - reference_date)) / 7
    )

  return(df_weekly)
}
