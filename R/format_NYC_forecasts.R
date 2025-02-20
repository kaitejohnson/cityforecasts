format_nyc_forecasts <- function(df_weekly) {
  df_weekly_quantiled <- df_weekly |>
    forecasttools::trajectories_to_quantiles(
      timepoint_cols = c("target_end_date"),
      value_col = "count",
      id_cols = c(
        "location", "reference_date",
        "horizon", "obs_data"
      )
    ) |>
    mutate(
      output_type = "quantile",
      target = "ILI ED visits",
      location = ifelse(location == "Citywide", "NYC", location)
    ) |>
    rename(
      output_type_id = quantile_level,
      value = quantile_value
    ) |>
    select(
      reference_date, location, horizon, obs_data,
      target, target_end_date,
      output_type, output_type_id, value
    )

  return(df_weekly_quantiled)
}
