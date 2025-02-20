format_forecasts <- function(df_weekly,
                             pred_type = "count",
                             target = "ILI ED visits") {
  df_weekly_quantiled <- df_weekly |>
    forecasttools::trajectories_to_quantiles(
      timepoint_cols = c("target_end_date"),
      value_col = {{ pred_type }},
      id_cols = c(
        "location", "reference_date",
        "horizon", "obs_data"
      )
    ) |>
    mutate(
      output_type = "quantile",
      target = {{ target }},
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

  if (pred_type == "pct") {
    df_weekly_quantiled <- df_weekly_quantiled |>
      mutate(
        value = 100 * value,
        obs_data = 100 * obs_data
      )
  }

  return(df_weekly_quantiled)
}
