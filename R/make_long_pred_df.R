make_long_pred_df <- function(forecast_obj,
                              model_data,
                              pred_type = "count",
                              timestep = "day") {
  for (i in seq_along(unique(model_data$location))) {
    matrix_preds <- forecast_obj$forecasts[[i]]
    matrix_hindcasts <- forecast_obj$hindcasts[[i]]
    df_i <- as.data.frame(matrix_preds) |>
      mutate(draw = seq(from = 1, to = nrow(matrix_preds)))
    colnames(df_i) <- c(as.character(
      seq(
        from = 1, to = ncol(matrix_preds),
        by = 1
      )
    ), "draw")
    df_i <- df_i |>
      pivot_longer(!draw,
        names_to = "t",
        values_to = {{ pred_type }}
      ) |>
      mutate(
        location = unique(model_data$location)[i],
        t = as.integer(t) + ncol(matrix_hindcasts),
        period = "forecast"
      )

    dfhind_i <- as.data.frame(matrix_hindcasts) |>
      mutate(draw = seq(
        from = 1,
        to = nrow(matrix_hindcasts)
      ))
    colnames(dfhind_i) <- c(as.character(
      seq(
        from = 1, to = ncol(matrix_hindcasts),
        by = 1
      )
    ), "draw")
    dfhind_i <- dfhind_i |>
      pivot_longer(!draw,
        names_to = "t",
        values_to = {{ pred_type }}
      ) |>
      mutate(
        location = unique(model_data$location)[i],
        t = as.integer(t),
        period = "hindcast"
      )
    if (i == 1) {
      df <- df_i
      df_hind <- dfhind_i
    } else {
      df <- bind_rows(df, df_i)
      df_hind <- bind_rows(df_hind, dfhind_i)
    }
  }

  dfall <- bind_rows(df, df_hind)
  dfall <- dfall |>
    left_join(data.frame(
      t = 1:max(dfall$t),
      date = seq(
        from = min(model_data$date),
        to = min(model_data$date) + max(dfall$t) - 1,
        by = {{ timestep }}
      )
    ), by = "t") |>
    left_join(
      model_data |>
        rename(obs_data = {{ pred_type }}) |>
        select(date, obs_data, location),
      by = c("date", "location")
    )

  return(dfall)
}
