##### Rscript that will generate hub formatted forecasts ####
library(argparser, quietly = TRUE)
library(RcppTOML)
library(dplyr)
library(mvgam)
library(cmdstanr)
library(lubridate)
library(purrr)

list.files(file.path("R"), full.names = TRUE) |>
  walk(source)

# Command Line Version --------------------------------------------------------
parsed_args <- arg_parser("Preprocess the data for a config") |>
  add_argument("config", help = "Path to TOML config file") |>
  add_argument("config_index",
    help = "index of entry in config to use",
    type = "integer"
  ) |>
  parse_args()

if (!is.na(parsed_args$config)) {
  # config <- parseTOML("input/example_config.toml") #nolint
  config <- parseTOML(parsed_args$config)
} else {
  message("File specified in config filepath does not exist")
}

# Large loop running around all the indices to get the files needed for
# all model fits
for (i in seq_along(config$regions_to_fit)) {
  index <- i
  ### Load in data and model ---------------------------------------------
  # forecast data
  load(file.path(
    config$input_data_path[index], config$forecast_date,
    glue::glue("{config$data_filename[index]}_forecast.rda")
  ))
  # model data
  load(file.path(
    config$input_data_path[index], config$forecast_date,
    glue::glue("{config$data_filename[index]}.rda")
  ))
  # Various filepaths to save extract the model and save the output
  fp_mod <- file.path(
    "output",
    "model",
    config$output_data_path,
    config$forecast_date,
    config$data_filename[index]
  )
  fp_figs <- file.path(
    "output",
    "figures",
    config$output_data_path,
    config$forecast_date,
    config$data_filename[index]
  )
  load(file.path(
    fp_mod,
    glue::glue("{config$model_filename[index]}.rda")
  ))


  # Use tidybayes to get a long tiy dataframe of draws
  forecast_obj <- forecast(ar_mod, newdata = forecast_data, type = "response")
  dfall <- make_long_pred_df(
    forecast_obj = forecast_obj,
    model_data = model_data,
    pred_type = config$pred_type[index],
    timestep = config$timestep_data[index]
  )

  sampled_draws <- sample(1:max(dfall$draw), 100)

  df_recent <- dfall |> filter(
    date >= ymd(config$forecast_date) - days(90)
  )

  plot_draws <- ggplot(df_recent |> filter(
    draw %in% c(sampled_draws)
  )) +
    geom_line(
      aes(
        x = date, y = count,
        group = draw,
        color = period
      ),
      alpha = 0.2,
      show.legend = FALSE
    ) +
    geom_line(aes(
      x = date,
      y = obs_data
    )) +
    coord_cartesian(xlim = ) +
    facet_wrap(~location, scales = "free_y") +
    theme_bw() +
    xlab("") +
    ylab("Incident ED visits due to ILI")
  ggsave(
    plot = plot_draws,
    filename = file.path(fp_figs, "forecast_draws.png")
  )

  if (config$timestep_data[index] != "week") {
    df_weekly <- daily_to_epiweekly_data(df_recent, config$forecast_date)
    if (!all(df_weekly$n_days_data[df_weekly$horizon >= 0 && horizon < 5] == 7)) { # nolint
      cli::cli_abort(
        message = "Not all weeks contain 7 days of data"
      )
    }
  } else {
    df_weekly <- dfall
  }

  if (config$regions_to_fit[index] == "NYC") {
    df_weekly_quantiled <- format_nyc_forecasts(df_weekly)
  } else {
    df_weekly_quantiled <- format_weekly_pct_data(df_weekly)
  }

  df_quantiles_wide <- df_weekly_quantiled |>
    filter(
      target_end_date >= (reference_date - weeks(1)),
      output_type_id %in% c(0.5, 0.025, 0.975, 0.25, 0.75)
    ) |>
    tidyr::pivot_wider(
      id_cols = c("location", "target_end_date"),
      names_from = "output_type_id"
    )



  plot_quantiles <-
    ggplot() +
    geom_line(
      data = df_weekly_quantiled |>
        filter(
          target_end_date >= reference_date - weeks(10),
          target_end_date < reference_date
        ),
      aes(x = target_end_date, y = obs_weekly_sum),
      linetype = "dashed"
    ) +
    geom_point(
      data = df_weekly_quantiled |>
        filter(
          target_end_date >= reference_date - weeks(10),
          target_end_date < reference_date
        ),
      aes(x = target_end_date, y = obs_weekly_sum)
    ) +
    facet_wrap(~location, scales = "free_y") +
    geom_line(
      data = df_quantiles_wide,
      aes(x = target_end_date, y = `0.5`)
    ) +
    geom_ribbon(
      data = df_quantiles_wide,
      aes(
        x = target_end_date,
        ymin = `0.25`,
        ymax = `0.75`
      ),
      alpha = 0.2
    ) +
    geom_ribbon(
      data = df_quantiles_wide,
      aes(
        x = target_end_date,
        ymin = `0.025`,
        ymax = `0.975`
      ),
      alpha = 0.2
    ) +
    xlab("") +
    ylab("ILI ED visits") +
    ggtitle("Dynamic GAM forecasts") +
    theme_bw()

  ggsave(
    plot = plot_quantiles,
    filename = file.path(fp_figs, "quantiled_forecasts.png")
  )

  #### Save the csvs--------------------------------------------------

  fp_forecasts <- file.path(
    "output",
    config$filepath_forecasts,
    config$forecast_date
  )
  fs::dir_create(fp_forecasts, recurse = TRUE)

  write.csv(
    df_weekly_quantiled,
    file.path(
      fp_forecasts,
      glue::glue("{config$forecast_date}-{config$team_name}-{model_name}.csv") # nolint
    )
  )
}
