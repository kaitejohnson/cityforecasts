#### R script to pre-process forecasting data from a config #####
library(argparser, quietly = TRUE)
library(RcppTOML)
library(dplyr)
library(lubridate)
library(tidyr)
library(ggplot2)



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
  index <- parsed_args$config_index
} else {
  message("File specified in config filepath does not exist")
}


# Use the config file to preprocess the data for the first index
# Will want to put these in functions but for now just do everything
# in a script

raw_data <- read.csv(config$data_url[index])

# Set up filepaths for saving data and figs -------------------------------
fp_figs <- file.path(
  "output",
  "figures",
  config$output_data_path,
  config$forecast_date,
  config$data_filename[index]
)
model_data_filename <- config$data_filename[index]
fp_data <- file.path(
  config$input_data_path[index],
  config$forecast_date
)
fs::dir_create(fp_figs, recurse = TRUE)
fs::dir_create(fp_data, recurse = TRUE)

# Case when data is retrospective (format is different for the `as_of` data

if (isTRUE(config$retrospective[index]) && config$targets[index] == "ILI ED visits") { # nolint
  data_formatted <- raw_data |>
    mutate(
      date = as.Date(Date, format = "%m/%d/%Y") + years(2000),
      count = as.integer(X),
      series = as.factor(Dim1Value),
      # Eventually we will want to scale these (compute z scores) but leave
      # as is for now
      year = year(date),
      week = week(date),
      day_of_week = wday(date)
    ) |>
    filter(Dim2Value == "All age groups") |>
    rename(location = Dim1Value) |>
    mutate(
      year = year - min(year) + 1,
      time = as.integer(date - min(date) + 1)
    ) |> # rescale year
    select(time, date, count, series, location, year, week, day_of_week)
} else {
  data_formatted <- raw_data |>
    mutate(
      date = ymd(target_end_date),
      count = observation,
      series = as.factor(location),
      year = year(date),
      week = week(date),
      day_of_week = wday(date)
    ) |>
    mutate(
      year = year - min(year) + 1,
      time = as.integer(date - min(date) + 1)
    ) |> # rescale year
    select(time, date, count, series, location, year, week, day_of_week)
}


if (isTRUE(config$exclude_COVID[index])) {
  data_formatted <- data_formatted |>
    filter(!date %in% seq(
      from = ymd(config$data_exclusion_period[1]),
      to = ymd(config$data_exclusion_period[2]),
      by = "day"
    ))
  message("Training data has been filtered to exclude COVID years")
}

plot_raw_data <-
  ggplot(data_formatted) +
  geom_line(aes(x = time, y = count)) +
  facet_wrap(~series, scales = "free_y") +
  theme_bw()
ggsave(
  filename = file.path(fp_figs, "raw_data.png"),
  plot = plot_raw_data
)

# Preprocessing needed for the daily data
if (config$targets[index] == "ILI ED visits" && config$regions_to_fit[index] == "NYC") { # nolint
  model_data <- data_formatted |>
    group_by(series, location) |>
    complete(date = seq(min(date), max(date), by = "day")) |>
    ungroup() |>
    # Remake the predictors to fill in the missing ones
    mutate(
      year = year(date),
      week = week(date),
      day_of_week = wday(date)
    ) |>
    mutate(
      year = year - min(year) + 1,
      time = as.integer(date - min(date) + 1)
    ) |>
    select(time, date, count, series, location, year, week, day_of_week)

  # Create daily forecast data to pass into mvgam
  forecast_data <- model_data |>
    group_by(series, location) |>
    tidyr::complete(date = seq(
      from = max(date) + days(1),
      to = max(date) +
        days(config$forecast_horizon[index]),
      by = "days"
    )) |>
    ungroup() |>
    mutate(
      year = year(date),
      week = week(date),
      day_of_week = wday(date)
    ) |>
    mutate(
      year = year - min(year) + 1,
      time = as.integer(date - min(date) + 1)
    ) |>
    filter(date > max(model_data$date))
} else { # Preprocessing for weekly data
}


# Save data -----------------------------------------------------------------
save(model_data,
  file = file.path(
    fp_data,
    glue::glue("{model_data_filename}.rda")
  )
)

save(forecast_data,
  file = file.path(
    fp_data,
    glue::glue("{model_data_filename}_forecast.rda")
  )
)
