#### R script to define and run the models #####
library(argparser, quietly = TRUE)
library(RcppTOML)
library(dplyr)
library(mvgam)
library(cmdstanr)
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



# Load the data---------------------------------------------------------
model_data_filename <- config$data_filename[index]
fp_data <- file.path(
  config$input_data_path[index],
  config$forecast_date
)
load(file.path(
  fp_data,
  glue::glue("{model_data_filename}.rda")
))
load(file.path(
  fp_data,
  glue::glue("{model_data_filename}_forecast.rda")
))
## Specify other filepaths for saving ---------------------------------
fp_figs <- file.path(
  "output",
  "figures",
  config$forecast_date,
  config$data_filename[index]
)
fp_mod <- fp <- file.path(
  "output",
  "model",
  config$forecast_date,
  config$data_filename[index]
)
fp_summary <- fp <- file.path(
  "output",
  "summary",
  config$forecast_date,
  config$data_filename[index]
)
fs::dir_create(fp, recurse = TRUE)
fs::dir_create(fp_mod, recurse = TRUE)
fs::dir_create(fp_summary, recurse = TRUE)


# Data source (count vs ED pct) dictates the observation model and whether time
# is indexed in days or weeks
# NYC daily count data
if (config$targets[index] == "ILI ED visits" && config$regions_to_fit[index] == "NYC" && config$timestep_data[index] == "day") { # nolint
  ##### Dynamical GAM with independent autoregression #####

  # y_{l,t} \sim Poisson(exp(x_{l,t})) \\
  # x_{l,t} \sim Normal(\mu_{l,t} + \delta_{l} x_{l,t-1},  \sigma_{process})\\
  # \mu_{l,t} = \beta_l + f_{global,t}(week) + f_{l,t}(week) + f_{global,t}(wday) \\ #nolint
  # \beta_l \sim Normal(\beta_{global}, \sigma_{count}) \\
  # \beta_{global} \sim Normal(log(avgcount), 1) \\
  # \sigma_{count} \sim exp(0.33) \\
  # \delta_l \sim Normal(0.5, 0.25) \\
  # \sigma \sim exp(1) \\

  prior_mode <- log(mean(model_data$obs_data, na.rm = TRUE))
  message("Prior for intercept (avg ED visits): ", exp(prior_mode))
  message("Recommend using ", round(prior_mode, 2), " as prior")

  ar_mod <- mvgam(
    # Observation formula, empty to only consider the Gamma observation process
    formula = obs_data ~ -1,

    # Process model formula that includes regional intercepts
    trend_formula = ~
      # Hierarchical intercepts capture variation in average count
      s(trend, bs = "re") +
        # Hierarchical effects of year(shared smooth)
        # s(year, k = 3) +
        # # Borough level deviations
        # s(year, trend, bs = "sz", k = 3) -1 +

        # Hierarchical effects of week(shared smooth)
        s(week, k = 12) +
        # Borough level deviations
        s(week, trend, bs = "sz", k = 12) - 1 +

        # Shared smooth of day of week
        s(day_of_week, k = 3),
    trend_model = "AR1",
    # Adjust the priors
    priors = c(
      prior(normal(6.69, 1), # data based prior
        class = mu_raw_trend
      ),
      prior(exponential(0.33), class = sigma_raw_trend),
      prior(exponential(1), class = sigma),
      prior(normal(0.5, 0.25), class = ar1, lb = -1, ub = 1)
    ),
    data = model_data,
    newdata = forecast_data,
    backend = "cmdstanr",
    family = poisson()
  )
  # NYC weekly count data
} else if (config$targets[index] == "ILI ED visits" && config$regions_to_fit[index] == "NYC" && config$timestep_data[index] == "week") { # nolint
  ##### Dynamical GAM with independent autoregression #####

  # y_{l,t} \sim Poisson(exp(x_{l,t})) \\
  # x_{l,t} \sim Normal(\mu_{l,t} + \delta_{l} x_{l,t-1},  \sigma_{process})\\
  # \mu_{l,t} = \beta_l + f_{global,t}(week) + f_{l,t}(week) \\
  # \beta_l \sim Normal(\beta_{global}, \sigma_{count}) \\
  # \beta_{global} \sim Normal(log(avgcount), 1) \\
  # \sigma_{count} \sim exp(0.33) \\
  # \delta_l \sim Normal(0.5, 0.25) \\
  # \sigma \sim exp(1) \\

  prior_mode <- log(mean(model_data$obs_data, na.rm = TRUE))
  message("Prior for intercept (avg ED visits): ", exp(prior_mode))
  message("Recommend using ", round(prior_mode, 2), " as prior")

  ar_mod <- mvgam(
    # Observation formula, empty to only consider the Gamma observation process
    formula = obs_data ~ -1,

    # Process model formula that includes regional intercepts
    trend_formula = ~
      # Hierarchical intercepts capture variation in average count
      s(trend, bs = "re") +
        # Hierarchical effects of year(shared smooth)
        # s(year, k = 3) +
        # # Borough level deviations
        # s(year, trend, bs = "sz", k = 3) -1 +

        # Hierarchical effects of week(shared smooth)
        s(week, k = 12) +
        # Borough level deviations
        s(week, trend, bs = "sz", k = 12) - 1,
    trend_model = "AR1",
    # Adjust the priors
    priors = c(
      prior(normal(6.69, 1), # data based prior
        class = mu_raw_trend
      ),
      prior(exponential(0.33), class = sigma_raw_trend),
      prior(exponential(1), class = sigma),
      prior(normal(0.5, 0.25), class = ar1, lb = -1, ub = 1)
    ),
    data = model_data,
    newdata = forecast_data,
    backend = "cmdstanr",
    family = poisson()
  )
} else if (config$targets[index] == "flu ED visits pct") {
  ##### Dynamical GAM with independent autoregression #####

  # y_{l,t} \sim Poisson(exp(x_{l,t})) \\
  # x_{l,t} \sim Normal(\mu_{l,t} + \delta_{l} x_{l,t-1},  \sigma_{process})\\
  # \mu_{l,t} = \beta_l + f_{global,t}(week) + f_{l,t}(week) \\
  # \beta_l \sim Normal(\beta_{global}, \sigma_{count}) \\
  # \beta_{global} \sim Normal(logit(avgprop), 1) \\
  # \sigma_{count} \sim exp(0.33) \\
  # \delta_l \sim Normal(0.5, 0.25) \\
  # \sigma \sim exp(1) \\

  ar_mod <- mvgam(
    # Observation formula, empty to only consider the Gamma observation process
    formula = obs_data ~ -1,

    # Process model formula that includes regional intercepts
    trend_formula = ~
      # Hierarchical intercepts capture variation in average count
      s(trend, bs = "re") +
        # Hierarchical effects of year(shared smooth)
        # s(year, k = 3) +
        # # Borough level deviations
        # s(year, trend, bs = "sz", k = 3) -1 +

        # Hierarchical effects of week(shared smooth)
        s(week, k = 12) +
        # Borough level deviations
        s(week, trend, bs = "sz", k = 12) - 1,
    trend_model = "AR1",
    # Adjust the priors
    priors = c(
      prior(normal(-4.5, 1),
        class = mu_raw_trend
      ),
      prior(exponential(0.33), class = sigma_raw_trend),
      prior(exponential(1), class = sigma),
      prior(normal(0.5, 0.25), class = ar1, lb = -1, ub = 1)
    ),
    data = model_data,
    newdata = forecast_data,
    backend = "cmdstanr",
    family = betar()
  )
}

#### Make a bunch of plots and save them ####----------------------------


summary <- summary(ar_mod)
save(ar_mod, file = file.path(
  fp_mod,
  glue::glue("{config$model_filename[index]}.rda")
))
save(summary, file = file.path(
  fp_summary,
  glue::glue("summary_{config$model_filename[index]}.rda")
))
week_coeffs <- plot_predictions(ar_mod,
  condition = c("week", "series"),
  points = 0.5, conf_level = 0.5
) +
  labs(y = "Counts", x = "week")
ggsave(
  plot = week_coeffs,
  filename = file.path(fp_figs, "week_coeffs.png")
)

if (config$targets[index] == "ILI ED visits" && config$regions_to_fit[index] == "NYC") { # nolint
  day_of_week <- plot_predictions(ar_mod,
    condition =
      c("week", "day_of_week", "series"),
    points = 0.5, conf_level = 0.5
  ) +
    labs(y = "Counts", x = "week")
  ggsave(
    plot = day_of_week,
    filename = file.path(fp_figs, "day_of_week.png")
  )
}


conditional_effects(ar_mod)

trace_sigma <- mcmc_plot(ar_mod,
  variable = "sigma",
  regex = TRUE,
  type = "trace"
)
ggsave(
  plot = trace_sigma,
  filename = file.path(fp_figs, "trace_sigma.png")
)

trace_ar_coeff <- mcmc_plot(ar_mod,
  variable = "ar1",
  regex = TRUE,
  type = "areas"
)
ggsave(
  plot = trace_ar_coeff,
  filename = file.path(fp_figs, "trace_ar_coeff.png")
)

slopes <- plot_slopes(ar_mod,
  variable = "week",
  condition = c("series", "series"),
  type = "link"
) +
  theme(legend.position = "none") +
  labs(y = "Log(counts)", x = "Location")
ggsave(
  plot = slopes,
  filename = file.path(fp_figs, "slopes.png")
)

# Hierarchical trend effects
trends <- plot(ar_mod, type = "smooths", trend_effects = TRUE)

# Hierarchical intercepts
intercepts <- plot(ar_mod, type = "re", trend_effects = TRUE)


example_forecast <- plot(ar_mod, type = "forecast", series = 1)
