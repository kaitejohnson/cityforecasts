#### R script to define and run the models #####
library(argparser, quietly = TRUE)
library(RcppTOML)
library(dplyr)
library(mvgam)
library(cmdstanr)

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



# Load the data
load(file.path(
  config$input_data_path[index], config$forecast_date,
  glue::glue("{config$data_filename[index]}.rda")
))
load(file.path(
  config$input_data_path[index], config$forecast_date,
  glue::glue("{config$data_filename[index]}_forecast.rda")
))


# Data source dictates the observation model and whether time is indexed
# in days or weeks
if (config$targets[index] == "ILI ED visits" && config$regions_to_fit[index] == "NYC") { # nolint
  ##### Dynamical GAM with independent autoregression #####

  # y_{l,t} \sim Poisson(exp(x_{l,t})) \\
  # x_{l,t} \sim Normal(\mu_{l,t} + \delta_{l} x_{l,t-1},  \sigma_{process})\\
  # \mu_{l,t} = \beta_l + f_{global,t}(week) + f_{l,t}(week) + f_{global,t}(wday) \\ #nolint
  # \beta_l \sim Normal(\beta_{global}, \sigma_{count}) \\
  # \beta_{global} \sim Normal(log(avgcount), 1) \\
  # \sigma_{count} \sim exp(0.33) \\
  # \delta_l \sim Normal(0.5, 0.25) \\
  # \sigma \sim exp(1) \\

  ar_mod <- mvgam(
    # Observation formula, empty to only consider the Gamma observation process
    formula = count ~ -1,

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
      prior(normal(log(mean(model_data$count, na.rm = TRUE)), 1),
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
  # \beta_{global} \sim Normal(log(avgcount), 1) \\
  # \sigma_{count} \sim exp(0.33) \\
  # \delta_l \sim Normal(0.5, 0.25) \\
  # \sigma \sim exp(1) \\

  ar_mod <- mvgam(
    # Observation formula, empty to only consider the Gamma observation process
    formula = prop_visits ~ -1,

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
      prior(normal(mean(qlogis(model_data$prop_visits, na.rm = TRUE)), 1),
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

#### Make a bunch of plots and save them ####
fp <- file.path(
  "output",
  "figures",
  config$output_data_path,
  config$forecast_date,
  config$data_filename[index]
)
fs::dir_create(fp, recurse = TRUE)
summary <- summary(ar_mod)
save(ar_mod, file = glue::glue(fp, "ar_mod.rda"))
save(summary, file = glue::glue(fp, "summary_ar_mod.rda"))
week_coeffs <- plot_predictions(ar_mod,
  condition = c("week", "series"),
  points = 0.5, conf_level = 0.5
) +
  labs(y = "Counts", x = "week")
ggsave(
  week_coeffs,
  file.path(fp, "week_coeffs.png")
)

if (config$targets[index] == "ILI ED visits" && config$regions_to_fit[index] == "NYC") { # nolint
  day_of_week <- plot_predictions(ar_mod,
    condition =
      c("week", "day_of_week", "series"),
    points = 0.5, conf_level = 0.5
  ) +
    labs(y = "Counts", x = "week")
  ggsave(
    day_of_week,
    file.path(fp, "day_of_week.png")
  )
}


conditional_effects(ar_mod)

trace_sigma <- mcmc_plot(ar_mod,
  variable = "sigma",
  regex = TRUE,
  type = "trace"
)
ggsave(
  trace_sigma,
  file.path(fp, "trace_sigma.png")
)

trace_ar_coeff <- mcmc_plot(ar_mod,
  variable = "ar1",
  regex = TRUE,
  type = "areas"
)
ggsave(
  trace_ar_coeff,
  file.path(fp, "trace_ar_coeff.png")
)

slopes <- plot_slopes(ar_mod,
  variable = "week",
  condition = c("series", "series"),
  type = "link"
) +
  theme(legend.position = "none") +
  labs(y = "Log(counts)", x = "Location")
ggsave(
  slopes,
  file.path(fp, "slopes.png")
)

# Hierarchical trend effects
trends <- plot(ar_mod, type = "smooths", trend_effects = TRUE)
ggsave(
  trends,
  file.path(fp, "trends.png")
)
# Hierarchical intercepts
intercepts <- plot(ar_mod, type = "re", trend_effects = TRUE)
ggsave(
  intercepts,
  file.path(fp, "intercepts.png")
)
