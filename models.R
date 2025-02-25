#### R script to define and run the models #####
library(argparser, quietly = TRUE)
library(RcppTOML)
library(dplyr)
library(mvgam)
library(cmdstanr)
library(ggplot2)
library(gratia)
library(marginaleffects)

# Command Line Version --------------------------------------------------------
parsed_args <- arg_parser("Preprocess the data for a config") |>
  add_argument("config", help = "Path to TOML config file") |>
  add_argument("config_index",
    help = "index of entry in config to use",
    type = "integer"
  ) |>
  parse_args()

if (!is.na(parsed_args$config)) {
  # config <- parseTOML("input/example_config_weekly.toml") #nolint
  config <- parseTOML(parsed_args$config)
  index <- parsed_args$config_index
} else {
  message("File specified in config filepath does not exist")
}



# Load the data---------------------------------------------------------
model_data_filename <- config$data_filename[index]
fp_data <- file.path(
  config$input_data_path,
  config$forecast_date
)
load(file.path(
  fp_data,
  glue::glue("{config$data_filename[index]}.rda")
))
load(file.path(
  fp_data,
  glue::glue("{config$data_filename[index]}_forecast.rda")
))
## Specify other filepaths for saving ---------------------------------
fp_figs <- file.path(
  "output",
  "figures",
  config$forecast_date,
  config$data_filename[index]
)
fp_mod <- file.path(
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
  ##### Dynamical GAM with vector autoregression #####

  # y_{l,t} \sim Poisson(exp(x_{l,t})) \\
  # x_{l,t} \sim MVNormal(\mu_{l,t} + A X{l,t-1},  \Sigma)\\
  # \mu_{l,t} = \beta_{l,season} + f_{global,t}(weekofyear) + f_{l,t}(weekofyear) + f_{global,t}(wday) \\ #nolint
  # \beta_{l,season} \sim Normal(\beta_l, \sigma_{count}) \\
  # \beta_{l} \sim Normal(\beta_{global}, \sigma_{count}) \\
  # \beta_{global} \sim Normal(log(avgcount), 1) \\
  # \sigma_{count} \sim exp(0.33) \\
  # A \in P(\mathbb{R}) \\
  # P \sim Normal(0, 0.5) T[-1,1] \\
  # \Sigma = \sigma \times C \times \sigma \\
  # \sigma \sim Beta(3,3) \\
  # C \sim LKJcorr(2) \\


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
        s(trend, season, bs = "re") +
        # Hierarchical effects of year(shared smooth)
        # s(year, k = 3) +
        # # Borough level deviations
        # s(year, trend, bs = "sz", k = 3) -1 +

        # Hierarchical effects of seasonality (not time varying for now)
        s(week, k = 12, bs = "cc") +
        # Location  level deviations
        s(week, k = 12, bs = "cc", by = trend) - 1,

    # Shared smooth of day of week
    s(day_of_week, k = 3),
    knots = list(week = c(1, 52)),
    trend_model = VAR(cor = TRUE),
    # Adjust the priors
    priors = c(
      prior(normal(6.69, 1), # data based prior
        class = mu_raw_trend
      ),
      prior(exponential(0.33), class = sigma_raw_trend),
      prior(beta(3, 3), class = sigma, lb = 0.2, ub = 1)
    ),
    data = model_data,
    newdata = forecast_data,
    backend = "cmdstanr",
    family = poisson()
  )
  # NYC weekly count data
} else if (config$targets[index] == "ILI ED visits" && config$regions_to_fit[index] == "NYC" && config$timestep_data[index] == "week") { # nolint
  ##### Dynamical GAM with vector autoregression #####

  # y_{l,t} \sim Poisson(exp(x_{l,t})) \\
  # x_{l,t} \sim MVNormal(\mu_{l,t} + A X{l,t-1},  \Sigma)\\
  # \mu_{l,t} = \beta_{l,season} + f_{global,t}(weekofyear) + f_{l,t}(weekofyear) \\ #nolint
  # \beta_{l,season} \sim Normal(\beta_l, \sigma_{count}) \\
  # \beta_{l} \sim Normal(\beta_{global}, \sigma_{count}) \\
  # \beta_{global} \sim Normal(log(avgcount), 1) \\
  # \sigma_{count} \sim exp(0.33) \\
  # A \in P(\mathbb{R}) \\
  # P \sim Normal(0, 0.5) T[-1,1] \\
  # \Sigma = \sigma \times C \times \sigma \\
  # \sigma \sim Beta(3,3) \\
  # C \sim LKJcorr(2) \\

  # Gets the default priors
  def_priors <- get_mvgam_priors(
    formula = obs_data ~ -1,
    trend_formula = ~
      s(trend, bs = "re") +
        s(trend, season, bs = "re") +
        s(week, k = 12, bs = "cc") +
        s(week, k = 12, bs = "cc", by = trend) - 1,
    knots = list(week = c(1, 52)),
    trend_model = VAR(cor = TRUE),
    data = model_data,
    family = poisson()
  )


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
        s(trend, season, bs = "re") +
        # Hierarchical effects of year(shared smooth)
        # s(year, k = 3) +
        # # Borough level deviations
        # s(year, trend, bs = "sz", k = 3) -1 +

        # Hierarchical effects of seasonality (not time varying for now)
        s(week, k = 12, bs = "cc") +
        # Location level deviations
        s(week, k = 12, bs = "cc", by = trend) - 1,
    knots = list(week = c(1, 52)),
    trend_model = VAR(cor = TRUE),
    # Adjust the priors
    priors = c(
      prior(normal(6.69, 1), # data based prior
        class = mu_raw_trend
      ),
      prior(exponential(0.33), class = sigma_raw_trend),
      prior(beta(3, 3), class = sigma, lb = 0.2, ub = 1)
    ),
    data = model_data,
    newdata = forecast_data,
    backend = "cmdstanr",
    family = poisson(),
    adapt_delta = 0.95
  )
} else if (config$targets[index] == "flu ED visits pct") {
  ##### Dynamical GAM with vector autoregression #####
  # p_{l,t} = y_{l,t} \times 100 \\
  # y_{l,t} \sim Beta (z_{l,t}, \phi) \\
  # logit(z_{l,t}) = x_{l,t} \\
  # x_{l,t} \sim MVNormal(\mu_{l,t} + A X{l,t-1},  \Sigma)\\
  # \mu_{l,t} = \beta_{l,season} + f_{global,t}(weekofyear) + f_{l,t}(weekofyear) \\ #nolint
  # \beta_{l,season} \sim Normal(\beta_l, \sigma_{count}) \\
  # \beta_{l} \sim Normal(\beta_{global}, \sigma_{count}) \\
  # \beta_{global} \sim Normal(log(avgprop), 1) \\
  # \sigma_{count} \sim exp(0.33) \\
  # A \in P(\mathbb{R}) \\
  # P \sim Normal(0, 0.5) T[-1,1] \\
  # \Sigma = \sigma \times C \times \sigma \\
  # \sigma \sim Beta(3,3) \\
  # C \sim LKJcorr(2) \\


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

        # Hierarchical effects of seasonality (not time varying for now)
        s(week, k = 12, bs = "cc") +
        # Location level deviations
        s(week, k = 12, bs = "cc", by = trend) - 1,
    knots = list(week = c(1, 52)),
    trend_model = "AR1",
    # Adjust the priors
    priors = c(
      prior(normal(-4.5, 1),
        class = mu_raw_trend
      ),
      prior(exponential(0.33), class = sigma_raw_trend),
      prior(exponential(1), class = sigma),
      prior(normal(0.5, 0.25), class = ar1, lb = 0, ub = 1)
    ),
    data = model_data,
    newdata = forecast_data,
    backend = "cmdstanr",
    family = betar()
  )
}

#### Make a bunch of plots and save them ####----------------------------


summary <- summary(ar_mod)
code <- stancode(ar_mod)
save(ar_mod, file = file.path(
  fp_mod,
  glue::glue("{config$model_filename[index]}.rda")
))
writeLines(code, file.path(
  fp_mod,
  glue::glue("{config$model_filename[index]}.stan")
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

plot_predictions(ar_mod, condition = "week", type = "link")

if (config$targets[index] == "ILI ED visits" && config$regions_to_fit[index] == "NYC" && config$timestep_data[index] == "day") { # nolint
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

# Make a plot of the matrix of autoregulation coefficients
n_locs <- length(unique(model_data$location))
a_pars <- matrix(NA,
  nrow = n_locs,
  ncol = n_locs
)
for (i in 1:n_locs) {
  for (j in 1:n_locs) {
    a_pars[i, j] <- paste0("A[", i, ",", j, "]")
  }
}
plot_a_matrix <- mcmc_plot(ar_mod,
  variable = as.vector(t(a_pars)),
  type = "hist"
) +
  geom_vline(
    xintercept = 0,
    col = "white",
    linewidth = 2
  ) +
  geom_vline(
    xintercept = 0,
    linewidth = 1
  )
ggsave(
  plot = plot_a_matrix,
  filename = file.path(fp_figs, "a_matrix.png")
)

sigma_pars <- matrix(NA, nrow = n_locs, ncol = n_locs)
for (i in 1:n_locs) {
  for (j in 1:n_locs) {
    sigma_pars[i, j] <- paste0("Sigma[", i, ",", j, "]")
  }
}
plot_sigma <- mcmc_plot(ar_mod,
  variable = as.vector(t(sigma_pars)),
  type = "hist"
) +
  geom_vline(
    xintercept = 0,
    col = "white",
    linewidth = 2
  ) +
  geom_vline(
    xintercept = 0,
    linewidth = 1
  )
ggsave(
  plot = plot_sigma,
  filename = file.path(fp_figs, "Sigma.png")
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
