# cityforecasts
R scripts to generate city-level forecasts

## Summary
This repository contains R code to generate city-level forecasts for submission to the [flu-metrocast hub](https://github.com/reichlab/flu-metrocast).
All models are exploratory/preliminary, though we will regularly update this document to describe the latest mathematical model used in the submission.

All outputs submitted to the Hub will be archived in this repository, along with additional model metadata (such as the model definition associated with a submission and details on any additional data sources used or decisions made in the submission process).
If significant changes to the model are made during submission, we will rename the model in the submission file.

Initially, we plan to fit the data from each state independently, using hierarchical partial pooling to jointly fit the cities within a state.
This initially includes producing forecast of:
|Forecast target| Location |
|-----------|-----------|
| ED visits due to ILI | New York City (5 boroughs, unknown, citywide) |
| Percent of ED visits due to flu | Texas (5 metro areas) |

We plan to use the same latent model structure for both forecast targets, modifying the observation model for count data (NYC) vs proportion data (Texas).

## Workflow
Because all data is available publicly, the forecasts generated should be completely reproducible from the specified configuration file.
We start by using the [`mvgam`](https://github.com/nicholasjclark/mvgam) package, which is a an R package that leverages both [`mgcv`](https://cran.r-project.org/web/packages/mgcv/index.html) and [`brms`](https://paulbuerkner.com/brms/) formula interface to fit Bayesian Dynamic Generalized Additive Models (GAMs).
These packages use metaprogramming to produce Stan files, and we also include the Stan code generated by the package.

To produce forecasts each week we follow the following workflow:


1. Modify the configuration file in `input/config.toml`
2. In the command line, run ` Rscript preprocess_data.R input/config.toml {index}` where index is used to track the individual model runs, which in this case, also have different pre-processing due to being from different data sources. 
3. Next run ` Rscript models.R input/config.toml {index}`
4. Lastly run `Rscript postprocess_forecasts.R input/{forecast_date}/config.toml`
5. This will populate the `output/cityforecasts/{forecast_date}` folder with a csv file formatted following the Hub submission guidelines.

Eventually, steps 2-4 will be automated with the Github Action `.git/workflows/generate_forecasts` and set on a schedule to run after 12 pm CST, corresponding to the time that the `target_data` is updated on the Hub.

## Model definition

The below describes the preliminary model used:
### Observation model
For the forecasts of counts due to ED visits, we assume a Poisson observation process

$$
y_{l,t} \sim Poisson(exp(x_{l,t}))
$$

For the forecasts of the percent of ED visits due to flu, we assume a Beta observation process on the proportion of ED visits due to flu:

```math
\begin{align}
p_{l,t} = y_{l,t} \times 100 \\
y_{l,t} \sim Beta (z_{l,t}, \phi) \\
logit(z_{l,t}) = x_{l,t}
\end{align}
```

### Latent state-space model: Dynamic hierachical GAM with vector autoregression
We model latent admissions with a hierarchical GAM component to capture shared seasonality and weekday effects and a multivariate vector autoregressive component to capture trends in the dynamics within and between each location.
Initially, we will use only the weekly data, so $t$ will be indexed in weeks.

```math
\begin{align}
x_{l,t} \sim Normal(\mu_{l,t} + A X_{l,t-1},  \Sigma)\\
\mu_{l,t} = \beta_{l,season} + f_{global,t}(weekofyear) + f_{l,t}(weekofyear) \\
\beta_{l,season} \sim Normal(\beta_l, \sigma_{count}) \\
\beta_{l} \sim Normal(\beta_{global}, \sigma_{count}) \\
\sigma_{count} \sim exp(0.33) \\
A \in P(\mathbb{R})\\
P \sim Normal(0, 0.5) T[-1,1] \\
\Sigma = \sigma \times C \times \sigma \\
\sigma \sim Beta(3,3) \\
C \sim LKJcorr(2) \\
\end{align}
```


For the NYC data, we have count data and so $\beta_{global}$ represents the intercept on the count scale, we place a prior on it using the mean observed count across the historical data:

$$
\beta_{global} \sim Normal(log(\frac{\sum_{l=1}^L \sum_{t=1}^T y_{l,t}}{N_{obs}}), 1) \\
$$

where $N_obs$ is the number of observations of $y_{l,t}$.

For the TX data, $\beta_{global}$ represents the intercept as a proportion, so we use:

$$
\beta_{global} \sim Normal(logit(\frac{\sum_{l=1}^L \sum_{t=1}^T y_{l,t}}{N_{obs}}), 1) \\
$$
