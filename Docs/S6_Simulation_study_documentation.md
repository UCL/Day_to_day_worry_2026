Supplementary Material: Monte Carlo Simulation: Q2 Parameter Recovery
and Statistical Power
================
Diego Vitali, Raffaele Buono, Tao Bi, Tony Pickering, Chris Eccleston,
Edmund Keogh, Nadia Bianchi-Berthouze, Amanda C de C Williams
2026-02-07

- [Summary & Logic](#summary--logic)
- [Methods](#methods)
  - [Model Fitted to Each Simulated
    Dataset](#model-fitted-to-each-simulated-dataset)
  - [Metrics](#metrics)
  - [Data-Generating Function](#data-generating-function)
  - [Single-Replication Function](#single-replication-function)
  - [Main Simulation Loop](#main-simulation-loop)
    - [Parameter Recovery (Bias)](#parameter-recovery-bias)
    - [Standard Error Calibration](#standard-error-calibration)
    - [Full Results Table](#full-results-table)
  - [power simulation chart](#power-simulation-chart)
    - [Power for the Learning
      Parameter](#power-for-the-learning-parameter)

# Summary & Logic

This supplementary material reports Monte Carlo simulations used to
validate the lagged-level model specification employed in Q2, which
tested whether participants adjust their future interference predictions
in response to prior worry–outcome errors. The worry–outcome error was
defined as the difference between morning worry and evening actual
interference (worry - actual), where positive values indicate over-worry
and negative values indicate under-worry.

1.  PARAMETER RECOVERY: Does the lagged-level specification return
    unbiased estimates of the AR(1) persistence parameter $\phi$ and the
    learning parameter $\alpha$?

2.  STATISTICAL POWER: Given our actual design (N=18, T=28 days), what
    is the probability of detecting learning effects of the sizes
    represented in our simulation grid, in order to contextualise
    whether the observed learning coefficient
    ($\hat\alpha \approx -0.08$) is more consistent with a genuine small
    effect or a chance finding.

Following the ADEMP framework for simulation studies ([Morris, White &
Crowther 2019)](https://pmc.ncbi.nlm.nih.gov/articles/PMC6492164/), we
specify aims, data-generating process (DGP), estimands, methods, and
performance measures.

# Methods

Each simulated dataset replicates our study design: N = 18 participants
each contributing T = 28 daily assessments, consistent with our actual
EMA protocol.

The data-generating process (DGP) implements a first-order
autoregressive model (AR(1)) with error-driven updating:

$$W_t = \phi \cdot W_{t-1} + \alpha \cdot E_{t-1} + u_i + \varepsilon_t$$

where:

- $W_t$ is morning worry about pain interference on day $t$.
- $E_{t-1} = W_{t-1} - A_{t-1}$ is the prediction error from the
  previous day (the difference between anticipated and experienced
  interference).
- $u_i\sim\mathcal{N}(0,\sigma_u^2)$ is a participant-level random
  intercept capturing stable individual differences
- $\varepsilon_t \sim \mathcal{N}(0, \sigma_\varepsilon^2)$ is day-level
  noise.

The estimands of interest are the conditional regression coefficients
$\phi$ (worry persistence) and $\alpha$ (error-driven learning rate),
both estimated from $$E[W_t\mid W_{t-1},E_{t-1},u_i]$$ Actual
interference $A_t$ is generated as daily white noise. This is a
deliberate simplifying assumption: the simulation targets recovery of
$\phi$ and $\alpha$ specifically, and introducing autocorrelated
interference would conflate the learning signal with interference
carry-over, obscuring the validation purpose of the simulation. All
distributional parameters and variance components were derived from the
first imputed dataset and are reported in the table below.

| Parameter   | Value | Source                             |
|:------------|:-----:|:-----------------------------------|
| W init mean | 3.95  | grand mean worry                   |
| W init SD   | 1.41  | mean within-person SD worry        |
| A mean      | 4.55  | grand mean interference            |
| A SD        | 1.68  | mean within-person SD interference |
| sigma_u     | 1.66  | SD of participant means            |
| sigma_eps   | 1.43  | residual SD from fitted Q2 model   |

The DGP is intentionally parsimonious: it includes only the core
autoregressive and learning components rather than the full set of
covariates included in the main analysis (3-day pain trend and
deviations, morning mood and stress, age).

Simulations were run across a factorial grid of:

- **AR(1) persistence** ($\phi$): 0.2, 0.4, 0.6 — spanning weak to
  moderate carry-over of worry across days.
- **Learning rate** ($\alpha$): -0.30, -0.20, -0.10, 0 — spanning
  clinically plausible negative updating to the null hypothesis of no
  learning; negative values indicate that overestimates on day t-1
  reduce worry on day t

This factorial grid span empirically plausible values of $\phi$ and
$\alpha$ so that we can explore a range of realistic scenarios rather
than a single parameter [(Morris, White & Crowther
2019)](https://pmc.ncbi.nlm.nih.gov/articles/PMC6492164/).

This yields 12 conditions. Each condition was replicated 1,000 times
(nsim = 1,000), giving 12,000 fitted models in total.

## Model Fitted to Each Simulated Dataset

Each simulated dataset was fitted with the pre-registered lagged-level
specification:

$$W_t = \beta_0 + \phi \cdot W_{t-1} + \alpha \cdot E_{t-1} + \beta_1 \cdot \text{day} + u_i + \varepsilon_t$$

estimated via restricted maximum likelihood (REML) using `lmer` (lme4
package; [Bates et al.,
2015](https://www.jstatsoft.org/article/view/v067i01)), with
Satterthwaite degrees of freedom for p-values (lmerTest; [Kuznetsova et
al., 2017](https://www.jstatsoft.org/article/view/v082i13)). A random
intercept for participant captures $u_i$; no random slopes were
included, matching our main analysis specification.

## Metrics

For each condition we computed:

- **Bias**: $\text{mean}(\hat{\phi}) - \phi_{\text{true}}$ and
  $\text{mean}(\hat{\alpha}) - \alpha_{\text{true}}$, quantifying
  systematic under- or over-estimation
- **Empirical SE**: SD of the 1,000 estimates, representing the true
  sampling variability
- **Model-based SE**: mean of the 1,000 model-reported standard
  errors,representing how well the model’s uncertainty estimates track
  true variability
- **SE ratio**: model-based SE / empirical SE; values near 1.0 indicate
  well-calibrated uncertainty; values \< 1.0 indicate underestimation
  (anti-conservative); values \> 1.0 indicate overestimation
  (conservative)
- **Power**: proportion of replications where p \< 0.05 for each
  parameter; for $\alpha = 0$ this gives the empirical Type I error rate

## Data-Generating Function

The simulate_panel() function generates a single longitudinal dataset
under the DGP described above. Each participant’s worry trajectory is
initialised at a random draw from $\mathcal{N}(3.95, 1.41)$ and evolved
according to the AR(1) + learning equation; actual interference ($A_t$)
is drawn each day from $\mathcal{N}(4.55, 1.68)$ with $\sigma_u=1.66$
and $\sigma_\varepsilon=1.43$.

Lagged variables are created using `ave()` within participants to
respect the panel structure.

``` r
simulate_panel <- function(N, T, phi, alpha, sigma_u, sigma_eps) {
  ids <- 1:N
  u   <- rnorm(N, 0, sigma_u)   # participant random intercepts
  out <- vector("list", N)

  for (i in ids) {
    W <- numeric(T)   # worry
    A <- numeric(T)   # actual interference
    E <- numeric(T)   # prediction error (W - A)

    # Initialise day 1
    W[1] <- rnorm(1, 3.95, 1.41) #  grand mean and mean within-person SD of morning worry;
    A[1] <- rnorm(1, 4.55, 1.68) #  grand mean and mean within-person SD of evening interference;
    E[1] <- W[1] - A[1]

    for (t in 2:T) {
      # Actual interference: white noise (no autocorrelation)
      A[t] <- rnorm(1, 4.55, 1.68)

      # Worry: AR(1) persistence + error-driven learning + person effect + noise
      W[t] <- phi * W[t-1] + alpha * E[t-1] + u[i] + rnorm(1, 0, sigma_eps)
      E[t] <- W[t] - A[t]
    }

    out[[i]] <- data.frame(
      alias            = i,
      day              = 1:T,
      day_count        = 1:T,
      interf_worry     = W,
      actual           = A,
      worry_outcome_error = E
    )
  }

  panel <- do.call(rbind, out)
  panel <- panel[order(panel$alias, panel$day), ]

  # Lagged variables created within-participant (lag-1)
  panel$lagged_interf_worry <- ave(
    panel$interf_worry, panel$alias,
    FUN = function(x) c(NA, x[-length(x)])
  )
  panel$lagged_worry_outcome_error <- ave(
    panel$worry_outcome_error, panel$alias,
    FUN = function(x) c(NA, x[-length(x)])
  )

  panel
}
```

## Single-Replication Function

`run_once()` calls `simulate_panel()`, fits the lagged-level model, and
returns the key estimates and p-values for $\phi$ and $\alpha$. The day
1 observations are dropped (NA lagged predictors) before fitting.

``` r
run_once <- function(N, T, phi_true, alpha_true, sigma_u, sigma_eps) {
  dat     <- simulate_panel(N, T, phi_true, alpha_true, sigma_u, sigma_eps)
  dat_fit <- subset(dat, !is.na(lagged_interf_worry))

  m <- lmer(
    interf_worry ~ lagged_interf_worry + lagged_worry_outcome_error +
      day_count + (1 | alias),
    data = dat_fit,
    REML = TRUE
  )

  coef_summary <- summary(m)$coefficients

  c(
    phi_hat   = coef_summary["lagged_interf_worry",        "Estimate"],
    alpha_hat = coef_summary["lagged_worry_outcome_error", "Estimate"],
    phi_se    = coef_summary["lagged_interf_worry",        "Std. Error"],
    alpha_se  = coef_summary["lagged_worry_outcome_error", "Std. Error"],
    phi_p     = coef_summary["lagged_interf_worry",        "Pr(>|t|)"],
    alpha_p   = coef_summary["lagged_worry_outcome_error", "Pr(>|t|)"]
  )
}
```

## Main Simulation Loop

The simulation iterates over the 12 conditions in `validation_grid`.

Within each condition, `replicate()` calls `run_once()` 1,000 times and
the results are summarised into bias, SE calibration, and power metrics.

`set.seed(42)` is set before the loop to ensure reproducibility.

``` r
set.seed(42)

nsim <- 1000

validation_grid <- expand.grid(
  phi_true   = c(0.2, 0.4, 0.6),
  alpha_true = c(-0.3, -0.2, -0.1, 0),
  stringsAsFactors = FALSE
)

results <- list()

for (i in 1:nrow(validation_grid)) {
  cat(sprintf(
    "Running condition %d/%d: phi=%.2f, alpha=%.2f\n",
    i, nrow(validation_grid),
    validation_grid$phi_true[i],
    validation_grid$alpha_true[i]
  ))

  res <- t(replicate(
    nsim,
    run_once(
      N          = 18, # sample size
      T          = 28, # number of days
      phi_true   = validation_grid$phi_true[i],
      alpha_true = validation_grid$alpha_true[i],
      sigma_u    = 1.66, # is the SD of participant-level means of morning worry;
      sigma_eps  = 1.43  # is the residual SD from the fitted Q2 model.  
    )
  ))

  results[[i]] <- data.frame(
    phi_true   = validation_grid$phi_true[i],
    alpha_true = validation_grid$alpha_true[i],

    # Bias: mean estimate minus true value
    phi_bias   = mean(res[, "phi_hat"])   - validation_grid$phi_true[i],
    alpha_bias = mean(res[, "alpha_hat"]) - validation_grid$alpha_true[i],

    # Empirical SE: SD across replications (true sampling variability)
    phi_se_empirical   = sd(res[, "phi_hat"]),
    alpha_se_empirical = sd(res[, "alpha_hat"]),

    # Model-based SE: mean of model-reported SEs
    phi_se_model   = mean(res[, "phi_se"]),
    alpha_se_model = mean(res[, "alpha_se"]),

    # Power / Type I error rate
    phi_power   = mean(res[, "phi_p"]   < 0.05),
    alpha_power = mean(res[, "alpha_p"] < 0.05)
  )
}

validation_results <- do.call(rbind, results)
#saveRDS(validation_results,"c:/Documents and Settings/d.vitali/Desktop/Github/CRIISP-WP6/Data/Sim_validation_results.rds")
```

One convergence warning was observed across 12,000 model fits (gradient
= 0.003, marginally above tolerance = 0.002)

### Parameter Recovery (Bias)

Bias was negligible across all 12 conditions. The maximum absolute bias
was 0.011 for $\phi$ and 0.007 for $\alpha$, representing at most 5.5%
of the smallest true $\phi$ value (0.2) and 2.3% of the smallest
non-zero true $\alpha$ value (0.30).

``` r
cat("Maximum absolute bias (phi)  :", 
    round(max(abs(validation_results$phi_bias)), 5), "\n")
```

    ## Maximum absolute bias (phi)  : 0.01129

``` r
cat("Maximum absolute bias (alpha):", 
    round(max(abs(validation_results$alpha_bias)), 5), "\n")
```

    ## Maximum absolute bias (alpha): 0.00645

This confirms that the lagged-level model recovers both the AR(1)
persistence parameter and the learning parameter without systematic
distortion, across the full range of persistence and learning rates
considered plausible for our population.

### Standard Error Calibration

The ratio of model-based SEs to empirical SEs quantifies whether the
model’s reported uncertainty matches the true sampling variability.
Ratios near 1.0 indicate accurate calibration; ratios above 1.0 indicate
conservative uncertainty estimates (acceptable); ratios below 1.0 would
indicate anti-conservative SEs (problematic for inference).

``` r
cat("Mean SE ratio model/empirical (phi)  :", 
    round(mean(validation_results$phi_se_model / 
               validation_results$phi_se_empirical), 3), "\n")
```

    ## Mean SE ratio model/empirical (phi)  : 0.988

``` r
cat("Mean SE ratio model/empirical (alpha):", 
    round(mean(validation_results$alpha_se_model / 
               validation_results$alpha_se_empirical), 3), "\n")
```

    ## Mean SE ratio model/empirical (alpha): 0.998

SE ratios near 1.0 confirm that the model-reported standard errors are
well calibrated relative to the true sampling variability, with a slight
conservative tendency. This indicates that hypothesis tests and
confidence intervals derived from the fitted models are reliable.

### Full Results Table

The table below reports all 12 conditions. `phi_bias` and `alpha_bias`
are the primary validation metrics; `phi_power` and `alpha_power` give
the proportion of replications achieving p \< 0.05 (power when the true
effect != 0; Type I error rate when the true effect = 0).

| phi (true) | alpha (true) | Bias phi | Bias alpha | SE phi (empirical) | SE alpha (empirical) | SE phi (model) | SE alpha (model) | Power phi | Power alpha |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.2 | -0.3 | -0.0050 | 0.0021 | 0.0821 | 0.0693 | 0.0820 | 0.0693 | 0.660 | 0.989 |
| 0.4 | -0.3 | -0.0064 | 0.0014 | 0.0820 | 0.0702 | 0.0822 | 0.0695 | 0.999 | 0.989 |
| 0.6 | -0.3 | -0.0102 | 0.0050 | 0.0827 | 0.0686 | 0.0814 | 0.0695 | 1.000 | 0.989 |
| 0.2 | -0.2 | -0.0044 | -0.0005 | 0.0820 | 0.0696 | 0.0821 | 0.0694 | 0.659 | 0.821 |
| 0.4 | -0.2 | -0.0113 | 0.0065 | 0.0805 | 0.0680 | 0.0817 | 0.0694 | 1.000 | 0.791 |
| 0.6 | -0.2 | -0.0042 | -0.0038 | 0.0839 | 0.0710 | 0.0805 | 0.0694 | 1.000 | 0.821 |
| 0.2 | -0.1 | -0.0063 | 0.0023 | 0.0845 | 0.0722 | 0.0817 | 0.0693 | 0.658 | 0.288 |
| 0.4 | -0.1 | -0.0092 | 0.0011 | 0.0809 | 0.0678 | 0.0811 | 0.0694 | 0.996 | 0.285 |
| 0.6 | -0.1 | -0.0088 | 0.0029 | 0.0807 | 0.0682 | 0.0792 | 0.0693 | 1.000 | 0.296 |
| 0.2 | 0.0 | -0.0065 | 0.0012 | 0.0820 | 0.0696 | 0.0812 | 0.0694 | 0.678 | 0.054 |
| 0.4 | 0.0 | -0.0099 | 0.0001 | 0.0806 | 0.0681 | 0.0798 | 0.0692 | 1.000 | 0.043 |
| 0.6 | 0.0 | -0.0084 | -0.0016 | 0.0811 | 0.0718 | 0.0775 | 0.0693 | 1.000 | 0.065 |

Monte Carlo simulation results: bias, SE calibration, and power across
12 conditions (N = 18, T = 28, nsim = 1,000, seed = 42).

## power simulation chart

![](C:/Users/d.vitali/Desktop/Github/EMA_Pain_Supl_2026/Docs/S6_Simulation_study_documentation_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

### Power for the Learning Parameter

The power results for $\alpha$ deserve specific attention, as they
directly contextualise the near-significant learning coefficient (
$\hat\alpha = -0.08$ ) observed in our main analysis. For the smallest
non-zero learning effect in our grid ( $\alpha = -0.10$ ), power ranged
from 28–30% across persistence levels — substantially below conventional
thresholds. Power increased to 79–82% for $\alpha = -0.20$, indicating
that our design is adequate for detecting moderate learning effects but
underpowered for small ones.
