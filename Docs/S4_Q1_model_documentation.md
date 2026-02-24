Supplementary Material: Q1 coupling model
================
Diego Vitali, Raffaele Buono, Tao Bi, Tony Pickering, Chris Eccleston,
Edmund Keogh, Nadia Bianchi-Berthouze, Amanda C de C Williams
2026-02-07

# Summary & Logic

This document presents the statistical models testing Question 1:
Worry-outcome coupling. We examine whether morning worry about
interference predicts actual evening interference on the same day,
controlling for pain context and concurrent affective states. The
analysis uses multiply imputed data (m=28 imputations) with linear
mixed-effects models accounting for repeated daily assessments nested
within participants.

All models use Rubin’s pooling rules to combine estimates across
imputations, thereby accounting for missing data uncertainty in
parameter estimates, standard errors, and significance tests.

# Q1: LMM for pred./Exp. coupling

The worry-outcome coupling model tests whether morning predictions about
interference have predictive validity for actual evening interference.

## DAG

The directed acyclic graph (DAG) shows our assumed causal structure.
Morning worry (exposure) and pain context both directly predict evening
interference (outcome), while mood and stress act as common causes of
both worry and interference. Age and study day are included as
background covariates. The random intercept structure is not shown in
this DAG but it is included in the model allowing these relationships to
vary around person-specific baselines.

![](C:/Users/d.vitali/Desktop/Github/EMA_Pain_Supl_2026/Docs/S4_Q1_model_documentation_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

The baseline model in the DAG tests morning worry about interference
(interf_worry) as a significant predictor of evening interference, with
a positive coefficient indicating that higher morning predictions
correspond to higher evening interference ratings.

## Model Specifications and Comparisons

We fit three primary models and several sensitivity analyses to assess
the robustness of the worry-interference coupling estimate.

- Model 0 (Baseline): Tests worry-interference coupling controlling for
  pain trend and current affective states
- Model 0a (Trend alternative): Same as model 0 but substitutes trend
  measures for current affective states

two models that extend model 0:

- Model 1 (Pain deviation): Extends model 0 adding acute pain
  fluctuation as predictor
- Model 2 (Interaction): Extends model 1 testing whether
  worry-interference coupling varies by acute pain

and an additional sensitivity analysis for the autocorrelation structure
of evening interference: \* Model 1_autocorr: Controls for
autoregressive structure of evening interference

All models are pooled using Rubin´s Rules. RR are designed to pool
parameter estimates, such as mean differences, regression coefficients,
standard errors and to derive confidence intervals and p-values.
\[<https://bookdown.org/mwheymans/bookmi/rubins-rules.html>\]

This means that by pooling parameter estimates we get the pooled means
of the same parameter estimated across all imputed datasets.

This is helpful to reflect the uncertainty around the imputation we have
done in our parameter estimates. Using mice::pool() we get appropriate
standard errors estimates that reflect both within-imputation variance
(the average of the standard errors from each dataset), and
between-imputation variance (the variability in parameter estimates
across the imputed datasets). Degrees of freedom are also adjusted to
account for the number of imputations and the for the proportion of
missing data; P-values, confidence intervals, and test statistics use
the pooled estimates and variances so they appropriately reflect the
extra uncertainty resulting from missing data imputation.
\[<https://bookdown.org/mwheymans/bookmi/rubins-rules.html#significance-testing>\]

``` r
# Fit Q1 models testing worry-outcome coupling across 28 imputed datasets.
# Model 0: Baseline with pain trend and current affective states
# Model 0a: Alternative using trend variables for stress/mood instead of current values

model_Coupling_pooled_0 <- with(imp_1_with_trends, 
                    lmer(evening_interference ~ interf_worry + 
                        pain_trend_3d +
                        morning_mood + 
                        morning_stress + 
                        day_count + 
                        Age_c +
                        (1 | alias)))
model_Coupling_pooled_0a <- with(imp_1_with_trends, 
                    lmer(evening_interference ~ interf_worry + 
                        pain_trend_3d +
                        stress_trend_3d + 
                        mood_trend_3d + 
                        day_count + 
                        Age_c +
                        (1 | alias)))
```

We compare Model 0 (current mood/stress) against Model 0a (trend
mood/stress) using information criteria to determine which specification
better captures the affective contribution to interference.

``` r
AIC_model_Coupling_pooled_0 <- mean(sapply(model_Coupling_pooled_0$analyses, AIC))
AIC_model_Coupling_pooled_0a <- mean(sapply(model_Coupling_pooled_0a$analyses, AIC))

BIC_model_Coupling_pooled_0 <- mean(sapply(model_Coupling_pooled_0$analyses, BIC))
BIC_model_Coupling_pooled_0a <- mean(sapply(model_Coupling_pooled_0a$analyses, BIC))


# Compare model fit using information criteria averaged across imputations
# Lower AIC/BIC indicates better balance of fit and parsimony
# Comparing model 0 (current mood/stress) vs 0a (trend mood/stress)
rbind(
  AIC_Coupling_pooled_mean_0 = AIC_model_Coupling_pooled_0,
  AIC_Coupling_pooled_mean_0a = AIC_model_Coupling_pooled_0a,
  BIC_Coupling_pooled_mean_0 = BIC_model_Coupling_pooled_0,
  BIC_Coupling_pooled_mean_0a = BIC_model_Coupling_pooled_0a
) %>% kable
```

|                             |          |
|:----------------------------|---------:|
| AIC_Coupling_pooled_mean_0  | 1977.222 |
| AIC_Coupling_pooled_mean_0a | 1988.272 |
| BIC_Coupling_pooled_mean_0  | 2015.172 |
| BIC_Coupling_pooled_mean_0a | 2026.221 |

Then we test the extended models (model 1 and 2) to examine whether
adding acute pain fluctuations (deviation from trend) and their
interaction with worry improves model fit beyond the baseline
specification.

The D1 test implements a Wald test for multiply imputed data, combining
test statistics across imputations using Rubin’s rules. Significant
p-values (\<0.05) indicate that the added predictor(s) significantly
improve model fit.

``` r
# pain_deviation_from_trend = morning_pain - pain_trend_3d,
model_Coupling_pooled_1 <- with(imp_1_with_trends, 
                    lmer(evening_interference ~ interf_worry + 
                        pain_deviation_from_trend + 
                        pain_trend_3d +
                        morning_mood + 
                        morning_stress + 
                        day_count + 
                        Age_c +
                        (1 | alias)))

model_Coupling_pooled_2 <- with(imp_1_with_trends, 
                    lmer(evening_interference ~ interf_worry * pain_deviation_from_trend + 
                        pain_trend_3d +
                        morning_mood + 
                        morning_stress + 
                        day_count + 
                        Age_c +
                        (1 | alias)))


# Compare nested models using D1 statistic (Wald test for multiply imputed data)
# D1 tests whether adding pain_deviation_from_trend improves fit beyond model_0
# Null hypothesis: additional predictor has zero effect
D1(model_Coupling_pooled_1, model_Coupling_pooled_0)
```

    ##    test statistic df1      df2 dfcom   p.value       riv
    ##  1 ~~ 2  1.415163   1 321.3422   491 0.2350795 0.1703686

``` r
D1(model_Coupling_pooled_2,model_Coupling_pooled_1)
```

    ##    test statistic df1      df2 dfcom   p.value       riv
    ##  1 ~~ 2 0.7912472   1 250.0268   490 0.3745783 0.2465724

### Sensitivity Analysis for Lagged Effect Of Interference

We test whether the worry-interference coupling remains robust when
controlling for autoregressive structure in evening interference (i.e.,
whether yesterday’s interference predicts today’s interference beyond
current predictors).

``` r
# sensitivity analysis controlling for autocorrelation
model_Coupling_pooled_autocorr <- with(imp_1_with_trends, 
                    lmer(evening_interference ~ interf_worry + 
                        lagged_evening_interference +
                        pain_trend_3d +
                        morning_mood + 
                        morning_stress + 
                        day_count + 
                        Age_c +
                        (1 | alias)))
Q1_Sens_A_pooled_summary <- summary(pool(model_Coupling_pooled_autocorr), conf.int = TRUE) 

Q1_Sens_A_pooled_summary <-Q1_Sens_A_pooled_summary %>%mutate(
  estimate = round(estimate, 3),
  std.error = round(std.error, 3),
  statistic = round(statistic, 3),
  df = round(df, 2),    # degrees of freedom often rounded to 2 decimals
  p.value = round(p.value, 3),
  `2.5 %` = round(`2.5 %`, 3),
  `97.5 %` = round(`97.5 %`, 3)
)

kableExtra::kable(Q1_Sens_A_pooled_summary[1:8])
```

| term | estimate | std.error | statistic | df | p.value | 2.5 % | 97.5 % |
|:---|---:|---:|---:|---:|---:|---:|---:|
| (Intercept) | 3.202 | 0.607 | 5.275 | 319.96 | 0.000 | 2.008 | 4.397 |
| interf_worry | 0.376 | 0.058 | 6.499 | 195.20 | 0.000 | 0.262 | 0.491 |
| lagged_evening_interference | 0.099 | 0.049 | 2.028 | 279.32 | 0.043 | 0.003 | 0.195 |
| pain_trend_3d | 0.250 | 0.078 | 3.195 | 376.59 | 0.002 | 0.096 | 0.404 |
| morning_mood | -0.249 | 0.069 | -3.604 | 295.65 | 0.000 | -0.385 | -0.113 |
| morning_stress | -0.040 | 0.049 | -0.822 | 368.73 | 0.411 | -0.136 | 0.056 |
| day_count | -0.005 | 0.009 | -0.593 | 378.21 | 0.553 | -0.024 | 0.013 |
| Age_c | 0.046 | 0.028 | 1.626 | 478.03 | 0.105 | -0.009 | 0.101 |

## Multicollinearity (pooled)

We calculate Variance Inflation Factors using
[check_collinearity()](https://www.rdocumentation.org/packages/performance/versions/0.12.4/topics/check_collinearity)
from the
[performance](https://www.rdocumentation.org/packages/performance/versions/0.12.4)
quantify how much predictor standard errors are inflated due to
correlations with other predictors. VIF \< 5 indicates acceptable
collinearity; VIF \> 10 suggests serious multicollinearity requiring
model revision.

``` r
# Check variance inflation factors (VIF) to detect multicollinearity
# VIF > 5 suggests problematic correlation between predictors
# Only checking first imputed dataset (pattern consistent across imputations)
check_collinearity(model_Coupling_pooled_0$analyses[[1]])
```

    ## # Check for Multicollinearity
    ## 
    ## Low Correlation
    ## 
    ##            Term  VIF    VIF 95% CI adj. VIF Tolerance Tolerance 95% CI
    ##    interf_worry 1.10 [1.04,  1.28]     1.05      0.91     [0.78, 0.96]
    ##   pain_trend_3d 1.03 [1.00,  1.74]     1.01      0.97     [0.57, 1.00]
    ##    morning_mood 1.18 [1.10,  1.35]     1.09      0.84     [0.74, 0.91]
    ##  morning_stress 1.24 [1.14,  1.41]     1.11      0.81     [0.71, 0.88]
    ##       day_count 1.01 [1.00, 16.52]     1.01      0.99     [0.06, 1.00]
    ##           Age_c 1.07 [1.02,  1.28]     1.04      0.93     [0.78, 0.98]

## Pooled Parameter Estimates

``` r
# Get the pooled summary with confidence intervals
Q1_pooled_summary <- summary(pool(model_Coupling_pooled_0), conf.int = TRUE) # The pool() function combines the estimates from m repeated complete data analyses

# Round relevant columns to 3 decimal places (adjust as needed)
Q1_pooled_summary <-Q1_pooled_summary %>%mutate(
  estimate = round(estimate, 3),
  std.error = round(std.error, 3),
  statistic = round(statistic, 3),
  df = round(df, 2),    # degrees of freedom often rounded to 2 decimals
  p.value = round(p.value, 3),
  `2.5 %` = round(`2.5 %`, 3),
  `97.5 %` = round(`97.5 %`, 3)
)

# Print rounded table
# pain_deviation_from_trend = morning_pain - pain_trend_3d,
kableExtra::kable(Q1_pooled_summary[1:8])
```

| term           | estimate | std.error | statistic |     df | p.value |  2.5 % | 97.5 % |
|:---------------|---------:|----------:|----------:|-------:|--------:|-------:|-------:|
| (Intercept)    |    3.369 |     0.612 |     5.504 | 339.20 |   0.000 |  2.165 |  4.573 |
| interf_worry   |    0.392 |     0.058 |     6.813 | 199.24 |   0.000 |  0.279 |  0.505 |
| pain_trend_3d  |    0.308 |     0.074 |     4.171 | 345.38 |   0.000 |  0.163 |  0.454 |
| morning_mood   |   -0.249 |     0.069 |    -3.589 | 305.00 |   0.000 | -0.385 | -0.112 |
| morning_stress |   -0.041 |     0.049 |    -0.837 | 382.04 |   0.403 | -0.137 |  0.055 |
| day_count      |   -0.006 |     0.009 |    -0.677 | 360.46 |   0.499 | -0.025 |  0.012 |
| Age_c          |    0.048 |     0.030 |     1.595 | 478.82 |   0.111 | -0.011 |  0.106 |

## Variance explained (pooled)

We decompose the total variance in evening interference into components
explained by fixed effects (predictors), random effects (between-person
differences), and residual (unexplained within-person fluctuation).

``` r
# Decompose variance explained using Nakagawa's R² for mixed models.
# R²_marginal = variance explained by fixed effects (predictors)
# R²_conditional = total variance explained (fixed + random effects)
r2_list <- lapply(model_Coupling_pooled_0$analyses, function(fit) {
  r2 <- r2_nakagawa(fit)
  c(R2_marginal = r2$R2_marginal, 
    R2_conditional = r2$R2_conditional)
})

# Pool across imputations
r2_pooled <- colMeans(do.call(rbind, r2_list))

# Extract values properly (remove names)
r2_marginal <- as.numeric(r2_pooled[1])
r2_conditional <- as.numeric(r2_pooled[2])
r2_random <- r2_conditional - r2_marginal
r2_residual <- 1 - r2_conditional

# Create visualization data
variance_components <- data.frame(
  Component = c("Fixed Effects (Predictors)", 
                "Random Effects (Between-Person)", 
                "Residual (Unexplained)"),
  Variance = c(r2_marginal, r2_random, r2_residual) * 100
)

variance_components %>% kable
```

| Component                       | Variance |
|:--------------------------------|---------:|
| Fixed Effects (Predictors)      | 27.38005 |
| Random Effects (Between-Person) | 20.15473 |
| Residual (Unexplained)          | 52.46521 |

### Additional Diagnostics

We examine additional model fit diagnostics using one representative
imputed dataset (imputation 1).

``` r
# Extract first imputed dataset for diagnostic checking
model_data<-complete(imp_1_with_trends, 1)

# Clean data: drop NAs from lagging (matches lmer's internal dropping)
model_data_clean <- model_data %>% 
  filter(!is.na(lagged_evening_interference))

mod_pain_coupling_unc <- lmer(evening_interference ~ 1+ (1 | alias), data = model_data,REML = FALSE)



# model with 3d pain trend
mod_pain_coupling_0_one_sample <-  lmer(
  evening_interference ~ interf_worry + 
                         pain_trend_3d +
                         morning_mood + 
                         morning_stress +
                         day_count +
                         Age_c +
                        (1 | alias),  
  data = model_data_clean,
  REML = FALSE  # Important for model comparison

)
```

check_model() from the performance package provides the visual check of
various model assumptions: \* normality of residuals \* normality of
random effects \* linear relationship \* homogeneity of variance \*
multicollinearity

``` r
check_model(mod_pain_coupling_0_one_sample, check = c("all"))
```

![](C:/Users/d.vitali/Desktop/Github/EMA_Pain_Supl_2026/Docs/S4_Q1_model_documentation_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

#### Autocorrelation

- ACF plot aid visual examination of autocorrelation (lines should stay
  within bands)
- Box-Ljung test Model residuals approximate white noise (should not be
  significant).

``` r
residuals <- resid(mod_pain_coupling_0_one_sample)
############## AUTOCORRELATION diagnostics #####################
# Check model residuals' autocorrelation structure against white noise
# ACF estimates the autocovariance or autocorrelation function
# ACF is useful for identifying the order of the MA (Moving Average) component
acf(residuals, main = "ACF of Model Residuals (Imputation 1)")
```

![](C:/Users/d.vitali/Desktop/Github/EMA_Pain_Supl_2026/Docs/S4_Q1_model_documentation_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

``` r
# PACF estimates the partial autocorrelations. Same as acf but pacf controls for 
# for the influence of intermediate lags. this way focuses only on the effect of immediate lags
#pacf(residuals, main = "PACF of Model Residuals")

#
Box.test(residuals, lag = 10, type = "Ljung-Box")
```

    ## 
    ##  Box-Ljung test
    ## 
    ## data:  residuals
    ## X-squared = 14.844, df = 10, p-value = 0.1379

``` r
## plot residuals vs day_count to see autocorrelation trends within participants
ggplot(model_data_clean, aes(x = day_count, y = residuals, group = alias)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = F) +
  labs(x = "Day Count", y = "Residuals", title = "Residuals vs Day Count") +
  theme_minimal()
```

    ## `geom_smooth()` using formula = 'y ~ x'

![](C:/Users/d.vitali/Desktop/Github/EMA_Pain_Supl_2026/Docs/S4_Q1_model_documentation_files/figure-gfm/unnamed-chunk-12-2.png)<!-- -->

``` r
############## END AUTOCORRELATION #####################

# VIF Thresholds
#    VIF < 5: No multicollinearity concern
#    VIF 5-10: Moderate multicollinearity (often acceptable)
#    VIF > 10: Severe multicollinearity (problematic)

# VIFs (~6-7) are within acceptable bounds and much lower than the severe threshold.
# Some collinearity is expected for Interaction Models because interaction terms often 
# correlate with their component terms and deviation scores often correlate with each 
# other when calculated from related time series.
```

## Coefficient Visualization

We visualize the pooled coefficient estimates from Model 0 (baseline) to
compare the relative magnitude of predictors. Error bars represent 95%
confidence intervals.

``` r
pooled_summary <- summary(pool(model_Coupling_pooled_0), conf.int = TRUE)

# Prepare data from your pooled results
coef_data <- pooled_summary %>%
  filter(term != "(Intercept)") %>%
  mutate(
    # Create cleaner labels
    term_label = case_when(
      term == "interf_worry" ~ "Worry about\nInterference",
      term == "pain_trend_3d" ~ "Pain Trend\n(3-day)",
      term == "morning_mood" ~ "Morning\nMood",
      term == "morning_stress" ~ "Morning\nStress",
      term == "day_count" ~ "Study\nDay",
      term == "Age_c" ~ "Age (centred)",
      TRUE ~ term
    ),
    # Flag significance
    significant = ifelse(p.value < 0.05, "sig", "ns"),
    # Reorder by absolute effect size
    term_label = reorder(term_label, abs(estimate))
  )

# Plot
ggplot(coef_data, aes(x = estimate, y = term_label, color = significant)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", size = 0.8) +
  geom_errorbarh(aes(xmin = `2.5 %`, xmax = `97.5 %`), 
                 width = 0.2, size = 1) +
  geom_point(size = 4) +
  scale_color_manual(values = c("sig" = "#2C7BB6", "ns" = "#D3D3D3"),
                     guide = "none") +
  labs(x = "Coefficient (b) with 95% CI", 
       y = "",
       title = "Worry-Outcome Coupling: Daily Predictors of Interference",
       caption = "Morning worry and contextual factors predicting evening pain interference (pooled estimates, m = 28)")+ 
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none",
    panel.grid.major.y = element_blank(),
    plot.title = element_text(hjust = 2, face = "bold"),
    plot.caption = element_text(hjust = -2, size = 9, color = "gray40")
  )
```

![](C:/Users/d.vitali/Desktop/Github/EMA_Pain_Supl_2026/Docs/S4_Q1_model_documentation_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->
