Supplementary Material: MAR assumptions
================
Diego Vitali, Raffaele Buono, Tao Bi, Tony Pickering, Chris Eccleston,
Edmund Keogh, Nadia Bianchi-Berthouze, Amanda C de C Williams
2026-02-07

# SUMMARY & LOGIC:

This code tests whether missing data meets Missing At Random (MAR)
assumptions required for multiple imputation.

Strategy:

1.  Quantify overall missingness rates
2.  Test for systematic patterns:
    - Survey fatigue (missingness ~ study day)
    - Person effects (missingness ~ participant)
    - Correlated missingness (variables missing together)
    - Predictive relationships (observed predictors ~ missingness)\*

MAR requires: missingness can depend on observed data but not unobserved
values.

\*Tests use logistic regression to check if observed morning variables
predict evening outcome missingness, controlling for person and time
effects.

# Missingness Rate

![](C:/Users/d.vitali/Desktop/Github/EMA_Pain_Supl_2026/Docs/S1_MAR_explore_documentation_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->![](C:/Users/d.vitali/Desktop/Github/EMA_Pain_Supl_2026/Docs/S1_MAR_explore_documentation_files/figure-gfm/unnamed-chunk-2-2.png)<!-- -->

| n_total | n_complete | n_incomplete | pct_incomplete | pct_complete |
|--------:|-----------:|-------------:|---------------:|-------------:|
|     501 |        399 |          102 |       20.35928 |     79.64072 |

# Summary Table

``` r
# Variable-level missingness rates
missingness_summary <- data_set_b4mice %>%
  select(all_of(model_vars)) %>% 
  filter(day_count > 0) %>%
  summarise(across(everything(), ~mean(is.na(.)) * 100)) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "pct_missing") %>%
  arrange(desc(pct_missing))

missingness_summary %>% kable()
```

| variable             | pct_missing |
|:---------------------|------------:|
| evening_interference |   12.175649 |
| interf_worry         |   10.978044 |
| morning_mood         |    5.588822 |
| morning_pain         |    5.588822 |
| morning_stress       |    5.588822 |
| alias                |    0.000000 |
| daily_avg_stress     |    0.000000 |
| daily_avg_pain       |    0.000000 |
| daily_avg_mood       |    0.000000 |
| Age                  |    0.000000 |
| day_count            |    0.000000 |

## MAR Assumption Tests

### Missingness ~ time

Here we explore temporal dropout patterns. The expectation under MCAR is
that missingness rates will be stable over time. We examine this chart
looking for:

1.  Upward trends in % missing as study day increases (survey fatigue)
2.  Systematic differences between early vs late study periods

Under MCAR, we would expect to see no clear trends in either case, with
only random day-to-day fluctuation in missingness rates.

``` r
# Test 1: Survey fatigue (missingness increases over time?)
missingness_by_day <- data_set_b4mice %>%
  select(all_of(model_vars)) %>% 
  filter(day_count > 0) %>%
  group_by(day_count) %>%
  summarise(
    pct_missing_worry = mean(is.na(interf_worry)) * 100,
    pct_missing_interference = mean(is.na(evening_interference)) * 100,
    pct_missing_pain = mean(is.na(morning_pain)) * 100
  )

# Visualize temporal trends
ggplot(missingness_by_day, aes(x = day_count)) +
  geom_line(aes(y = pct_missing_worry, color = "Worry")) +
  geom_line(aes(y = pct_missing_interference, color = "Interference")) +
  geom_line(aes(y = pct_missing_pain, color = "Pain")) +
  labs(x = "Study Day", y = "% Missing", 
       title = "Missingness Over Time (Survey Fatigue Test)")
```

![](C:/Users/d.vitali/Desktop/Github/EMA_Pain_Supl_2026/Docs/S1_MAR_explore_documentation_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

### Missingness ~ Participant?

This code examines person-level heterogeneity in missingness patterns to
assess whether certain participants systematically provide less data.

Scatter plot shows relationship between worry and interference
missingness:

- Points on diagonal = participants miss both equally (whole survey
  sessions skipped)
- Points off diagonal = differential missingness (one variable more
  likely missing than the other)

under MCAR we would expect a random scatter around diagonal and no
systematic pattern of morning worry vs evening interference prevalence.

In the next steps glm() will test whether the levels of worry detected
in the morning predict the report of interference being missingness on
the same evening.

``` r
# Test 2: Person-level differences in compliance
missingness_by_person <- data_set_b4mice %>%
  select(all_of(model_vars)) %>% 
  filter(day_count > 0) %>%
  group_by(alias) %>%
  summarise(
    n_days = n(),
    pct_missing_worry = mean(is.na(interf_worry)) * 100,
    pct_missing_interference = mean(is.na(evening_interference)) * 100,
    pct_any_missing = mean(is.na(interf_worry) | is.na(evening_interference)) * 100
  ) %>%
  arrange(desc(pct_any_missing))

missingness_by_person %>% kable()
```

| alias | n_days | pct_missing_worry | pct_missing_interference | pct_any_missing |
|------:|-------:|------------------:|-------------------------:|----------------:|
|    15 |     28 |         32.142857 |                25.000000 |       46.428571 |
|    22 |     28 |         14.285714 |                28.571429 |       42.857143 |
|    18 |     27 |         11.111111 |                29.629630 |       40.740741 |
|    10 |     16 |         31.250000 |                25.000000 |       37.500000 |
|    14 |     29 |         17.241379 |                17.241379 |       31.034483 |
|    21 |     29 |         20.689655 |                13.793103 |       31.034483 |
|    23 |     29 |         24.137931 |                 3.448276 |       24.137931 |
|    19 |     28 |         10.714286 |                10.714286 |       21.428571 |
|     5 |     29 |         10.344828 |                13.793103 |       20.689655 |
|    11 |     29 |          3.448276 |                10.344828 |       13.793103 |
|    12 |     29 |          3.448276 |                10.344828 |       13.793103 |
|    25 |     29 |         13.793103 |                 3.448276 |       13.793103 |
|    24 |     28 |          3.571429 |                10.714286 |       10.714286 |
|    26 |     28 |          3.571429 |                 7.142857 |       10.714286 |
|     9 |     29 |          0.000000 |                 6.896552 |        6.896552 |
|    13 |     28 |          3.571429 |                 3.571429 |        3.571429 |
|    17 |     29 |          0.000000 |                 3.448276 |        3.448276 |
|    20 |     29 |          3.448276 |                 3.448276 |        3.448276 |

``` r
# Correlation between worry and interference missingness
ggplot(missingness_by_person, aes(x = pct_missing_worry, 
                                   y = pct_missing_interference)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_text(aes(label = alias), vjust = -0.5, size = 3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  labs(x = "% Missing Worry",
       y = "% Missing Interference",
       title = "Relationship Between Missingness Types (Person-Level)") +
  theme_minimal()
```

![](C:/Users/d.vitali/Desktop/Github/EMA_Pain_Supl_2026/Docs/S1_MAR_explore_documentation_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

### Missingness Co-occurrence?

This UpSet plot shows which combinations of variables tend to be missing
together. Understanding co-occurrence patterns helps distinguish between
item-level skipping (participants selectively skip certain questions)
versus survey-level dropout (participants miss entire survey sessions).

Note: Morning mood, pain, and stress are always asked in the same
morning survey, so these three variables will always be missing together
when the morning survey is skipped. Similarly, evening interference and
morning worry come from different survey sessions (evening and morning
respectively), so their co-occurrence patterns reflect between-survey
rather than within-survey missingness.

``` r
# Test 3: Which variables tend to be missing together?
# UpSet plot shows combinations of missing variables
gg_miss_upset(data_set_b4mice %>%
              select(all_of(model_vars)) %>%
              filter(day_count > 0))
```

![](C:/Users/d.vitali/Desktop/Github/EMA_Pain_Supl_2026/Docs/S1_MAR_explore_documentation_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

### Outcome Missing ~ Predictors

This section tests the core MAR assumption: whether observed variables
predict missingness patterns. We create binary indicators for whether
each key variable is missing, then use logistic regression to test if
observed morning states predict the probability of missing data.

**MAR vs MCAR distinction**: - MCAR requires missingness to be
completely independent of all variables (observed and unobserved) - MAR
allows missingness to depend on observed data, but not on the unobserved
values themselves

This means that if data is MAR missingness depends on observed data
(lagged/same-day observables), but not on unobserved/missing values
conditional on those observables.

``` r
# Create binary missingness indicators for logistic regression
Missngness_glm_DF <- data_set_b4mice %>%
              select(all_of(model_vars)) %>%
              filter(day_count > 0) %>%
              mutate(interference_missing = is.na(evening_interference),
                     worry_missing = is.na(interf_worry))
```

#### Interference Missing ~ Observed Predictors

\*Test:**Does morning worry predict evening interference being
missing?**

This tests whether morning worry (net of pain/mood/stress) predicts
evening interference missingness.

- Significant worry → supports MAR (high morning worry → skip evening →
  imputation conditions on observed worry and other morning states).
- Non-significant → supports MCAR.

``` r
# Test 4a: Do observed morning states predict evening interference missingness?
# MAR allows this; MCAR would require no relationship
missingness_model_1 <- glm(
  interference_missing ~ morning_pain + morning_stress + 
    morning_mood + interf_worry +day_count ,
  family = binomial,
  data = Missngness_glm_DF
)

summary(missingness_model_1)
```

    ## 
    ## Call:
    ## glm(formula = interference_missing ~ morning_pain + morning_stress + 
    ##     morning_mood + interf_worry + day_count, family = binomial, 
    ##     data = Missngness_glm_DF)
    ## 
    ## Coefficients:
    ##                 Estimate Std. Error z value Pr(>|z|)  
    ## (Intercept)    -0.757158   0.714497  -1.060   0.2893  
    ## morning_pain   -0.107025   0.070758  -1.513   0.1304  
    ## morning_stress -0.014053   0.080726  -0.174   0.8618  
    ## morning_mood   -0.176996   0.093328  -1.896   0.0579 .
    ## interf_worry    0.008462   0.084272   0.100   0.9200  
    ## day_count      -0.001811   0.018959  -0.096   0.9239  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 300.38  on 445  degrees of freedom
    ## Residual deviance: 295.65  on 440  degrees of freedom
    ##   (55 observations deleted due to missingness)
    ## AIC: 307.65
    ## 
    ## Number of Fisher Scoring iterations: 5

#### Worry Missing ~ Observed Predictors

**Do observed lagged states predict morning worry data being missing?**
This tests whether previous states predict morning survey non-completion
(which includes worry).

- If lagged worry/pain/mood/stress predict missingness → supports for
  MAR imputation can succesfully condition on them).
- If lagged worry/pain/mood/stress do not predict missingness → supports
  for MCAR.

``` r
# Test 4b: Do observed morning states predict worry missingness?
Missingness_model_2b <- glm(
  worry_missing ~ lag(interf_worry) + 
    lag(evening_interference) + 
    lag(daily_avg_pain) + lag(daily_avg_stress) + lag(daily_avg_mood) + 
    day_count,
  family = binomial,
  data = Missngness_glm_DF
)


summary(Missingness_model_2b)
```

    ## 
    ## Call:
    ## glm(formula = worry_missing ~ lag(interf_worry) + lag(evening_interference) + 
    ##     lag(daily_avg_pain) + lag(daily_avg_stress) + lag(daily_avg_mood) + 
    ##     day_count, family = binomial, data = Missngness_glm_DF)
    ## 
    ## Coefficients:
    ##                            Estimate Std. Error z value Pr(>|z|)   
    ## (Intercept)               -2.730791   1.011477  -2.700  0.00694 **
    ## lag(interf_worry)          0.302126   0.114042   2.649  0.00807 **
    ## lag(evening_interference) -0.066760   0.112616  -0.593  0.55331   
    ## lag(daily_avg_pain)       -0.235426   0.125667  -1.873  0.06101 . 
    ## lag(daily_avg_stress)     -0.007991   0.113240  -0.071  0.94374   
    ## lag(daily_avg_mood)       -0.057545   0.128280  -0.449  0.65373   
    ## day_count                  0.037016   0.023261   1.591  0.11154   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 227.52  on 397  degrees of freedom
    ## Residual deviance: 213.86  on 391  degrees of freedom
    ##   (103 observations deleted due to missingness)
    ## AIC: 227.86
    ## 
    ## Number of Fisher Scoring iterations: 6

Lagged prior-day worry significantly predicted higher morning worry
missingness (OR=1.35, p=.008), while prior pain marginally predicted
lower missingness (OR=0.79, p=.061). This supports MAR: imputation
conditions on these observed lagged states to adjust for systematic
patterns

## MAR-based sensitivity analysis

\[Text copied from Q2 supplementary to complement the exploration,
please see Q2 supplementary for more details\]

Our MAR diagnostics revealed that lagged worry (t-1) significantly
predicted current worry missingness ($\beta=0.30$, p=.008), suggesting
imputation models should ideally condition on this predictor. To assess
robustness to this MAR pattern, we compared Q2 model estimates under two
imputation strategies:

- **Model A**: Two-level PMM without lagged worry predictor (primary
  specification)
- **Model B**: Two-level PMM including lagged worry predictor
  (MAR-informed)

| term                       | est_A |  se_A |  lwr_A | upr_A |     | est_B |  se_B |  lwr_B | upr_B |
|:---------------------------|------:|------:|-------:|------:|:---:|------:|------:|-------:|------:|
| (Intercept)                |  1.37 | 0.862 | -0.325 | 3.065 |     |  1.35 | 0.935 | -0.487 | 3.188 |
| lagged_interf_worry        |  0.28 | 0.058 |  0.167 | 0.394 |     |  0.26 | 0.061 |  0.143 | 0.382 |
| lagged_worry_outcome_error | -0.08 | 0.044 | -0.169 | 0.005 |     | -0.08 | 0.046 | -0.172 | 0.009 |
| lagged_pain_trend_3d       | -0.08 | 0.076 | -0.233 | 0.064 |     | -0.11 | 0.081 | -0.269 | 0.048 |
| lagged_mood_trend_3d       |  0.18 | 0.107 | -0.031 | 0.390 |     |  0.19 | 0.118 | -0.040 | 0.425 |
| lagged_stress_trend_3d     |  0.19 | 0.071 |  0.051 | 0.329 |     |  0.20 | 0.076 |  0.047 | 0.345 |
| day_count                  |  0.01 | 0.009 | -0.009 | 0.024 |     |  0.00 | 0.010 | -0.015 | 0.022 |
| Age_c                      | -0.01 | 0.033 | -0.072 | 0.057 |     | -0.01 | 0.031 | -0.067 | 0.053 |

The pooled fixed‑effect estimates from the analysis model were highly
similar across the two imputation strategies (See table above). For all
parameters, point estimates, standard errors, and confidence intervals
overlapped closely. In particular, the lagged morning interference
coefficient remained positive and of comparable magnitude in both models
(Model A: 0.28, 95% CI \[0.17, 0.39\]; Model B: 0.26, 95% CI \[0.14,
0.38\]), and the pattern of results for the lagged symptom trends, day
count, and age was essentially unchanged.The maximum difference in point
estimates across any parameter was 0.03 (lagged pain trend),

representing \<5% relative change. Considering these findings, our
longitudinal inferences are appear robust to reasonable variation in the
imputation model for morning worry nd we therefore selected the more
parsimonious and numerically stable imputation model without the lagged
term as our primary specification.
