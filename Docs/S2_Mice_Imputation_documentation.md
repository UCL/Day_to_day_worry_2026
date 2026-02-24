Supplementary Material: Imputation Model
================
Diego Vitali, Raffaele Buono, Tao Bi, Tony Pickering, Chris Eccleston,
Edmund Keogh, Nadia Bianchi-Berthouze, Amanda C de C Williams
2026-02-07

# Summary & Logic

This document specifies the multiple imputation models used to handle
missing data in our daily diary study of pain interference predictions.
We employed two-level predictive mean matching (2l.pmm) via the mice
package in R, which preserves the multilevel structure of repeated daily
assessments nested within participants [Wijesuriya et al
2025](https://doi.org/10.1002/sim.10274).

Our imputation strategy was informed by Missing At Random (MAR)
diagnostics presented in the companion supplementary material
(MAR_explore_documentation), which identified systematic patterns in
missingness that we address through careful predictor selection.

We present two imputation specifications:

1.  Primary model (imp_1): Standard specification based on theoretical
    relationships and general imputation practice
2.  Sensitivity model (imp_2): Adds lagged predictors to address
    specific MAR patterns identified in diagnostics

This leads to two predictor matrices differing **only** by the inclusion
of `lagged_interf_worry` as a predictor of morning worry. We call the
reference predictor matrix “pred” and the alternative predictor matrix
with addition of `lagged_interf_worry` named “pred_b”.

# Imputation

The reference MICE model runs with “pred” and attain model `imp_1`.

## imputation 1

``` r
set.seed(200720)
imp_1 = mice(data = data_set_b4mice, 
            m = 28, 
            predictorMatrix = pred, 
            method = meth1, 
            maxit = 12)

#saveRDS(imp_1, "c:/Documents and Settings/d.vitali/Desktop/Github/CRIISP-WP6/Data/imputation_m28_12_05_25.rds")
```

The alternative reference MICE model runs with “pred_b” and attain model
`imp_1`. All other variables use identical predictor specifications as
the primary model.

By including lagged worry as a predictor of current worry, we condition
the imputation model on a variable that explains missingness patterns,
potentially reducing bias under MAR assumptions. The comparison between
models assesses how robust are our downstream findings (Q2) to this
specification choice.

## imputation 2

``` r
pred_b<-pred

pred_b["interf_worry", c("day_count", 
                       "lagged_interf_worry",
                       "morning_pain", "morning_mood", "morning_stress", 
                       "daily_avg_interference", "Age")] <- 1


diff_matrix <- (pred_b != pred) 

diff_pos <- which(diff_matrix, arr.ind = TRUE)


# Show how many differences
n_diff <- nrow(diff_pos)

# Show positions (with dimnames)
rownames_diff <- rownames(pred)[diff_pos[,"row"]]
colnames_diff <- colnames(pred)[diff_pos[,"col"]]

data.frame(
  row = rownames_diff,
  col = colnames_diff,
  pred_a = pred[diff_pos],
  pred_b = pred_b[diff_pos]
)


imp_2 = mice(data = data_set_b4mice, 
            m = 28, 
            predictorMatrix = pred_b, 
            method = meth1, 
            maxit = 12)

#saveRDS(imp_2, "c:/Documents and Settings/d.vitali/Desktop/Github/CRIISP-WP6/Data/imputation_m28_2026_02_09.rds")
```

## IMP_1 Convergence plots

We assessed convergence via trace plots of means/standard deviations
across 28 chains and iterations, plus density plots of imputed
vs. observed values [van Buuren & Groothuis-Oudshoorn
2011](https://doi.org/10.18637/jss.v045.i03). Acceptable convergence
requires:

- random chain mixing: We assess multiple MCMC chains to make sure that
  they explore the parameter space freely, that they interweave without
  separation or sticking, and we check whether the sampler gets
  “trapped” in local modes
- parameter stabilization: We expectet that means and standard
  deviations of imputed values plateau after initial burn-in iterations.
  This progressive “convergence” shows that the chains have forgotten
  initial values and reach sowewhat stationary posterior distributions.
- plausible imputed distributions: Densityplots should confirm that
  imputed values (red) have similar shape, spread, and range to observed
  data (blue). This means that the original features of the data like
  (e.g) the skewness in pain interference is preserved in the imputation
  datasets.

To extend this diagnostic we load the imputation with extendeded 20
iteration
![](C:/Users/d.vitali/Desktop/Github/EMA_Pain_Supl_2026/Docs/S2_Mice_Imputation_documentation_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->![](C:/Users/d.vitali/Desktop/Github/EMA_Pain_Supl_2026/Docs/S2_Mice_Imputation_documentation_files/figure-gfm/unnamed-chunk-5-2.png)<!-- -->

Trace plots show random intermingling with no trends/divergence, and
stabilization by iterations 4-6

In addition we examine ([van Buuren & Groothuis-Oudshoorn
2011](https://doi.org/10.18637/jss.v045.i03)):

- densityplots: it is useful to show differences in kernel density
  estimates of the imputed and observed data
- bwplot: is useful to study the distributions in datasets. Gives a
  side-by-side box-and-whisker plots for the observed and synthetic data

``` r
densityplot(imp_1_ext, ~ interf_worry + evening_interference + morning_pain + morning_mood + morning_stress)
```

![](C:/Users/d.vitali/Desktop/Github/EMA_Pain_Supl_2026/Docs/S2_Mice_Imputation_documentation_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

Density plots showed imputed values (red) overlapping observed data and
preserving skewness in pain interference

``` r
bwplot(imp_1, interf_worry + evening_interference + morning_pain + morning_mood + morning_stress ~ .imp, pch = 20, cex = 1.2)
```

![](C:/Users/d.vitali/Desktop/Github/EMA_Pain_Supl_2026/Docs/S2_Mice_Imputation_documentation_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

bwplot across all 28 imputations confirmed that imputed values fell
within plausible ranges for both interf_worry and evening_interference
(0-10 scale). Imputed values showed appropriate dispersion across the
scale without clustering at extremes or producing impossible values. The
evident overlap of red and blue boxplots indicates that imputed values
matched the distributional characteristics of observed data.

### Additional MAR Explore

From ([van Buuren & Groothuis-Oudshoorn
2011](https://doi.org/10.18637/jss.v045.i03))

Bondarenko and Raghunathan
[2016](https://stefvanbuuren.name/fimd/references.html#ref-BONDARENKO2016)
proposed a more refined diagnostic tool that aims to compare the
distributions of observed and imputed data conditional on the
missingness probability. The idea is that under MAR the conditional
distributions should be similar if the assumed model for creating
multiple imputations has a good fit.

If MAR assumption holds:

- Imputed values (filled circles) distributed similarly to observed
  values (open circles) across all propensity score levels
- No systematic separation between observed and imputed

If MAR does not hold: \* Imputed values cluster at different BMI levels
than observed values with similar missingness probability \* Gap between
filled and open circles suggests unobserved confounding

``` r
# Predict probability that record has missing outcome (evening_interference)
# using observed morning predictors
fit <- with(imp_1, 
            glm(ici(imp_1) ~ interf_worry + morning_pain + 
                             morning_mood + morning_stress + day_count,
                family = binomial))

# Calculate propensity scores
ps <- rep(rowMeans(sapply(fit$analyses, fitted.values)), imp_1$m + 1)

# Plot evening_interference vs propensity score
xyplot(imp_1, evening_interference ~ ps | as.factor(.imp),
       xlab = "P(incomplete | observed morning variables)",
       ylab = "Evening Interference",
       pch = c(1, 19),  # 1=observed, 19=imputed
       col = mdc(1:2))
```

![](C:/Users/d.vitali/Desktop/Github/EMA_Pain_Supl_2026/Docs/S2_Mice_Imputation_documentation_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

## IMP_2 convergence

Imp_2 differs from imp_1 by exactly ONE predictor (lagged_interf_worry
for interf_worry). This is a narrow sensitivity check for later models,
it is not a fundamentally different model and it does not require the
same convergence testing.
