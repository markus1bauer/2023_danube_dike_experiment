Analysis of Bauer et al. (unpublished) Field experiment: <br> Favourable
Conservation Status (FCS)
================
<b>Markus Bauer\*</b> <br>
<b>2022-12-05</b>

- <a href="#preparation" id="toc-preparation">Preparation</a>
- <a href="#statistics" id="toc-statistics">Statistics</a>
  - <a href="#data-exploration" id="toc-data-exploration">Data
    exploration</a>
    - <a href="#graphs-of-raw-data" id="toc-graphs-of-raw-data">Graphs of raw
      data</a>
    - <a href="#outliers-zero-inflation-transformations"
      id="toc-outliers-zero-inflation-transformations">Outliers,
      zero-inflation, transformations?</a>
  - <a href="#models" id="toc-models">Models</a>
    - <a href="#priors" id="toc-priors">Priors</a>
  - <a href="#model-check" id="toc-model-check">Model check</a>
    - <a href="#dharma" id="toc-dharma">DHARMa</a>
    - <a href="#preparation-1" id="toc-preparation-1">Preparation</a>
    - <a href="#sampling-efficency-and-effectiveness"
      id="toc-sampling-efficency-and-effectiveness">Sampling efficency and
      effectiveness</a>
    - <a href="#mcmc-diagnostics" id="toc-mcmc-diagnostics">MCMC
      diagnostics</a>
    - <a href="#posterior-predictive-check"
      id="toc-posterior-predictive-check">Posterior predictive check</a>
    - <a href="#autocorrelation-check"
      id="toc-autocorrelation-check">Autocorrelation check</a>
  - <a href="#output-of-choosen-model"
    id="toc-output-of-choosen-model">Output of choosen model</a>
    - <a href="#model-output" id="toc-model-output">Model output</a>
    - <a href="#effect-sizes" id="toc-effect-sizes">Effect sizes</a>
- <a href="#session-info" id="toc-session-info">Session info</a>

Technichal University of Munich, TUM School of Life Sciences, Chair of
Restoration Ecology, Emil-Ramann-Straße 6, 85354 Freising, Germany

\* Corresponding author: <markus1.bauer@tum.de> <br> ORCiD ID:
[0000-0001-5372-4174](https://orcid.org/0000-0001-5372-4174) <br>
[Google
Scholar](https://scholar.google.de/citations?user=oHhmOkkAAAAJ&hl=de&oi=ao)
<br> GitHub: [markus1bauer](https://github.com/markus1bauer)

Favourable Conservation Status (FSC) sensu Helm et al. (2015) Divers
Distrib [DOI: 10.1111/ddi.12285](https://doi.org/10.1111/ddi.12285)

# Preparation

#### Packages

``` r
library(here)
library(tidyverse)
library(ggbeeswarm)
library(patchwork)
library(brms)
library(DHARMa)
library(DHARMa.helpers)
library(bayesplot)
library(loo)
library(tidybayes)
library(emmeans)
```

#### Load data

``` r
sites <- read_csv(
  here("data", "processed", "data_processed_sites.csv"),
  col_names = TRUE, na = c("na", "NA", ""), col_types =
    cols(
      .default = "?",
      plot = "f",
      site = "f",
      sand_ratio = "f",
      substrate_depth = col_factor(levels = c("30", "15")),
      target_type = col_factor(levels = c(
        "hay_meadow", "dry_grassland"
        )),
      seed_density = "f",
      exposition = col_factor(levels = c("north", "south")),
      survey_year = "c"
      )
  ) %>%
  ### Exclude data of seed mixtures
  filter(survey_year != "seeded") %>%
  mutate(
    survey_year_fct = factor(survey_year),
    botanist_year = str_c(survey_year, botanist, sep = " "),
    botanist_year = factor(botanist_year),
    n = fcs_target,
    id = factor(id)
    ) %>%
  select(
    id, plot, site, exposition, sand_ratio, substrate_depth, target_type,
    seed_density, survey_year_fct, survey_year, botanist_year, n
    )
```

# Statistics

## Data exploration

### Graphs of raw data

![](model_fcs_check_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->![](model_fcs_check_files/figure-gfm/unnamed-chunk-4-2.png)<!-- -->![](model_fcs_check_files/figure-gfm/unnamed-chunk-4-3.png)<!-- -->![](model_fcs_check_files/figure-gfm/unnamed-chunk-4-4.png)<!-- -->![](model_fcs_check_files/figure-gfm/unnamed-chunk-4-5.png)<!-- -->

### Outliers, zero-inflation, transformations?

    ## # A tibble: 12 × 3
    ##    exposition site      n
    ##    <fct>      <fct> <int>
    ##  1 north      1        96
    ##  2 north      2        96
    ##  3 north      3        96
    ##  4 north      4        96
    ##  5 north      5        96
    ##  6 north      6        96
    ##  7 south      1        96
    ##  8 south      2        96
    ##  9 south      3        96
    ## 10 south      4        96
    ## 11 south      5        96
    ## 12 south      6        96

![](model_fcs_check_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->![](model_fcs_check_files/figure-gfm/unnamed-chunk-5-2.png)<!-- -->![](model_fcs_check_files/figure-gfm/unnamed-chunk-5-3.png)<!-- -->![](model_fcs_check_files/figure-gfm/unnamed-chunk-5-4.png)<!-- -->

## Models

``` r
load(file = here("outputs", "models", "model_fcs_2.Rdata"))
load(file = here("outputs", "models", "model_fcs_2_new.Rdata"))
## Registered S3 methods overwritten by 'adegraphics':
##   method         from
##   biplot.dudi    ade4
##   kplot.foucart  ade4
##   kplot.mcoa     ade4
##   kplot.mfa      ade4
##   kplot.pta      ade4
##   kplot.sepan    ade4
##   kplot.statis   ade4
##   scatter.coa    ade4
##   scatter.dudi   ade4
##   scatter.nipals ade4
##   scatter.pco    ade4
##   score.acm      ade4
##   score.mix      ade4
##   score.pca      ade4
##   screeplot.dudi ade4
#load(file = here("outputs", "models", "model_fcs_full.Rdata"))
m_1 <- m2
m_2 <- m2_new
#m_2 <- m_full
```

``` r
m_1$formula
## n ~ sand_ratio * target_type * exposition * survey_year_fct + substrate_depth + seed_density + substrate_depth:exposition + seed_density:exposition + substrate_depth:survey_year_fct + seed_density:survey_year_fct + botanist_year + (1 | site/plot)
m_2$formula
## n ~ sand_ratio * target_type * exposition * survey_year_fct + substrate_depth + seed_density + substrate_depth:exposition + seed_density:exposition + substrate_depth:survey_year_fct + seed_density:survey_year_fct + botanist_year + (1 | site/plot)
```

``` r
m_1$family
## 
## Family: gaussian 
## Link function: identity
m_2$family
## 
## Family: gaussian 
## Link function: identity
```

### Priors

#### Possible prior distributions

``` r
ggplot(data = data.frame(x = c(-1, 1)), aes(x = x)) +
  stat_function(
    fun = dnorm, n = 101, args = list(mean = 0.1, sd = 1)
    ) +
  expand_limits(y = 0) + ggtitle("Normal distribution")
```

![](model_fcs_check_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

``` r
ggplot(data = data.frame(x = c(-1, 1)), aes(x = x)) +
  stat_function(
    fun = dcauchy, n = 101, args = list(location = 0, scale = 1)
    ) +
  expand_limits(y = 0) + ggtitle("Cauchy distribution")
```

![](model_fcs_check_files/figure-gfm/unnamed-chunk-9-2.png)<!-- -->

``` r
ggplot(data.frame(x = c(-1, 1)), aes(x = x)) +
  stat_function(
    fun = dstudent_t, args = list(df = 3, mu = 0, sigma = 2.5)
    ) +
  expand_limits(y = 0) + ggtitle(expression(Student~italic(t)*"-distribution"))
```

![](model_fcs_check_files/figure-gfm/unnamed-chunk-9-3.png)<!-- -->

#### Prior summary

``` r
prior_summary(m_1, all = FALSE)
```

    ## # A tibble: 9 × 10
    ##   prior                class    coef  group resp  dpar  nlpar lb    ub    source
    ##   <chr>                <chr>    <chr> <chr> <chr> <chr> <chr> <chr> <chr> <chr> 
    ## 1 normal(0, 1)         b        ""    ""    ""    ""    ""    ""    ""    user  
    ## 2 normal(0.1, 1)       b        "san… ""    ""    ""    ""    ""    ""    user  
    ## 3 normal(0.2, 1)       b        "san… ""    ""    ""    ""    ""    ""    user  
    ## 4 normal(0.1, 1)       b        "sur… ""    ""    ""    ""    ""    ""    user  
    ## 5 normal(0.2, 1)       b        "sur… ""    ""    ""    ""    ""    ""    user  
    ## 6 normal(0.3, 1)       b        "sur… ""    ""    ""    ""    ""    ""    user  
    ## 7 normal(0, 1)         Interce… ""    ""    ""    ""    ""    ""    ""    user  
    ## 8 student_t(3, 0, 2.5) sd       ""    ""    ""    ""    ""    "0"   ""    defau…
    ## 9 cauchy(0, 1)         sigma    ""    ""    ""    ""    ""    "0"   ""    user

Conditional <i>R</i>² values

``` r
bayes_R2(m_1, probs = c(0.05, 0.5, 0.95),
         re_formula =  ~ (1 | site/plot) + (1 | botanist_year)) 
##     Estimate   Est.Error        Q5       Q50       Q95
## R2 0.8462899 0.005097922 0.8375052 0.8465433 0.8542504
bayes_R2(m_2, probs = c(0.05, 0.5, 0.95),
         re_formula =  ~ (1 | site/plot) + (1 | botanist_year))
##    Estimate   Est.Error        Q5       Q50       Q95
## R2 0.846296 0.005195245 0.8373229 0.8465799 0.8544651
```

Marginal <i>R</i>² values

``` r
bayes_R2(m_1, probs = c(0.05, 0.5, 0.95),
         re_formula = 1 ~ 1)
##     Estimate   Est.Error        Q5       Q50       Q95
## R2 0.8080202 0.004573187 0.8002889 0.8082438 0.8151444
bayes_R2(m_2, probs = c(0.05, 0.5, 0.95),
         re_formula = 1 ~ 1)
##     Estimate   Est.Error       Q5       Q50       Q95
## R2 0.8082861 0.004599353 0.800521 0.8085158 0.8154934
```

## Model check

### DHARMa

``` r
DHARMa.helpers::dh_check_brms(m_1, integer = TRUE)
```

![](model_fcs_check_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

``` r
DHARMa.helpers::dh_check_brms(m_2, integer = TRUE)
```

![](model_fcs_check_files/figure-gfm/unnamed-chunk-13-2.png)<!-- -->

### Preparation

``` r
posterior1 <- m_1 %>%
  posterior::as_draws() %>%
  posterior::subset_draws(
    variable = c(
      "b_sand_ratio25",
      "b_sand_ratio50",
      "b_substrate_depth15",
      "b_target_typedry_grassland",
      "b_seed_density8",
      "b_expositionsouth",
      "b_survey_year_fct2019",
      "b_survey_year_fct2020",
      "b_survey_year_fct2021",
      "sd_site__Intercept",
      "sd_site:plot__Intercept",
      "sigma"
    )
  )
posterior2 <- m_2 %>%
  posterior::as_draws() %>%
  posterior::subset_draws(
    variable = c(
      "b_sand_ratio25",
      "b_sand_ratio50",
      "b_substrate_depth15",
      "b_target_typedry_grassland",
      "b_seed_density8",
      "b_expositionsouth",
      "b_survey_year_fct2019",
      "b_survey_year_fct2020",
      "b_survey_year_fct2021",
      "sd_site__Intercept",
      "sd_site:plot__Intercept",
      "sigma"
    )
  )
hmc_diagnostics1 <- nuts_params(m_1)
hmc_diagnostics2 <- nuts_params(m_2)
y <- sites$n
yrep1 <- posterior_predict(m_1, draws = 500)
yrep2 <- posterior_predict(m_2, draws = 500)
loo1 <- loo(m_1, save_psis = TRUE, moment_match = FALSE)
```

    ## Warning: Found 4 observations with a pareto_k > 0.7 in model 'm_1'. It is
    ## recommended to set 'moment_match = TRUE' in order to perform moment matching for
    ## problematic observations.

``` r
loo2 <- loo(m_2, save_psis = TRUE, moment_match = FALSE)
```

    ## Warning: Found 1 observations with a pareto_k > 0.7 in model 'm_2'. It is
    ## recommended to set 'moment_match = TRUE' in order to perform moment matching for
    ## problematic observations.

``` r
draws1 <- m_1 %>%
  posterior::as_draws() %>%
  posterior::summarize_draws() %>%
  filter(str_starts(variable, "b_"))
draws2 <- m_2 %>%
  posterior::as_draws() %>%
  posterior::summarize_draws() %>%
  filter(str_starts(variable, "b_"))
```

### Sampling efficency and effectiveness

#### Rhat

``` r
mcmc_rhat(draws1$rhat)
```

![](model_fcs_check_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

``` r
mcmc_rhat(draws2$rhat)
```

![](model_fcs_check_files/figure-gfm/unnamed-chunk-15-2.png)<!-- -->

#### Effective sampling size (ESS)

``` r
mcmc_neff(neff_ratio(m_1))
```

![](model_fcs_check_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

``` r
mcmc_neff(neff_ratio(m_2))
```

![](model_fcs_check_files/figure-gfm/unnamed-chunk-16-2.png)<!-- -->

### MCMC diagnostics

``` r
mcmc_trace(posterior1, np = hmc_diagnostics1)
```

    ## No divergences to plot.

![](model_fcs_check_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

``` r
mcmc_trace(posterior2, np = hmc_diagnostics2)
```

    ## No divergences to plot.

![](model_fcs_check_files/figure-gfm/unnamed-chunk-17-2.png)<!-- -->

``` r
mcmc_pairs(m_1, off_diag_args = list(size = 1.2),
           pars = c(
             "b_sand_ratio25", "b_sand_ratio50", "b_substrate_depth15",
             "b_target_typedry_grassland", "b_seed_density8",
             "b_expositionsouth", "sigma"
           ))
```

![](model_fcs_check_files/figure-gfm/unnamed-chunk-17-3.png)<!-- -->

``` r
mcmc_pairs(m_2, off_diag_args = list(size = 1.2),
           pars = c(
             "b_sand_ratio25", "b_sand_ratio50", "b_substrate_depth15",
             "b_target_typedry_grassland", "b_seed_density8",
             "b_expositionsouth", "sigma"
           ))
```

![](model_fcs_check_files/figure-gfm/unnamed-chunk-17-4.png)<!-- -->

``` r
mcmc_parcoord(posterior1, np = hmc_diagnostics1)
```

![](model_fcs_check_files/figure-gfm/unnamed-chunk-17-5.png)<!-- -->

``` r
mcmc_parcoord(posterior2, np = hmc_diagnostics2)
```

![](model_fcs_check_files/figure-gfm/unnamed-chunk-17-6.png)<!-- -->

### Posterior predictive check

#### Kernel density

``` r
p1 <- ppc_dens_overlay(y, yrep1[1:50, ])
p2 <- ppc_dens_overlay(y, yrep2[1:50, ])
p1 / p2
```

![](model_fcs_check_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

``` r
ppc_dens_overlay_grouped(y, yrep1[1:50, ], group = sites$site)
```

![](model_fcs_check_files/figure-gfm/unnamed-chunk-18-2.png)<!-- -->

``` r
ppc_dens_overlay_grouped(y, yrep2[1:50, ], group = sites$site)
```

![](model_fcs_check_files/figure-gfm/unnamed-chunk-18-3.png)<!-- -->

``` r
p1 <- ppc_dens_overlay_grouped(y, yrep1[1:50, ], group = sites$exposition)
p2 <- ppc_dens_overlay_grouped(y, yrep2[1:50, ], group = sites$exposition)
p1 / p2
```

![](model_fcs_check_files/figure-gfm/unnamed-chunk-18-4.png)<!-- -->

``` r
ppc_dens_overlay_grouped(y, yrep1[1:50, ], group = sites$survey_year_fct)
```

![](model_fcs_check_files/figure-gfm/unnamed-chunk-18-5.png)<!-- -->

``` r
ppc_dens_overlay_grouped(y, yrep2[1:50, ], group = sites$survey_year_fct)
```

![](model_fcs_check_files/figure-gfm/unnamed-chunk-18-6.png)<!-- -->

``` r
p1 <- ppc_dens_overlay_grouped(y, yrep1[1:50, ], group = sites$target_type)
p2 <- ppc_dens_overlay_grouped(y, yrep2[1:50, ], group = sites$target_type)
p1 / p2
```

![](model_fcs_check_files/figure-gfm/unnamed-chunk-18-7.png)<!-- -->

``` r
p1 <- ppc_dens_overlay_grouped(y, yrep1[1:50, ], group = sites$seed_density)
p2 <- ppc_dens_overlay_grouped(y, yrep2[1:50, ], group = sites$seed_density)
p1 / p2
```

![](model_fcs_check_files/figure-gfm/unnamed-chunk-18-8.png)<!-- -->

``` r
p1 <- ppc_dens_overlay_grouped(y, yrep1[1:50, ], group = sites$sand_ratio)
p2 <- ppc_dens_overlay_grouped(y, yrep2[1:50, ], group = sites$sand_ratio)
p1 / p2
```

![](model_fcs_check_files/figure-gfm/unnamed-chunk-18-9.png)<!-- -->

``` r
p1 <- ppc_dens_overlay_grouped(y, yrep1[1:50, ], group = sites$substrate_depth)
p2 <- ppc_dens_overlay_grouped(y, yrep2[1:50, ], group = sites$substrate_depth)
p1 / p2
```

![](model_fcs_check_files/figure-gfm/unnamed-chunk-18-10.png)<!-- -->

#### Histograms of statistics skew

``` r
p1 <- ppc_stat(y, yrep1, binwidth = 0.001)
p2 <- ppc_stat(y, yrep2, binwidth = 0.001)
p1 / p2
```

![](model_fcs_check_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

``` r
ppc_stat_grouped(y, yrep1, group = sites$site, binwidth = 0.001)
```

![](model_fcs_check_files/figure-gfm/unnamed-chunk-19-2.png)<!-- -->

``` r
ppc_stat_grouped(y, yrep2, group = sites$site, binwidth = 0.001)
```

![](model_fcs_check_files/figure-gfm/unnamed-chunk-19-3.png)<!-- -->

``` r
p1 <- ppc_stat_grouped(y, yrep1, group = sites$exposition, binwidth = 0.001)
p2 <- ppc_stat_grouped(y, yrep2, group = sites$exposition, binwidth = 0.001)
p1 / p2
```

![](model_fcs_check_files/figure-gfm/unnamed-chunk-19-4.png)<!-- -->

``` r
ppc_stat_grouped(y, yrep1, group = sites$survey_year_fct, binwidth = 0.001)
```

![](model_fcs_check_files/figure-gfm/unnamed-chunk-19-5.png)<!-- -->

``` r
ppc_stat_grouped(y, yrep2, group = sites$survey_year_fct, binwidth = 0.001)
```

![](model_fcs_check_files/figure-gfm/unnamed-chunk-19-6.png)<!-- -->

``` r
p1 <- ppc_stat_grouped(y, yrep1, group = sites$target_type, binwidth = 0.001)
p2 <- ppc_stat_grouped(y, yrep2, group = sites$target_type, binwidth = 0.001)
p1 / p2
```

![](model_fcs_check_files/figure-gfm/unnamed-chunk-19-7.png)<!-- -->

``` r
p1 <- ppc_stat_grouped(y, yrep1, group = sites$seed_density, binwidth = 0.001)
p2 <- ppc_stat_grouped(y, yrep2, group = sites$seed_density, binwidth = 0.001)
p1 / p2
```

![](model_fcs_check_files/figure-gfm/unnamed-chunk-19-8.png)<!-- -->

``` r
p1 <- ppc_stat_grouped(y, yrep1, group = sites$sand_ratio, binwidth = 0.001)
p2 <- ppc_stat_grouped(y, yrep2, group = sites$sand_ratio, binwidth = 0.001)
p1 / p2
```

![](model_fcs_check_files/figure-gfm/unnamed-chunk-19-9.png)<!-- -->

``` r
p1 <- ppc_stat_grouped(y, yrep1, group = sites$substrate_depth, binwidth = 0.001)
p2 <- ppc_stat_grouped(y, yrep2, group = sites$substrate_depth, binwidth = 0.001)
p1 / p2
```

![](model_fcs_check_files/figure-gfm/unnamed-chunk-19-10.png)<!-- -->

#### LOO (Leave one out)

``` r
loo1
```

    ## 
    ## Computed from 10000 by 1152 log-likelihood matrix
    ## 
    ##          Estimate   SE
    ## elpd_loo   -464.7 30.4
    ## p_loo       182.4  8.8
    ## looic       929.4 60.8
    ## ------
    ## Monte Carlo SE of elpd_loo is NA.
    ## 
    ## Pareto k diagnostic values:
    ##                          Count Pct.    Min. n_eff
    ## (-Inf, 0.5]   (good)     1137  98.7%   762       
    ##  (0.5, 0.7]   (ok)         11   1.0%   451       
    ##    (0.7, 1]   (bad)         4   0.3%   168       
    ##    (1, Inf)   (very bad)    0   0.0%   <NA>      
    ## See help('pareto-k-diagnostic') for details.

``` r
loo2
```

    ## 
    ## Computed from 10000 by 1152 log-likelihood matrix
    ## 
    ##          Estimate   SE
    ## elpd_loo   -466.4 30.3
    ## p_loo       184.6  8.9
    ## looic       932.7 60.6
    ## ------
    ## Monte Carlo SE of elpd_loo is NA.
    ## 
    ## Pareto k diagnostic values:
    ##                          Count Pct.    Min. n_eff
    ## (-Inf, 0.5]   (good)     1132  98.3%   789       
    ##  (0.5, 0.7]   (ok)         19   1.6%   204       
    ##    (0.7, 1]   (bad)         1   0.1%   399       
    ##    (1, Inf)   (very bad)    0   0.0%   <NA>      
    ## See help('pareto-k-diagnostic') for details.

``` r
plot(loo1)
```

![](model_fcs_check_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

``` r
plot(loo2)
```

![](model_fcs_check_files/figure-gfm/unnamed-chunk-20-2.png)<!-- -->

Leave one out probability integral transform

``` r
p1 <- ppc_loo_pit_overlay(y, yrep1, lw = weights(loo1$psis_object))
```

    ## NOTE: The kernel density estimate assumes continuous observations and is not optimal for discrete observations.

``` r
p2 <- ppc_loo_pit_overlay(y, yrep2, lw = weights(loo2$psis_object))
```

    ## NOTE: The kernel density estimate assumes continuous observations and is not optimal for discrete observations.

``` r
p1 / p2
```

![](model_fcs_check_files/figure-gfm/unnamed-chunk-21-1.png)<!-- -->

### Autocorrelation check

``` r
mcmc_acf(posterior1, lags = 10)
```

![](model_fcs_check_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->

``` r
mcmc_acf(posterior2, lags = 10)
```

![](model_fcs_check_files/figure-gfm/unnamed-chunk-22-2.png)<!-- -->

## Output of choosen model

### Model output

Conditional and marignal <i>R</i>²

``` r
bayes_R2(m_1, probs = c(0.05, 0.5, 0.95),
         re_formula =  ~ (1 | site/plot) + (1 | botanist_year)) 
##     Estimate   Est.Error        Q5       Q50       Q95
## R2 0.8462899 0.005097922 0.8375052 0.8465433 0.8542504
bayes_R2(m_1, probs = c(0.05, 0.5, 0.95),
         re_formula = 1 ~ 1)
##     Estimate   Est.Error        Q5       Q50       Q95
## R2 0.8080202 0.004573187 0.8002889 0.8082438 0.8151444
```

Posteriors of chosen model

``` r
draws1
```

    ## # A tibble: 65 × 10
    ##    variable           mean   median     sd    mad      q5      q95  rhat ess_b…¹
    ##    <chr>             <dbl>    <dbl>  <dbl>  <dbl>   <dbl>    <dbl> <dbl>   <dbl>
    ##  1 b_Intercept    -0.678   -0.678   0.0907 0.0895 -0.829  -5.30e-1  1.00   6245.
    ##  2 b_sand_ratio25  0.151    0.151   0.0984 0.0976 -0.0100  3.15e-1  1.00   5437.
    ##  3 b_sand_ratio50  0.134    0.135   0.0984 0.0989 -0.0287  2.93e-1  1.00   5723.
    ##  4 b_target_type…  0.00175  0.00156 0.0965 0.0965 -0.157   1.62e-1  1.00   4965.
    ##  5 b_expositions… -0.654   -0.653   0.105  0.105  -0.827  -4.84e-1  1.00   5204.
    ##  6 b_survey_year…  0.816    0.809   0.605  0.613  -0.179   1.81e+0  1.00   9633.
    ##  7 b_survey_year…  1.02     1.02    0.588  0.591   0.0429  1.99e+0  1.00   9892.
    ##  8 b_survey_year…  0.918    0.924   0.713  0.712  -0.245   2.09e+0  1.00   9129.
    ##  9 b_substrate_d… -0.0828  -0.0827  0.0509 0.0508 -0.166   9.93e-4  1.00   9159.
    ## 10 b_seed_densit…  0.108    0.107   0.0510 0.0511  0.0260  1.93e-1  1.00   9632.
    ## # … with 55 more rows, 1 more variable: ess_tail <dbl>, and abbreviated
    ## #   variable name ¹​ess_bulk

``` r
mcmc_intervals(
  posterior1,
  prob = 0.66,
  prob_outer = 0.95,
  point_est = "mean"
  ) +
  theme_classic()
```

![](model_fcs_check_files/figure-gfm/unnamed-chunk-24-1.png)<!-- -->

Posteriors of second model:

``` r
mcmc_intervals(
  posterior2,
  prob = 0.66,
  prob_outer = 0.95,
  point_est = "mean"
  ) +
  theme_classic()
```

![](model_fcs_check_files/figure-gfm/unnamed-chunk-25-1.png)<!-- -->

### Effect sizes

Just, to get exact values if necessary, which is not possible from the
figure

``` r
(emm <- emmeans(m_1, revpairwise ~ target_type + sand_ratio |
                  exposition | survey_year_fct, type = "response"))
```

    ## NOTE: A nesting structure was detected in the fitted model:
    ##     botanist_year %in% survey_year_fct

    ## $emmeans
    ## exposition = north, survey_year_fct = 2018:
    ##  target_type   sand_ratio  emmean lower.HPD upper.HPD
    ##  hay_meadow    0          -0.7290   -0.8868   -0.5731
    ##  dry_grassland 0          -0.7273   -0.8840   -0.5611
    ##  hay_meadow    25         -0.5787   -0.7408   -0.4108
    ##  dry_grassland 25         -0.7217   -0.8912   -0.5552
    ##  hay_meadow    50         -0.5959   -0.7597   -0.4303
    ##  dry_grassland 50         -0.5495   -0.7196   -0.3925
    ## 
    ## exposition = south, survey_year_fct = 2018:
    ##  target_type   sand_ratio  emmean lower.HPD upper.HPD
    ##  hay_meadow    0          -1.3404   -1.4986   -1.1774
    ##  dry_grassland 0          -1.1769   -1.3353   -1.0030
    ##  hay_meadow    25         -1.6058   -1.7681   -1.4367
    ##  dry_grassland 25         -1.2188   -1.3857   -1.0564
    ##  hay_meadow    50         -1.8548   -2.0186   -1.6886
    ##  dry_grassland 50         -1.4954   -1.6554   -1.3245
    ## 
    ## exposition = north, survey_year_fct = 2019:
    ##  target_type   sand_ratio  emmean lower.HPD upper.HPD
    ##  hay_meadow    0           0.4075   -0.5417    1.3434
    ##  dry_grassland 0           0.3170   -0.6856    1.2259
    ##  hay_meadow    25          0.2483   -0.6907    1.2065
    ##  dry_grassland 25          0.1906   -0.7119    1.1870
    ##  hay_meadow    50          0.4629   -0.5037    1.4055
    ##  dry_grassland 50          0.4923   -0.4315    1.4693
    ## 
    ## exposition = south, survey_year_fct = 2019:
    ##  target_type   sand_ratio  emmean lower.HPD upper.HPD
    ##  hay_meadow    0           0.0389   -0.4234    0.5642
    ##  dry_grassland 0          -0.0696   -0.5627    0.4331
    ##  hay_meadow    25          0.0218   -0.4808    0.5211
    ##  dry_grassland 25          0.0476   -0.4392    0.5619
    ##  hay_meadow    50         -0.0336   -0.5526    0.4453
    ##  dry_grassland 50          0.1666   -0.3359    0.6720
    ## 
    ## exposition = north, survey_year_fct = 2020:
    ##  target_type   sand_ratio  emmean lower.HPD upper.HPD
    ##  hay_meadow    0           0.7959    0.6340    0.9628
    ##  dry_grassland 0           0.8591    0.6973    1.0254
    ##  hay_meadow    25          0.8038    0.6407    0.9665
    ##  dry_grassland 25          0.7718    0.6096    0.9408
    ##  hay_meadow    50          0.8953    0.7303    1.0632
    ##  dry_grassland 50          0.8681    0.6947    1.0275
    ## 
    ## exposition = south, survey_year_fct = 2020:
    ##  target_type   sand_ratio  emmean lower.HPD upper.HPD
    ##  hay_meadow    0           0.0713   -0.0962    0.2342
    ##  dry_grassland 0          -0.0788   -0.2402    0.0886
    ##  hay_meadow    25          0.1481   -0.0186    0.3139
    ##  dry_grassland 25          0.2571    0.0928    0.4243
    ##  hay_meadow    50          0.0527   -0.1132    0.2220
    ##  dry_grassland 50          0.3006    0.1283    0.4590
    ## 
    ## exposition = north, survey_year_fct = 2021:
    ##  target_type   sand_ratio  emmean lower.HPD upper.HPD
    ##  hay_meadow    0           0.8984    0.7339    1.0606
    ##  dry_grassland 0           0.9189    0.7570    1.0913
    ##  hay_meadow    25          0.9454    0.7783    1.1122
    ##  dry_grassland 25          0.9064    0.7524    1.0817
    ##  hay_meadow    50          0.9930    0.8242    1.1565
    ##  dry_grassland 50          1.0168    0.8594    1.1944
    ## 
    ## exposition = south, survey_year_fct = 2021:
    ##  target_type   sand_ratio  emmean lower.HPD upper.HPD
    ##  hay_meadow    0           0.2731    0.1074    0.4381
    ##  dry_grassland 0           0.1349   -0.0323    0.2938
    ##  hay_meadow    25          0.1473   -0.0164    0.3133
    ##  dry_grassland 25          0.2968    0.1301    0.4669
    ##  hay_meadow    50          0.1474   -0.0227    0.3092
    ##  dry_grassland 50          0.2352    0.0650    0.3960
    ## 
    ## Results are averaged over the levels of: substrate_depth, seed_density, botanist_year 
    ## Point estimate displayed: median 
    ## HPD interval probability: 0.95 
    ## 
    ## $contrasts
    ## exposition = north, survey_year_fct = 2018:
    ##  contrast                                                estimate lower.HPD
    ##  dry_grassland sand_ratio0 - hay_meadow sand_ratio0       0.00156  -0.18759
    ##  hay_meadow sand_ratio25 - hay_meadow sand_ratio0         0.15103  -0.04974
    ##  hay_meadow sand_ratio25 - dry_grassland sand_ratio0      0.14973  -0.05107
    ##  dry_grassland sand_ratio25 - hay_meadow sand_ratio0      0.00610  -0.19370
    ##  dry_grassland sand_ratio25 - dry_grassland sand_ratio0   0.00696  -0.19400
    ##  dry_grassland sand_ratio25 - hay_meadow sand_ratio25    -0.14231  -0.34600
    ##  hay_meadow sand_ratio50 - hay_meadow sand_ratio0         0.13489  -0.05256
    ##  hay_meadow sand_ratio50 - dry_grassland sand_ratio0      0.13170  -0.07231
    ##  hay_meadow sand_ratio50 - hay_meadow sand_ratio25       -0.01896  -0.22546
    ##  hay_meadow sand_ratio50 - dry_grassland sand_ratio25     0.12563  -0.08387
    ##  dry_grassland sand_ratio50 - hay_meadow sand_ratio0      0.18044  -0.02247
    ##  dry_grassland sand_ratio50 - dry_grassland sand_ratio0   0.17921  -0.02489
    ##  dry_grassland sand_ratio50 - hay_meadow sand_ratio25     0.02790  -0.17391
    ##  dry_grassland sand_ratio50 - dry_grassland sand_ratio25  0.17328  -0.01965
    ##  dry_grassland sand_ratio50 - hay_meadow sand_ratio50     0.04653  -0.15562
    ##  upper.HPD
    ##    0.19174
    ##    0.33923
    ##    0.35425
    ##    0.20785
    ##    0.21009
    ##    0.05573
    ##    0.32858
    ##    0.32923
    ##    0.17980
    ##    0.32781
    ##    0.37533
    ##    0.37594
    ##    0.22973
    ##    0.39144
    ##    0.24050
    ## 
    ## exposition = south, survey_year_fct = 2018:
    ##  contrast                                                estimate lower.HPD
    ##  dry_grassland sand_ratio0 - hay_meadow sand_ratio0       0.16372  -0.03078
    ##  hay_meadow sand_ratio25 - hay_meadow sand_ratio0        -0.26527  -0.46785
    ##  hay_meadow sand_ratio25 - dry_grassland sand_ratio0     -0.42899  -0.64831
    ##  dry_grassland sand_ratio25 - hay_meadow sand_ratio0      0.12357  -0.08157
    ##  dry_grassland sand_ratio25 - dry_grassland sand_ratio0  -0.04126  -0.23835
    ##  dry_grassland sand_ratio25 - hay_meadow sand_ratio25     0.38887   0.18221
    ##  hay_meadow sand_ratio50 - hay_meadow sand_ratio0        -0.51265  -0.71349
    ##  hay_meadow sand_ratio50 - dry_grassland sand_ratio0     -0.67564  -0.89026
    ##  hay_meadow sand_ratio50 - hay_meadow sand_ratio25       -0.24759  -0.45714
    ##  hay_meadow sand_ratio50 - dry_grassland sand_ratio25    -0.63537  -0.84118
    ##  dry_grassland sand_ratio50 - hay_meadow sand_ratio0     -0.15300  -0.36233
    ##  dry_grassland sand_ratio50 - dry_grassland sand_ratio0  -0.31802  -0.52420
    ##  dry_grassland sand_ratio50 - hay_meadow sand_ratio25     0.11218  -0.08835
    ##  dry_grassland sand_ratio50 - dry_grassland sand_ratio25 -0.27747  -0.47659
    ##  dry_grassland sand_ratio50 - hay_meadow sand_ratio50     0.35848   0.15930
    ##  upper.HPD
    ##    0.36013
    ##   -0.07129
    ##   -0.23966
    ##    0.32597
    ##    0.16298
    ##    0.58304
    ##   -0.31229
    ##   -0.47703
    ##   -0.04942
    ##   -0.43468
    ##    0.04246
    ##   -0.11952
    ##    0.31751
    ##   -0.07067
    ##    0.56780
    ## 
    ## exposition = north, survey_year_fct = 2019:
    ##  contrast                                                estimate lower.HPD
    ##  dry_grassland sand_ratio0 - hay_meadow sand_ratio0      -0.08630  -0.28454
    ##  hay_meadow sand_ratio25 - hay_meadow sand_ratio0        -0.15685  -0.35684
    ##  hay_meadow sand_ratio25 - dry_grassland sand_ratio0     -0.07197  -0.27471
    ##  dry_grassland sand_ratio25 - hay_meadow sand_ratio0     -0.21454  -0.42081
    ##  dry_grassland sand_ratio25 - dry_grassland sand_ratio0  -0.12829  -0.33037
    ##  dry_grassland sand_ratio25 - hay_meadow sand_ratio25    -0.05711  -0.25306
    ##  hay_meadow sand_ratio50 - hay_meadow sand_ratio0         0.05874  -0.13366
    ##  hay_meadow sand_ratio50 - dry_grassland sand_ratio0      0.14402  -0.05534
    ##  hay_meadow sand_ratio50 - hay_meadow sand_ratio25        0.21587   0.01219
    ##  hay_meadow sand_ratio50 - dry_grassland sand_ratio25     0.27044   0.07634
    ##  dry_grassland sand_ratio50 - hay_meadow sand_ratio0      0.08998  -0.12384
    ##  dry_grassland sand_ratio50 - dry_grassland sand_ratio0   0.17359  -0.03818
    ##  dry_grassland sand_ratio50 - hay_meadow sand_ratio25     0.24440   0.03974
    ##  dry_grassland sand_ratio50 - dry_grassland sand_ratio25  0.30034   0.10156
    ##  dry_grassland sand_ratio50 - hay_meadow sand_ratio50     0.02804  -0.18238
    ##  upper.HPD
    ##    0.11876
    ##    0.04649
    ##    0.12934
    ##   -0.01586
    ##    0.07373
    ##    0.15207
    ##    0.27011
    ##    0.35070
    ##    0.42447
    ##    0.48701
    ##    0.28116
    ##    0.37134
    ##    0.44825
    ##    0.50846
    ##    0.22671
    ## 
    ## exposition = south, survey_year_fct = 2019:
    ##  contrast                                                estimate lower.HPD
    ##  dry_grassland sand_ratio0 - hay_meadow sand_ratio0      -0.10796  -0.31659
    ##  hay_meadow sand_ratio25 - hay_meadow sand_ratio0        -0.01864  -0.22612
    ##  hay_meadow sand_ratio25 - dry_grassland sand_ratio0      0.08799  -0.11185
    ##  dry_grassland sand_ratio25 - hay_meadow sand_ratio0      0.00995  -0.18877
    ##  dry_grassland sand_ratio25 - dry_grassland sand_ratio0   0.11731  -0.09400
    ##  dry_grassland sand_ratio25 - hay_meadow sand_ratio25     0.02722  -0.18343
    ##  hay_meadow sand_ratio50 - hay_meadow sand_ratio0        -0.07259  -0.28830
    ##  hay_meadow sand_ratio50 - dry_grassland sand_ratio0      0.03490  -0.17624
    ##  hay_meadow sand_ratio50 - hay_meadow sand_ratio25       -0.05565  -0.25361
    ##  hay_meadow sand_ratio50 - dry_grassland sand_ratio25    -0.08132  -0.28580
    ##  dry_grassland sand_ratio50 - hay_meadow sand_ratio0      0.12710  -0.08169
    ##  dry_grassland sand_ratio50 - dry_grassland sand_ratio0   0.23396   0.02294
    ##  dry_grassland sand_ratio50 - hay_meadow sand_ratio25     0.14588  -0.06147
    ##  dry_grassland sand_ratio50 - dry_grassland sand_ratio25  0.11832  -0.08278
    ##  dry_grassland sand_ratio50 - hay_meadow sand_ratio50     0.20031  -0.00424
    ##  upper.HPD
    ##    0.09394
    ##    0.18667
    ##    0.30328
    ##    0.22236
    ##    0.32239
    ##    0.23074
    ##    0.12692
    ##    0.24289
    ##    0.15259
    ##    0.12770
    ##    0.33152
    ##    0.43819
    ##    0.34900
    ##    0.32380
    ##    0.40310
    ## 
    ## exposition = north, survey_year_fct = 2020:
    ##  contrast                                                estimate lower.HPD
    ##  dry_grassland sand_ratio0 - hay_meadow sand_ratio0       0.06340  -0.15044
    ##  hay_meadow sand_ratio25 - hay_meadow sand_ratio0         0.00738  -0.19495
    ##  hay_meadow sand_ratio25 - dry_grassland sand_ratio0     -0.05473  -0.26954
    ##  dry_grassland sand_ratio25 - hay_meadow sand_ratio0     -0.02509  -0.22636
    ##  dry_grassland sand_ratio25 - dry_grassland sand_ratio0  -0.08752  -0.28585
    ##  dry_grassland sand_ratio25 - hay_meadow sand_ratio25    -0.03105  -0.23704
    ##  hay_meadow sand_ratio50 - hay_meadow sand_ratio0         0.09772  -0.10305
    ##  hay_meadow sand_ratio50 - dry_grassland sand_ratio0      0.03582  -0.16188
    ##  hay_meadow sand_ratio50 - hay_meadow sand_ratio25        0.09160  -0.11379
    ##  hay_meadow sand_ratio50 - dry_grassland sand_ratio25     0.12367  -0.08754
    ##  dry_grassland sand_ratio50 - hay_meadow sand_ratio0      0.07134  -0.14240
    ##  dry_grassland sand_ratio50 - dry_grassland sand_ratio0   0.00871  -0.19718
    ##  dry_grassland sand_ratio50 - hay_meadow sand_ratio25     0.06465  -0.14336
    ##  dry_grassland sand_ratio50 - dry_grassland sand_ratio25  0.09550  -0.10780
    ##  dry_grassland sand_ratio50 - hay_meadow sand_ratio50    -0.02659  -0.23377
    ##  upper.HPD
    ##    0.25670
    ##    0.21314
    ##    0.14330
    ##    0.18662
    ##    0.12088
    ##    0.17506
    ##    0.30587
    ##    0.25513
    ##    0.30487
    ##    0.32654
    ##    0.27532
    ##    0.21450
    ##    0.27007
    ##    0.30207
    ##    0.18387
    ## 
    ## exposition = south, survey_year_fct = 2020:
    ##  contrast                                                estimate lower.HPD
    ##  dry_grassland sand_ratio0 - hay_meadow sand_ratio0      -0.15059  -0.35942
    ##  hay_meadow sand_ratio25 - hay_meadow sand_ratio0         0.07761  -0.12210
    ##  hay_meadow sand_ratio25 - dry_grassland sand_ratio0      0.22740   0.02189
    ##  dry_grassland sand_ratio25 - hay_meadow sand_ratio0      0.18738  -0.01010
    ##  dry_grassland sand_ratio25 - dry_grassland sand_ratio0   0.33521   0.12899
    ##  dry_grassland sand_ratio25 - hay_meadow sand_ratio25     0.10899  -0.09972
    ##  hay_meadow sand_ratio50 - hay_meadow sand_ratio0        -0.01677  -0.21317
    ##  hay_meadow sand_ratio50 - dry_grassland sand_ratio0      0.13197  -0.07153
    ##  hay_meadow sand_ratio50 - hay_meadow sand_ratio25       -0.09498  -0.30938
    ##  hay_meadow sand_ratio50 - dry_grassland sand_ratio25    -0.20348  -0.40997
    ##  dry_grassland sand_ratio50 - hay_meadow sand_ratio0      0.23086   0.01671
    ##  dry_grassland sand_ratio50 - dry_grassland sand_ratio0   0.37945   0.17780
    ##  dry_grassland sand_ratio50 - hay_meadow sand_ratio25     0.15234  -0.05383
    ##  dry_grassland sand_ratio50 - dry_grassland sand_ratio25  0.04465  -0.16330
    ##  dry_grassland sand_ratio50 - hay_meadow sand_ratio50     0.24695   0.03419
    ##  upper.HPD
    ##    0.05282
    ##    0.29061
    ##    0.42930
    ##    0.40345
    ##    0.53794
    ##    0.30791
    ##    0.19485
    ##    0.33635
    ##    0.10680
    ##    0.00579
    ##    0.42967
    ##    0.59009
    ##    0.35893
    ##    0.24253
    ##    0.44499
    ## 
    ## exposition = north, survey_year_fct = 2021:
    ##  contrast                                                estimate lower.HPD
    ##  dry_grassland sand_ratio0 - hay_meadow sand_ratio0       0.01901  -0.18317
    ##  hay_meadow sand_ratio25 - hay_meadow sand_ratio0         0.04649  -0.15695
    ##  hay_meadow sand_ratio25 - dry_grassland sand_ratio0      0.02634  -0.18902
    ##  dry_grassland sand_ratio25 - hay_meadow sand_ratio0      0.00791  -0.20005
    ##  dry_grassland sand_ratio25 - dry_grassland sand_ratio0  -0.01170  -0.22242
    ##  dry_grassland sand_ratio25 - hay_meadow sand_ratio25    -0.03769  -0.24242
    ##  hay_meadow sand_ratio50 - hay_meadow sand_ratio0         0.09373  -0.11496
    ##  hay_meadow sand_ratio50 - dry_grassland sand_ratio0      0.07385  -0.13928
    ##  hay_meadow sand_ratio50 - hay_meadow sand_ratio25        0.04840  -0.16008
    ##  hay_meadow sand_ratio50 - dry_grassland sand_ratio25     0.08625  -0.12310
    ##  dry_grassland sand_ratio50 - hay_meadow sand_ratio0      0.11938  -0.09004
    ##  dry_grassland sand_ratio50 - dry_grassland sand_ratio0   0.09818  -0.11644
    ##  dry_grassland sand_ratio50 - hay_meadow sand_ratio25     0.07167  -0.13382
    ##  dry_grassland sand_ratio50 - dry_grassland sand_ratio25  0.11096  -0.09715
    ##  dry_grassland sand_ratio50 - hay_meadow sand_ratio50     0.02409  -0.17957
    ##  upper.HPD
    ##    0.22040
    ##    0.25313
    ##    0.22791
    ##    0.20508
    ##    0.18940
    ##    0.16570
    ##    0.29260
    ##    0.27638
    ##    0.25306
    ##    0.28153
    ##    0.32210
    ##    0.29980
    ##    0.27953
    ##    0.31476
    ##    0.22724
    ## 
    ## exposition = south, survey_year_fct = 2021:
    ##  contrast                                                estimate lower.HPD
    ##  dry_grassland sand_ratio0 - hay_meadow sand_ratio0      -0.13864  -0.34880
    ##  hay_meadow sand_ratio25 - hay_meadow sand_ratio0        -0.12418  -0.32892
    ##  hay_meadow sand_ratio25 - dry_grassland sand_ratio0      0.01302  -0.19304
    ##  dry_grassland sand_ratio25 - hay_meadow sand_ratio0      0.02459  -0.18159
    ##  dry_grassland sand_ratio25 - dry_grassland sand_ratio0   0.16314  -0.05069
    ##  dry_grassland sand_ratio25 - hay_meadow sand_ratio25     0.14883  -0.05209
    ##  hay_meadow sand_ratio50 - hay_meadow sand_ratio0        -0.12612  -0.33492
    ##  hay_meadow sand_ratio50 - dry_grassland sand_ratio0      0.01109  -0.18506
    ##  hay_meadow sand_ratio50 - hay_meadow sand_ratio25       -0.00109  -0.20457
    ##  hay_meadow sand_ratio50 - dry_grassland sand_ratio25    -0.14980  -0.34956
    ##  dry_grassland sand_ratio50 - hay_meadow sand_ratio0     -0.03753  -0.24174
    ##  dry_grassland sand_ratio50 - dry_grassland sand_ratio0   0.10121  -0.10433
    ##  dry_grassland sand_ratio50 - hay_meadow sand_ratio25     0.08764  -0.12168
    ##  dry_grassland sand_ratio50 - dry_grassland sand_ratio25 -0.06264  -0.27265
    ##  dry_grassland sand_ratio50 - hay_meadow sand_ratio50     0.08874  -0.12266
    ##  upper.HPD
    ##    0.06201
    ##    0.08418
    ##    0.21917
    ##    0.23263
    ##    0.36493
    ##    0.36520
    ##    0.08002
    ##    0.23119
    ##    0.20822
    ##    0.06412
    ##    0.16779
    ##    0.30248
    ##    0.28976
    ##    0.14428
    ##    0.29339
    ## 
    ## Results are averaged over the levels of: substrate_depth, seed_density, botanist_year 
    ## Point estimate displayed: median 
    ## HPD interval probability: 0.95

# Session info

    ## R version 4.2.2 (2022-10-31 ucrt)
    ## Platform: x86_64-w64-mingw32/x64 (64-bit)
    ## Running under: Windows 10 x64 (build 22621)
    ## 
    ## Matrix products: default
    ## 
    ## locale:
    ## [1] LC_COLLATE=German_Germany.utf8  LC_CTYPE=German_Germany.utf8   
    ## [3] LC_MONETARY=German_Germany.utf8 LC_NUMERIC=C                   
    ## [5] LC_TIME=German_Germany.utf8    
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] emmeans_1.8.2             tidybayes_3.0.2          
    ##  [3] loo_2.5.1                 bayesplot_1.10.0         
    ##  [5] DHARMa.helpers_0.0.0.9000 DHARMa_0.4.6             
    ##  [7] brms_2.18.0               Rcpp_1.0.9               
    ##  [9] patchwork_1.1.2           ggbeeswarm_0.6.0         
    ## [11] forcats_0.5.2             stringr_1.4.1            
    ## [13] dplyr_1.0.10              purrr_0.3.5              
    ## [15] readr_2.1.3               tidyr_1.2.1              
    ## [17] tibble_3.1.8              ggplot2_3.4.0            
    ## [19] tidyverse_1.3.2           here_1.0.1               
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] utf8_1.2.2           tidyselect_1.2.0     lme4_1.1-31         
    ##   [4] htmlwidgets_1.5.4    grid_4.2.2           RNeXML_2.4.8        
    ##   [7] munsell_0.5.0        codetools_0.2-18     interp_1.1-3        
    ##  [10] DT_0.26              miniUI_0.1.1.1       withr_2.5.0         
    ##  [13] Brobdingnag_1.2-9    colorspace_2.0-3     qgam_1.3.4          
    ##  [16] uuid_1.1-0           highr_0.9            knitr_1.41          
    ##  [19] rstudioapi_0.14      stats4_4.2.2         labeling_0.4.2      
    ##  [22] rstan_2.26.13        bit64_4.0.5          farver_2.1.1        
    ##  [25] gap.datasets_0.0.5   bridgesampling_1.1-2 rprojroot_2.0.3     
    ##  [28] coda_0.19-4          vctrs_0.5.1          generics_0.1.3      
    ##  [31] xfun_0.35            timechange_0.1.1     adegenet_2.1.8      
    ##  [34] doParallel_1.0.17    R6_2.5.1             markdown_1.4        
    ##  [37] assertthat_0.2.1     promises_1.2.0.1     scales_1.2.1        
    ##  [40] vroom_1.6.0          googlesheets4_1.0.1  phylobase_0.8.10    
    ##  [43] beeswarm_0.4.0       gtable_0.3.1         processx_3.8.0      
    ##  [46] rlang_1.0.6          splines_4.2.2        gargle_1.2.1        
    ##  [49] broom_1.0.1          checkmate_2.1.0      inline_0.3.19       
    ##  [52] yaml_2.3.6           reshape2_1.4.4       abind_1.4-5         
    ##  [55] modelr_0.1.10        threejs_0.3.3        crosstalk_1.2.0     
    ##  [58] backports_1.4.1      httpuv_1.6.6         tensorA_0.36.2      
    ##  [61] tools_4.2.2          ellipsis_0.3.2       posterior_1.3.1     
    ##  [64] RColorBrewer_1.1-3   plyr_1.8.8           progress_1.2.2      
    ##  [67] base64enc_0.1-3      ps_1.7.2             prettyunits_1.1.1   
    ##  [70] deldir_1.0-6         zoo_1.8-11           haven_2.5.1         
    ##  [73] cluster_2.1.4        fs_1.5.2             magrittr_2.0.3      
    ##  [76] ggdist_3.2.0         colourpicker_1.2.0   reprex_2.0.2        
    ##  [79] googledrive_2.0.0    mvtnorm_1.1-3        matrixStats_0.63.0  
    ##  [82] hms_1.1.2            shinyjs_2.1.0        mime_0.12           
    ##  [85] evaluate_0.18        arrayhelpers_1.1-0   xtable_1.8-4        
    ##  [88] XML_3.99-0.12        shinystan_2.6.0      jpeg_0.1-10         
    ##  [91] readxl_1.4.1         gridExtra_2.3        rstantools_2.2.0    
    ##  [94] compiler_4.2.2       KernSmooth_2.23-20   V8_4.2.2            
    ##  [97] crayon_1.5.2         minqa_1.2.5          StanHeaders_2.26.13 
    ## [100] htmltools_0.5.3      mgcv_1.8-41          later_1.3.0         
    ## [103] tzdb_0.3.0           RcppParallel_5.1.5   lubridate_1.9.0     
    ## [106] DBI_1.1.3            dbplyr_2.2.1         MASS_7.3-58.1       
    ## [109] boot_1.3-28          Matrix_1.5-3         ade4_1.7-20         
    ## [112] permute_0.9-7        cli_3.4.1            adegraphics_1.0-16  
    ## [115] parallel_4.2.2       igraph_1.3.5         pkgconfig_2.0.3     
    ## [118] rncl_0.8.6           sp_1.5-1             foreach_1.5.2       
    ## [121] xml2_1.3.3           svUnit_1.0.6         dygraphs_1.1.1.6    
    ## [124] vipor_0.4.5          estimability_1.4.1   rvest_1.0.3         
    ## [127] distributional_0.3.1 callr_3.7.3          digest_0.6.30       
    ## [130] vegan_2.6-4          rmarkdown_2.18       cellranger_1.1.0    
    ## [133] gap_1.3-1            curl_4.3.3           shiny_1.7.3         
    ## [136] gtools_3.9.4         nloptr_2.0.3         lifecycle_1.0.3     
    ## [139] nlme_3.1-160         jsonlite_1.8.3       seqinr_4.2-23       
    ## [142] fansi_1.0.3          pillar_1.8.1         lattice_0.20-45     
    ## [145] fastmap_1.1.0        httr_1.4.4           pkgbuild_1.4.0      
    ## [148] glue_1.6.2           xts_0.12.2           iterators_1.0.14    
    ## [151] png_0.1-8            shinythemes_1.2.0    bit_4.0.5           
    ## [154] stringi_1.7.8        latticeExtra_0.6-30  ape_5.6-2