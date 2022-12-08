Analysis of Bauer et al. (unpublished) Field experiment: <br>
Persistence
================
<b>Markus Bauer</b> <br>
<b>2022-12-07</b>

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
    - <a href="#load-models-barg-6g" id="toc-load-models-barg-6g">Load models
      (BARG 6.G)</a>
    - <a href="#model-specifications-barg-1d"
      id="toc-model-specifications-barg-1d">Model specifications (BARG
      1.D)</a>
    - <a href="#priors-barg-1d" id="toc-priors-barg-1d">Priors (BARG 1.D)</a>
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
  - <a href="#model-comparison" id="toc-model-comparison">Model
    comparison</a>
    - <a href="#conditional-r²-values"
      id="toc-conditional-r²-values">Conditional <i>R</i>² values</a>
    - <a href="#marginal-r²-values" id="toc-marginal-r²-values">Marginal
      <i>R</i>² values</a>
    - <a href="#bayes-factor-barg-3c" id="toc-bayes-factor-barg-3c">Bayes
      factor (BARG 3.C)</a>
  - <a href="#posterior-distributions-barg-3b"
    id="toc-posterior-distributions-barg-3b">Posterior distributions (BARG
    3.B)</a>
    - <a href="#forest-plot" id="toc-forest-plot">Forest plot</a>
    - <a href="#effect-sizes" id="toc-effect-sizes">Effect sizes</a>
- <a href="#session-info-barg-2a6a6b"
  id="toc-session-info-barg-2a6a6b">Session info (BARG 2.A/6.A/6.B)</a>

<b>Markus Bauer</b>

Technichal University of Munich, TUM School of Life Sciences, Chair of
Restoration Ecology, Emil-Ramann-Straße 6, 85354 Freising, Germany

<markus1.bauer@tum.de>

ORCiD ID: [0000-0001-5372-4174](https://orcid.org/0000-0001-5372-4174)
<br> [Google
Scholar](https://scholar.google.de/citations?user=oHhmOkkAAAAJ&hl=de&oi=ao)
<br> GitHub: [markus1bauer](https://github.com/markus1bauer)

To compare different models, you only have to change the models in
section ‘Load models’

# Preparation

Favourable Conservation Status (FSC) sensu Helm et al. (2015) Divers
Distrib [DOI: 10.1111/ddi.12285](https://doi.org/10.1111/ddi.12285)

Analysis motivated by Applestein et al. (2021) Restor Ecol [DOI:
10.1111/rec.13596](https://doi.org/10.1111/rec.13596)

Analysis guided by <br> <b>BARG</b> (Bayesian Analysis Reporting
Guidelines): Kruschke (2021) Nat Hum Behav [DOI:
10.1038/s41562-021-01177-7](https://doi.org/10.1038/s41562-021-01177-7)
<br> Model check: Gabry et al. (2019) J R Stat Soc A Stat [DOI:
10.1111/rssa.12378](https://doi.org/10.1111/rssa.12378) <br> Priors:
Lemoine (2019) Oikos [DOI:
10.1111/oik.05985](https://doi.org/10.1111/oik.05985) <br> Model check:
Deapoli & Schoot (2017) Psychol Methods [DOI:
10.1037/met0000065](https://doi.org/10.1037/met0000065)

#### Packages (BARG 2.A)

``` r
library(here)
library(tidyverse)
library(ggbeeswarm)
library(patchwork)
library(brms)
library(DHARMa.helpers)
library(bayesplot)
library(loo)
library(tidybayes)
library(emmeans)
```

#### Load data

``` r
sites <- read_csv(
  here("data", "processed", "data_processed_sites_temporal.csv"),
  col_names = TRUE, na = c("na", "NA", ""),
  col_types = cols(
    .default = "?",
    plot = "f",
    site = "f",
    sand_ratio = "f",
    substrate_depth = col_factor(levels = c("30", "15")),
    target_type = col_factor(levels = c("hay_meadow", "dry_grassland")),
    seed_density = "f",
    exposition = col_factor(levels = c("north", "south")),
    survey_year = "c"
  )
  ) %>%
  ### Exclude data of seed mixtures
  filter(presabu == "presence") %>%
  mutate(
    survey_year_fct = factor(survey_year),
    botanist_year = str_c(survey_year, botanist, exposition, sep = " "),
    botanist_year = factor(botanist_year),
    id = factor(id),
    n = persistence / 100
  ) %>%
  select(
    id, plot, site, exposition, sand_ratio, substrate_depth, target_type,
    seed_density, survey_year_fct, survey_year, botanist_year, n
  )
```

# Statistics

## Data exploration

### Graphs of raw data

![](model_persistence_check_files/figure-gfm/data-exploration-1.png)<!-- -->![](model_persistence_check_files/figure-gfm/data-exploration-2.png)<!-- -->![](model_persistence_check_files/figure-gfm/data-exploration-3.png)<!-- -->![](model_persistence_check_files/figure-gfm/data-exploration-4.png)<!-- -->![](model_persistence_check_files/figure-gfm/data-exploration-5.png)<!-- -->

### Outliers, zero-inflation, transformations?

    ## # A tibble: 12 × 3
    ## # Groups:   exposition [2]
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

![](model_persistence_check_files/figure-gfm/outliers-1.png)<!-- -->![](model_persistence_check_files/figure-gfm/outliers-2.png)<!-- -->![](model_persistence_check_files/figure-gfm/outliers-3.png)<!-- -->![](model_persistence_check_files/figure-gfm/outliers-4.png)<!-- -->

## Models

### Load models (BARG 6.G)

Only here you have to modify the script to compare other models

``` r
load(file = here("outputs", "models", "model_persistence_2.Rdata"))
#load(file = here("outputs", "models", "model_persistence_full.Rdata"))
# BARG 1.E
load(file = here("outputs", "models", "model_persistence_2_prior.Rdata"))
# BARG 5.A/B/C
#load(file = here("outputs", "models", "model_persistence_2_flat.Rdata"))
m_1 <- m2
#m_2 <- m_full
m_2 <- m2_prior
```

### Model specifications (BARG 1.D)

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

Amount of chains for MCMC

``` r
m_1$fit@sim$chains
## [1] 4
m_2$fit@sim$chains
## [1] 4
```

Total amount of iterations for MCMC

``` r
m_1$fit@sim$iter
## [1] 20000
m_2$fit@sim$iter
## [1] 10000
```

Amount of iterations before burn-in

``` r
m_1$fit@sim$warmup
## [1] 10000
m_2$fit@sim$warmup
## [1] 5000
```

Thinning rate

``` r
m_1$fit@sim$thin
## [1] 2
m_2$fit@sim$thin
## [1] 2
```

### Priors (BARG 1.D)

#### Possible prior distributions

``` r
get_prior(n ~ target_type + exposition + sand_ratio + survey_year_fct +
            seed_density + substrate_depth +
            (1 | site/plot) + (1 | botanist_year),
          data = sites)
```

    ##                   prior     class                     coef         group resp
    ##                  (flat)         b                                            
    ##                  (flat)         b          expositionsouth                   
    ##                  (flat)         b             sand_ratio25                   
    ##                  (flat)         b             sand_ratio50                   
    ##                  (flat)         b            seed_density8                   
    ##                  (flat)         b        substrate_depth15                   
    ##                  (flat)         b      survey_year_fct2019                   
    ##                  (flat)         b      survey_year_fct2020                   
    ##                  (flat)         b      survey_year_fct2021                   
    ##                  (flat)         b target_typedry_grassland                   
    ##  student_t(3, 0.6, 2.5) Intercept                                            
    ##    student_t(3, 0, 2.5)        sd                                            
    ##    student_t(3, 0, 2.5)        sd                          botanist_year     
    ##    student_t(3, 0, 2.5)        sd                Intercept botanist_year     
    ##    student_t(3, 0, 2.5)        sd                                   site     
    ##    student_t(3, 0, 2.5)        sd                Intercept          site     
    ##    student_t(3, 0, 2.5)        sd                              site:plot     
    ##    student_t(3, 0, 2.5)        sd                Intercept     site:plot     
    ##    student_t(3, 0, 2.5)     sigma                                            
    ##  dpar nlpar lb ub       source
    ##                        default
    ##                   (vectorized)
    ##                   (vectorized)
    ##                   (vectorized)
    ##                   (vectorized)
    ##                   (vectorized)
    ##                   (vectorized)
    ##                   (vectorized)
    ##                   (vectorized)
    ##                   (vectorized)
    ##                        default
    ##              0         default
    ##              0    (vectorized)
    ##              0    (vectorized)
    ##              0    (vectorized)
    ##              0    (vectorized)
    ##              0    (vectorized)
    ##              0    (vectorized)
    ##              0         default

``` r
ggplot(data = data.frame(x = c(0, 1)), aes(x = x)) +
  stat_function(fun = dbeta, n = 101, args = list(shape1 = 3.5, shape2 = 2.2)) +
  expand_limits(y = 0) +
  ggtitle("Beta distribution for Intercept")
```

![](model_persistence_check_files/figure-gfm/possible-priors-1.png)<!-- -->

``` r
ggplot(data = data.frame(x = c(-.4, .4)), aes(x = x)) +
  stat_function(fun = dnorm, n = 101, args = list(mean = 0, sd = .2)) +
  expand_limits(y = 0) +
  ggtitle("Normal distribution for treatments")
```

![](model_persistence_check_files/figure-gfm/possible-priors-2.png)<!-- -->

``` r
ggplot(data = data.frame(x = c(-.4, .4)), aes(x = x)) +
  stat_function(fun = dcauchy, n = 101, args = list(location = 0, scale = 1)) +
  expand_limits(y = 0) +
  ggtitle("Cauchy distribution") # See Lemoine 2019 https://doi.org/10.1111/oik.05985
```

![](model_persistence_check_files/figure-gfm/possible-priors-3.png)<!-- -->

``` r
ggplot(data.frame(x = c(-.4, .4)), aes(x = x)) +
  stat_function(fun = dstudent_t, args = list(df = 3, mu = 0, sigma = 2.5)) +
  expand_limits(y = 0) +
  ggtitle(expression(Student~italic(t)*"-distribution")) # Software standard
```

![](model_persistence_check_files/figure-gfm/possible-priors-4.png)<!-- -->

#### Prior summary (BARG 1.D)

``` r
brms::prior_summary(m_1, all = FALSE)
```

    ##                 prior     class                coef group resp dpar nlpar lb ub
    ##         normal(0, .2)         b                                           -1  1
    ##      normal(-.05, .2)         b     expositionsouth                       -1  1
    ##      normal(.025, .2)         b        sand_ratio25                       -1  1
    ##       normal(.05, .2)         b        sand_ratio50                       -1  1
    ##      normal(.025, .2)         b survey_year_fct2019                       -1  1
    ##       normal(.05, .2)         b survey_year_fct2020                       -1  1
    ##      normal(.075, .2)         b survey_year_fct2021                       -1  1
    ##        beta(3.5, 2.2) Intercept                                            0  1
    ##  student_t(3, 0, 2.5)        sd                                            0   
    ##         cauchy(0, .1)     sigma                                            0   
    ##   source
    ##     user
    ##     user
    ##     user
    ##     user
    ##     user
    ##     user
    ##     user
    ##     user
    ##  default
    ##     user

## Model check

### DHARMa

``` r
DHARMa.helpers::dh_check_brms(m_1, integer = TRUE)
```

![](model_persistence_check_files/figure-gfm/dharma-1.png)<!-- -->

``` r
DHARMa.helpers::dh_check_brms(m_2, integer = TRUE)
```

![](model_persistence_check_files/figure-gfm/dharma-2.png)<!-- -->

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
hmc_diagnostics1 <- brms::nuts_params(m_1)
hmc_diagnostics2 <- brms::nuts_params(m_2)
y <- sites$n
yrep1 <- brms::posterior_predict(m_1, draws = 500)
yrep2 <- brms::posterior_predict(m_2, draws = 500)
loo1 <- brms::loo(m_1, save_psis = TRUE, moment_match = FALSE)
loo2 <- brms::loo(m_2, save_psis = TRUE, moment_match = FALSE)
```

    ## Warning: Found 1152 observations with a pareto_k > 0.7 in model 'm_2'. It is
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

#### Rhat (BARG 2.B)

``` r
mcmc_rhat(draws1$rhat)
```

![](model_persistence_check_files/figure-gfm/rhat-1.png)<!-- -->

``` r
mcmc_rhat(draws2$rhat)
```

![](model_persistence_check_files/figure-gfm/rhat-2.png)<!-- -->

#### Effective sampling size (ESS) (BARG 2.C)

``` r
mcmc_neff(neff_ratio(m_1))
```

![](model_persistence_check_files/figure-gfm/ess-1.png)<!-- -->

``` r
mcmc_neff(neff_ratio(m_2))
```

![](model_persistence_check_files/figure-gfm/ess-2.png)<!-- -->

### MCMC diagnostics

#### Trace plots

``` r
mcmc_trace(posterior1, np = hmc_diagnostics1, facet_args = list(ncol = 2))
```

    ## No divergences to plot.

![](model_persistence_check_files/figure-gfm/mcmc-trace-1.png)<!-- -->

``` r
mcmc_trace(posterior2, np = hmc_diagnostics2, facet_args = list(ncol = 2))
```

    ## No divergences to plot.

![](model_persistence_check_files/figure-gfm/mcmc-trace-2.png)<!-- -->

#### Pairs plot

``` r
mcmc_pairs(m_1, off_diag_args = list(size = 1.2),
           pars = c(
             "b_sand_ratio25", "b_sand_ratio50", "b_substrate_depth15",
             "b_target_typedry_grassland", "b_seed_density8",
             "b_expositionsouth", "sigma"
           ))
```

![](model_persistence_check_files/figure-gfm/mcmc-pairs-1.png)<!-- -->

``` r
mcmc_pairs(m_2, off_diag_args = list(size = 1.2),
           pars = c(
             "b_sand_ratio25", "b_sand_ratio50", "b_substrate_depth15",
             "b_target_typedry_grassland", "b_seed_density8",
             "b_expositionsouth", "sigma"
           ))
```

![](model_persistence_check_files/figure-gfm/mcmc-pairs-2.png)<!-- -->

#### Parallel coordinate plot

``` r
mcmc_parcoord(posterior1, np = hmc_diagnostics1) +
  theme(axis.text.x = element_text(angle = 45))
```

![](model_persistence_check_files/figure-gfm/mcmc-parcoord-1.png)<!-- -->

``` r
mcmc_parcoord(posterior2, np = hmc_diagnostics2) +
  theme(axis.text.x = element_text(angle = 45))
```

![](model_persistence_check_files/figure-gfm/mcmc-parcoord-2.png)<!-- -->

### Posterior predictive check

#### Kernel density

``` r
p1 <- ppc_dens_overlay(y, yrep1[1:50, ])
p2 <- ppc_dens_overlay(y, yrep2[1:50, ])
p1 / p2
```

![](model_persistence_check_files/figure-gfm/kernel-density-1.png)<!-- -->

``` r
ppc_dens_overlay_grouped(y, yrep1[1:50, ], group = sites$site)
```

![](model_persistence_check_files/figure-gfm/kernel-density-2.png)<!-- -->

``` r
ppc_dens_overlay_grouped(y, yrep2[1:50, ], group = sites$site)
```

![](model_persistence_check_files/figure-gfm/kernel-density-3.png)<!-- -->

``` r
p1 <- ppc_dens_overlay_grouped(y, yrep1[1:50, ], group = sites$exposition)
p2 <- ppc_dens_overlay_grouped(y, yrep2[1:50, ], group = sites$exposition)
p1 / p2
```

![](model_persistence_check_files/figure-gfm/kernel-density-4.png)<!-- -->

``` r
ppc_dens_overlay_grouped(y, yrep1[1:50, ], group = sites$survey_year_fct)
```

![](model_persistence_check_files/figure-gfm/kernel-density-5.png)<!-- -->

``` r
ppc_dens_overlay_grouped(y, yrep2[1:50, ], group = sites$survey_year_fct)
```

![](model_persistence_check_files/figure-gfm/kernel-density-6.png)<!-- -->

``` r
p1 <- ppc_dens_overlay_grouped(y, yrep1[1:50, ], group = sites$target_type)
p2 <- ppc_dens_overlay_grouped(y, yrep2[1:50, ], group = sites$target_type)
p1 / p2
```

![](model_persistence_check_files/figure-gfm/kernel-density-7.png)<!-- -->

``` r
p1 <- ppc_dens_overlay_grouped(y, yrep1[1:50, ], group = sites$seed_density)
p2 <- ppc_dens_overlay_grouped(y, yrep2[1:50, ], group = sites$seed_density)
p1 / p2
```

![](model_persistence_check_files/figure-gfm/kernel-density-8.png)<!-- -->

``` r
p1 <- ppc_dens_overlay_grouped(y, yrep1[1:50, ], group = sites$sand_ratio)
p2 <- ppc_dens_overlay_grouped(y, yrep2[1:50, ], group = sites$sand_ratio)
p1 / p2
```

![](model_persistence_check_files/figure-gfm/kernel-density-9.png)<!-- -->

``` r
p1 <- ppc_dens_overlay_grouped(y, yrep1[1:50, ], group = sites$substrate_depth)
p2 <- ppc_dens_overlay_grouped(y, yrep2[1:50, ], group = sites$substrate_depth)
p1 / p2
```

![](model_persistence_check_files/figure-gfm/kernel-density-10.png)<!-- -->

#### Histograms of statistics skew

``` r
p1 <- ppc_stat(y, yrep1, binwidth = 0.001)
p2 <- ppc_stat(y, yrep2, binwidth = 0.001)
p1 / p2
```

![](model_persistence_check_files/figure-gfm/histograms-1.png)<!-- -->

``` r
ppc_stat_grouped(y, yrep1, group = sites$site, binwidth = 0.001)
```

![](model_persistence_check_files/figure-gfm/histograms-2.png)<!-- -->

``` r
ppc_stat_grouped(y, yrep2, group = sites$site, binwidth = 0.001)
```

![](model_persistence_check_files/figure-gfm/histograms-3.png)<!-- -->

``` r
p1 <- ppc_stat_grouped(y, yrep1, group = sites$exposition, binwidth = 0.001)
p2 <- ppc_stat_grouped(y, yrep2, group = sites$exposition, binwidth = 0.001)
p1 / p2
```

![](model_persistence_check_files/figure-gfm/histograms-4.png)<!-- -->

``` r
ppc_stat_grouped(y, yrep1, group = sites$survey_year_fct, binwidth = 0.001)
```

![](model_persistence_check_files/figure-gfm/histograms-5.png)<!-- -->

``` r
ppc_stat_grouped(y, yrep2, group = sites$survey_year_fct, binwidth = 0.001)
```

![](model_persistence_check_files/figure-gfm/histograms-6.png)<!-- -->

``` r
p1 <- ppc_stat_grouped(y, yrep1, group = sites$target_type, binwidth = 0.001)
p2 <- ppc_stat_grouped(y, yrep2, group = sites$target_type, binwidth = 0.001)
p1 / p2
```

![](model_persistence_check_files/figure-gfm/histograms-7.png)<!-- -->

``` r
p1 <- ppc_stat_grouped(y, yrep1, group = sites$seed_density, binwidth = 0.001)
p2 <- ppc_stat_grouped(y, yrep2, group = sites$seed_density, binwidth = 0.001)
p1 / p2
```

![](model_persistence_check_files/figure-gfm/histograms-8.png)<!-- -->

``` r
p1 <- ppc_stat_grouped(y, yrep1, group = sites$sand_ratio, binwidth = 0.001)
p2 <- ppc_stat_grouped(y, yrep2, group = sites$sand_ratio, binwidth = 0.001)
p1 / p2
```

![](model_persistence_check_files/figure-gfm/histograms-9.png)<!-- -->

``` r
p1 <- ppc_stat_grouped(y, yrep1, group = sites$substrate_depth, binwidth = 0.001)
p2 <- ppc_stat_grouped(y, yrep2, group = sites$substrate_depth, binwidth = 0.001)
p1 / p2
```

![](model_persistence_check_files/figure-gfm/histograms-10.png)<!-- -->

#### LOO (Leave one out)

``` r
loo1
```

    ## 
    ## Computed from 20000 by 1152 log-likelihood matrix
    ## 
    ##          Estimate   SE
    ## elpd_loo   1449.6 25.7
    ## p_loo       199.8  8.1
    ## looic     -2899.1 51.4
    ## ------
    ## Monte Carlo SE of elpd_loo is 0.2.
    ## 
    ## Pareto k diagnostic values:
    ##                          Count Pct.    Min. n_eff
    ## (-Inf, 0.5]   (good)     1135  98.5%   1080      
    ##  (0.5, 0.7]   (ok)         17   1.5%   592       
    ##    (0.7, 1]   (bad)         0   0.0%   <NA>      
    ##    (1, Inf)   (very bad)    0   0.0%   <NA>      
    ## 
    ## All Pareto k estimates are ok (k < 0.7).
    ## See help('pareto-k-diagnostic') for details.

``` r
loo2
```

    ## 
    ## Computed from 10000 by 1152 log-likelihood matrix
    ## 
    ##               Estimate           SE
    ## elpd_loo -1.196748e+16 3.876921e+14
    ## p_loo     1.196748e+16 3.876921e+14
    ## looic     2.393496e+16 7.753842e+14
    ## ------
    ## Monte Carlo SE of elpd_loo is NA.
    ## 
    ## Pareto k diagnostic values:
    ##                          Count Pct.    Min. n_eff
    ## (-Inf, 0.5]   (good)        0    0.0%  <NA>      
    ##  (0.5, 0.7]   (ok)          0    0.0%  <NA>      
    ##    (0.7, 1]   (bad)         0    0.0%  <NA>      
    ##    (1, Inf)   (very bad) 1152  100.0%  1         
    ## See help('pareto-k-diagnostic') for details.

``` r
#plot(loo1)
#plot(loo2)
```

Leave one out probability integral transform

``` r
#p1 <- bayesplot::ppc_loo_pit_overlay(y, yrep1, lw = weights(loo1$psis_object))
#p2 <- bayesplot::ppc_loo_pit_overlay(y, yrep2, lw = weights(loo2$psis_object))
#p1 / p2
```

### Autocorrelation check

``` r
mcmc_acf(posterior1, lags = 10)
```

![](model_persistence_check_files/figure-gfm/autocorrelation-1.png)<!-- -->

``` r
mcmc_acf(posterior2, lags = 10)
```

![](model_persistence_check_files/figure-gfm/autocorrelation-2.png)<!-- -->

## Model comparison

### Conditional <i>R</i>² values

``` r
brms::bayes_R2(m_1, probs = c(0.05, 0.5, 0.95),
         re_formula =  ~ (1 | site/plot)) 
##     Estimate   Est.Error        Q5       Q50       Q95
## R2 0.8943888 0.003371453 0.8885608 0.8945292 0.8997099
brms::bayes_R2(m_2, probs = c(0.05, 0.5, 0.95),
         re_formula =  ~ (1 | site/plot))
##     Estimate   Est.Error        Q5      Q50       Q95
## R2 0.4985877 0.006852914 0.4897801 0.499612 0.5038451
```

### Marginal <i>R</i>² values

``` r
brms::bayes_R2(m_1, probs = c(0.05, 0.5, 0.95),
         re_formula = 1 ~ 1)
##     Estimate   Est.Error        Q5       Q50       Q95
## R2 0.8587836 0.003533339 0.8526202 0.8589903 0.8642242
brms::bayes_R2(m_2, probs = c(0.05, 0.5, 0.95),
         re_formula = 1 ~ 1)
##    Estimate  Est.Error      Q5      Q50       Q95
## R2 0.479069 0.04335869 0.41154 0.477332 0.5529388
```

### Bayes factor (BARG 3.C)

``` r
brms::bayes_factor(m_1, m_2)
```

    ## Iteration: 1
    ## Iteration: 2
    ## Iteration: 3
    ## Iteration: 4
    ## Iteration: 5
    ## Iteration: 1
    ## Iteration: 2
    ## Iteration: 3
    ## Iteration: 4
    ## Iteration: 5
    ## Iteration: 6
    ## Iteration: 7
    ## Iteration: 8
    ## Iteration: 9
    ## Iteration: 10
    ## Iteration: 11

    ## Estimated Bayes factor in favor of m_1 over m_2:    Inf

## Posterior distributions (BARG 3.B)

### Forest plot

Chosen model:

``` r
mcmc_intervals(
  posterior1,
  prob = 0.66,
  prob_outer = 0.95,
  point_est = "mean"
  ) +
  theme_classic()
```

![](model_persistence_check_files/figure-gfm/posteriors-chosen-model-1.png)<!-- -->

Second model:

``` r
mcmc_intervals(
  posterior2,
  prob = 0.66,
  prob_outer = 0.95,
  point_est = "mean"
  ) +
  theme_classic()
```

![](model_persistence_check_files/figure-gfm/posterors-2nd-model-1.png)<!-- -->

### Effect sizes

Effect sizes of chosen model just to get exact values of means etc. if
necessary.

``` r
(emm <- emmeans(m_1, revpairwise ~ target_type + sand_ratio |
                  exposition | survey_year_fct, type = "response"))
```

    ## NOTE: A nesting structure was detected in the fitted model:
    ##     botanist_year %in% (exposition*survey_year_fct)

    ## $emmeans
    ## exposition = north, survey_year_fct = 2018:
    ##  target_type   sand_ratio emmean lower.HPD upper.HPD
    ##  hay_meadow    0           0.475     0.443     0.509
    ##  dry_grassland 0           0.474     0.440     0.508
    ##  hay_meadow    25          0.628     0.594     0.661
    ##  dry_grassland 25          0.560     0.525     0.593
    ##  hay_meadow    50          0.589     0.555     0.623
    ##  dry_grassland 50          0.555     0.522     0.590
    ## 
    ## exposition = south, survey_year_fct = 2018:
    ##  target_type   sand_ratio emmean lower.HPD upper.HPD
    ##  hay_meadow    0           0.306     0.272     0.340
    ##  dry_grassland 0           0.329     0.296     0.364
    ##  hay_meadow    25          0.311     0.277     0.346
    ##  dry_grassland 25          0.368     0.334     0.402
    ##  hay_meadow    50          0.281     0.247     0.315
    ##  dry_grassland 50          0.318     0.283     0.351
    ## 
    ## exposition = north, survey_year_fct = 2019:
    ##  target_type   sand_ratio emmean lower.HPD upper.HPD
    ##  hay_meadow    0           0.806     0.771     0.839
    ##  dry_grassland 0           0.763     0.728     0.797
    ##  hay_meadow    25          0.837     0.803     0.871
    ##  dry_grassland 25          0.755     0.721     0.789
    ##  hay_meadow    50          0.788     0.754     0.823
    ##  dry_grassland 50          0.761     0.727     0.795
    ## 
    ## exposition = south, survey_year_fct = 2019:
    ##  target_type   sand_ratio emmean lower.HPD upper.HPD
    ##  hay_meadow    0           0.457     0.420     0.491
    ##  dry_grassland 0           0.434     0.399     0.468
    ##  hay_meadow    25          0.460     0.425     0.496
    ##  dry_grassland 25          0.450     0.415     0.485
    ##  hay_meadow    50          0.421     0.387     0.459
    ##  dry_grassland 50          0.419     0.383     0.454
    ## 
    ## exposition = north, survey_year_fct = 2020:
    ##  target_type   sand_ratio emmean lower.HPD upper.HPD
    ##  hay_meadow    0           0.826     0.792     0.861
    ##  dry_grassland 0           0.796     0.762     0.830
    ##  hay_meadow    25          0.852     0.818     0.885
    ##  dry_grassland 25          0.805     0.771     0.839
    ##  hay_meadow    50          0.841     0.807     0.875
    ##  dry_grassland 50          0.800     0.766     0.834
    ## 
    ## exposition = south, survey_year_fct = 2020:
    ##  target_type   sand_ratio emmean lower.HPD upper.HPD
    ##  hay_meadow    0           0.517     0.483     0.552
    ##  dry_grassland 0           0.519     0.484     0.552
    ##  hay_meadow    25          0.543     0.509     0.576
    ##  dry_grassland 25          0.549     0.514     0.581
    ##  hay_meadow    50          0.501     0.468     0.536
    ##  dry_grassland 50          0.555     0.520     0.588
    ## 
    ## exposition = north, survey_year_fct = 2021:
    ##  target_type   sand_ratio emmean lower.HPD upper.HPD
    ##  hay_meadow    0           0.771     0.737     0.805
    ##  dry_grassland 0           0.756     0.722     0.790
    ##  hay_meadow    25          0.819     0.785     0.852
    ##  dry_grassland 25          0.785     0.752     0.819
    ##  hay_meadow    50          0.821     0.788     0.856
    ##  dry_grassland 50          0.807     0.774     0.842
    ## 
    ## exposition = south, survey_year_fct = 2021:
    ##  target_type   sand_ratio emmean lower.HPD upper.HPD
    ##  hay_meadow    0           0.528     0.494     0.562
    ##  dry_grassland 0           0.495     0.462     0.529
    ##  hay_meadow    25          0.486     0.452     0.520
    ##  dry_grassland 25          0.508     0.473     0.542
    ##  hay_meadow    50          0.512     0.479     0.547
    ##  dry_grassland 50          0.518     0.483     0.552
    ## 
    ## Results are averaged over the levels of: substrate_depth, seed_density, botanist_year 
    ## Point estimate displayed: median 
    ## HPD interval probability: 0.95 
    ## 
    ## $contrasts
    ## exposition = north, survey_year_fct = 2018:
    ##  contrast                                                estimate lower.HPD
    ##  dry_grassland sand_ratio0 - hay_meadow sand_ratio0      -0.00121 -0.038258
    ##  hay_meadow sand_ratio25 - hay_meadow sand_ratio0         0.15254  0.115522
    ##  hay_meadow sand_ratio25 - dry_grassland sand_ratio0      0.15367  0.114096
    ##  dry_grassland sand_ratio25 - hay_meadow sand_ratio0      0.08478  0.044659
    ##  dry_grassland sand_ratio25 - dry_grassland sand_ratio0   0.08578  0.047530
    ##  dry_grassland sand_ratio25 - hay_meadow sand_ratio25    -0.06773 -0.105646
    ##  hay_meadow sand_ratio50 - hay_meadow sand_ratio0         0.11393  0.075743
    ##  hay_meadow sand_ratio50 - dry_grassland sand_ratio0      0.11530  0.076369
    ##  hay_meadow sand_ratio50 - hay_meadow sand_ratio25       -0.03839 -0.077614
    ##  hay_meadow sand_ratio50 - dry_grassland sand_ratio25     0.02940 -0.011089
    ##  dry_grassland sand_ratio50 - hay_meadow sand_ratio0      0.08034  0.041490
    ##  dry_grassland sand_ratio50 - dry_grassland sand_ratio0   0.08136  0.043704
    ##  dry_grassland sand_ratio50 - hay_meadow sand_ratio25    -0.07241 -0.111836
    ##  dry_grassland sand_ratio50 - dry_grassland sand_ratio25 -0.00431 -0.045486
    ##  dry_grassland sand_ratio50 - hay_meadow sand_ratio50    -0.03381 -0.073280
    ##  upper.HPD
    ##   3.62e-02
    ##   1.90e-01
    ##   1.92e-01
    ##   1.23e-01
    ##   1.25e-01
    ##  -2.81e-02
    ##   1.50e-01
    ##   1.54e-01
    ##   3.89e-05
    ##   6.83e-02
    ##   1.20e-01
    ##   1.20e-01
    ##  -3.29e-02
    ##   3.40e-02
    ##   3.86e-03
    ## 
    ## exposition = south, survey_year_fct = 2018:
    ##  contrast                                                estimate lower.HPD
    ##  dry_grassland sand_ratio0 - hay_meadow sand_ratio0       0.02308 -0.014449
    ##  hay_meadow sand_ratio25 - hay_meadow sand_ratio0         0.00480 -0.032905
    ##  hay_meadow sand_ratio25 - dry_grassland sand_ratio0     -0.01843 -0.058093
    ##  dry_grassland sand_ratio25 - hay_meadow sand_ratio0      0.06218  0.023127
    ##  dry_grassland sand_ratio25 - dry_grassland sand_ratio0   0.03904  0.000809
    ##  dry_grassland sand_ratio25 - hay_meadow sand_ratio25     0.05749  0.019390
    ##  hay_meadow sand_ratio50 - hay_meadow sand_ratio0        -0.02554 -0.065793
    ##  hay_meadow sand_ratio50 - dry_grassland sand_ratio0     -0.04868 -0.088035
    ##  hay_meadow sand_ratio50 - hay_meadow sand_ratio25       -0.03037 -0.070416
    ##  hay_meadow sand_ratio50 - dry_grassland sand_ratio25    -0.08762 -0.128828
    ##  dry_grassland sand_ratio50 - hay_meadow sand_ratio0      0.01153 -0.027108
    ##  dry_grassland sand_ratio50 - dry_grassland sand_ratio0  -0.01139 -0.051381
    ##  dry_grassland sand_ratio50 - hay_meadow sand_ratio25     0.00686 -0.033259
    ##  dry_grassland sand_ratio50 - dry_grassland sand_ratio25 -0.05066 -0.090843
    ##  dry_grassland sand_ratio50 - hay_meadow sand_ratio50     0.03706 -0.002314
    ##  upper.HPD
    ##   6.21e-02
    ##   4.49e-02
    ##   2.21e-02
    ##   1.03e-01
    ##   7.83e-02
    ##   9.78e-02
    ##   1.22e-02
    ##  -9.20e-03
    ##   8.80e-03
    ##  -4.85e-02
    ##   5.18e-02
    ##   2.69e-02
    ##   4.68e-02
    ##  -1.12e-02
    ##   7.65e-02
    ## 
    ## exposition = north, survey_year_fct = 2019:
    ##  contrast                                                estimate lower.HPD
    ##  dry_grassland sand_ratio0 - hay_meadow sand_ratio0      -0.04306 -0.081595
    ##  hay_meadow sand_ratio25 - hay_meadow sand_ratio0         0.03161 -0.007356
    ##  hay_meadow sand_ratio25 - dry_grassland sand_ratio0      0.07477  0.034810
    ##  dry_grassland sand_ratio25 - hay_meadow sand_ratio0     -0.05066 -0.089682
    ##  dry_grassland sand_ratio25 - dry_grassland sand_ratio0  -0.00730 -0.048709
    ##  dry_grassland sand_ratio25 - hay_meadow sand_ratio25    -0.08197 -0.121885
    ##  hay_meadow sand_ratio50 - hay_meadow sand_ratio0        -0.01754 -0.056575
    ##  hay_meadow sand_ratio50 - dry_grassland sand_ratio0      0.02566 -0.012420
    ##  hay_meadow sand_ratio50 - hay_meadow sand_ratio25       -0.04913 -0.088726
    ##  hay_meadow sand_ratio50 - dry_grassland sand_ratio25     0.03308 -0.006338
    ##  dry_grassland sand_ratio50 - hay_meadow sand_ratio0     -0.04482 -0.085015
    ##  dry_grassland sand_ratio50 - dry_grassland sand_ratio0  -0.00125 -0.041111
    ##  dry_grassland sand_ratio50 - hay_meadow sand_ratio25    -0.07605 -0.114999
    ##  dry_grassland sand_ratio50 - dry_grassland sand_ratio25  0.00596 -0.033863
    ##  dry_grassland sand_ratio50 - hay_meadow sand_ratio50    -0.02682 -0.067930
    ##  upper.HPD
    ##  -2.46e-03
    ##   7.12e-02
    ##   1.15e-01
    ##  -9.95e-03
    ##   3.10e-02
    ##  -4.27e-02
    ##   2.27e-02
    ##   6.68e-02
    ##  -9.23e-03
    ##   7.38e-02
    ##  -4.59e-03
    ##   3.81e-02
    ##  -3.41e-02
    ##   4.65e-02
    ##   1.18e-02
    ## 
    ## exposition = south, survey_year_fct = 2019:
    ##  contrast                                                estimate lower.HPD
    ##  dry_grassland sand_ratio0 - hay_meadow sand_ratio0      -0.02331 -0.063724
    ##  hay_meadow sand_ratio25 - hay_meadow sand_ratio0         0.00337 -0.035719
    ##  hay_meadow sand_ratio25 - dry_grassland sand_ratio0      0.02684 -0.012626
    ##  dry_grassland sand_ratio25 - hay_meadow sand_ratio0     -0.00697 -0.047639
    ##  dry_grassland sand_ratio25 - dry_grassland sand_ratio0   0.01647 -0.023301
    ##  dry_grassland sand_ratio25 - hay_meadow sand_ratio25    -0.01026 -0.049240
    ##  hay_meadow sand_ratio50 - hay_meadow sand_ratio0        -0.03599 -0.076210
    ##  hay_meadow sand_ratio50 - dry_grassland sand_ratio0     -0.01261 -0.053031
    ##  hay_meadow sand_ratio50 - hay_meadow sand_ratio25       -0.03972 -0.078522
    ##  hay_meadow sand_ratio50 - dry_grassland sand_ratio25    -0.02928 -0.069563
    ##  dry_grassland sand_ratio50 - hay_meadow sand_ratio0     -0.03771 -0.077168
    ##  dry_grassland sand_ratio50 - dry_grassland sand_ratio0  -0.01428 -0.053714
    ##  dry_grassland sand_ratio50 - hay_meadow sand_ratio25    -0.04113 -0.081395
    ##  dry_grassland sand_ratio50 - dry_grassland sand_ratio25 -0.03057 -0.070582
    ##  dry_grassland sand_ratio50 - hay_meadow sand_ratio50    -0.00156 -0.041375
    ##  upper.HPD
    ##   1.55e-02
    ##   4.43e-02
    ##   6.79e-02
    ##   3.16e-02
    ##   5.68e-02
    ##   3.07e-02
    ##   3.17e-03
    ##   2.74e-02
    ##   1.46e-03
    ##   1.02e-02
    ##   2.53e-03
    ##   2.64e-02
    ##  -1.74e-03
    ##   9.39e-03
    ##   3.81e-02
    ## 
    ## exposition = north, survey_year_fct = 2020:
    ##  contrast                                                estimate lower.HPD
    ##  dry_grassland sand_ratio0 - hay_meadow sand_ratio0      -0.02979 -0.069624
    ##  hay_meadow sand_ratio25 - hay_meadow sand_ratio0         0.02541 -0.012756
    ##  hay_meadow sand_ratio25 - dry_grassland sand_ratio0      0.05533  0.017420
    ##  dry_grassland sand_ratio25 - hay_meadow sand_ratio0     -0.02100 -0.059846
    ##  dry_grassland sand_ratio25 - dry_grassland sand_ratio0   0.00884 -0.031559
    ##  dry_grassland sand_ratio25 - hay_meadow sand_ratio25    -0.04647 -0.085380
    ##  hay_meadow sand_ratio50 - hay_meadow sand_ratio0         0.01491 -0.023821
    ##  hay_meadow sand_ratio50 - dry_grassland sand_ratio0      0.04471  0.005488
    ##  hay_meadow sand_ratio50 - hay_meadow sand_ratio25       -0.01050 -0.051450
    ##  hay_meadow sand_ratio50 - dry_grassland sand_ratio25     0.03596 -0.003826
    ##  dry_grassland sand_ratio50 - hay_meadow sand_ratio0     -0.02587 -0.066463
    ##  dry_grassland sand_ratio50 - dry_grassland sand_ratio0   0.00400 -0.035031
    ##  dry_grassland sand_ratio50 - hay_meadow sand_ratio25    -0.05160 -0.092796
    ##  dry_grassland sand_ratio50 - dry_grassland sand_ratio25 -0.00500 -0.044645
    ##  dry_grassland sand_ratio50 - hay_meadow sand_ratio50    -0.04069 -0.081325
    ##  upper.HPD
    ##   9.08e-03
    ##   6.59e-02
    ##   9.58e-02
    ##   1.89e-02
    ##   4.73e-02
    ##  -6.74e-03
    ##   5.43e-02
    ##   8.50e-02
    ##   2.84e-02
    ##   7.54e-02
    ##   1.35e-02
    ##   4.39e-02
    ##  -1.36e-02
    ##   3.53e-02
    ##  -1.81e-03
    ## 
    ## exposition = south, survey_year_fct = 2020:
    ##  contrast                                                estimate lower.HPD
    ##  dry_grassland sand_ratio0 - hay_meadow sand_ratio0       0.00200 -0.038127
    ##  hay_meadow sand_ratio25 - hay_meadow sand_ratio0         0.02511 -0.014004
    ##  hay_meadow sand_ratio25 - dry_grassland sand_ratio0      0.02314 -0.016467
    ##  dry_grassland sand_ratio25 - hay_meadow sand_ratio0      0.03114 -0.009554
    ##  dry_grassland sand_ratio25 - dry_grassland sand_ratio0   0.02911 -0.009612
    ##  dry_grassland sand_ratio25 - hay_meadow sand_ratio25     0.00587 -0.033561
    ##  hay_meadow sand_ratio50 - hay_meadow sand_ratio0        -0.01620 -0.056104
    ##  hay_meadow sand_ratio50 - dry_grassland sand_ratio0     -0.01798 -0.057032
    ##  hay_meadow sand_ratio50 - hay_meadow sand_ratio25       -0.04116 -0.080969
    ##  hay_meadow sand_ratio50 - dry_grassland sand_ratio25    -0.04718 -0.086195
    ##  dry_grassland sand_ratio50 - hay_meadow sand_ratio0      0.03713 -0.001576
    ##  dry_grassland sand_ratio50 - dry_grassland sand_ratio0   0.03508 -0.004652
    ##  dry_grassland sand_ratio50 - hay_meadow sand_ratio25     0.01190 -0.029110
    ##  dry_grassland sand_ratio50 - dry_grassland sand_ratio25  0.00608 -0.033030
    ##  dry_grassland sand_ratio50 - hay_meadow sand_ratio50     0.05325  0.013076
    ##  upper.HPD
    ##   4.11e-02
    ##   6.54e-02
    ##   6.27e-02
    ##   7.00e-02
    ##   6.94e-02
    ##   4.60e-02
    ##   2.33e-02
    ##   2.25e-02
    ##  -2.14e-03
    ##  -7.02e-03
    ##   7.74e-02
    ##   7.58e-02
    ##   5.05e-02
    ##   4.64e-02
    ##   9.26e-02
    ## 
    ## exposition = north, survey_year_fct = 2021:
    ##  contrast                                                estimate lower.HPD
    ##  dry_grassland sand_ratio0 - hay_meadow sand_ratio0      -0.01528 -0.055466
    ##  hay_meadow sand_ratio25 - hay_meadow sand_ratio0         0.04824  0.008960
    ##  hay_meadow sand_ratio25 - dry_grassland sand_ratio0      0.06354  0.022199
    ##  dry_grassland sand_ratio25 - hay_meadow sand_ratio0      0.01378 -0.027434
    ##  dry_grassland sand_ratio25 - dry_grassland sand_ratio0   0.02930 -0.011748
    ##  dry_grassland sand_ratio25 - hay_meadow sand_ratio25    -0.03420 -0.073114
    ##  hay_meadow sand_ratio50 - hay_meadow sand_ratio0         0.04985  0.010403
    ##  hay_meadow sand_ratio50 - dry_grassland sand_ratio0      0.06531  0.024960
    ##  hay_meadow sand_ratio50 - hay_meadow sand_ratio25        0.00168 -0.037895
    ##  hay_meadow sand_ratio50 - dry_grassland sand_ratio25     0.03600 -0.003188
    ##  dry_grassland sand_ratio50 - hay_meadow sand_ratio0      0.03620 -0.004649
    ##  dry_grassland sand_ratio50 - dry_grassland sand_ratio0   0.05148  0.013131
    ##  dry_grassland sand_ratio50 - hay_meadow sand_ratio25    -0.01190 -0.051967
    ##  dry_grassland sand_ratio50 - dry_grassland sand_ratio25  0.02235 -0.016554
    ##  dry_grassland sand_ratio50 - hay_meadow sand_ratio50    -0.01339 -0.052404
    ##  upper.HPD
    ##   2.35e-02
    ##   8.69e-02
    ##   1.01e-01
    ##   5.16e-02
    ##   6.82e-02
    ##   5.31e-03
    ##   8.92e-02
    ##   1.04e-01
    ##   4.03e-02
    ##   7.64e-02
    ##   7.53e-02
    ##   9.25e-02
    ##   2.77e-02
    ##   6.31e-02
    ##   2.61e-02
    ## 
    ## exposition = south, survey_year_fct = 2021:
    ##  contrast                                                estimate lower.HPD
    ##  dry_grassland sand_ratio0 - hay_meadow sand_ratio0      -0.03228 -0.070576
    ##  hay_meadow sand_ratio25 - hay_meadow sand_ratio0        -0.04210 -0.081571
    ##  hay_meadow sand_ratio25 - dry_grassland sand_ratio0     -0.00970 -0.047786
    ##  dry_grassland sand_ratio25 - hay_meadow sand_ratio0     -0.01966 -0.059556
    ##  dry_grassland sand_ratio25 - dry_grassland sand_ratio0   0.01265 -0.026055
    ##  dry_grassland sand_ratio25 - hay_meadow sand_ratio25     0.02243 -0.018149
    ##  hay_meadow sand_ratio50 - hay_meadow sand_ratio0        -0.01562 -0.053881
    ##  hay_meadow sand_ratio50 - dry_grassland sand_ratio0      0.01673 -0.022084
    ##  hay_meadow sand_ratio50 - hay_meadow sand_ratio25        0.02665 -0.011583
    ##  hay_meadow sand_ratio50 - dry_grassland sand_ratio25     0.00402 -0.036416
    ##  dry_grassland sand_ratio50 - hay_meadow sand_ratio0     -0.00932 -0.049133
    ##  dry_grassland sand_ratio50 - dry_grassland sand_ratio0   0.02294 -0.017448
    ##  dry_grassland sand_ratio50 - hay_meadow sand_ratio25     0.03267 -0.008367
    ##  dry_grassland sand_ratio50 - dry_grassland sand_ratio25  0.01013 -0.031146
    ##  dry_grassland sand_ratio50 - hay_meadow sand_ratio50     0.00604 -0.034518
    ##  upper.HPD
    ##   8.39e-03
    ##  -2.56e-03
    ##   3.06e-02
    ##   2.02e-02
    ##   5.32e-02
    ##   6.15e-02
    ##   2.42e-02
    ##   5.67e-02
    ##   6.76e-02
    ##   4.28e-02
    ##   3.02e-02
    ##   6.19e-02
    ##   7.16e-02
    ##   4.93e-02
    ##   4.51e-02
    ## 
    ## Results are averaged over the levels of: substrate_depth, seed_density, botanist_year 
    ## Point estimate displayed: median 
    ## HPD interval probability: 0.95

# Session info (BARG 2.A/6.A/6.B)

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
    ##  [1] emmeans_1.8.3             tidybayes_3.0.2          
    ##  [3] loo_2.5.1                 bayesplot_1.10.0         
    ##  [5] DHARMa.helpers_0.0.0.9000 brms_2.18.0              
    ##  [7] Rcpp_1.0.9                patchwork_1.1.2          
    ##  [9] ggbeeswarm_0.6.0          forcats_0.5.2            
    ## [11] stringr_1.5.0             dplyr_1.0.10             
    ## [13] purrr_0.3.5               readr_2.1.3              
    ## [15] tidyr_1.2.1               tibble_3.1.8             
    ## [17] ggplot2_3.4.0             tidyverse_1.3.2          
    ## [19] here_1.0.1               
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
    ##  [70] deldir_1.0-6         zoo_1.8-11           cluster_2.1.4       
    ##  [73] haven_2.5.1          fs_1.5.2             magrittr_2.0.3      
    ##  [76] ggdist_3.2.0         colourpicker_1.2.0   reprex_2.0.2        
    ##  [79] googledrive_2.0.0    mvtnorm_1.1-3        matrixStats_0.63.0  
    ##  [82] hms_1.1.2            shinyjs_2.1.0        mime_0.12           
    ##  [85] evaluate_0.18        arrayhelpers_1.1-0   xtable_1.8-4        
    ##  [88] XML_3.99-0.13        shinystan_2.6.0      jpeg_0.1-10         
    ##  [91] readxl_1.4.1         gridExtra_2.3        rstantools_2.2.0    
    ##  [94] compiler_4.2.2       KernSmooth_2.23-20   V8_4.2.2            
    ##  [97] crayon_1.5.2         minqa_1.2.5          StanHeaders_2.26.13 
    ## [100] htmltools_0.5.3      mgcv_1.8-41          later_1.3.0         
    ## [103] tzdb_0.3.0           RcppParallel_5.1.5   lubridate_1.9.0     
    ## [106] DBI_1.1.3            dbplyr_2.2.1         MASS_7.3-58.1       
    ## [109] boot_1.3-28          Matrix_1.5-3         ade4_1.7-20         
    ## [112] permute_0.9-7        cli_3.4.1            adegraphics_1.0-16  
    ## [115] DHARMa_0.4.6         parallel_4.2.2       igraph_1.3.5        
    ## [118] pkgconfig_2.0.3      rncl_0.8.6           sp_1.5-1            
    ## [121] foreach_1.5.2        xml2_1.3.3           svUnit_1.0.6        
    ## [124] dygraphs_1.1.1.6     vipor_0.4.5          estimability_1.4.1  
    ## [127] rvest_1.0.3          distributional_0.3.1 callr_3.7.3         
    ## [130] digest_0.6.30        vegan_2.6-4          rmarkdown_2.18      
    ## [133] cellranger_1.1.0     gap_1.3-1            curl_4.3.3          
    ## [136] shiny_1.7.3          gtools_3.9.4         nloptr_2.0.3        
    ## [139] lifecycle_1.0.3      nlme_3.1-160         jsonlite_1.8.4      
    ## [142] seqinr_4.2-23        fansi_1.0.3          pillar_1.8.1        
    ## [145] lattice_0.20-45      fastmap_1.1.0        httr_1.4.4          
    ## [148] pkgbuild_1.4.0       glue_1.6.2           xts_0.12.2          
    ## [151] iterators_1.0.14     png_0.1-8            shinythemes_1.2.0   
    ## [154] bit_4.0.5            stringi_1.7.8        latticeExtra_0.6-30 
    ## [157] ape_5.6-2
