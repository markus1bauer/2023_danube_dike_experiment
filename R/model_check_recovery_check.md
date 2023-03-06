Analysis of Bauer et al. (2023) bioRxiv: <br> Recovery completeness
================
<b>Markus Bauer</b> <br>
<b>2023-03-03</b>

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
    - <a href="#preparation-for-analysis"
      id="toc-preparation-for-analysis">Preparation for analysis</a>
    - <a href="#priors-barg-1de" id="toc-priors-barg-1de">Priors (BARG
      1.D/E)</a>
  - <a href="#model-check" id="toc-model-check">Model check</a>
    - <a href="#check-computation-mcmc-diagnostics-barg-2bc"
      id="toc-check-computation-mcmc-diagnostics-barg-2bc">Check computation
      (MCMC diagnostics, BARG 2.B/C)</a>
    - <a href="#autocorrelation-check"
      id="toc-autocorrelation-check">Autocorrelation check</a>
    - <a href="#posterior-predictive-check-barg-3a"
      id="toc-posterior-predictive-check-barg-3a">Posterior predictive check
      (BARG 3.A)</a>
    - <a href="#dharma" id="toc-dharma">DHARMa</a>
  - <a href="#model-comparison" id="toc-model-comparison">Model
    comparison</a>
    - <a href="#conditional-r2-values"
      id="toc-conditional-r2-values">Conditional <i>R</i><sup>2</sup>
      values</a>
    - <a href="#marginal-r2-values" id="toc-marginal-r2-values">Marginal
      <i>R</i><sup>2</sup> values</a>
    - <a href="#bayes-factor-barg-3c" id="toc-bayes-factor-barg-3c">Bayes
      factor (BARG 3.C)</a>
  - <a href="#posterior-distributions-barg-3b"
    id="toc-posterior-distributions-barg-3b">Posterior distributions (BARG
    3.B)</a>
    - <a href="#forest-plot-barg-3b5b" id="toc-forest-plot-barg-3b5b">Forest
      plot (BARG 3.B/5.B)</a>
    - <a href="#effect-sizes" id="toc-effect-sizes">Effect sizes</a>
- <a href="#session-info-barg-2a6a6b"
  id="toc-session-info-barg-2a6a6b">Session info (BARG 2.A/6.A/6.B)</a>

<br/> <br/> <b>Markus Bauer</b>

Technichal University of Munich, TUM School of Life Sciences, Chair of
Restoration Ecology, Emil-Ramann-Straße 6, 85354 Freising, Germany

<markus1.bauer@tum.de>

ORCiD ID: [0000-0001-5372-4174](https://orcid.org/0000-0001-5372-4174)
<br> [Google
Scholar](https://scholar.google.de/citations?user=oHhmOkkAAAAJ&hl=de&oi=ao)
<br> GitHub: [markus1bauer](https://github.com/markus1bauer)

To compare different models, you only have to change the models in
section `Load models`

# Preparation

Recovery completeness sensu Rydgren et al. (2019) J Appl Ecol [DOI:
10.1111/1365-2664.13254](https://doi.org/10.1111/1365-2664.13254)

Bayesian analysis motivated by Applestein et al. (2021) Restor Ecol
[DOI: 10.1111/rec.13596](https://doi.org/10.1111/rec.13596)

Analysis guided by <br> <b>BARG</b> (Bayesian Analysis Reporting
Guidelines): Kruschke (2021) Nat Hum Behav [DOI:
10.1038/s41562-021-01177-7](https://doi.org/10.1038/s41562-021-01177-7)
<br> Model check: Gabry et al. (2019) J R Stat Soc A Stat [DOI:
10.1111/rssa.12378](https://doi.org/10.1111/rssa.12378) <br> Priors:
Lemoine (2019) Oikos [DOI:
10.1111/oik.05985](https://doi.org/10.1111/oik.05985) <br> Model check:
Deapoli & Schoot (2017) Psychol Methods [DOI:
10.1037/met0000065](https://doi.org/10.1037/met0000065)

#### Packages

``` r
library(here)
library(tidyverse)
library(ggbeeswarm)
library(patchwork)
library(brms)
library(DHARMa.helpers)
library(bayesplot)
library(loo)
library(emmeans)
```

#### Load data

``` r
sites <- read_csv(
  here("data", "processed", "data_processed_sites_nmds.csv"),
  col_names = TRUE, na = c("na", "NA", ""),
  col_types = cols(
    .default = "?",
    id = "f",
    plot = "f",
    site = "f",
    exposition = col_factor(levels = c("north", "south", "other")),
    sand_ratio = "f",
    substrate_depth = col_factor(levels = c("30", "15")),
    target_type = col_factor(levels = c(
      "hay_meadow", "dry_grassland", "other"
      )),
    seed_density = "f"
    )
  ) %>%
  filter(reference == "2018" | reference == "2019" | reference == "2020" |
           reference == "2021") %>%
  mutate(
    survey_year_fct = factor(survey_year),
    survey_year = as.numeric(survey_year),
    botanist_year = str_c(survey_year, botanist, exposition, sep = " "),
    botanist_year = factor(botanist_year),
    n = recovery_completeness
    )
sites
```

    ## # A tibble: 1,152 × 26
    ##    id         plot  site  esy   refere…¹ expos…² sand_…³ subst…⁴ targe…⁵ seed_…⁶
    ##    <fct>      <fct> <fct> <chr> <chr>    <fct>   <fct>   <fct>   <fct>   <fct>  
    ##  1 L1_01_2018 L1_01 1     V37   2018     south   0       15      hay_me… 4      
    ##  2 L1_01_2019 L1_01 1     V37   2019     south   0       15      hay_me… 4      
    ##  3 L1_01_2020 L1_01 1     V37   2020     south   0       15      hay_me… 4      
    ##  4 L1_01_2021 L1_01 1     R22   2021     south   0       15      hay_me… 4      
    ##  5 L1_02_2018 L1_02 1     V37   2018     south   0       15      dry_gr… 4      
    ##  6 L1_02_2019 L1_02 1     V37   2019     south   0       15      dry_gr… 4      
    ##  7 L1_02_2020 L1_02 1     V37   2020     south   0       15      dry_gr… 4      
    ##  8 L1_02_2021 L1_02 1     R22   2021     south   0       15      dry_gr… 4      
    ##  9 L1_03_2018 L1_03 1     V37   2018     south   0       15      dry_gr… 8      
    ## 10 L1_03_2019 L1_03 1     R     2019     south   0       15      dry_gr… 8      
    ## # … with 1,142 more rows, 16 more variables: survey_year <dbl>,
    ## #   longitude <dbl>, latitude <dbl>, elevation <dbl>, plot_size <dbl>,
    ## #   botanist <chr>, givd_id <chr>, source <chr>, NMDS1 <dbl>, NMDS2 <dbl>,
    ## #   mean_reference <dbl>, sd_reference <dbl>, recovery_completeness <dbl>,
    ## #   survey_year_fct <fct>, botanist_year <fct>, n <dbl>, and abbreviated
    ## #   variable names ¹​reference, ²​exposition, ³​sand_ratio, ⁴​substrate_depth,
    ## #   ⁵​target_type, ⁶​seed_density

# Statistics

## Data exploration

### Graphs of raw data

![](model_check_recovery_check_files/figure-gfm/data-exploration-1.png)<!-- -->![](model_check_recovery_check_files/figure-gfm/data-exploration-2.png)<!-- -->![](model_check_recovery_check_files/figure-gfm/data-exploration-3.png)<!-- -->![](model_check_recovery_check_files/figure-gfm/data-exploration-4.png)<!-- -->![](model_check_recovery_check_files/figure-gfm/data-exploration-5.png)<!-- -->

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

![](model_check_recovery_check_files/figure-gfm/outliers-1.png)<!-- -->

## Models

### Load models (BARG 6.G)

Only here you have to modify the script to compare other models

``` r
load(file = here("outputs", "models", "model_recovery_2.Rdata"))
load(file = here("outputs", "models", "model_recovery_full.Rdata"))
load(file = here("outputs", "models", "model_recovery_2_prior.Rdata"))
# BARG 5.A/B/C
load(file = here("outputs", "models", "model_recovery_2_flat.Rdata"))
m_1 <- m2
m_2 <- m_full
m_prior <- m2_prior
m_flat <- m2_flat
```

### Model specifications (BARG 1.D)

``` r
m_1$formula
## n ~ (target_type + exposition + sand_ratio + survey_year_fct)^4 + substrate_depth + seed_density + sand_ratio:substrate_depth + substrate_depth:exposition + seed_density:exposition + (1 | site/plot) + (1 | botanist_year)
m_2$formula
## n ~ sand_ratio * target_type * exposition * survey_year_fct + substrate_depth * seed_density + substrate_depth:exposition + seed_density:exposition + substrate_depth:survey_year_fct + seed_density:survey_year_fct + substrate_depth:exposition:survey_year_fct + seed_density:exposition:survey_year_fct + botanist_year + (1 | site/plot)
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

Markov Chain Monte Carlo (MCMC) method used with No-U-Turn Sampler
(NUTS) which is an extension of Hamiltonian Monte Carlo (HMC).

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
## [1] 10000
m_2$fit@sim$iter
## [1] 10000
```

Amount of iterations before burn-in

``` r
m_1$fit@sim$warmup
## [1] 5000
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

### Preparation for analysis

``` r
# Chose variables
variables <- c(
      "Intercept",
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
# Subset draws
posterior1 <- m_1 %>%
  posterior::as_draws() %>%
  posterior::subset_draws(variable = variables)
posterior2 <- m_2 %>%
  posterior::as_draws() %>%
  posterior::subset_draws(variable = variables)
posterior_prior <- m_prior %>%
  posterior::as_draws() %>%
  posterior::subset_draws(variable = variables)
posterior_flat <- m_flat %>%
  posterior::as_draws() %>%
  posterior::subset_draws(variable = variables)
# Long format of draws
hmc_diagnostics1 <- m_1 %>% brms::nuts_params()
hmc_diagnostics2 <- m_2 %>% brms::nuts_params()
y <- sites$n
# Posterior predictive distribution
yrep1 <- m_1 %>% brms::posterior_predict(draws = 500)
yrep2 <- m_2 %>% brms::posterior_predict(draws = 500)
yrep_prior <- m_prior %>% brms::posterior_predict(draws = 500)
# Summary statistics
draws1 <- m_1 %>%
  posterior::as_draws() %>%
  posterior::summarize_draws() %>%
  filter(str_starts(variable, "b_"))
draws2 <- m_2 %>%
  posterior::as_draws() %>%
  posterior::summarize_draws() %>%
  filter(str_starts(variable, "b_"))
```

### Priors (BARG 1.D/E)

#### Possible prior distributions

``` r
brms::get_prior(
  n ~ target_type + sand_ratio + seed_density + substrate_depth +
    exposition + survey_year_fct + botanist_year + (1 | site / plot),
  data = sites
  ) %>%
  select(prior, class, coef, group, source)
```

    ##                    prior     class                                coef
    ##                   (flat)         b                                    
    ##                   (flat)         b    botanist_year2018JakobHubersouth
    ##                   (flat)         b    botanist_year2018SimonReithnorth
    ##                   (flat)         b    botanist_year2018SimonReithsouth
    ##                   (flat)         b    botanist_year2019JakobHubernorth
    ##                   (flat)         b  botanist_year2019LindaWegglersouth
    ##                   (flat)         b   botanist_year2019MarkusBauersouth
    ##                   (flat)         b    botanist_year2020JakobHubernorth
    ##                   (flat)         b    botanist_year2020JakobHubersouth
    ##                   (flat)         b botanist_year2020KatharinaBecknorth
    ##                   (flat)         b botanist_year2020KatharinaBecksouth
    ##                   (flat)         b   botanist_year2021MarkusBauernorth
    ##                   (flat)         b   botanist_year2021MarkusBauersouth
    ##                   (flat)         b                     expositionsouth
    ##                   (flat)         b                        sand_ratio25
    ##                   (flat)         b                        sand_ratio50
    ##                   (flat)         b                       seed_density8
    ##                   (flat)         b                   substrate_depth15
    ##                   (flat)         b                 survey_year_fct2019
    ##                   (flat)         b                 survey_year_fct2020
    ##                   (flat)         b                 survey_year_fct2021
    ##                   (flat)         b            target_typedry_grassland
    ##  student_t(3, -0.7, 2.5) Intercept                                    
    ##     student_t(3, 0, 2.5)        sd                                    
    ##                   (flat)        sd                                    
    ##                   (flat)        sd                           Intercept
    ##                   (flat)        sd                                    
    ##                   (flat)        sd                           Intercept
    ##     student_t(3, 0, 2.5)     sigma                                    
    ##      group  source
    ##            default
    ##            default
    ##            default
    ##            default
    ##            default
    ##            default
    ##            default
    ##            default
    ##            default
    ##            default
    ##            default
    ##            default
    ##            default
    ##            default
    ##            default
    ##            default
    ##            default
    ##            default
    ##            default
    ##            default
    ##            default
    ##            default
    ##            default
    ##            default
    ##       site default
    ##       site default
    ##  site:plot default
    ##  site:plot default
    ##            default

``` r
plot1 <- ggplot(data = data.frame(x = c(-2, 0)), aes(x = x)) +
  stat_function(fun = dnorm, n = 101, args = list(mean = -1, sd = 1)) +
  expand_limits(y = 0) +
  ggtitle("Normal distribution for Intecept") + theme_mb()
data <- data.frame(x = c(-2, 2))
plot2 <- ggplot(data, aes(x = x)) +
  stat_function(fun = dnorm, n = 101, args = list(mean = 0, sd = .8)) +
  expand_limits(y = 0) +
  ggtitle("Normal distribution for treatments") + theme_mb()
# See Lemoine 2019 https://doi.org/10.1111/oik.05985
plot3 <- ggplot(data, aes(x = x)) +
  stat_function(fun = dcauchy, n = 101, args = list(location = 0, scale = 1)) +
  expand_limits(y = 0) +
  ggtitle("Cauchy distribution") + theme_mb()
# Software standard
plot4 <- ggplot(data, aes(x = x)) +
  stat_function(fun = dstudent_t, args = list(df = 3, mu = 0, sigma = 2.5)) +
  expand_limits(y = 0) +
  ggtitle(expression(Student ~ italic(t) * "-distribution")) + theme_mb()
(plot1 + plot2) / (plot3 + plot4)
```

![](model_check_recovery_check_files/figure-gfm/possible-priors-1.png)<!-- -->

#### Prior summary (BARG 1.D)

``` r
m_1 %>%
  brms::prior_summary(all = FALSE) %>%
  select(prior, class, coef, group, source)
```

    ##                    prior     class                coef group  source
    ##             normal(0, 2)         b                              user
    ##          normal(-0.1, 2)         b     expositionsouth          user
    ##           normal(0.1, 2)         b        sand_ratio25          user
    ##           normal(0.2, 2)         b        sand_ratio50          user
    ##           normal(0.1, 2)         b survey_year_fct2019          user
    ##           normal(0.2, 2)         b survey_year_fct2020          user
    ##           normal(0.3, 2)         b survey_year_fct2021          user
    ##  student_t(3, -0.8, 2.5) Intercept                           default
    ##     student_t(3, 0, 2.5)        sd                           default
    ##             cauchy(0, 1)     sigma                              user

#### Prior predictive check (BARG 1.E)

``` r
bayesplot::ppc_stat(y, yrep_prior[1:500, ], binwidth = 0.5) +
  coord_cartesian(xlim = c(-4, 4)) + bayesplot::theme_default()
```

![](model_check_recovery_check_files/figure-gfm/prior-predictive-check-1.png)<!-- -->

``` r
ppc_stat_grouped(
  y, yrep_prior[1:500, ], group = sites$site, binwidth = 0.5
  ) +
  coord_cartesian(xlim = c(-4, 4)) + bayesplot::theme_default()
```

![](model_check_recovery_check_files/figure-gfm/prior-predictive-check-2.png)<!-- -->

``` r
ppc_stat_grouped(
  y, yrep_prior[1:500, ], group = sites$exposition, binwidth = 0.5
  ) +  coord_cartesian(xlim = c(-4, 4)) + bayesplot::theme_default()
```

![](model_check_recovery_check_files/figure-gfm/prior-predictive-check-3.png)<!-- -->

``` r
ppc_stat_grouped(
  y, yrep_prior[1:500, ], group = sites$survey_year_fct, binwidth = 0.5
  ) +
  coord_cartesian(xlim = c(-4, 4)) + bayesplot::theme_default()
```

![](model_check_recovery_check_files/figure-gfm/prior-predictive-check-4.png)<!-- -->

``` r
ppc_stat_grouped(
  y, yrep_prior[1:500, ], group = sites$target_type, binwidth = 0.5
  ) +
  coord_cartesian(xlim = c(-4, 4)) + bayesplot::theme_default()
```

![](model_check_recovery_check_files/figure-gfm/prior-predictive-check-5.png)<!-- -->

``` r
ppc_stat_grouped(
  y, yrep_prior[1:500, ], group = sites$seed_density, binwidth = 0.5
  ) +
  coord_cartesian(xlim = c(-4, 4)) + bayesplot::theme_default()
```

![](model_check_recovery_check_files/figure-gfm/prior-predictive-check-6.png)<!-- -->

``` r
ppc_stat_grouped(
  y, yrep_prior[1:500, ], group = sites$sand_ratio, binwidth = 0.5
  ) +
  coord_cartesian(xlim = c(-4, 4)) + bayesplot::theme_default()
```

![](model_check_recovery_check_files/figure-gfm/prior-predictive-check-7.png)<!-- -->

``` r
ppc_stat_grouped(
  y, yrep_prior[1:500, ], group = sites$substrate_depth, binwidth = 0.5
  ) +
  coord_cartesian(xlim = c(-4, 4)) + bayesplot::theme_default()
```

![](model_check_recovery_check_files/figure-gfm/prior-predictive-check-8.png)<!-- -->

## Model check

### Check computation (MCMC diagnostics, BARG 2.B/C)

#### Trace plots

``` r
bayesplot::mcmc_trace(
  posterior1, np = hmc_diagnostics1, facet_args = list(ncol = 2)
  ) + bayesplot::theme_default()
```

![](model_check_recovery_check_files/figure-gfm/mcmc-trace-1.png)<!-- -->

``` r
bayesplot::mcmc_trace(
  posterior2, np = hmc_diagnostics2, facet_args = list(ncol = 2)
  ) + bayesplot::theme_default()
```

    ## No divergences to plot.

![](model_check_recovery_check_files/figure-gfm/mcmc-trace-2.png)<!-- -->

#### Sampling efficency: R-hat (BARG 2.B)

``` r
m_1 %>% brms::rhat() %>% bayesplot::mcmc_rhat() + bayesplot::theme_default() + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
```

![](model_check_recovery_check_files/figure-gfm/rhat-1.png)<!-- -->

``` r
m_2 %>% brms::rhat() %>% bayesplot::mcmc_rhat() + bayesplot::theme_default() + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
```

![](model_check_recovery_check_files/figure-gfm/rhat-2.png)<!-- -->

#### Sampling effectiveness: Effective sampling size (ESS) (BARG 2.C)

``` r
m_1 %>% brms::neff_ratio() %>% bayesplot::mcmc_neff() + bayesplot::theme_default() + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
```

![](model_check_recovery_check_files/figure-gfm/ess-1.png)<!-- -->

``` r
m_2 %>% brms::neff_ratio() %>% bayesplot::mcmc_neff() + bayesplot::theme_default() + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
```

![](model_check_recovery_check_files/figure-gfm/ess-2.png)<!-- -->

#### Pairs plot

``` r
m_1 %>% bayesplot::mcmc_pairs(
  off_diag_args = list(size = 1.2),
  pars = c(
    "b_sand_ratio25", "b_sand_ratio50", "b_substrate_depth15",
    "b_target_typedry_grassland", "b_seed_density8",
    "b_expositionsouth", "sigma"
           )
)
```

![](model_check_recovery_check_files/figure-gfm/mcmc-pairs-1.png)<!-- -->

``` r
m_2 %>% bayesplot::mcmc_pairs(
  off_diag_args = list(size = 1.2),
  pars = c(
    "b_sand_ratio25", "b_sand_ratio50", "b_substrate_depth15",
    "b_target_typedry_grassland", "b_seed_density8",
    "b_expositionsouth", "sigma"
           )
)
```

![](model_check_recovery_check_files/figure-gfm/mcmc-pairs-2.png)<!-- -->

#### Parallel coordinate plot

``` r
posterior1 %>% bayesplot::mcmc_parcoord(np = hmc_diagnostics1) +
  bayesplot::theme_default() + theme(axis.text.x = element_text(angle = 45))
```

![](model_check_recovery_check_files/figure-gfm/mcmc-parcoord-1.png)<!-- -->

``` r
posterior2 %>% bayesplot::mcmc_parcoord(np = hmc_diagnostics2) + 
  bayesplot::theme_default() + theme(axis.text.x = element_text(angle = 45))
```

![](model_check_recovery_check_files/figure-gfm/mcmc-parcoord-2.png)<!-- -->

### Autocorrelation check

``` r
posterior1 %>% bayesplot::mcmc_acf(lags = 10) + bayesplot::theme_default()
```

![](model_check_recovery_check_files/figure-gfm/autocorrelation-1.png)<!-- -->

``` r
posterior2 %>% bayesplot::mcmc_acf(lags = 10) + bayesplot::theme_default()
```

![](model_check_recovery_check_files/figure-gfm/autocorrelation-2.png)<!-- -->

### Posterior predictive check (BARG 3.A)

#### Kernel density

``` r
p1 <- bayesplot::ppc_dens_overlay(y, yrep1[1:50, ]) + bayesplot::theme_default()
p2 <- bayesplot::ppc_dens_overlay(y, yrep2[1:50, ]) + bayesplot::theme_default()
p1 / p2
```

![](model_check_recovery_check_files/figure-gfm/kernel-density-1.png)<!-- -->

``` r
p1 <- ppc_dens_overlay_grouped(y, yrep1[1:50, ], group = sites$site) + bayesplot::theme_default()
p2 <- ppc_dens_overlay_grouped(y, yrep2[1:50, ], group = sites$site) + bayesplot::theme_default()
p1 / p2
```

![](model_check_recovery_check_files/figure-gfm/kernel-density-2.png)<!-- -->

``` r
p1 <- ppc_dens_overlay_grouped(y, yrep1[1:50, ], group = sites$exposition) + bayesplot::theme_default()
p2 <- ppc_dens_overlay_grouped(y, yrep2[1:50, ], group = sites$exposition) + bayesplot::theme_default()
p1 / p2
```

![](model_check_recovery_check_files/figure-gfm/kernel-density-3.png)<!-- -->

``` r
p1 <- ppc_dens_overlay_grouped(y, yrep1[1:50, ], group = sites$survey_year_fct) + bayesplot::theme_default()
p2 <- ppc_dens_overlay_grouped(y, yrep2[1:50, ], group = sites$survey_year_fct) + bayesplot::theme_default()
p1 / p2
```

![](model_check_recovery_check_files/figure-gfm/kernel-density-4.png)<!-- -->

``` r
p1 <- ppc_dens_overlay_grouped(y, yrep1[1:50, ], group = sites$target_type) + bayesplot::theme_default()
p2 <- ppc_dens_overlay_grouped(y, yrep2[1:50, ], group = sites$target_type) + bayesplot::theme_default()
p1 / p2
```

![](model_check_recovery_check_files/figure-gfm/kernel-density-5.png)<!-- -->

``` r
p1 <- ppc_dens_overlay_grouped(y, yrep1[1:50, ], group = sites$seed_density) + bayesplot::theme_default()
p2 <- ppc_dens_overlay_grouped(y, yrep2[1:50, ], group = sites$seed_density) + bayesplot::theme_default()
p1 / p2
```

![](model_check_recovery_check_files/figure-gfm/kernel-density-6.png)<!-- -->

``` r
p1 <- ppc_dens_overlay_grouped(y, yrep1[1:50, ], group = sites$sand_ratio) + bayesplot::theme_default()
p2 <- ppc_dens_overlay_grouped(y, yrep2[1:50, ], group = sites$sand_ratio) + bayesplot::theme_default()
p1 / p2
```

![](model_check_recovery_check_files/figure-gfm/kernel-density-7.png)<!-- -->

``` r
p1 <- ppc_dens_overlay_grouped(y, yrep1[1:50, ], group = sites$substrate_depth) + bayesplot::theme_default()
p2 <- ppc_dens_overlay_grouped(y, yrep2[1:50, ], group = sites$substrate_depth) + bayesplot::theme_default()
p1 / p2
```

![](model_check_recovery_check_files/figure-gfm/kernel-density-8.png)<!-- -->

#### Histograms of statistics skew

``` r
p1 <- bayesplot::ppc_stat(y, yrep1, binwidth = 0.001) + bayesplot::theme_default()
p2 <- bayesplot::ppc_stat(y, yrep2, binwidth = 0.001) + bayesplot::theme_default()
p1 / p2
```

![](model_check_recovery_check_files/figure-gfm/histograms-1.png)<!-- -->

``` r
p1 <- ppc_stat_grouped(y, yrep1, group = sites$site, binwidth = 0.001) + bayesplot::theme_default()
p2 <- ppc_stat_grouped(y, yrep2, group = sites$site, binwidth = 0.001) + bayesplot::theme_default()
p1 / p2
```

![](model_check_recovery_check_files/figure-gfm/histograms-2.png)<!-- -->

``` r
p1 <- ppc_stat_grouped(y, yrep1, group = sites$exposition, binwidth = 0.001) + bayesplot::theme_default()
p2 <- ppc_stat_grouped(y, yrep2, group = sites$exposition, binwidth = 0.001) + bayesplot::theme_default()
p1 / p2
```

![](model_check_recovery_check_files/figure-gfm/histograms-3.png)<!-- -->

``` r
p1 <- ppc_stat_grouped(y, yrep1, group = sites$survey_year_fct, binwidth = 0.001) + bayesplot::theme_default()
p2 <- ppc_stat_grouped(y, yrep2, group = sites$survey_year_fct, binwidth = 0.001) + bayesplot::theme_default()
p1 / p2
```

![](model_check_recovery_check_files/figure-gfm/histograms-4.png)<!-- -->

``` r
p1 <- ppc_stat_grouped(y, yrep1, group = sites$target_type, binwidth = 0.001) + bayesplot::theme_default()
p2 <- ppc_stat_grouped(y, yrep2, group = sites$target_type, binwidth = 0.001) + bayesplot::theme_default()
p1 / p2
```

![](model_check_recovery_check_files/figure-gfm/histograms-5.png)<!-- -->

``` r
p1 <- ppc_stat_grouped(y, yrep1, group = sites$seed_density, binwidth = 0.001) + bayesplot::theme_default()
p2 <- ppc_stat_grouped(y, yrep2, group = sites$seed_density, binwidth = 0.001) + bayesplot::theme_default()
p1 / p2
```

![](model_check_recovery_check_files/figure-gfm/histograms-6.png)<!-- -->

``` r
p1 <- ppc_stat_grouped(y, yrep1, group = sites$sand_ratio, binwidth = 0.001) + bayesplot::theme_default()
p2 <- ppc_stat_grouped(y, yrep2, group = sites$sand_ratio, binwidth = 0.001) + bayesplot::theme_default()
p1 / p2
```

![](model_check_recovery_check_files/figure-gfm/histograms-7.png)<!-- -->

``` r
p1 <- ppc_stat_grouped(y, yrep1, group = sites$substrate_depth, binwidth = 0.001) + bayesplot::theme_default()
p2 <- ppc_stat_grouped(y, yrep2, group = sites$substrate_depth, binwidth = 0.001) + bayesplot::theme_default()
p1 / p2
```

![](model_check_recovery_check_files/figure-gfm/histograms-8.png)<!-- -->

#### LOO cross-validation (Leave one out)

``` r
# Leave-one-out cross validation based on posterior likelihood
loo1 <- m_1 %>% brms::loo(save_psis = TRUE, moment_match = TRUE)
```

    ## Warning: Some Pareto k diagnostic values are slightly high. See help('pareto-k-diagnostic') for details.

``` r
loo2 <- m_2 %>% brms::loo(save_psis = TRUE, moment_match = TRUE)
```

    ## Warning: Some Pareto k diagnostic values are slightly high. See help('pareto-k-diagnostic') for details.

``` r
loo1
```

    ## 
    ## Computed from 10000 by 1152 log-likelihood matrix
    ## 
    ##          Estimate   SE
    ## elpd_loo    690.1 26.0
    ## p_loo       147.2  6.2
    ## looic     -1380.2 52.1
    ## ------
    ## Monte Carlo SE of elpd_loo is 0.2.
    ## 
    ## Pareto k diagnostic values:
    ##                          Count Pct.    Min. n_eff
    ## (-Inf, 0.5]   (good)     1149  99.7%   516       
    ##  (0.5, 0.7]   (ok)          3   0.3%   754       
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
    ##          Estimate   SE
    ## elpd_loo    686.9 25.5
    ## p_loo       168.3  6.9
    ## looic     -1373.8 51.0
    ## ------
    ## Monte Carlo SE of elpd_loo is 0.2.
    ## 
    ## Pareto k diagnostic values:
    ##                          Count Pct.    Min. n_eff
    ## (-Inf, 0.5]   (good)     1143  99.2%   1097      
    ##  (0.5, 0.7]   (ok)          9   0.8%   661       
    ##    (0.7, 1]   (bad)         0   0.0%   <NA>      
    ##    (1, Inf)   (very bad)    0   0.0%   <NA>      
    ## 
    ## All Pareto k estimates are ok (k < 0.7).
    ## See help('pareto-k-diagnostic') for details.

``` r
plot(loo1)
```

![](model_check_recovery_check_files/figure-gfm/loo-1.png)<!-- -->

``` r
plot(loo2)
```

![](model_check_recovery_check_files/figure-gfm/loo-2.png)<!-- -->

Leave one out probability integral transform

``` r
p1 <- bayesplot::ppc_loo_pit_overlay(y, yrep1, lw = weights(loo1$psis_object)) + bayesplot::theme_default()
```

    ## NOTE: The kernel density estimate assumes continuous observations and is not optimal for discrete observations.

``` r
p2 <- bayesplot::ppc_loo_pit_overlay(y, yrep2, lw = weights(loo2$psis_object)) + bayesplot::theme_default()
```

    ## NOTE: The kernel density estimate assumes continuous observations and is not optimal for discrete observations.

``` r
p1 / p2
```

![](model_check_recovery_check_files/figure-gfm/loo-pit-1.png)<!-- -->

### DHARMa

``` r
m_1 %>% DHARMa.helpers::dh_check_brms(integer = TRUE)
```

![](model_check_recovery_check_files/figure-gfm/dharma-1.png)<!-- -->

``` r
m_2 %>% DHARMa.helpers::dh_check_brms(integer = TRUE)
```

![](model_check_recovery_check_files/figure-gfm/dharma-2.png)<!-- -->

## Model comparison

### Conditional <i>R</i><sup>2</sup> values

``` r
m_1 %>% brms::bayes_R2(
  probs = c(0.05, 0.5, 0.95), re_formula =  ~ (1 | site / plot)
  )
##     Estimate  Est.Error       Q5       Q50       Q95
## R2 0.9068055 0.01980042 0.871471 0.9128796 0.9205551
m_2 %>% brms::bayes_R2(
  probs = c(0.05, 0.5, 0.95), re_formula =  ~ (1 | site / plot)
  )
##    Estimate   Est.Error        Q5       Q50       Q95
## R2 0.918869 0.002526252 0.9145194 0.9189666 0.9227846
```

### Marginal <i>R</i><sup>2</sup> values

``` r
m_1 %>% brms::bayes_R2(
  probs = c(0.05, 0.5, 0.95), re_formula = 1 ~ 1
  )
##     Estimate  Est.Error        Q5       Q50       Q95
## R2 0.8936716 0.01985293 0.8587515 0.8998706 0.9073303
m_2 %>% brms::bayes_R2(
  probs = c(0.05, 0.5, 0.95), re_formula = 1 ~ 1
  )
##     Estimate   Est.Error        Q5       Q50       Q95
## R2 0.9040615 0.001862663 0.9008053 0.9041456 0.9069507
```

### Bayes factor (BARG 3.C)

``` r
bayes_factor <- brms::bayes_factor(m_1, m_2)
```

``` r
bayes_factor
```

    ## Estimated Bayes factor in favor of m_1 over m_2: 6860.54455

## Posterior distributions (BARG 3.B)

### Forest plot (BARG 3.B/5.B)

``` r
combined_models <- bind_rows(
  bayesplot::mcmc_intervals_data(
    posterior1, prob = 0.66, prob_outer = 0.95
    ) %>%
    mutate(model = "m_1"),
  bayesplot::mcmc_intervals_data(
    posterior2, prob = 0.66, prob_outer = 0.95
    ) %>%
    mutate(model = "m_2"),
  bayesplot::mcmc_intervals_data(
    posterior_flat, prob = 0.66, prob_outer = 0.95
    ) %>%
    mutate(model = "m_flat"),
  bayesplot::mcmc_intervals_data(
    posterior_prior, prob = 0.66, prob_outer = 0.95
    ) %>%
    mutate(model = "m_prior")
  )

pos <- position_nudge(
  y = if_else(
    combined_models$model == "m_2", -.2, if_else(
      combined_models$model == "m_flat", -.4, if_else(
        combined_models$model == "m_prior", -.6, 0
        )
      )
    )
  )

ggplot(
  data = combined_models,
  aes(x = m, y = forcats::fct_rev(factor(parameter)), color = model)
  ) +
  geom_vline(xintercept = 0, color = "grey") +
  geom_linerange(aes(xmin = l, xmax = h), position = pos, linewidth = 2) +
  geom_linerange(aes(xmin = ll, xmax = hh), position = pos) +
  geom_point(position = pos, color = "black") +
  coord_cartesian(xlim = c(-.9, .9)) +
  bayesplot::theme_default() +
  ggtitle("Posterior distbributions (mean, CI66, CI95)")
```

![](model_check_recovery_check_files/figure-gfm/posteriors-1.png)<!-- -->

### Effect sizes

Effect sizes of chosen model just to get exact values of means etc. if
necessary.

``` r
posterior1 %>%
  posterior::summarize_draws() %>%
  knitr::kable()
```

| variable                   |    mean |  median |      sd |     mad |       q5 |     q95 | rhat | ess_bulk | ess_tail |
|:---------------------------|--------:|--------:|--------:|--------:|---------:|--------:|-----:|---------:|---------:|
| Intercept                  |  -0.820 |  -0.820 |  0.0311 |  0.0240 |   -0.869 |  -0.768 | 1.00 |   4033\. |   1212\. |
| b_sand_ratio25             |   0.159 |   0.159 |  0.0395 |  0.0400 |   0.0943 |   0.224 | 1.00 |   2234\. |   4702\. |
| b_sand_ratio50             |   0.101 |   0.101 |  0.0396 |  0.0392 |   0.0366 |   0.167 | 1.00 |   2255\. |   4651\. |
| b_substrate_depth15        | -0.0445 | -0.0445 |  0.0181 |  0.0181 |  -0.0741 | -0.0148 | 1.00 |   5739\. |   7770\. |
| b_target_typedry_grassland |  -0.170 |  -0.170 |  0.0383 |  0.0383 |   -0.233 |  -0.107 | 1.00 |   1832\. |   3442\. |
| b_seed_density8            |  0.0162 |  0.0159 |  0.0130 |  0.0128 | -0.00499 |  0.0380 | 1.00 |   5072\. |   1454\. |
| b_expositionsouth          |  -0.395 |  -0.395 |  0.0407 |  0.0409 |   -0.462 |  -0.327 | 1.00 |   2054\. |   4641\. |
| b_survey_year_fct2019      |   0.774 |   0.776 |  0.0869 |  0.0677 |    0.638 |   0.907 | 1.00 |   2101\. |   1219\. |
| b_survey_year_fct2020      |   0.869 |   0.869 |  0.0727 |  0.0605 |    0.754 |   0.980 | 1.00 |   3011\. |   4251\. |
| b_survey_year_fct2021      |   0.866 |   0.866 |  0.0826 |  0.0687 |    0.736 |   0.995 | 1.00 |   3285\. |   4386\. |
| sd_site\_\_Intercept       |  0.0370 |  0.0320 |  0.0214 |  0.0132 |   0.0167 |  0.0722 | 1.00 |   4593\. |   5780\. |
| sd_site:plot\_\_Intercept  |  0.0447 |  0.0449 | 0.00636 | 0.00627 |   0.0340 |  0.0549 | 1.00 |   4435\. |   5609\. |
| sigma                      |   0.124 |   0.124 | 0.00303 | 0.00301 |    0.119 |   0.129 | 1.00 |   6100\. |   7907\. |

``` r
emm <- m_1 %>% emmeans(
  revpairwise ~ target_type + sand_ratio | exposition | survey_year_fct,
  type = "response"
  )
emm$emmeans
```

    ## exposition = north, survey_year_fct = 2018:
    ##  target_type   sand_ratio emmean lower.HPD upper.HPD
    ##  hay_meadow    0          -1.162    -1.269    -1.039
    ##  dry_grassland 0          -1.331    -1.452    -1.221
    ##  hay_meadow    25         -0.971    -1.080    -0.847
    ##  dry_grassland 25         -1.233    -1.354    -1.120
    ##  hay_meadow    50         -1.025    -1.138    -0.909
    ##  dry_grassland 50         -1.215    -1.326    -1.091
    ## 
    ## exposition = south, survey_year_fct = 2018:
    ##  target_type   sand_ratio emmean lower.HPD upper.HPD
    ##  hay_meadow    0          -1.545    -1.658    -1.429
    ##  dry_grassland 0          -1.602    -1.724    -1.488
    ##  hay_meadow    25         -1.548    -1.658    -1.426
    ##  dry_grassland 25         -1.550    -1.666    -1.430
    ##  hay_meadow    50         -1.625    -1.740    -1.507
    ##  dry_grassland 50         -1.654    -1.777    -1.544
    ## 
    ## exposition = north, survey_year_fct = 2019:
    ##  target_type   sand_ratio emmean lower.HPD upper.HPD
    ##  hay_meadow    0          -0.387    -0.543    -0.247
    ##  dry_grassland 0          -0.575    -0.725    -0.429
    ##  hay_meadow    25         -0.440    -0.589    -0.303
    ##  dry_grassland 25         -0.651    -0.815    -0.521
    ##  hay_meadow    50         -0.391    -0.555    -0.254
    ##  dry_grassland 50         -0.516    -0.668    -0.377
    ## 
    ## exposition = south, survey_year_fct = 2019:
    ##  target_type   sand_ratio emmean lower.HPD upper.HPD
    ##  hay_meadow    0          -1.030    -1.150    -0.912
    ##  dry_grassland 0          -1.160    -1.275    -1.044
    ##  hay_meadow    25         -0.993    -1.113    -0.882
    ##  dry_grassland 25         -1.050    -1.168    -0.934
    ##  hay_meadow    50         -1.046    -1.167    -0.933
    ##  dry_grassland 50         -1.054    -1.175    -0.941
    ## 
    ## exposition = north, survey_year_fct = 2020:
    ##  target_type   sand_ratio emmean lower.HPD upper.HPD
    ##  hay_meadow    0          -0.293    -0.401    -0.174
    ##  dry_grassland 0          -0.437    -0.549    -0.323
    ##  hay_meadow    25         -0.327    -0.436    -0.211
    ##  dry_grassland 25         -0.471    -0.578    -0.352
    ##  hay_meadow    50         -0.295    -0.412    -0.186
    ##  dry_grassland 50         -0.427    -0.538    -0.316
    ## 
    ## exposition = south, survey_year_fct = 2020:
    ##  target_type   sand_ratio emmean lower.HPD upper.HPD
    ##  hay_meadow    0          -0.798    -0.903    -0.671
    ##  dry_grassland 0          -0.885    -0.996    -0.768
    ##  hay_meadow    25         -0.745    -0.853    -0.630
    ##  dry_grassland 25         -0.791    -0.909    -0.678
    ##  hay_meadow    50         -0.799    -0.910    -0.683
    ##  dry_grassland 50         -0.797    -0.912    -0.686
    ## 
    ## exposition = north, survey_year_fct = 2021:
    ##  target_type   sand_ratio emmean lower.HPD upper.HPD
    ##  hay_meadow    0          -0.297    -0.444    -0.150
    ##  dry_grassland 0          -0.397    -0.552    -0.250
    ##  hay_meadow    25         -0.295    -0.443    -0.147
    ##  dry_grassland 25         -0.430    -0.579    -0.284
    ##  hay_meadow    50         -0.275    -0.416    -0.114
    ##  dry_grassland 50         -0.361    -0.506    -0.209
    ## 
    ## exposition = south, survey_year_fct = 2021:
    ##  target_type   sand_ratio emmean lower.HPD upper.HPD
    ##  hay_meadow    0          -0.680    -0.829    -0.530
    ##  dry_grassland 0          -0.839    -0.997    -0.696
    ##  hay_meadow    25         -0.718    -0.873    -0.573
    ##  dry_grassland 25         -0.768    -0.920    -0.618
    ##  hay_meadow    50         -0.703    -0.852    -0.549
    ##  dry_grassland 50         -0.760    -0.905    -0.610
    ## 
    ## Results are averaged over the levels of: substrate_depth, seed_density 
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
    ##  [1] emmeans_1.8.4-1           loo_2.5.1                
    ##  [3] bayesplot_1.10.0          DHARMa.helpers_0.0.0.9000
    ##  [5] brms_2.18.0               Rcpp_1.0.10              
    ##  [7] patchwork_1.1.2           ggbeeswarm_0.7.1         
    ##  [9] lubridate_1.9.2           forcats_1.0.0            
    ## [11] stringr_1.5.0             dplyr_1.1.0              
    ## [13] purrr_1.0.1               readr_2.1.4              
    ## [15] tidyr_1.3.0               tibble_3.1.8             
    ## [17] ggplot2_3.4.1             tidyverse_2.0.0          
    ## [19] here_1.0.1               
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] uuid_1.1-0           backports_1.4.1      plyr_1.8.8          
    ##   [4] igraph_1.4.1         sp_1.6-0             splines_4.2.2       
    ##   [7] crosstalk_1.2.0      gap.datasets_0.0.5   rncl_0.8.7          
    ##  [10] rstantools_2.2.0     inline_0.3.19        digest_0.6.31       
    ##  [13] foreach_1.5.2        htmltools_0.5.4      fansi_1.0.4         
    ##  [16] magrittr_2.0.3       checkmate_2.1.0      DHARMa_0.4.6        
    ##  [19] doParallel_1.0.17    cluster_2.1.4        tzdb_0.3.0          
    ##  [22] RcppParallel_5.1.7   matrixStats_0.63.0   vroom_1.6.1         
    ##  [25] xts_0.13.0           timechange_0.2.0     prettyunits_1.1.1   
    ##  [28] jpeg_0.1-10          colorspace_2.1-0     xfun_0.37           
    ##  [31] callr_3.7.3          crayon_1.5.2         jsonlite_1.8.4      
    ##  [34] lme4_1.1-31          phylobase_0.8.10     iterators_1.0.14    
    ##  [37] ape_5.7              zoo_1.8-11           glue_1.6.2          
    ##  [40] gtable_0.3.1         seqinr_4.2-23        V8_4.2.2            
    ##  [43] distributional_0.3.1 pkgbuild_1.4.0       rstan_2.26.13       
    ##  [46] adegraphics_1.0-17   abind_1.4-5          scales_1.2.1        
    ##  [49] mvtnorm_1.1-3        miniUI_0.1.1.1       progress_1.2.2      
    ##  [52] xtable_1.8-4         diffobj_0.3.5        bit_4.0.5           
    ##  [55] stats4_4.2.2         StanHeaders_2.26.13  DT_0.27             
    ##  [58] httr_1.4.5           htmlwidgets_1.6.1    threejs_0.3.3       
    ##  [61] RColorBrewer_1.1-3   posterior_1.4.0      ellipsis_0.3.2      
    ##  [64] XML_3.99-0.13        pkgconfig_2.0.3      farver_2.1.1        
    ##  [67] qgam_1.3.4           deldir_1.0-6         utf8_1.2.3          
    ##  [70] tidyselect_1.2.0     labeling_0.4.2       rlang_1.0.6         
    ##  [73] reshape2_1.4.4       later_1.3.0          munsell_0.5.0       
    ##  [76] tools_4.2.2          cli_3.6.0            generics_0.1.3      
    ##  [79] ade4_1.7-22          evaluate_0.20        fastmap_1.1.1       
    ##  [82] yaml_2.3.7           processx_3.8.0       knitr_1.42          
    ##  [85] bit64_4.0.5          nlme_3.1-162         mime_0.12           
    ##  [88] adegenet_2.1.10      xml2_1.3.3           gap_1.5-1           
    ##  [91] compiler_4.2.2       shinythemes_1.2.0    rstudioapi_0.14     
    ##  [94] beeswarm_0.4.0       curl_5.0.0           png_0.1-8           
    ##  [97] RNeXML_2.4.11        stringi_1.7.12       highr_0.10          
    ## [100] ps_1.7.2             Brobdingnag_1.2-9    lattice_0.20-45     
    ## [103] Matrix_1.5-3         nloptr_2.0.3         markdown_1.5        
    ## [106] permute_0.9-7        vegan_2.6-4          shinyjs_2.1.0       
    ## [109] tensorA_0.36.2       vctrs_0.5.2          pillar_1.8.1        
    ## [112] lifecycle_1.0.3      bridgesampling_1.1-2 estimability_1.4.1  
    ## [115] httpuv_1.6.9         R6_2.5.1             latticeExtra_0.6-30 
    ## [118] promises_1.2.0.1     KernSmooth_2.23-20   gridExtra_2.3       
    ## [121] vipor_0.4.5          codetools_0.2-18     boot_1.3-28         
    ## [124] colourpicker_1.2.0   MASS_7.3-58.2        gtools_3.9.4        
    ## [127] rprojroot_2.0.3      withr_2.5.0          shinystan_2.6.0     
    ## [130] mgcv_1.8-41          parallel_4.2.2       hms_1.1.2           
    ## [133] grid_4.2.2           coda_0.19-4          minqa_1.2.5         
    ## [136] rmarkdown_2.20       shiny_1.7.4          base64enc_0.1-3     
    ## [139] dygraphs_1.1.1.6     interp_1.1-3
