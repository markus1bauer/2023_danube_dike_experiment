Analysis of Bauer et al. (unpublished) Field experiment: <br> Recovery
completeness
================
<b>Markus Bauer</b> <br>
<b>2022-12-16</b>

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
    - <a href="#dharma" id="toc-dharma">DHARMa</a>
    - <a href="#sampling-efficency-and-effectiveness-barg-2bc"
      id="toc-sampling-efficency-and-effectiveness-barg-2bc">Sampling
      efficency and effectiveness (BARG 2.B/C)</a>
    - <a href="#mcmc-diagnostics" id="toc-mcmc-diagnostics">MCMC
      diagnostics</a>
    - <a href="#posterior-predictive-check-barg-3a"
      id="toc-posterior-predictive-check-barg-3a">Posterior predictive check
      (BARG 3.A)</a>
    - <a href="#autocorrelation-check"
      id="toc-autocorrelation-check">Autocorrelation check</a>
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
section ‘Load models’

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
library(tidybayes)
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
    n = recovery_time
    )
sites
```

    ## # A tibble: 1,152 × 24
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
    ## # … with 1,142 more rows, 14 more variables: survey_year <dbl>,
    ## #   longitude <dbl>, latitude <dbl>, elevation <dbl>, plot_size <dbl>,
    ## #   botanist <chr>, NMDS1 <dbl>, NMDS2 <dbl>, mean_reference <dbl>,
    ## #   sd_reference <dbl>, recovery_time <dbl>, survey_year_fct <fct>,
    ## #   botanist_year <fct>, n <dbl>, and abbreviated variable names ¹​reference,
    ## #   ²​exposition, ³​sand_ratio, ⁴​substrate_depth, ⁵​target_type, ⁶​seed_density

# Statistics

## Data exploration

### Graphs of raw data

![](model_recovery_check_files/figure-gfm/data-exploration-1.png)<!-- -->![](model_recovery_check_files/figure-gfm/data-exploration-2.png)<!-- -->![](model_recovery_check_files/figure-gfm/data-exploration-3.png)<!-- -->![](model_recovery_check_files/figure-gfm/data-exploration-4.png)<!-- -->![](model_recovery_check_files/figure-gfm/data-exploration-5.png)<!-- -->

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

![](model_recovery_check_files/figure-gfm/outliers-1.png)<!-- -->![](model_recovery_check_files/figure-gfm/outliers-2.png)<!-- -->![](model_recovery_check_files/figure-gfm/outliers-3.png)<!-- -->![](model_recovery_check_files/figure-gfm/outliers-4.png)<!-- -->

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

NUTS sampler is used.

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
# R hat
rhat1 <- rhat(m_1)
rhat2 <- rhat(m_2)
# NEFF ratio
neff1 <- neff_ratio(m_1)
neff2 <- neff_ratio(m_2)
# Long format of draws
hmc_diagnostics1 <- brms::nuts_params(m_1)
hmc_diagnostics2 <- brms::nuts_params(m_2)
y <- sites$n
# Posterior predictive distribution
yrep1 <- brms::posterior_predict(m_1, draws = 500)
yrep2 <- brms::posterior_predict(m_2, draws = 500)
yrep_prior <- brms::posterior_predict(m_prior, draws = 500)
# Leave-one-out cross validation based on posterior likelihood
loo1 <- brms::loo(m_1, save_psis = TRUE, moment_match = FALSE)
loo2 <- brms::loo(m_2, save_psis = TRUE, moment_match = FALSE)
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
get_prior(n ~ target_type + exposition + sand_ratio + survey_year_fct +
            seed_density + substrate_depth +
            (1 | site/plot) + (1|botanist_year),
          data = sites)
```

    ##                    prior     class                     coef         group resp
    ##                   (flat)         b                                            
    ##                   (flat)         b          expositionsouth                   
    ##                   (flat)         b             sand_ratio25                   
    ##                   (flat)         b             sand_ratio50                   
    ##                   (flat)         b            seed_density8                   
    ##                   (flat)         b        substrate_depth15                   
    ##                   (flat)         b      survey_year_fct2019                   
    ##                   (flat)         b      survey_year_fct2020                   
    ##                   (flat)         b      survey_year_fct2021                   
    ##                   (flat)         b target_typedry_grassland                   
    ##  student_t(3, -0.8, 2.5) Intercept                                            
    ##     student_t(3, 0, 2.5)        sd                                            
    ##     student_t(3, 0, 2.5)        sd                          botanist_year     
    ##     student_t(3, 0, 2.5)        sd                Intercept botanist_year     
    ##     student_t(3, 0, 2.5)        sd                                   site     
    ##     student_t(3, 0, 2.5)        sd                Intercept          site     
    ##     student_t(3, 0, 2.5)        sd                              site:plot     
    ##     student_t(3, 0, 2.5)        sd                Intercept     site:plot     
    ##     student_t(3, 0, 2.5)     sigma                                            
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
ggplot(data = data.frame(x = c(-2, 0)), aes(x = x)) +
  stat_function(fun = dnorm, n = 101, args = list(mean = -1, sd = 1)) +
  expand_limits(y = 0) +
  ggtitle("Normal distribution for Intecept")
```

![](model_recovery_check_files/figure-gfm/possible-priors-1.png)<!-- -->

``` r
data <- data.frame(x = c(-2, 2))
ggplot(data, aes(x = x)) +
  stat_function(fun = dnorm, n = 101, args = list(mean = 0, sd = .8)) +
  expand_limits(y = 0) +
  ggtitle("Normal distribution for treatments")
```

![](model_recovery_check_files/figure-gfm/possible-priors-2.png)<!-- -->

``` r
ggplot(data, aes(x = x)) +
  stat_function(fun = dcauchy, n = 101, args = list(location = 0, scale = 1)) +
  expand_limits(y = 0) +
  ggtitle("Cauchy distribution") # See Lemoine 2019 https://doi.org/10.1111/oik.05985
```

![](model_recovery_check_files/figure-gfm/possible-priors-3.png)<!-- -->

``` r
ggplot(data, aes(x = x)) +
  stat_function(fun = dstudent_t, args = list(df = 3, mu = 0, sigma = 2.5)) +
  expand_limits(y = 0) +
  ggtitle(expression(Student~italic(t)*"-distribution")) # Software standard
```

![](model_recovery_check_files/figure-gfm/possible-priors-4.png)<!-- -->

#### Prior summary (BARG 1.D)

``` r
brms::prior_summary(m_1, all = FALSE)
```

    ##                    prior     class                coef group resp dpar nlpar lb
    ##             normal(0, 2)         b                                             
    ##          normal(-0.1, 2)         b     expositionsouth                         
    ##           normal(0.1, 2)         b        sand_ratio25                         
    ##           normal(0.2, 2)         b        sand_ratio50                         
    ##           normal(0.1, 2)         b survey_year_fct2019                         
    ##           normal(0.2, 2)         b survey_year_fct2020                         
    ##           normal(0.3, 2)         b survey_year_fct2021                         
    ##  student_t(3, -0.8, 2.5) Intercept                                             
    ##     student_t(3, 0, 2.5)        sd                                            0
    ##             cauchy(0, 1)     sigma                                            0
    ##  ub  source
    ##        user
    ##        user
    ##        user
    ##        user
    ##        user
    ##        user
    ##        user
    ##     default
    ##     default
    ##        user

#### Prior predictive check (BARG 1.E)

``` r
ppd_stat(yrep_prior[1:500, ], binwidth = 0.2) +
  coord_cartesian(xlim = c(-2, 2))
```

![](model_recovery_check_files/figure-gfm/prior-predictive-check-1.png)<!-- -->

``` r
ppd_stat_grouped(yrep_prior[1:500, ], group = sites$site,
                 binwidth = 0.2) +
  coord_cartesian(xlim = c(-2, 2))
```

![](model_recovery_check_files/figure-gfm/prior-predictive-check-2.png)<!-- -->

``` r
ppd_stat_grouped(yrep_prior[1:500, ], group = sites$exposition,
                 binwidth = 0.2) +
  coord_cartesian(xlim = c(-2, 2))
```

![](model_recovery_check_files/figure-gfm/prior-predictive-check-3.png)<!-- -->

``` r
ppd_stat_grouped(yrep_prior[1:500, ], group = sites$survey_year_fct,
                 binwidth = 0.2) +
  coord_cartesian(xlim = c(-2, 2))
```

![](model_recovery_check_files/figure-gfm/prior-predictive-check-4.png)<!-- -->

``` r
ppd_stat_grouped(yrep_prior[1:500, ], group = sites$target_type,
                 binwidth = 0.2) +
  coord_cartesian(xlim = c(-2, 2))
```

![](model_recovery_check_files/figure-gfm/prior-predictive-check-5.png)<!-- -->

``` r
ppd_stat_grouped(yrep_prior[1:500, ], group = sites$seed_density,
                 binwidth = 0.1) +
  coord_cartesian(xlim = c(-2, 2))
```

![](model_recovery_check_files/figure-gfm/prior-predictive-check-6.png)<!-- -->

``` r
ppd_stat_grouped(yrep_prior[1:500, ], group = sites$sand_ratio,
                 binwidth = 0.2) +
  coord_cartesian(xlim = c(-2, 2))
```

![](model_recovery_check_files/figure-gfm/prior-predictive-check-7.png)<!-- -->

``` r
ppd_stat_grouped(yrep_prior[1:500, ], group = sites$substrate_depth,
                 binwidth = 0.2) +
  coord_cartesian(xlim = c(-2, 2))
```

![](model_recovery_check_files/figure-gfm/prior-predictive-check-8.png)<!-- -->

## Model check

### DHARMa

``` r
DHARMa.helpers::dh_check_brms(m_1, integer = TRUE)
```

![](model_recovery_check_files/figure-gfm/dharma-1.png)<!-- -->

``` r
DHARMa.helpers::dh_check_brms(m_2, integer = TRUE)
```

![](model_recovery_check_files/figure-gfm/dharma-2.png)<!-- -->

### Sampling efficency and effectiveness (BARG 2.B/C)

#### Rhat (BARG 2.B)

``` r
mcmc_rhat(rhat1)
```

![](model_recovery_check_files/figure-gfm/rhat-1.png)<!-- -->

``` r
mcmc_rhat(rhat2)
```

![](model_recovery_check_files/figure-gfm/rhat-2.png)<!-- -->

#### Effective sampling size (ESS) (BARG 2.C)

``` r
mcmc_neff(neff1)
```

![](model_recovery_check_files/figure-gfm/ess-1.png)<!-- -->

``` r
mcmc_neff(neff2)
```

![](model_recovery_check_files/figure-gfm/ess-2.png)<!-- -->

### MCMC diagnostics

#### Trace plots

``` r
mcmc_trace(posterior1, np = hmc_diagnostics1, facet_args = list(ncol = 2))
```

![](model_recovery_check_files/figure-gfm/mcmc-trace-1.png)<!-- -->

``` r
mcmc_trace(posterior2, np = hmc_diagnostics2, facet_args = list(ncol = 2))
```

    ## No divergences to plot.

![](model_recovery_check_files/figure-gfm/mcmc-trace-2.png)<!-- -->

#### Pairs plot

``` r
mcmc_pairs(m_1, off_diag_args = list(size = 1.2),
           pars = c(
             "b_sand_ratio25", "b_sand_ratio50", "b_substrate_depth15",
             "b_target_typedry_grassland", "b_seed_density8",
             "b_expositionsouth", "sigma"
           ))
```

![](model_recovery_check_files/figure-gfm/mcmc-pairs-1.png)<!-- -->

``` r
mcmc_pairs(m_2, off_diag_args = list(size = 1.2),
           pars = c(
             "b_sand_ratio25", "b_sand_ratio50", "b_substrate_depth15",
             "b_target_typedry_grassland", "b_seed_density8",
             "b_expositionsouth", "sigma"
           ))
```

![](model_recovery_check_files/figure-gfm/mcmc-pairs-2.png)<!-- -->

#### Parallel coordinate plot

``` r
mcmc_parcoord(posterior1, np = hmc_diagnostics1) +
  theme(axis.text.x = element_text(angle = 45))
```

![](model_recovery_check_files/figure-gfm/mcmc-parcoord-1.png)<!-- -->

``` r
mcmc_parcoord(posterior2, np = hmc_diagnostics2) +
  theme(axis.text.x = element_text(angle = 45))
```

![](model_recovery_check_files/figure-gfm/mcmc-parcoord-2.png)<!-- -->

### Posterior predictive check (BARG 3.A)

#### Kernel density

``` r
p1 <- ppc_dens_overlay(y, yrep1[1:50, ])
p2 <- ppc_dens_overlay(y, yrep2[1:50, ])
p1 / p2
```

![](model_recovery_check_files/figure-gfm/kernel-density-1.png)<!-- -->

``` r
ppc_dens_overlay_grouped(y, yrep1[1:50, ], group = sites$site)
```

![](model_recovery_check_files/figure-gfm/kernel-density-2.png)<!-- -->

``` r
ppc_dens_overlay_grouped(y, yrep2[1:50, ], group = sites$site)
```

![](model_recovery_check_files/figure-gfm/kernel-density-3.png)<!-- -->

``` r
p1 <- ppc_dens_overlay_grouped(y, yrep1[1:50, ], group = sites$exposition)
p2 <- ppc_dens_overlay_grouped(y, yrep2[1:50, ], group = sites$exposition)
p1 / p2
```

![](model_recovery_check_files/figure-gfm/kernel-density-4.png)<!-- -->

``` r
ppc_dens_overlay_grouped(y, yrep1[1:50, ], group = sites$survey_year_fct)
```

![](model_recovery_check_files/figure-gfm/kernel-density-5.png)<!-- -->

``` r
ppc_dens_overlay_grouped(y, yrep2[1:50, ], group = sites$survey_year_fct)
```

![](model_recovery_check_files/figure-gfm/kernel-density-6.png)<!-- -->

``` r
p1 <- ppc_dens_overlay_grouped(y, yrep1[1:50, ], group = sites$target_type)
p2 <- ppc_dens_overlay_grouped(y, yrep2[1:50, ], group = sites$target_type)
p1 / p2
```

![](model_recovery_check_files/figure-gfm/kernel-density-7.png)<!-- -->

``` r
p1 <- ppc_dens_overlay_grouped(y, yrep1[1:50, ], group = sites$seed_density)
p2 <- ppc_dens_overlay_grouped(y, yrep2[1:50, ], group = sites$seed_density)
p1 / p2
```

![](model_recovery_check_files/figure-gfm/kernel-density-8.png)<!-- -->

``` r
p1 <- ppc_dens_overlay_grouped(y, yrep1[1:50, ], group = sites$sand_ratio)
p2 <- ppc_dens_overlay_grouped(y, yrep2[1:50, ], group = sites$sand_ratio)
p1 / p2
```

![](model_recovery_check_files/figure-gfm/kernel-density-9.png)<!-- -->

``` r
p1 <- ppc_dens_overlay_grouped(y, yrep1[1:50, ], group = sites$substrate_depth)
p2 <- ppc_dens_overlay_grouped(y, yrep2[1:50, ], group = sites$substrate_depth)
p1 / p2
```

![](model_recovery_check_files/figure-gfm/kernel-density-10.png)<!-- -->

#### Histograms of statistics skew

``` r
p1 <- ppc_stat(y, yrep1, binwidth = 0.001)
p2 <- ppc_stat(y, yrep2, binwidth = 0.001)
p1 / p2
```

![](model_recovery_check_files/figure-gfm/histograms-1.png)<!-- -->

``` r
ppc_stat_grouped(y, yrep1, group = sites$site, binwidth = 0.001)
```

![](model_recovery_check_files/figure-gfm/histograms-2.png)<!-- -->

``` r
ppc_stat_grouped(y, yrep2, group = sites$site, binwidth = 0.001)
```

![](model_recovery_check_files/figure-gfm/histograms-3.png)<!-- -->

``` r
p1 <- ppc_stat_grouped(y, yrep1, group = sites$exposition, binwidth = 0.001)
p2 <- ppc_stat_grouped(y, yrep2, group = sites$exposition, binwidth = 0.001)
p1 / p2
```

![](model_recovery_check_files/figure-gfm/histograms-4.png)<!-- -->

``` r
ppc_stat_grouped(y, yrep1, group = sites$survey_year_fct, binwidth = 0.001)
```

![](model_recovery_check_files/figure-gfm/histograms-5.png)<!-- -->

``` r
ppc_stat_grouped(y, yrep2, group = sites$survey_year_fct, binwidth = 0.001)
```

![](model_recovery_check_files/figure-gfm/histograms-6.png)<!-- -->

``` r
p1 <- ppc_stat_grouped(y, yrep1, group = sites$target_type, binwidth = 0.001)
p2 <- ppc_stat_grouped(y, yrep2, group = sites$target_type, binwidth = 0.001)
p1 / p2
```

![](model_recovery_check_files/figure-gfm/histograms-7.png)<!-- -->

``` r
p1 <- ppc_stat_grouped(y, yrep1, group = sites$seed_density, binwidth = 0.001)
p2 <- ppc_stat_grouped(y, yrep2, group = sites$seed_density, binwidth = 0.001)
p1 / p2
```

![](model_recovery_check_files/figure-gfm/histograms-8.png)<!-- -->

``` r
p1 <- ppc_stat_grouped(y, yrep1, group = sites$sand_ratio, binwidth = 0.001)
p2 <- ppc_stat_grouped(y, yrep2, group = sites$sand_ratio, binwidth = 0.001)
p1 / p2
```

![](model_recovery_check_files/figure-gfm/histograms-9.png)<!-- -->

``` r
p1 <- ppc_stat_grouped(y, yrep1, group = sites$substrate_depth, binwidth = 0.001)
p2 <- ppc_stat_grouped(y, yrep2, group = sites$substrate_depth, binwidth = 0.001)
p1 / p2
```

![](model_recovery_check_files/figure-gfm/histograms-10.png)<!-- -->

#### LOO (Leave one out)

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

![](model_recovery_check_files/figure-gfm/loo-1.png)<!-- -->

``` r
plot(loo2)
```

![](model_recovery_check_files/figure-gfm/loo-2.png)<!-- -->

Leave one out probability integral transform

``` r
p1 <- bayesplot::ppc_loo_pit_overlay(y, yrep1, lw = weights(loo1$psis_object))
```

    ## NOTE: The kernel density estimate assumes continuous observations and is not optimal for discrete observations.

``` r
p2 <- bayesplot::ppc_loo_pit_overlay(y, yrep2, lw = weights(loo2$psis_object))
```

    ## NOTE: The kernel density estimate assumes continuous observations and is not optimal for discrete observations.

``` r
p1 / p2
```

![](model_recovery_check_files/figure-gfm/loo-pit-1.png)<!-- -->

### Autocorrelation check

``` r
mcmc_acf(posterior1, lags = 10)
```

![](model_recovery_check_files/figure-gfm/autocorrelation-1.png)<!-- -->

``` r
mcmc_acf(posterior2, lags = 10)
```

![](model_recovery_check_files/figure-gfm/autocorrelation-2.png)<!-- -->

## Model comparison

### Conditional <i>R</i><sup>2</sup> values

``` r
brms::bayes_R2(m_1, probs = c(0.05, 0.5, 0.95),
         re_formula =  ~ (1 | site/plot)) 
##     Estimate  Est.Error       Q5       Q50       Q95
## R2 0.9068055 0.01980042 0.871471 0.9128796 0.9205551
brms::bayes_R2(m_2, probs = c(0.05, 0.5, 0.95),
         re_formula =  ~ (1 | site/plot))
##    Estimate   Est.Error        Q5       Q50       Q95
## R2 0.918869 0.002526252 0.9145194 0.9189666 0.9227846
```

### Marginal <i>R</i><sup>2</sup> values

``` r
brms::bayes_R2(m_1, probs = c(0.05, 0.5, 0.95),
         re_formula = 1 ~ 1)
##     Estimate  Est.Error        Q5       Q50       Q95
## R2 0.8936716 0.01985293 0.8587515 0.8998706 0.9073303
brms::bayes_R2(m_2, probs = c(0.05, 0.5, 0.95),
         re_formula = 1 ~ 1)
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

    ## Estimated Bayes factor in favor of m_1 over m_2: 8678.32629

## Posterior distributions (BARG 3.B)

### Forest plot (BARG 3.B/5.B)

``` r
combined <- bind_rows(
  bayesplot::mcmc_intervals_data(posterior1, prob = 0.66, prob_outer = 0.95) %>%
    mutate(model = "m_1"),
  bayesplot::mcmc_intervals_data(posterior2, prob = 0.66, prob_outer = 0.95) %>%
    mutate(model = "m_2"),
  bayesplot::mcmc_intervals_data(posterior_flat, prob = 0.66, prob_outer = 0.95) %>%
    mutate(model = "m_flat"),
  bayesplot::mcmc_intervals_data(posterior_prior, prob = 0.66, prob_outer = 0.95) %>%
    mutate(model = "m_prior")
  )

pos <- position_nudge(
  y = if_else(
    combined$model == "m_2", -.2, if_else(
      combined$model == "m_flat", -.4, if_else(
        combined$model == "m_prior", -.6, 0
        )
      )
    )
  )

ggplot(data = combined, aes(x = m, y = forcats::fct_rev(factor(parameter)), color = model)) + 
  geom_vline(xintercept = 0, color = "grey") +
  geom_linerange(aes(xmin = l, xmax = h), position = pos, linewidth = 2) +
  geom_linerange(aes(xmin = ll, xmax = hh), position = pos) +
  geom_point(position = pos, color = "black") +
  coord_cartesian(xlim = c(-1, .6)) +
  bayesplot::theme_default() +
  ggtitle("Posterior distbributions (mean, CI66, CI95)")
```

![](model_recovery_check_files/figure-gfm/posteriors-1.png)<!-- -->

### Effect sizes

Effect sizes of chosen model just to get exact values of means etc. if
necessary.

``` r
posterior1 %>%
  posterior::summarize_draws() %>%
  knitr::kable()
```

| variable                   |       mean |     median |        sd |       mad |         q5 |        q95 |     rhat | ess_bulk | ess_tail |
|:---------------------------|-----------:|-----------:|----------:|----------:|-----------:|-----------:|---------:|---------:|---------:|
| Intercept                  | -0.8196026 | -0.8196327 | 0.0311115 | 0.0240296 | -0.8690680 | -0.7684832 | 1.001980 | 4033.348 | 1212.149 |
| b_sand_ratio25             |  0.1587360 |  0.1587951 | 0.0395010 | 0.0399671 |  0.0942963 |  0.2236053 | 1.001866 | 2234.361 | 4702.160 |
| b_sand_ratio50             |  0.1014113 |  0.1014346 | 0.0395569 | 0.0391919 |  0.0366146 |  0.1665276 | 1.001278 | 2255.234 | 4651.450 |
| b_substrate_depth15        | -0.0444619 | -0.0445411 | 0.0181157 | 0.0181248 | -0.0741270 | -0.0148276 | 1.000122 | 5739.387 | 7770.151 |
| b_target_typedry_grassland | -0.1698148 | -0.1698754 | 0.0382637 | 0.0382804 | -0.2327421 | -0.1066600 | 1.001402 | 1832.157 | 3441.552 |
| b_seed_density8            |  0.0162078 |  0.0159075 | 0.0129580 | 0.0127678 | -0.0049885 |  0.0380383 | 1.001446 | 5072.016 | 1454.450 |
| b_expositionsouth          | -0.3951320 | -0.3951642 | 0.0407208 | 0.0409170 | -0.4617760 | -0.3273771 | 1.001064 | 2054.457 | 4641.283 |
| b_survey_year_fct2019      |  0.7736543 |  0.7757078 | 0.0869140 | 0.0677217 |  0.6377376 |  0.9070336 | 1.001554 | 2101.248 | 1219.468 |
| b_survey_year_fct2020      |  0.8687929 |  0.8692599 | 0.0726742 | 0.0604726 |  0.7540026 |  0.9803451 | 1.000148 | 3011.371 | 4250.676 |
| b_survey_year_fct2021      |  0.8658432 |  0.8660651 | 0.0825965 | 0.0687188 |  0.7355024 |  0.9949243 | 1.000802 | 3285.086 | 4386.183 |
| sd_site\_\_Intercept       |  0.0369950 |  0.0320495 | 0.0213797 | 0.0131801 |  0.0167116 |  0.0722090 | 1.000214 | 4592.954 | 5780.205 |
| sd_site:plot\_\_Intercept  |  0.0447483 |  0.0448517 | 0.0063593 | 0.0062661 |  0.0340245 |  0.0548786 | 1.000165 | 4435.302 | 5608.777 |
| sigma                      |  0.1238085 |  0.1237160 | 0.0030336 | 0.0030058 |  0.1189475 |  0.1289529 | 1.000460 | 6099.910 | 7907.036 |

``` r
emmeans(m_1, revpairwise ~ target_type + sand_ratio | exposition |
          survey_year_fct, type = "response")
```

    ## $emmeans
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
    ## 
    ## $contrasts
    ## exposition = north, survey_year_fct = 2018:
    ##  contrast                                                 estimate lower.HPD
    ##  dry_grassland sand_ratio0 - hay_meadow sand_ratio0      -0.169875  -0.24443
    ##  hay_meadow sand_ratio25 - hay_meadow sand_ratio0         0.190674   0.11703
    ##  hay_meadow sand_ratio25 - dry_grassland sand_ratio0      0.360227   0.28729
    ##  dry_grassland sand_ratio25 - hay_meadow sand_ratio0     -0.070645  -0.14617
    ##  dry_grassland sand_ratio25 - dry_grassland sand_ratio0   0.099216   0.02942
    ##  dry_grassland sand_ratio25 - hay_meadow sand_ratio25    -0.260961  -0.33639
    ##  hay_meadow sand_ratio50 - hay_meadow sand_ratio0         0.137755   0.06256
    ##  hay_meadow sand_ratio50 - dry_grassland sand_ratio0      0.307784   0.23342
    ##  hay_meadow sand_ratio50 - hay_meadow sand_ratio25       -0.053515  -0.12722
    ##  hay_meadow sand_ratio50 - dry_grassland sand_ratio25     0.207421   0.13175
    ##  dry_grassland sand_ratio50 - hay_meadow sand_ratio0     -0.052177  -0.12312
    ##  dry_grassland sand_ratio50 - dry_grassland sand_ratio0   0.117938   0.04228
    ##  dry_grassland sand_ratio50 - hay_meadow sand_ratio25    -0.243053  -0.31939
    ##  dry_grassland sand_ratio50 - dry_grassland sand_ratio25  0.017632  -0.05840
    ##  dry_grassland sand_ratio50 - hay_meadow sand_ratio50    -0.189666  -0.26331
    ##  upper.HPD
    ##  -9.41e-02
    ##   2.65e-01
    ##   4.34e-01
    ##   3.76e-03
    ##   1.77e-01
    ##  -1.89e-01
    ##   2.10e-01
    ##   3.80e-01
    ##   2.13e-02
    ##   2.81e-01
    ##   2.26e-02
    ##   1.90e-01
    ##  -1.72e-01
    ##   9.08e-02
    ##  -1.17e-01
    ## 
    ## exposition = south, survey_year_fct = 2018:
    ##  contrast                                                 estimate lower.HPD
    ##  dry_grassland sand_ratio0 - hay_meadow sand_ratio0      -0.057559  -0.13053
    ##  hay_meadow sand_ratio25 - hay_meadow sand_ratio0        -0.002871  -0.07531
    ##  hay_meadow sand_ratio25 - dry_grassland sand_ratio0      0.053575  -0.01896
    ##  dry_grassland sand_ratio25 - hay_meadow sand_ratio0     -0.005285  -0.07995
    ##  dry_grassland sand_ratio25 - dry_grassland sand_ratio0   0.051624  -0.02307
    ##  dry_grassland sand_ratio25 - hay_meadow sand_ratio25    -0.002800  -0.07523
    ##  hay_meadow sand_ratio50 - hay_meadow sand_ratio0        -0.079314  -0.15453
    ##  hay_meadow sand_ratio50 - dry_grassland sand_ratio0     -0.022389  -0.09880
    ##  hay_meadow sand_ratio50 - hay_meadow sand_ratio25       -0.076849  -0.15038
    ##  hay_meadow sand_ratio50 - dry_grassland sand_ratio25    -0.074946  -0.14948
    ##  dry_grassland sand_ratio50 - hay_meadow sand_ratio0     -0.108714  -0.18194
    ##  dry_grassland sand_ratio50 - dry_grassland sand_ratio0  -0.052308  -0.12135
    ##  dry_grassland sand_ratio50 - hay_meadow sand_ratio25    -0.106064  -0.17970
    ##  dry_grassland sand_ratio50 - dry_grassland sand_ratio25 -0.103592  -0.18019
    ##  dry_grassland sand_ratio50 - hay_meadow sand_ratio50    -0.030196  -0.10688
    ##  upper.HPD
    ##   1.55e-02
    ##   7.31e-02
    ##   1.31e-01
    ##   6.90e-02
    ##   1.24e-01
    ##   7.28e-02
    ##  -4.67e-03
    ##   4.98e-02
    ##  -1.68e-03
    ##  -6.08e-05
    ##  -3.39e-02
    ##   2.65e-02
    ##  -3.13e-02
    ##  -3.12e-02
    ##   4.31e-02
    ## 
    ## exposition = north, survey_year_fct = 2019:
    ##  contrast                                                 estimate lower.HPD
    ##  dry_grassland sand_ratio0 - hay_meadow sand_ratio0      -0.188396  -0.26377
    ##  hay_meadow sand_ratio25 - hay_meadow sand_ratio0        -0.051501  -0.12241
    ##  hay_meadow sand_ratio25 - dry_grassland sand_ratio0      0.136564   0.06075
    ##  dry_grassland sand_ratio25 - hay_meadow sand_ratio0     -0.264180  -0.33667
    ##  dry_grassland sand_ratio25 - dry_grassland sand_ratio0  -0.075420  -0.15155
    ##  dry_grassland sand_ratio25 - hay_meadow sand_ratio25    -0.212557  -0.28627
    ##  hay_meadow sand_ratio50 - hay_meadow sand_ratio0        -0.005320  -0.08041
    ##  hay_meadow sand_ratio50 - dry_grassland sand_ratio0      0.183276   0.11254
    ##  hay_meadow sand_ratio50 - hay_meadow sand_ratio25        0.046646  -0.02621
    ##  hay_meadow sand_ratio50 - dry_grassland sand_ratio25     0.259317   0.18290
    ##  dry_grassland sand_ratio50 - hay_meadow sand_ratio0     -0.128753  -0.20134
    ##  dry_grassland sand_ratio50 - dry_grassland sand_ratio0   0.059631  -0.01373
    ##  dry_grassland sand_ratio50 - hay_meadow sand_ratio25    -0.077446  -0.14905
    ##  dry_grassland sand_ratio50 - dry_grassland sand_ratio25  0.135518   0.05923
    ##  dry_grassland sand_ratio50 - hay_meadow sand_ratio50    -0.123776  -0.20052
    ##  upper.HPD
    ##  -1.13e-01
    ##   2.30e-02
    ##   2.11e-01
    ##  -1.90e-01
    ##  -2.49e-03
    ##  -1.38e-01
    ##   6.86e-02
    ##   2.61e-01
    ##   1.23e-01
    ##   3.35e-01
    ##  -5.39e-02
    ##   1.36e-01
    ##  -3.10e-04
    ##   2.06e-01
    ##  -5.05e-02
    ## 
    ## exposition = south, survey_year_fct = 2019:
    ##  contrast                                                 estimate lower.HPD
    ##  dry_grassland sand_ratio0 - hay_meadow sand_ratio0      -0.130682  -0.20308
    ##  hay_meadow sand_ratio25 - hay_meadow sand_ratio0         0.036166  -0.04031
    ##  hay_meadow sand_ratio25 - dry_grassland sand_ratio0      0.166842   0.08933
    ##  dry_grassland sand_ratio25 - hay_meadow sand_ratio0     -0.020493  -0.09696
    ##  dry_grassland sand_ratio25 - dry_grassland sand_ratio0   0.109772   0.03252
    ##  dry_grassland sand_ratio25 - hay_meadow sand_ratio25    -0.057020  -0.13259
    ##  hay_meadow sand_ratio50 - hay_meadow sand_ratio0        -0.017084  -0.09019
    ##  hay_meadow sand_ratio50 - dry_grassland sand_ratio0      0.113885   0.03919
    ##  hay_meadow sand_ratio50 - hay_meadow sand_ratio25       -0.053135  -0.12495
    ##  hay_meadow sand_ratio50 - dry_grassland sand_ratio25     0.003988  -0.07188
    ##  dry_grassland sand_ratio50 - hay_meadow sand_ratio0     -0.024721  -0.09963
    ##  dry_grassland sand_ratio50 - dry_grassland sand_ratio0   0.105691   0.03170
    ##  dry_grassland sand_ratio50 - hay_meadow sand_ratio25    -0.061223  -0.13475
    ##  dry_grassland sand_ratio50 - dry_grassland sand_ratio25 -0.003438  -0.07861
    ##  dry_grassland sand_ratio50 - hay_meadow sand_ratio50    -0.007596  -0.08049
    ##  upper.HPD
    ##  -5.44e-02
    ##   1.07e-01
    ##   2.39e-01
    ##   5.31e-02
    ##   1.81e-01
    ##   1.41e-02
    ##   6.16e-02
    ##   1.89e-01
    ##   2.10e-02
    ##   7.63e-02
    ##   4.78e-02
    ##   1.83e-01
    ##   1.31e-02
    ##   6.99e-02
    ##   6.70e-02
    ## 
    ## exposition = north, survey_year_fct = 2020:
    ##  contrast                                                 estimate lower.HPD
    ##  dry_grassland sand_ratio0 - hay_meadow sand_ratio0      -0.144660  -0.21796
    ##  hay_meadow sand_ratio25 - hay_meadow sand_ratio0        -0.033937  -0.10866
    ##  hay_meadow sand_ratio25 - dry_grassland sand_ratio0      0.110991   0.03438
    ##  dry_grassland sand_ratio25 - hay_meadow sand_ratio0     -0.178207  -0.25095
    ##  dry_grassland sand_ratio25 - dry_grassland sand_ratio0  -0.033850  -0.11288
    ##  dry_grassland sand_ratio25 - hay_meadow sand_ratio25    -0.144366  -0.22141
    ##  hay_meadow sand_ratio50 - hay_meadow sand_ratio0        -0.001143  -0.07314
    ##  hay_meadow sand_ratio50 - dry_grassland sand_ratio0      0.142953   0.06981
    ##  hay_meadow sand_ratio50 - hay_meadow sand_ratio25        0.031601  -0.04061
    ##  hay_meadow sand_ratio50 - dry_grassland sand_ratio25     0.176617   0.10021
    ##  dry_grassland sand_ratio50 - hay_meadow sand_ratio0     -0.132782  -0.20676
    ##  dry_grassland sand_ratio50 - dry_grassland sand_ratio0   0.010791  -0.06320
    ##  dry_grassland sand_ratio50 - hay_meadow sand_ratio25    -0.099928  -0.17467
    ##  dry_grassland sand_ratio50 - dry_grassland sand_ratio25  0.044994  -0.03200
    ##  dry_grassland sand_ratio50 - hay_meadow sand_ratio50    -0.132136  -0.20551
    ##  upper.HPD
    ##  -6.94e-02
    ##   3.85e-02
    ##   1.85e-01
    ##  -1.03e-01
    ##   3.70e-02
    ##  -7.43e-02
    ##   7.59e-02
    ##   2.20e-01
    ##   1.07e-01
    ##   2.50e-01
    ##  -5.79e-02
    ##   8.59e-02
    ##  -2.76e-02
    ##   1.17e-01
    ##  -5.59e-02
    ## 
    ## exposition = south, survey_year_fct = 2020:
    ##  contrast                                                 estimate lower.HPD
    ##  dry_grassland sand_ratio0 - hay_meadow sand_ratio0      -0.087424  -0.16521
    ##  hay_meadow sand_ratio25 - hay_meadow sand_ratio0         0.052579  -0.02227
    ##  hay_meadow sand_ratio25 - dry_grassland sand_ratio0      0.140555   0.06556
    ##  dry_grassland sand_ratio25 - hay_meadow sand_ratio0      0.007091  -0.06412
    ##  dry_grassland sand_ratio25 - dry_grassland sand_ratio0   0.095172   0.02100
    ##  dry_grassland sand_ratio25 - hay_meadow sand_ratio25    -0.045689  -0.11784
    ##  hay_meadow sand_ratio50 - hay_meadow sand_ratio0        -0.001632  -0.07445
    ##  hay_meadow sand_ratio50 - dry_grassland sand_ratio0      0.085681   0.01308
    ##  hay_meadow sand_ratio50 - hay_meadow sand_ratio25       -0.054081  -0.12807
    ##  hay_meadow sand_ratio50 - dry_grassland sand_ratio25    -0.008493  -0.08268
    ##  dry_grassland sand_ratio50 - hay_meadow sand_ratio0      0.000601  -0.07313
    ##  dry_grassland sand_ratio50 - dry_grassland sand_ratio0   0.088477   0.01207
    ##  dry_grassland sand_ratio50 - hay_meadow sand_ratio25    -0.051974  -0.12444
    ##  dry_grassland sand_ratio50 - dry_grassland sand_ratio25 -0.006867  -0.07942
    ##  dry_grassland sand_ratio50 - hay_meadow sand_ratio50     0.001643  -0.07391
    ##  upper.HPD
    ##  -1.47e-02
    ##   1.24e-01
    ##   2.15e-01
    ##   8.45e-02
    ##   1.69e-01
    ##   2.91e-02
    ##   7.36e-02
    ##   1.63e-01
    ##   2.05e-02
    ##   6.43e-02
    ##   7.57e-02
    ##   1.61e-01
    ##   2.24e-02
    ##   6.88e-02
    ##   7.51e-02
    ## 
    ## exposition = north, survey_year_fct = 2021:
    ##  contrast                                                 estimate lower.HPD
    ##  dry_grassland sand_ratio0 - hay_meadow sand_ratio0      -0.100414  -0.17231
    ##  hay_meadow sand_ratio25 - hay_meadow sand_ratio0         0.001221  -0.07300
    ##  hay_meadow sand_ratio25 - dry_grassland sand_ratio0      0.102335   0.02827
    ##  dry_grassland sand_ratio25 - hay_meadow sand_ratio0     -0.134329  -0.20611
    ##  dry_grassland sand_ratio25 - dry_grassland sand_ratio0  -0.033611  -0.10512
    ##  dry_grassland sand_ratio25 - hay_meadow sand_ratio25    -0.135265  -0.20920
    ##  hay_meadow sand_ratio50 - hay_meadow sand_ratio0         0.020316  -0.05497
    ##  hay_meadow sand_ratio50 - dry_grassland sand_ratio0      0.121163   0.04462
    ##  hay_meadow sand_ratio50 - hay_meadow sand_ratio25        0.020189  -0.05506
    ##  hay_meadow sand_ratio50 - dry_grassland sand_ratio25     0.154889   0.07964
    ##  dry_grassland sand_ratio50 - hay_meadow sand_ratio0     -0.064182  -0.13431
    ##  dry_grassland sand_ratio50 - dry_grassland sand_ratio0   0.036396  -0.03876
    ##  dry_grassland sand_ratio50 - hay_meadow sand_ratio25    -0.064894  -0.14163
    ##  dry_grassland sand_ratio50 - dry_grassland sand_ratio25  0.069909  -0.00715
    ##  dry_grassland sand_ratio50 - hay_meadow sand_ratio50    -0.085441  -0.15882
    ##  upper.HPD
    ##  -2.35e-02
    ##   7.57e-02
    ##   1.74e-01
    ##  -5.73e-02
    ##   4.31e-02
    ##  -6.05e-02
    ##   9.24e-02
    ##   1.94e-01
    ##   9.22e-02
    ##   2.30e-01
    ##   1.45e-02
    ##   1.10e-01
    ##   6.31e-03
    ##   1.43e-01
    ##  -8.47e-03
    ## 
    ## exposition = south, survey_year_fct = 2021:
    ##  contrast                                                 estimate lower.HPD
    ##  dry_grassland sand_ratio0 - hay_meadow sand_ratio0      -0.158733  -0.23340
    ##  hay_meadow sand_ratio25 - hay_meadow sand_ratio0        -0.038170  -0.11447
    ##  hay_meadow sand_ratio25 - dry_grassland sand_ratio0      0.120111   0.04643
    ##  dry_grassland sand_ratio25 - hay_meadow sand_ratio0     -0.087189  -0.16084
    ##  dry_grassland sand_ratio25 - dry_grassland sand_ratio0   0.071751  -0.00088
    ##  dry_grassland sand_ratio25 - hay_meadow sand_ratio25    -0.048811  -0.12549
    ##  hay_meadow sand_ratio50 - hay_meadow sand_ratio0        -0.023782  -0.09768
    ##  hay_meadow sand_ratio50 - dry_grassland sand_ratio0      0.134623   0.05994
    ##  hay_meadow sand_ratio50 - hay_meadow sand_ratio25        0.014904  -0.06020
    ##  hay_meadow sand_ratio50 - dry_grassland sand_ratio25     0.063010  -0.00773
    ##  dry_grassland sand_ratio50 - hay_meadow sand_ratio0     -0.080285  -0.15092
    ##  dry_grassland sand_ratio50 - dry_grassland sand_ratio0   0.078128   0.00541
    ##  dry_grassland sand_ratio50 - hay_meadow sand_ratio25    -0.040884  -0.11775
    ##  dry_grassland sand_ratio50 - dry_grassland sand_ratio25  0.007716  -0.06802
    ##  dry_grassland sand_ratio50 - hay_meadow sand_ratio50    -0.055981  -0.13120
    ##  upper.HPD
    ##  -8.72e-02
    ##   3.42e-02
    ##   1.96e-01
    ##  -1.21e-02
    ##   1.45e-01
    ##   2.53e-02
    ##   4.90e-02
    ##   2.07e-01
    ##   8.97e-02
    ##   1.41e-01
    ##  -3.31e-03
    ##   1.55e-01
    ##   3.21e-02
    ##   8.20e-02
    ##   1.90e-02
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
    ##  [85] evaluate_0.19        arrayhelpers_1.1-0   xtable_1.8-4        
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
    ## [130] digest_0.6.30        vegan_2.6-4          rmarkdown_2.19      
    ## [133] cellranger_1.1.0     gap_1.3-1            curl_4.3.3          
    ## [136] shiny_1.7.3          gtools_3.9.4         nloptr_2.0.3        
    ## [139] lifecycle_1.0.3      nlme_3.1-160         jsonlite_1.8.4      
    ## [142] seqinr_4.2-23        fansi_1.0.3          pillar_1.8.1        
    ## [145] lattice_0.20-45      fastmap_1.1.0        httr_1.4.4          
    ## [148] pkgbuild_1.4.0       glue_1.6.2           xts_0.12.2          
    ## [151] diffobj_0.3.5        iterators_1.0.14     png_0.1-8           
    ## [154] shinythemes_1.2.0    bit_4.0.5            stringi_1.7.8       
    ## [157] latticeExtra_0.6-30  ape_5.6-2
