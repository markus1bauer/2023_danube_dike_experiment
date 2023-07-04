Analysis of Bauer et al. (2023) bioRxiv: <br> Favourable Conservation
Status (FCS)
================
<b>Markus Bauer</b> <br>
<b>2023-07-04</b>

- [Preparation](#preparation)
- [Statistics](#statistics)
  - [Data exploration](#data-exploration)
    - [Graphs of raw data](#graphs-of-raw-data)
    - [Outliers, zero-inflation,
      transformations?](#outliers-zero-inflation-transformations)
  - [Models](#models)
    - [Load models (BARG 6.G)](#load-models-barg-6g)
    - [Model specifications (BARG 1.D)](#model-specifications-barg-1d)
    - [Preparation for analysis](#preparation-for-analysis)
    - [Priors (BARG 1.D/E)](#priors-barg-1de)
  - [Model check](#model-check)
    - [Check computation (MCMC diagnostics, BARG
      2.B/C)](#check-computation-mcmc-diagnostics-barg-2bc)
    - [Posterior predictive check (BARG
      3.A)](#posterior-predictive-check-barg-3a)
    - [DHARMa](#dharma)
  - [Model comparison](#model-comparison)
    - [Conditional <i>R</i><sup>2</sup> values](#conditional-r2-values)
    - [Marginal <i>R</i><sup>2</sup> values](#marginal-r2-values)
    - [Bayes factor (BARG 3.C)](#bayes-factor-barg-3c)
  - [Posterior distributions (BARG
    3.B)](#posterior-distributions-barg-3b)
    - [Forest plot of effect sizes (BARG
      3.B/5.B)](#forest-plot-of-effect-sizes-barg-3b5b)
    - [Effect sizes](#effect-sizes)
- [Session info (BARG 2.A/6.A/6.B)](#session-info-barg-2a6a6b)

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

Favourable Conservation Status (FCS) sensu Helm et al. (2015) Divers
Distrib [DOI: 10.1111/ddi.12285](https://doi.org/10.1111/ddi.12285)

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
library(emmeans)
```

#### Load data

``` r
sites <- read_csv(
  here("data", "processed", "data_processed_sites_spatial.csv"),
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
    botanist_year = str_c(survey_year, botanist, exposition, sep = " "),
    botanist_year = factor(botanist_year),
    n = fcs_target,
    id = factor(id)
    ) %>%
  select(
    id, plot, site, exposition, sand_ratio, substrate_depth, target_type,
    seed_density, survey_year_fct, survey_year, botanist_year, n
    )
sites
```

    ## # A tibble: 1,152 × 12
    ##    id        plot  site  expos…¹ sand_…² subst…³ targe…⁴ seed_…⁵ surve…⁶ surve…⁷
    ##    <fct>     <fct> <fct> <fct>   <fct>   <fct>   <fct>   <fct>   <fct>   <chr>  
    ##  1 L1_01_20… L1_01 1     south   0       15      hay_me… 4       2018    2018   
    ##  2 L1_01_20… L1_01 1     south   0       15      hay_me… 4       2019    2019   
    ##  3 L1_01_20… L1_01 1     south   0       15      hay_me… 4       2020    2020   
    ##  4 L1_01_20… L1_01 1     south   0       15      hay_me… 4       2021    2021   
    ##  5 L1_02_20… L1_02 1     south   0       15      dry_gr… 4       2018    2018   
    ##  6 L1_02_20… L1_02 1     south   0       15      dry_gr… 4       2019    2019   
    ##  7 L1_02_20… L1_02 1     south   0       15      dry_gr… 4       2020    2020   
    ##  8 L1_02_20… L1_02 1     south   0       15      dry_gr… 4       2021    2021   
    ##  9 L1_03_20… L1_03 1     south   0       15      dry_gr… 8       2018    2018   
    ## 10 L1_03_20… L1_03 1     south   0       15      dry_gr… 8       2019    2019   
    ## # … with 1,142 more rows, 2 more variables: botanist_year <fct>, n <dbl>, and
    ## #   abbreviated variable names ¹​exposition, ²​sand_ratio, ³​substrate_depth,
    ## #   ⁴​target_type, ⁵​seed_density, ⁶​survey_year_fct, ⁷​survey_year

# Statistics

## Data exploration

### Graphs of raw data

![](model_check_fcs_files/figure-gfm/data-exploration-1.png)<!-- -->![](model_check_fcs_files/figure-gfm/data-exploration-2.png)<!-- -->![](model_check_fcs_files/figure-gfm/data-exploration-3.png)<!-- -->![](model_check_fcs_files/figure-gfm/data-exploration-4.png)<!-- -->![](model_check_fcs_files/figure-gfm/data-exploration-5.png)<!-- -->

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

![](model_check_fcs_files/figure-gfm/outliers-1.png)<!-- -->

## Models

### Load models (BARG 6.G)

Only here you have to modify the script to compare other models

``` r
load(file = here("outputs", "models", "model_fcs_2.Rdata"))
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
load(file = here("outputs", "models", "model_fcs_full.Rdata"))
load(file = here("outputs", "models", "model_fcs_2_prior.Rdata"))
# BARG 5.A/B/C
load(file = here("outputs", "models", "model_fcs_2_flat.Rdata"))
m_1 <- m2
m_2 <- m_full
m_prior <- m2_prior
m_flat <- m2_flat
```

### Model specifications (BARG 1.D)

``` r
m_1$formula
## n ~ sand_ratio * target_type * exposition * survey_year_fct + substrate_depth + seed_density + substrate_depth:exposition + seed_density:exposition + substrate_depth:survey_year_fct + seed_density:survey_year_fct + botanist_year + (1 | site/plot)
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

    ##                   prior     class                                coef     group
    ##                  (flat)         b                                              
    ##                  (flat)         b    botanist_year2018JakobHubersouth          
    ##                  (flat)         b    botanist_year2018SimonReithnorth          
    ##                  (flat)         b    botanist_year2018SimonReithsouth          
    ##                  (flat)         b    botanist_year2019JakobHubernorth          
    ##                  (flat)         b  botanist_year2019LindaWegglersouth          
    ##                  (flat)         b   botanist_year2019MarkusBauersouth          
    ##                  (flat)         b    botanist_year2020JakobHubernorth          
    ##                  (flat)         b    botanist_year2020JakobHubersouth          
    ##                  (flat)         b botanist_year2020KatharinaBecknorth          
    ##                  (flat)         b botanist_year2020KatharinaBecksouth          
    ##                  (flat)         b   botanist_year2021MarkusBauernorth          
    ##                  (flat)         b   botanist_year2021MarkusBauersouth          
    ##                  (flat)         b                     expositionsouth          
    ##                  (flat)         b                        sand_ratio25          
    ##                  (flat)         b                        sand_ratio50          
    ##                  (flat)         b                       seed_density8          
    ##                  (flat)         b                   substrate_depth15          
    ##                  (flat)         b                 survey_year_fct2019          
    ##                  (flat)         b                 survey_year_fct2020          
    ##                  (flat)         b                 survey_year_fct2021          
    ##                  (flat)         b            target_typedry_grassland          
    ##  student_t(3, 0.2, 2.5) Intercept                                              
    ##    student_t(3, 0, 2.5)        sd                                              
    ##                  (flat)        sd                                          site
    ##                  (flat)        sd                           Intercept      site
    ##                  (flat)        sd                                     site:plot
    ##                  (flat)        sd                           Intercept site:plot
    ##    student_t(3, 0, 2.5)     sigma                                              
    ##   source
    ##  default
    ##  default
    ##  default
    ##  default
    ##  default
    ##  default
    ##  default
    ##  default
    ##  default
    ##  default
    ##  default
    ##  default
    ##  default
    ##  default
    ##  default
    ##  default
    ##  default
    ##  default
    ##  default
    ##  default
    ##  default
    ##  default
    ##  default
    ##  default
    ##  default
    ##  default
    ##  default
    ##  default
    ##  default

``` r
data <- data.frame(x = c(-5, 5))
plot1 <- ggplot(data, aes(x = x)) +
  stat_function(fun = dnorm, n = 101, args = list(mean = 0, sd = 2)) +
  expand_limits(y = 0) +
  ggtitle("Normal distribution for Intercept") + theme_mb()
plot2 <- ggplot(data, aes(x = x)) +
  stat_function(fun = dnorm, n = 101, args = list(mean = 0.3, sd = 2)) +
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

![](model_check_fcs_files/figure-gfm/possible-priors-1.png)<!-- -->

#### Prior summary (BARG 1.D)

``` r
m_1 %>%
  brms::prior_summary(all = FALSE) %>%
  select(prior, class, coef, group, source)
```

    ##                 prior     class                coef group  source
    ##          normal(0, 2)         b                              user
    ##        normal(0.1, 2)         b        sand_ratio25          user
    ##        normal(0.2, 2)         b        sand_ratio50          user
    ##        normal(0.1, 2)         b survey_year_fct2019          user
    ##        normal(0.2, 2)         b survey_year_fct2020          user
    ##        normal(0.3, 2)         b survey_year_fct2021          user
    ##          normal(0, 2) Intercept                              user
    ##  student_t(3, 0, 2.5)        sd                           default
    ##          cauchy(0, 1)     sigma                              user

#### Prior predictive check (BARG 1.E)

The blue line (observed value) should be in the center of the
distribution of 500 simulated data sets.

``` r
bayesplot::ppc_stat(y, yrep_prior[1:500, ], binwidth = 0.5) +
  coord_cartesian(xlim = c(-4, 4)) + bayesplot::theme_default()
```

![](model_check_fcs_files/figure-gfm/prior-predictive-check-1.png)<!-- -->

``` r
ppc_stat_grouped(
  y, yrep_prior[1:500, ], group = sites$site, binwidth = 0.5
  ) +
  coord_cartesian(xlim = c(-4, 4)) + bayesplot::theme_default()
```

![](model_check_fcs_files/figure-gfm/prior-predictive-check-2.png)<!-- -->

``` r
ppc_stat_grouped(
  y, yrep_prior[1:500, ], group = sites$exposition, binwidth = 0.5
  ) +
  coord_cartesian(xlim = c(-4, 4)) + bayesplot::theme_default()
```

![](model_check_fcs_files/figure-gfm/prior-predictive-check-3.png)<!-- -->

``` r
ppc_stat_grouped(
  y, yrep_prior[1:500, ], group = sites$survey_year_fct, binwidth = 0.5
  ) +
  coord_cartesian(xlim = c(-4, 4)) + bayesplot::theme_default()
```

![](model_check_fcs_files/figure-gfm/prior-predictive-check-4.png)<!-- -->

``` r
ppc_stat_grouped(
  y, yrep_prior[1:500, ], group = sites$target_type, binwidth = 0.5
  ) +
  coord_cartesian(xlim = c(-4, 4)) + bayesplot::theme_default()
```

![](model_check_fcs_files/figure-gfm/prior-predictive-check-5.png)<!-- -->

``` r
ppc_stat_grouped(
  y, yrep_prior[1:500, ], group = sites$seed_density, binwidth = 0.5
  ) +
  coord_cartesian(xlim = c(-4, 4)) + bayesplot::theme_default()
```

![](model_check_fcs_files/figure-gfm/prior-predictive-check-6.png)<!-- -->

``` r
ppc_stat_grouped(
  y, yrep_prior[1:500, ], group = sites$sand_ratio, binwidth = 0.5
  ) +
  coord_cartesian(xlim = c(-4, 4)) + bayesplot::theme_default()
```

![](model_check_fcs_files/figure-gfm/prior-predictive-check-7.png)<!-- -->

``` r
ppc_stat_grouped(
  y, yrep_prior[1:500, ], group = sites$substrate_depth, binwidth = 0.5
  ) +
  coord_cartesian(xlim = c(-4, 4)) + bayesplot::theme_default()
```

![](model_check_fcs_files/figure-gfm/prior-predictive-check-8.png)<!-- -->

## Model check

### Check computation (MCMC diagnostics, BARG 2.B/C)

#### Trace plots

The traces of the four chains should overlap each other.

``` r
bayesplot::mcmc_trace(
  posterior1, np = hmc_diagnostics1, facet_args = list(ncol = 2)
  ) + bayesplot::theme_default()
```

    ## No divergences to plot.

![](model_check_fcs_files/figure-gfm/mcmc-trace-1.png)<!-- -->

``` r
bayesplot::mcmc_trace(
  posterior2, np = hmc_diagnostics2, facet_args = list(ncol = 2)
  ) + bayesplot::theme_default()
```

    ## No divergences to plot.

![](model_check_fcs_files/figure-gfm/mcmc-trace-2.png)<!-- -->

#### Sampling efficency: R-hat (BARG 2.B)

R-hat values should be lower than 1.1.

``` r
m_1 %>% brms::rhat() %>% bayesplot::mcmc_rhat() + bayesplot::theme_default() + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
```

![](model_check_fcs_files/figure-gfm/rhat-1.png)<!-- -->

``` r
m_2 %>% brms::rhat() %>% bayesplot::mcmc_rhat() + bayesplot::theme_default() + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
```

![](model_check_fcs_files/figure-gfm/rhat-2.png)<!-- -->

#### Sampling effectiveness: Effective sampling size (ESS) (BARG 2.C)

ESS should be greater than 0.1.

``` r
m_1 %>% brms::neff_ratio() %>% bayesplot::mcmc_neff() + bayesplot::theme_default() + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
```

![](model_check_fcs_files/figure-gfm/ess-1.png)<!-- -->

``` r
m_2 %>% brms::neff_ratio() %>% bayesplot::mcmc_neff() + bayesplot::theme_default() + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
```

![](model_check_fcs_files/figure-gfm/ess-2.png)<!-- -->

#### Pairs plot

If the technical problem of divergent transitions occurs, points would
be highlighted in green.

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

![](model_check_fcs_files/figure-gfm/mcmc-pairs-1.png)<!-- -->

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

![](model_check_fcs_files/figure-gfm/mcmc-pairs-2.png)<!-- -->

#### Parallel coordinate plot

If the technical problem of divergent transitions occurs, points would
be highlighted in green.

``` r
posterior1 %>% bayesplot::mcmc_parcoord(np = hmc_diagnostics1) + bayesplot::theme_default() + theme(axis.text.x = element_text(angle = 45))
```

![](model_check_fcs_files/figure-gfm/mcmc-parcoord-1.png)<!-- -->

``` r
posterior2 %>% bayesplot::mcmc_parcoord(np = hmc_diagnostics2) + bayesplot::theme_default() + theme(axis.text.x = element_text(angle = 45))
```

![](model_check_fcs_files/figure-gfm/mcmc-parcoord-2.png)<!-- -->

#### Autocorrelation check

Autocorrelation should go down to zero very fast.

``` r
posterior1 %>% bayesplot::mcmc_acf(lags = 10) + bayesplot::theme_default()
```

![](model_check_fcs_files/figure-gfm/autocorrelation-1.png)<!-- -->

``` r
posterior2 %>% bayesplot::mcmc_acf(lags = 10) + bayesplot::theme_default()
```

![](model_check_fcs_files/figure-gfm/autocorrelation-2.png)<!-- -->

### Posterior predictive check (BARG 3.A)

#### Kernel density

The kernel density estimate of the observed data (dark line) compared
with estimates for 50 simulated data sets drawn from the posterior
predictive distribution (thinner and lighter lines). The dark line
should be near the simulated curves.

``` r
p1 <- bayesplot::ppc_dens_overlay(y, yrep1[1:50, ]) + bayesplot::theme_default()
p2 <- bayesplot::ppc_dens_overlay(y, yrep2[1:50, ]) + bayesplot::theme_default()
p1 / p2
```

![](model_check_fcs_files/figure-gfm/kernel-density-1.png)<!-- -->

``` r
p1 <- ppc_dens_overlay_grouped(y, yrep1[1:50, ], group = sites$site) + bayesplot::theme_default()
p2 <- ppc_dens_overlay_grouped(y, yrep2[1:50, ], group = sites$site) + bayesplot::theme_default()
p1 / p2
```

![](model_check_fcs_files/figure-gfm/kernel-density-2.png)<!-- -->

``` r
p1 <- ppc_dens_overlay_grouped(y, yrep1[1:50, ], group = sites$exposition) + bayesplot::theme_default()
p2 <- ppc_dens_overlay_grouped(y, yrep2[1:50, ], group = sites$exposition) + bayesplot::theme_default()
p1 / p2
```

![](model_check_fcs_files/figure-gfm/kernel-density-3.png)<!-- -->

``` r
p1 <- ppc_dens_overlay_grouped(y, yrep1[1:50, ], group = sites$survey_year_fct) + bayesplot::theme_default()
p2 <- ppc_dens_overlay_grouped(y, yrep2[1:50, ], group = sites$survey_year_fct) + bayesplot::theme_default()
p1 / p2
```

![](model_check_fcs_files/figure-gfm/kernel-density-4.png)<!-- -->

``` r
p1 <- ppc_dens_overlay_grouped(y, yrep1[1:50, ], group = sites$target_type) + bayesplot::theme_default()
p2 <- ppc_dens_overlay_grouped(y, yrep2[1:50, ], group = sites$target_type) + bayesplot::theme_default()
p1 / p2
```

![](model_check_fcs_files/figure-gfm/kernel-density-5.png)<!-- -->

``` r
p1 <- ppc_dens_overlay_grouped(y, yrep1[1:50, ], group = sites$seed_density) + bayesplot::theme_default()
p2 <- ppc_dens_overlay_grouped(y, yrep2[1:50, ], group = sites$seed_density) + bayesplot::theme_default()
p1 / p2
```

![](model_check_fcs_files/figure-gfm/kernel-density-6.png)<!-- -->

``` r
p1 <- ppc_dens_overlay_grouped(y, yrep1[1:50, ], group = sites$sand_ratio) + bayesplot::theme_default()
p2 <- ppc_dens_overlay_grouped(y, yrep2[1:50, ], group = sites$sand_ratio) + bayesplot::theme_default()
p1 / p2
```

![](model_check_fcs_files/figure-gfm/kernel-density-7.png)<!-- -->

``` r
p1 <- ppc_dens_overlay_grouped(y, yrep1[1:50, ], group = sites$substrate_depth) + bayesplot::theme_default()
p2 <- ppc_dens_overlay_grouped(y, yrep2[1:50, ], group = sites$substrate_depth) + bayesplot::theme_default()
p1 / p2
```

![](model_check_fcs_files/figure-gfm/kernel-density-8.png)<!-- -->

#### Histograms of statistics skew

The dark line (computed from the observed data) should be in the center
of the posterior predictive distribution.

``` r
p1 <- bayesplot::ppc_stat(y, yrep1, binwidth = 0.001) + bayesplot::theme_default()
p2 <- bayesplot::ppc_stat(y, yrep2, binwidth = 0.001) + bayesplot::theme_default()
p1 / p2
```

![](model_check_fcs_files/figure-gfm/histogram-1.png)<!-- -->

``` r
p1 <- ppc_stat_grouped(y, yrep1, group = sites$site, binwidth = 0.001) + bayesplot::theme_default()
p2 <- ppc_stat_grouped(y, yrep2, group = sites$site, binwidth = 0.001) + bayesplot::theme_default()
p1 / p2
```

![](model_check_fcs_files/figure-gfm/histogram-2.png)<!-- -->

``` r
p1 <- ppc_stat_grouped(y, yrep1, group = sites$exposition, binwidth = 0.001) + bayesplot::theme_default()
p2 <- ppc_stat_grouped(y, yrep2, group = sites$exposition, binwidth = 0.001) + bayesplot::theme_default()
p1 / p2
```

![](model_check_fcs_files/figure-gfm/histogram-3.png)<!-- -->

``` r
p1 <- ppc_stat_grouped(y, yrep1, group = sites$survey_year_fct, binwidth = 0.001) + bayesplot::theme_default()
p2 <- ppc_stat_grouped(y, yrep2, group = sites$survey_year_fct, binwidth = 0.001) + bayesplot::theme_default()
p1 / p2
```

![](model_check_fcs_files/figure-gfm/histogram-4.png)<!-- -->

``` r
p1 <- ppc_stat_grouped(y, yrep1, group = sites$target_type, binwidth = 0.001) + bayesplot::theme_default()
p2 <- ppc_stat_grouped(y, yrep2, group = sites$target_type, binwidth = 0.001) + bayesplot::theme_default()
p1 / p2
```

![](model_check_fcs_files/figure-gfm/histogram-5.png)<!-- -->

``` r
p1 <- ppc_stat_grouped(y, yrep1, group = sites$seed_density, binwidth = 0.001) + bayesplot::theme_default()
p2 <- ppc_stat_grouped(y, yrep2, group = sites$seed_density, binwidth = 0.001) + bayesplot::theme_default()
p1 / p2
```

![](model_check_fcs_files/figure-gfm/histogram-6.png)<!-- -->

``` r
p1 <- ppc_stat_grouped(y, yrep1, group = sites$sand_ratio, binwidth = 0.001) + bayesplot::theme_default()
p2 <- ppc_stat_grouped(y, yrep2, group = sites$sand_ratio, binwidth = 0.001) + bayesplot::theme_default()
p1 / p2
```

![](model_check_fcs_files/figure-gfm/histogram-7.png)<!-- -->

``` r
p1 <- ppc_stat_grouped(y, yrep1, group = sites$substrate_depth, binwidth = 0.001) + bayesplot::theme_default()
p2 <- ppc_stat_grouped(y, yrep2, group = sites$substrate_depth, binwidth = 0.001) + bayesplot::theme_default()
p1 / p2
```

![](model_check_fcs_files/figure-gfm/histogram-8.png)<!-- -->

#### LOO cross-validation (Leave one out)

All Pareto-k estimates should be \< 0.7.

``` r
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
    ## elpd_loo   -466.8 30.3
    ## p_loo       183.9  8.6
    ## looic       933.5 60.6
    ## ------
    ## Monte Carlo SE of elpd_loo is 0.3.
    ## 
    ## Pareto k diagnostic values:
    ##                          Count Pct.    Min. n_eff
    ## (-Inf, 0.5]   (good)     1143  99.2%   701       
    ##  (0.5, 0.7]   (ok)          9   0.8%   338       
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
    ## elpd_loo   -473.3 30.3
    ## p_loo       188.3  8.8
    ## looic       946.6 60.7
    ## ------
    ## Monte Carlo SE of elpd_loo is 0.3.
    ## 
    ## Pareto k diagnostic values:
    ##                          Count Pct.    Min. n_eff
    ## (-Inf, 0.5]   (good)     1136  98.6%   985       
    ##  (0.5, 0.7]   (ok)         16   1.4%   345       
    ##    (0.7, 1]   (bad)         0   0.0%   <NA>      
    ##    (1, Inf)   (very bad)    0   0.0%   <NA>      
    ## 
    ## All Pareto k estimates are ok (k < 0.7).
    ## See help('pareto-k-diagnostic') for details.

``` r
plot(loo1)
```

![](model_check_fcs_files/figure-gfm/loo-1.png)<!-- -->

``` r
plot(loo2)
```

![](model_check_fcs_files/figure-gfm/loo-2.png)<!-- -->

Density of the computed leave one out probability integral transform
(LOO PIT; dark curve) compared with the posterior predictive
distribution (light curves). The LOO PIT curve should be within the
range of light curves.

``` r
p1 <- bayesplot::ppc_loo_pit_overlay(y, yrep1, lw = weights(loo1$psis_object)) + bayesplot::theme_default()
p2 <- bayesplot::ppc_loo_pit_overlay(y, yrep2, lw = weights(loo2$psis_object)) + bayesplot::theme_default()
p1 / p2
```

![](model_check_fcs_files/figure-gfm/loo-pit-1.png)<!-- -->

### DHARMa

``` r
m_1 %>% DHARMa.helpers::dh_check_brms(integer = TRUE)
```

![](model_check_fcs_files/figure-gfm/dharma-1.png)<!-- -->

``` r
m_2 %>% DHARMa.helpers::dh_check_brms(integer = TRUE)
```

![](model_check_fcs_files/figure-gfm/dharma-2.png)<!-- -->

## Model comparison

### Conditional <i>R</i><sup>2</sup> values

``` r
m_1 %>% brms::bayes_R2(
  probs = c(0.05, 0.5, 0.95), re_formula =  ~ (1 | site / plot)
  )
##     Estimate   Est.Error        Q5       Q50       Q95
## R2 0.8464731 0.005100259 0.8376967 0.8466772 0.8544419
m_2 %>% brms::bayes_R2(
  probs = c(0.05, 0.5, 0.95), re_formula =  ~ (1 | site / plot)
  )
##     Estimate  Est.Error        Q5       Q50       Q95
## R2 0.8458954 0.00515156 0.8369915 0.8461621 0.8539969
```

### Marginal <i>R</i><sup>2</sup> values

``` r
m_1 %>% brms::bayes_R2(
  probs = c(0.05, 0.5, 0.95), re_formula = 1 ~ 1
  )
##     Estimate   Est.Error        Q5       Q50       Q95
## R2 0.8084947 0.004562813 0.8005822 0.8087267 0.8155298
m_2 %>% brms::bayes_R2(
  probs = c(0.05, 0.5, 0.95), re_formula = 1 ~ 1
  )
##    Estimate   Est.Error        Q5       Q50       Q95
## R2 0.808271 0.004648837 0.8002222 0.8084726 0.8155058
```

### Bayes factor (BARG 3.C)

``` r
bayes_factor <- brms::bayes_factor(m_1, m_2)
```

``` r
bayes_factor
```

    ## Estimated Bayes factor in favor of m_1 over m_2: 576256234.77279

## Posterior distributions (BARG 3.B)

### Forest plot of effect sizes (BARG 3.B/5.B)

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
  coord_cartesian(xlim = c(-.65, 1)) +
  bayesplot::theme_default() +
  ggtitle("Posterior distbributions (mean, CI66, CI95)")
```

![](model_check_fcs_files/figure-gfm/posteriors-1.png)<!-- -->

### Effect sizes

Effect sizes of chosen model just to get exact values of means etc. if
necessary.

``` r
posterior1 %>%
  posterior::summarize_draws() %>%
  knitr::kable()
```

| variable                   |    mean |  median |      sd |     mad |       q5 |      q95 | rhat | ess_bulk | ess_tail |
|:---------------------------|--------:|--------:|--------:|--------:|---------:|---------:|-----:|---------:|---------:|
| Intercept                  |  0.0643 |  0.0642 |  0.0409 |  0.0319 |  0.00157 |    0.127 | 1.00 |   8070\. |   8255\. |
| b_sand_ratio25             |   0.184 |   0.183 |   0.103 |   0.104 |   0.0142 |    0.354 | 1.00 |   4956\. |   7020\. |
| b_sand_ratio50             |   0.161 |   0.162 |   0.104 |   0.104 | -0.00906 |    0.334 | 1.00 |   4749\. |   7063\. |
| b_substrate_depth15        | -0.0831 | -0.0829 |  0.0507 |  0.0504 |   -0.167 | 0.000106 | 1.00 |   9531\. |   8645\. |
| b_target_typedry_grassland |  0.0219 |  0.0233 |   0.102 |   0.101 |   -0.145 |    0.190 | 1.00 |   3964\. |   6678\. |
| b_seed_density8            |   0.110 |   0.110 |  0.0509 |  0.0510 |   0.0244 |    0.192 | 1.00 |   9792\. |   9571\. |
| b_expositionsouth          |  -0.390 |  -0.398 |   0.956 |   0.957 |    -1.96 |     1.16 | 1.00 |   8893\. |   8547\. |
| b_survey_year_fct2019      |   0.764 |   0.753 |    1.25 |    1.26 |    -1.29 |     2.82 | 1.00 |   9685\. |   9376\. |
| b_survey_year_fct2020      |    1.04 |    1.04 |    1.05 |    1.06 |   -0.701 |     2.76 | 1.00 |   9203\. |   9453\. |
| b_survey_year_fct2021      |    1.03 |    1.03 |    1.28 |    1.29 |    -1.09 |     3.15 | 1.00 |   9148\. |   9344\. |
| sd_site\_\_Intercept       |  0.0826 |  0.0718 |  0.0489 |  0.0334 |   0.0308 |    0.169 | 1.00 |   4787\. |   6153\. |
| sd_site:plot\_\_Intercept  |   0.158 |   0.158 |  0.0153 |  0.0152 |    0.133 |    0.183 | 1.00 |   5816\. |   8105\. |
| sigma                      |   0.331 |   0.331 | 0.00817 | 0.00826 |    0.318 |    0.344 | 1.00 |   6945\. |   8599\. |

``` r
emm <- m_1 %>% emmeans(
  revpairwise ~ target_type + sand_ratio | exposition | survey_year_fct,
  type = "response"
  )
emm$emmeans
```

    ## exposition = north, survey_year_fct = 2018:
    ##  target_type   sand_ratio  emmean lower.HPD upper.HPD
    ##  hay_meadow    0          -0.7502   -0.9068  -0.58484
    ##  dry_grassland 0          -0.7290   -0.8910  -0.56640
    ##  hay_meadow    25         -0.5663   -0.7279  -0.40143
    ##  dry_grassland 25         -0.7253   -0.8929  -0.56032
    ##  hay_meadow    50         -0.5894   -0.7529  -0.41891
    ##  dry_grassland 50         -0.5484   -0.7040  -0.38147
    ## 
    ## exposition = south, survey_year_fct = 2018:
    ##  target_type   sand_ratio  emmean lower.HPD upper.HPD
    ##  hay_meadow    0          -1.3221   -1.4897  -1.15812
    ##  dry_grassland 0          -1.1631   -1.3298  -1.00197
    ##  hay_meadow    25         -1.6033   -1.7691  -1.43755
    ##  dry_grassland 25         -1.2084   -1.3771  -1.04165
    ##  hay_meadow    50         -1.8508   -2.0273  -1.69469
    ##  dry_grassland 50         -1.4862   -1.6469  -1.31592
    ## 
    ## exposition = north, survey_year_fct = 2019:
    ##  target_type   sand_ratio  emmean lower.HPD upper.HPD
    ##  hay_meadow    0           0.6318    0.4620   0.79461
    ##  dry_grassland 0           0.5361    0.3732   0.69943
    ##  hay_meadow    25          0.4623    0.2990   0.63338
    ##  dry_grassland 25          0.4124    0.2397   0.57344
    ##  hay_meadow    50          0.6829    0.5252   0.85295
    ##  dry_grassland 50          0.7133    0.5539   0.87910
    ## 
    ## exposition = south, survey_year_fct = 2019:
    ##  target_type   sand_ratio  emmean lower.HPD upper.HPD
    ##  hay_meadow    0          -0.0691   -0.2323   0.10101
    ##  dry_grassland 0          -0.1774   -0.3374  -0.00573
    ##  hay_meadow    25         -0.0834   -0.2594   0.08456
    ##  dry_grassland 25         -0.0611   -0.2407   0.11113
    ##  hay_meadow    50         -0.1366   -0.3148   0.03230
    ##  dry_grassland 50          0.0598   -0.1165   0.23225
    ## 
    ## exposition = north, survey_year_fct = 2020:
    ##  target_type   sand_ratio  emmean lower.HPD upper.HPD
    ##  hay_meadow    0           0.8009    0.6462   0.97538
    ##  dry_grassland 0           0.8611    0.7016   1.02550
    ##  hay_meadow    25          0.8020    0.6413   0.97112
    ##  dry_grassland 25          0.7707    0.6012   0.93117
    ##  hay_meadow    50          0.8951    0.7276   1.05677
    ##  dry_grassland 50          0.8659    0.7017   1.02915
    ## 
    ## exposition = south, survey_year_fct = 2020:
    ##  target_type   sand_ratio  emmean lower.HPD upper.HPD
    ##  hay_meadow    0           0.0684   -0.0917   0.23502
    ##  dry_grassland 0          -0.0843   -0.2444   0.08771
    ##  hay_meadow    25          0.1484   -0.0129   0.31666
    ##  dry_grassland 25          0.2575    0.0931   0.42323
    ##  hay_meadow    50          0.0542   -0.1152   0.21635
    ##  dry_grassland 50          0.3038    0.1378   0.46850
    ## 
    ## exposition = north, survey_year_fct = 2021:
    ##  target_type   sand_ratio  emmean lower.HPD upper.HPD
    ##  hay_meadow    0           0.9003    0.7430   1.07066
    ##  dry_grassland 0           0.9205    0.7595   1.09095
    ##  hay_meadow    25          0.9448    0.7839   1.11503
    ##  dry_grassland 25          0.9055    0.7506   1.07630
    ##  hay_meadow    50          0.9900    0.8340   1.15768
    ##  dry_grassland 50          1.0142    0.8384   1.16759
    ## 
    ## exposition = south, survey_year_fct = 2021:
    ##  target_type   sand_ratio  emmean lower.HPD upper.HPD
    ##  hay_meadow    0           0.2723    0.1117   0.44244
    ##  dry_grassland 0           0.1319   -0.0334   0.29690
    ##  hay_meadow    25          0.1462   -0.0198   0.30584
    ##  dry_grassland 25          0.2967    0.1343   0.46849
    ##  hay_meadow    50          0.1496   -0.0195   0.30858
    ##  dry_grassland 50          0.2329    0.0620   0.39194
    ## 
    ## Results are averaged over the levels of: substrate_depth, seed_density, botanist_year 
    ## Point estimate displayed: median 
    ## HPD interval probability: 0.95

# Session info (BARG 2.A/6.A/6.B)

    ## R version 4.2.3 (2023-03-15 ucrt)
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
    ##   [4] igraph_1.4.1         sp_1.6-0             splines_4.2.3       
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
    ##  [52] xtable_1.8-4         bit_4.0.5            stats4_4.2.3        
    ##  [55] StanHeaders_2.26.13  DT_0.27              httr_1.4.5          
    ##  [58] htmlwidgets_1.6.1    threejs_0.3.3        RColorBrewer_1.1-3  
    ##  [61] posterior_1.4.0      ellipsis_0.3.2       XML_3.99-0.13       
    ##  [64] pkgconfig_2.0.3      farver_2.1.1         qgam_1.3.4          
    ##  [67] deldir_1.0-6         utf8_1.2.3           tidyselect_1.2.0    
    ##  [70] labeling_0.4.2       rlang_1.0.6          reshape2_1.4.4      
    ##  [73] later_1.3.0          munsell_0.5.0        tools_4.2.3         
    ##  [76] cli_3.6.0            generics_0.1.3       ade4_1.7-22         
    ##  [79] evaluate_0.20        fastmap_1.1.1        yaml_2.3.7          
    ##  [82] processx_3.8.0       knitr_1.42           bit64_4.0.5         
    ##  [85] nlme_3.1-162         mime_0.12            adegenet_2.1.10     
    ##  [88] xml2_1.3.3           gap_1.5-1            compiler_4.2.3      
    ##  [91] shinythemes_1.2.0    rstudioapi_0.14      beeswarm_0.4.0      
    ##  [94] curl_5.0.0           png_0.1-8            RNeXML_2.4.11       
    ##  [97] stringi_1.7.12       highr_0.10           ps_1.7.2            
    ## [100] Brobdingnag_1.2-9    lattice_0.20-45      Matrix_1.5-3        
    ## [103] nloptr_2.0.3         markdown_1.5         permute_0.9-7       
    ## [106] vegan_2.6-4          shinyjs_2.1.0        tensorA_0.36.2      
    ## [109] vctrs_0.5.2          pillar_1.8.1         lifecycle_1.0.3     
    ## [112] bridgesampling_1.1-2 estimability_1.4.1   httpuv_1.6.9        
    ## [115] R6_2.5.1             latticeExtra_0.6-30  promises_1.2.0.1    
    ## [118] KernSmooth_2.23-20   gridExtra_2.3        vipor_0.4.5         
    ## [121] codetools_0.2-19     boot_1.3-28.1        colourpicker_1.2.0  
    ## [124] MASS_7.3-58.2        gtools_3.9.4         rprojroot_2.0.3     
    ## [127] withr_2.5.0          shinystan_2.6.0      mgcv_1.8-41         
    ## [130] parallel_4.2.3       hms_1.1.2            grid_4.2.3          
    ## [133] coda_0.19-4          minqa_1.2.5          rmarkdown_2.20      
    ## [136] shiny_1.7.4          base64enc_0.1-3      dygraphs_1.1.1.6    
    ## [139] interp_1.1-3
