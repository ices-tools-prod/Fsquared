
# `Fsquared`: FLR toolset for calculating ICES’ F-based reference points

# Authors

- Iago Mosqueira (WMR) <iago.mosqueira@wur.nl>
- John Trochta (IMR) <john.tyler.trochta@hi.no>
- Ernesto Jardim (IPMA) <ernesto.jardim@ipma.pt>
- Arni Magnusson (SPC) <arnim@spc.int>
- Max Cardinale (SLU) <massimiliano.cardinale@slu.se>
- Dorleta Garcia (ICES) <dorleta.garcia@ices.dk>
- Colin Millar (ICES) <colin.millar@ices.dk>

# Installation

Packages required to run this analysis (`mse`, `msemodules`, `mseviz`,
`FLSRTMB` and their dependencies) can be installed from the relevant
repositories by executing the following code:

``` r
install.packages("TAF", repos="https://cloud.r-project.org")

pkgs <- TAF::deps()

install.packages(pkgs, repos=c(FLR="https://flr.r-universe.dev",
  CRAN="https://cloud.r-project.org"))
```

The minimum versions needed are stored in the `boot/SOFTWARE.bib` file.
A check can be carried out on whether the versions available in the
system are up to date to those by calling

``` r
library(TAF)
check.software()
```

## Progress

The `mse` package makes use of `progressr`, if available, to report on
the progress of simulation runs. Standard progress bars are used by
default as defined in the local `.Rprofile` file. To switch them off
please run

``` r
handlers(global=FALSE)
```

# Introduction

This tutorial demonstrates the setting up and running of the `Fsquared`
simulation modeling toolset for calculating ICES fishing mortality
reference points (Fp05 and Fmsy) using the
[mse](https://flr-project.org/mse) package. This toolset comprises R
scripts structured according to the TAF system. The code within these R
scripts was adapted from the toolset developed for ICES WKREBUILD2, for
which a tutorial is available that describes setting up, running, and
summarising results of simulations
([here](https://htmlpreview.github.io/?https://github.com/ices-tools-prod/WKREBUILD_toolset/blob/main/tutorial.html)).
We refer users to the WKREBUILD2 tutorial for further details not
covered in this tutorial.

# Workflow

At the heart of this toolset is a shortcut MSE framework that introduces
some level of estimation error from the management procedure stock
assessment model and, uses a short-term forecast to bridge the data and
management lags, as it is commonly done for advice in ICES.

`model_fsquared.R` sets up and runs the shortcut MSE, generally in these
steps:

1.  Set simulation parameters such as the number of iterations and
    projection years.
2.  Configure and/or condition parameters of the operating model (OM)
    and the management procedure (MP).
3.  Run closed loop simulations of the standard ICES advice rule
    parameterized by the MSY Btrigger over a range of Ftarget values.
4.  Derive equilibrium distributions of stock quantities (spawning
    biomass, catches, etc.) for each Ftarget and calculate ICES’ Fp05
    and Fmsy reference points.

- Run first with it = cores, so that the call to mps() runs in parallel

# Summary and default assumptions of tool

- Operating model (om)
  - Future weights, maturity, and exploitation pattern is the average of
    the most recent 5 years (EQSIM default)
- Observation error model (oem)
  - No explicit OEM (i.e. perfect observations of spawning biomass)
- Management procedure (mp)
  - Estimation model (est)
    - Default CV (normal distribution) of ‘estimated’ log(SB) is 0
      (generated from function `shortcut_devs()`); suggested starting
      CV=0.1
  - Harvest control rule (hcr)
    - Input ‘estimated’ SB from preceding year (i.e. data_lag=1) into
      ICES advice rule to obtain target F in the next year
      (i.e. management_lag=1)
    - Target F is directly applied on OM stock
    - Target F corresponds to fbar, with the fbar age range specified in
      the input `FLStock` object
  - Implementation system (isys)
    - Short-term forecast that uses geometric mean of recruitment from
      the 3 most years
    - Applies Target F from hcr with multiplicative errors on
      ‘estimated’ SSB
    - Multiplicative errors on target F are generated as normally
      distributed and autocorrelated log(F) values from
      `shortcut_devs()`; default arguments are normal CV=0.212 and
      autocorrelation coefficient=0.423
- Implementation error model (IEM)
  - No explicit IEM

# Prerequisites

- R (\>= 4.2 recommended)
- FLR packages: `mse`, `msemodules`, `FLSRTMB` (and their dependencies)
- Local helper scripts: `utilities_fsquared.R`
- Example assessment results/objects: `boot/data/sol274.rda`

# Example: North Sea plaice assessed with SAM

Below are condensed code snippets extracted from `model_fsquared.R`.

## Specifying simulation parameters

``` r
library(mse)    # Contains MSE functionality (shortcut or full-feedback)
library(msemodules)
library(FLSRTMB)# SR estimation and uncertainty simulation

# Local helpers (adjust paths as needed)
source("utilities_fsquared.R")

# Load stock assessment results (toy example, AAP model with increased F)
# 'run' should be an FLStock-like or assessment object used downstream.
load("boot/data/sol274.rda")

# Key simulation parameters
stkname <- "sol.27.4"
srmodels <- c("segreg")          
iy <- 2022           # Starting year of projection period                 
dy <- iy - 1         # Current data year (i.e. when observations last collected)
fy <- 2055           # Final year of projection period               
it <- 5              # Number of iterations (i.e. individual projections)
set.seed(987)        # Fixed for reproducibility during initial testing
conditioning_ny <- 5 # Number of years in tail of historical period used to project biology and exploitation         
bcv_sa <- 0.10       # CV of shortcut errors on log(SSB)
fcv_sa <- 0.10       # CV of shortcut errors on log(F)
fphi_sa <- 0.0       # Autocorrelation coefficient of shortcut errors on log(F)
fg_mp <- seq(0, 1.5, length=51) 
recyrs_mp <- -2      # Maximum lags indicating start of recent historical recruitment for which geometric mean is taken to use in short term forecast          
```

What these settings mean:

- `srmodels` specifies the functional forms to estimate and bootstrap
  from in OM.  
- `iy`, `dy`, `fy` define the initial year, last year of observations,
  and final year of the MSE projection period.  
- `it` controls how many stochastic replicates are run in the OM.
- `conditioning_ny` specifies number of years backward from `dy` to use
  for conditioning future weights, maturity, and exploitation
  patterns.  
- `bcv_sa` the CV of estimation uncertainty on SB in the shortcut
  estimator.  
- `fg_mp` is the grid of candidate F targets you will test.  
- `recyrs_mp` is the range of years (from current assessment year) to
  take geometric mean of recruitment; used in short-term forecast within
  `isys`.

## Conditioning stock-recruitment

The default approach is estimate SR parameters via bootstrapping and
select the model with the largest log-likelihood via `FLSRTMB`:

``` r
spryrs <- ac(2018:2022) # years to compute mean SPR0

# Estimate SR parameters by bootstrapping and model selection
srpars <- bootstrapSR(
  run,
  iters = it,
  spr0  = mean(spr0y(run)[, spryrs]),
  models = srmodels
)

# Generate future recruitment deviances (autocorrelated lognormal)
srdevs <- rlnormar1(
  sdlog = srpars$sigmaR,
  rho = srpars$rho,
  years = seq(dy, fy),
  bias.correct = FALSE
)               
```

**TODO: Double check this is accurate and specify further** Key
functions (briefly):  
- `spr0y(run)`: computes spawning-per-recruit (SPR0) by year from the
assessment object.  
- `bootstrapSR(...)`: fits specified SR models across bootstrap
resamples, returning chosen parameters (e.g., s, sigmaR, rho).  
- `rlnormar1(...)`: simulates autocorrelated lognormal deviates with
AR(1) structure given `sdlog` and `rho`.

Notes:  
- $\sigma_R$ controls recruitment variability; $\rho$ controls temporal
autocorrelation (AR(1)).  
- `spr0` is used by certain SR models for scaling.  
- **Importantly, `bias.correct` should always be explicitly specified,
and usually will be FALSE as it is in this example, but please read the
[How to specify
`bias.correct`](#how-to-specify-bias-correct-argument-in-rlnormar1-function)
section below.**

## Creating the OM object (`FLom` class)

Key functions:  
- `FLom(...)`: creates an Operating Model object combining
stock-specific biology (stock weights, maturity, etc.), reference
points, SR model(s) and parameters, and recruitment deviances.  
- `propagate(run, it)`: replicates the `run` object for `it` number of
iterations.  
- `fwdWindow(...)`: extends the year range of the `om` object to the end
of projection period (year `fy`), carrying forward the mean of the most
recent `conditioning_ny` years.

``` r
# 'refpts' should contain reference points such as Fmsy and Btrigger
# If 'refpts' is defined in your environment (from your assessment or SS3), use it here.
# In many workflows, refpts() are derived from a fitted SR and SPR analyses.

# Build the OM
om <- FLom(
  stock = propagate(run, it),  # replicate stock object across iterations
  refpts = refpts,             # must include Fmsy, Btrigger (etc.)
  model = "mixedsrr",          # SR model type used in the OM
  params = srpars,             # SR parameters
  deviances = srdevs,          # recruitment deviances
  name = stkname
)

# Extend OM into the future, conditioning on the mean of recent years
# The `fun` argument can be one of c("mean", "geomean", "sample")
om <- fwdWindow(om, end = fy, nsq = conditioning_ny, fun = "mean")
```

## Specifying Management Procedure (MP)

The MP is consists of three modules, specified as elements within the
`mpCtrl` class:  
- Estimation (`est`): a shortcut stock assessment (with SSB
uncertainty). - HCR (`hcr`): hockey-stick rule mapping SSB to a target F
with a trigger (Btrigger) and minimum (min) floor.  
- Implementation system (`isys`): TAC setting via a one-year short-term
forecast.

``` r
# General MSE arguments: timing and lags
mseargs <- list(
  iy = iy, fy = fy,      # start and end of projections
  data_lag = 1,          # 1-year lag between observation and availability
  management_lag = 1,    # 1-year lead from advice to implementation
  frq = 1                # annual frequency of advice
)

# Shortcut estimator uncertainty (normal deviances on log(SSB) and log(F))
# To check default arguments, type `msemodules:::shortcut_devs` in your console)
sdevs <- shortcut_devs(om, SSBcv=bcv_sa, Fcv=fcv_sa, Fphi=fphi_sa)

# Define the ICES-style advice rule
arule <- mpCtrl(list(
  # Estimation: shortcut SA with SSB deviances
  est = mseCtrl(
    method = shortcut.sa,
    args = list(SSBdevs = sdevs$SSB)
  ),
  # HCR: hockey-stick maps ssb -> fbar with lim, trigger, target, min
  hcr = mseCtrl(
    method = hockeystick.hcr,
    args = list(
      lim = 0,  # 0 means the origin; when testing ICES AR, this may equal Blim 
      trigger = refpts(om)$Btrigger,
      target = refpts(om)$Fmsy,
      min = 0,
      metric = "ssb",
      output = "fbar"
    )
  ),
  # Implementation system: TAC based on F
  isys = mseCtrl(
    method = tac.is,
    args = list(recyrs = recyrs_mp, Fdevs=sdevs$F)
  )
))
```

Key functions:  
- `shortcut_devs(om, SSBcv)`: simulates estimation uncertainty in SSB
(and optionally F).  
- `mpCtrl(list(...))`: bundles MP modules needed for the shortcut MSE:
estimation (est), HCR (hcr), and implementation system (isys) models.  
- `mseCtrl(method=..., args=list(...))`: defines one module of the MP by
defining a method (function) and arguments to this method.  
- `shortcut.sa`: internal function of `mse` that multiplies estimation
errors to SSB from the OM to obtain ‘estimated’ SSB.  
- `hockeystick.hcr`: internal function of `mse` that implements the ICES
standard hockey-stick HCR that maps SSB to an F target. The ICES
standard HCR sets F target to Fmsy if SSB ≥ trigger, and linearly
reduces F target down to lim if SSB \< trigger.  
- `tac.is`: internal function of `mse` that calculates a catch advice
(TAC) by projecting F target with ‘estimated’ SSB in the current year
and short-term recruitment assumptions.

## Running the F grid

We define and run simulations of a grid of candidate F targets and run
the MSE for each value.

``` r
# Candidate F targets to test
fg <- list(target = fg_mp)

# Run MSEs over the F target grid
fgrid <- mps(
  om,
  ctrl = arule,
  args = mseargs,
  hcr = fg,
  names = paste0("F", fg_mp)
)

# 'fgrid' is typically a list-like object with outputs per F target
```

Key functions:  
- `mps(om, ctrl, args, hcr, names)`: runs multiple MSEs varying
specified HCR parameters; here we vary the target F.

## Inspecting results (quick checks)

Below are minimal examples to extract and visualize a few quantities.
Adapt as needed to your object classes.

``` r
# Example: list available runs
names(fgrid)

# Plots summary
plot(om, fgrid)

# Computs average performance over the projection years (pys)
performance(fgrid) <- performance(fgrid, statistics=icestats["PBlim"], year=pys,
  type="arule")
```

## Troubleshooting and tips

- Missing reference points: ensure `refpts` includes at least `Fmsy` and
  `Btrigger`. If not, derive them from your assessment or SR fit (e.g.,
  MSY-based analysis).  
- SR model choice: if `segreg` is unsuitable, try other SR models
  supported in your workflow (e.g., Beverton–Holt) and adjust
  `srmodels`.  
- Reproducibility: keep `set.seed()` fixed while debugging; vary it when
  exploring uncertainty.

# Operational guidance

## First run a few iterations

`Fsquared` can take a long time to run, depending on the number of
cores. While many iterations are ultimately required to derive stable
performance metrics, we suggest to run the script with a small number of
iters (e.g. 50) and a small search grid. If the grid’s length is set to
be the number of cores all of them will run in parallel. Make sure the
script is running and the results seem right; don’t worry too much if
results are a bit off on this phase.

Once simulations run without error, then increase your number of iters
to a minimum of 250 (1000 is better, but be aware of the computing and
memory capacity you have) and the length of the grid to 50, or whichever
value comes closer to a multiple of the number of cores you are using
for efficient parallelization.

## Consider turning off parallelization when running first few iters

Parallelization is the process of distributing computing tasks across
multiple processors (also called cores) on your PC. The number of cores
over which iterations distributed is set up by default in the
`model_fsquared.R` at 5, as such:

``` r
cores <- 5

if(os.linux()) {
  plan(multicore, workers=cores)
} else {
  plan(multisession, workers=cores)
}
```

When first running simulations, you may also turn off parallelization
(i.e. run sequentially) in order to make debugging potential errors
easier. To turn off parallelization, use the following code (e.g. inside
`config.R` or executing directly in your console if you are running code
line by line):

``` r
plan(sequential)
```

## How to specify `bias.correct` argument in `rlnormar1()` function

This argument turns on and off the standard log-normal bias correction
that subtracts $0.5\sigma^2$ from the log of a (normal) mean. You need
to use `bias.correct=FALSE`, which is already the default, if your
stock-recruitment parameters are estimated without the log-normal bias
correction. **The stock-recruitment parameters output from
`bootstrapSR()` are NOT estimated with the log-normal bias correction
(so use `rlnormar1(bias.correct=FALSE)`).**

Additionally, earlier versions of `FLCore` from which `rlnormar1()` is
sourced had the default argument `bias.correct=TRUE`. Thus, you should
always explicitly set the `bias.correct`. You can always check the
default argument for `bias.correct` in your locally installed version of
`FLCore` by typing `View(rlnormar1)` and seeing what is
`bias.correct=...` in the function definition.

Below are additional details for users of the `srrTMB()` function, which
is internally called by `bootstrapSR()`. The function `srrTMB()` also
has a `bias.correct` argument (with default `bias.correct=TRUE`), **but
this is NOT the log-normal bias correction**. This argument turns on a
bias correction that is 1) for a uniform logistic prior (on Blim/B0) and
2) only applied when `model=segreg` and/or `model=segregDa`. So if you
are working directly with `srrTMB()` and use `bias.correct=TRUE` (as
well as `bias.correct=FALSE`), you still call
`rlnormar1(bias.correct=FALSE)` if simulating new recruitment from the
estimated parameters output by `srrTMB()` (e.g. `sigmaR` and `rho`).

## Obtaining or calculating estimation error for shortcut assessment

**- Describe how to extract SB error (CV, but also bias) from SAM**  
**- Describe how to extract SB error (CV, but also bias) from SS3**  
**- Describe how to extract errors in N-at-age (e.g. variance-covariance
matrix)**

## (Sensitivity) Testing multiple OMs

**- Should develop set of OMs, each of which has its own Fp05 and Fmsy
computed**

# Asking for help

**UPDATE**  
To report bugs or make suggestions regarding this example, please use
the **TODO**. If the problem is with code in any one particular package,
please use the corresponding issue page for its repository at the [FLR
github project](https://github.com/flr/).

**LICENCE** This document is licensed under a Creative Commons
Attribution-ShareAlike 4.0 International Public License
