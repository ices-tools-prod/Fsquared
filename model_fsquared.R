# model_fsquared.R - runs Fsquared search for Ftarget
# Fsquared/model_fsquared.R

#===============================================================================
# Fsquared - estimation of F reference points for the ICES advice rule
# Authors: participants of ICES MGWG workshop on the validation of new tool for refpts
#	  Iago Mosqueira (WMR) <iago.mosqueira@wur.nl>
#	  Ernesto Jardim (IPMA) <ernesto.jardim@ipma.pt>
#	  John Trochta (IMR) <john.tyler.trochta@hi.no>
#	  Arni Magnusson (SPC) <arnim@spc.int>
#	  Max Cardinale (SLU) <massimiliano.cardinale@slu.se>
#	  Dorleta Garcia (ICES) <dorleta.garcia@ices.dk>
#	  Colin Millar (ICES) <colin.millar@ices.dk>
#
# Distributed under the terms of the EUPL 1.2
#===============================================================================

# load libraries

library(TAF)
library(mse)
library(msemodules)
library(FLasher)
library(FLSRTMB)  # to estimate SR parameters and add estimation uncertainty

source("utilities_fsquared.R")

mkdir("model")

# LOAD stock assessment results, 'run' is output FLStock
# if you have your own assessment load it, for simplicity rename the output
# to "run", but can use other name, just make sure all calls to "run" 
# are updated
load("model/model.rda")

#===============================================================================
# SETUP
#===============================================================================

# Name of the stock
stkname <- name(run)
# TODO: Recruitment models to be used in the OM conditioning
srmodels <- c("segreg") # segreg, bevholt, ricker
# Initial year of projections
iy <- dims(run)$maxyear
# Years to be used to compute SPR0 for stock-recruitment model, last 5
spryrs <- seq(iy - 5, iy)
# Data year
dy <- iy - 1
# Final year
fy <- iy + 50
# Probability years
pys <- seq(fy - 5, fy)
# How many years from the past to condition the future
conditioning_ny <- 5
# CV for SSB to add uncertainty in the shortcut estimator
bcv_sa <- 0.1
# CV for F to add uncertainty in the shortcut estimator
fcv_sa <- 0.1
# Years for geometric mean in short term forecast
recyrs_mp <- -2
# TODO: Blim and Btrigger
Blim <- 788
Btrigger <- 1095
refpts <- FLPar(c(Blim = Blim, Btrigger = Btrigger))
# TODO: no. of cores to use in parallel, defauls to 2/3 of those in machine
cores <- round(availableCores() * 0.6)
# TODO: F search grid
fg_mp <- seq(0, 1.5, length=cores)
# Number of iterations (minimum of 50 for testing, 500 for final)
it <- 25
# it <- max(25, cores * 25)
# Random seed
set.seed(987)

# PARALLEL setup via doFuture
if(os.linux()) {
  plan(multicore, workers=cores)
} else {
  plan(multisession, workers=cores)
}

options(doFuture.rng.onMisuse="ignore")

#===============================================================================
# OM conditioning
#===============================================================================

# Stock-recruitment relationship(s)
# The file utilities_fsquared.R has code examples to condition the OM
# using stock recruitment parameters estimated by the stock assessment
# model, like SS3, SAM and a4a.

# BOOTSTRAP and SELECT model by largest logLik
srpars <- bootstrapSR(run, iters=it, spr0=mean(spr0y(run)[, ac(spryrs)]),
  models=srmodels)

# GENERATE future deviances: lognormal autocorrelated
srdevs <- rlnormar1(sdlog=srpars$sigmaR, rho=srpars$rho, years=seq(dy, fy),
  bias.correct=FALSE)

# BUILD FLom, OM FLR object
om <- FLom(stock=propagate(run, it), refpts=refpts, model="segreg",
  params=srpars, deviances=srdevs, name=stkname)

# TODO: SETUP om future: average of most recent years set by conditioning_ny
om <- fwdWindow(om, end=fy, nsq=conditioning_ny)

#===============================================================================
# diagnostic(s)
#===============================================================================

# TEST Blim with F=0 projection
f0 <- fwd(om, control=fwdControl(quant='fbar', value=0, year=seq(iy + 1, fy)))

# COMPUTE P(SB<Blim) in time
# when F=0 this probability should be 0, otherwise is a sign BLim is not 
# coherent with the stock-recruitment model and parameters 
performance(f0, statistics=icestats["PBlim"])[year %in% seq(iy, fy, by=5),
  .(PBlim=mean(data)), by=year]

# COMPUTE Inter-annual variability of biomass
# to check when biomass stabilizes and set number of years for projections
# leave about 10 years of stable biomass to compute metrics
# review MSE setup section to avoid projecting for longer than necessary
# rerun setup and objects if you changed the final year of projections
dt0 <- performance(f0, statistics=icestats["IACB"])[year %in% iy:fy,
  .(IACB=mean(data)), by=year]

plot(dt0, type="l", main="Inter-annual changes in biomass")

#===============================================================================
# MP
#===============================================================================

# SET intermediate year + start of runs, lags and frequency
mseargs <- list(iy=iy, fy=fy, data_lag=1, management_lag=1, frq=1) #link to setup

# SET shortcut estimator uncertainty: F and SSB deviances and auto-correlation
# Note your SSB deviances and auto-correlation have very little impact on the P(SB<Blim)
sdevs <- shortcut_devs(om, SSBcv=bcv_sa, Fcv=fcv_sa, Fphi=0)

# SETUP standard ICES advice rule
arule <- mpCtrl(

  # (est)imation method: shortcut.sa + SSB deviances
  est = mseCtrl(method=shortcut.sa,
    args=list(SSBdevs=sdevs$SSB)),

  # hcr: hockeystick (fbar ~ ssb | lim, trigger, target, min)
  hcr = mseCtrl(method=hockeystick.hcr,
    args=list(lim=0, trigger=refpts(om)$Btrigger, target=0.22, min=0,
    metric="ssb", output="fbar")),

  # (i)mplementation (sys)tem: tac.is (C ~ F)
  isys = mseCtrl(method=tac.is, args=list(recyrs=recyrs_mp, Fdevs=sdevs$F))
)

#===============================================================================
# Run simulations
#===============================================================================

# RUN over Ftarget grid
fgrid <- mps(om, ctrl=arule, args=mseargs, hcr=list(target=fg_mp),
  names=paste0("F", fg_mp))

# PLOT
plot(om, fgrid)

# COMPUTE average performance over pys
performance(fgrid) <- performance(fgrid, statistics=icestats["PBlim"], year=pys,
  type="arule")

# OR ... RUN over Ftarget grid and return only performance stats
# fgrid <- mps(om, ctrl=arule, args=mseargs, hcr=list(target=fg_mp),
#   names=paste("F" fg_mp), statistics=icestats, type="arule")

# FIND Ftarget that gives mean P(B < Blim) = 5%
tune <- tunebisect(om, control=arule, args=mseargs,
  tune=list(target=0.3 * c(0.5, 1.5)),
  statistic=icestats["PBlim"], prob=0.05, tol=0.005, years=pys)

# PLOT
plot(om, tune)

# COMPUTE performance
performance(tune) <- performance(tune, statistics=icestats, type="arule", run="tune")

# CHECK Ftarget value
args(control(tune, "hcr"))$target

# SAVE
save("fgrid", "om", "tune", file="model/fsquared.rda")

