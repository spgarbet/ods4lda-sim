  #############################################################################
 ##
## Perform multiple simulations from a set of known parameters of ods4lda
## Compares random, ACML, WL, and Imputation for int, slp and bivar
##  
.libPaths(c("~/R/rlib-3.4.0", .libPaths()))  # Append ACCRE local directory

library(lme4)
#library(ODS4LDA)
library(reshape)
library(rms)
library(nlme)
library(mitools)
library(MASS)
library(mvtnorm)
library(TwoPhaseReg)
source("functions.R") # instead of library(ODS4LDA)
source("generate-data.R")
source("impute.R")
source("setup.R")
source("two_phase.R")

options(width=200)

  #############################################################################
 ##
## Helper Functions For Run
progress   <- function(...)
{
  cat(date(), ' ')
  lapply(list(...), function(x) cat(x,'\n'))
}

# Coverage of fit
covered <- function(fit, truth)
{
  rng   <- 1:5
  ses   <- sqrt(diag(fit$covariance))[rng]
  lci   <- fit$coefficients[rng] - qnorm(.975)*ses
  uci   <- fit$coefficients[rng] + qnorm(.975)*ses

  (truth[rng] >= lci) & (truth[rng] <= uci)
}

# Coverage of a list of fits
coverage <- function(fits, truth)
{
  result <- t(sapply(fits, function(f) covered(f, truth)))
  colnames(result) <- c("cover.int", "cover.time", "cover.grp", "cover.conf", "cover.tg")
  as.data.frame(result)
}

# Difference of a list of fits from truth
diffs    <- function(fits, truth)
{
  result <- t(sapply(fits, function(f) f$coefficients[1:5] - truth[1:5]))
  colnames(result) <- c("diff.int", "diff.time", "diff.grp", "diff.conf", "diff.tg")
  as.data.frame(result)
}

# The estimates returned from list of fits
ests     <- function(fits)
{
  result <- t(sapply(fits, function(f) f$coefficients[1:5]))
  colnames(result) <- c("est.int", "est.time", "est.grp", "est.conf", "est.tg")
  as.data.frame(result)
}

# Variance of estimates from list of fits
vars     <- function(fits)
{
  result <- t(sapply(fits, function(f) diag(f$covariance)[1:5]))
  colnames(result) <- c("var.int", "var.time", "var.grp", "var.conf", "var.tg")
  as.data.frame(result)
}

  #############################################################################
 ##
## Main Simulation Loop
simulation <- function(run, count)
{
  progress("Generating random data from known parameters")
  conf.param       <- conf.param.sim[count,]
  ni               <- ni.sim[count,]
  inits            <- inits.sim[count,]
  NsPerStratumUniv <- NsPerStratumUniv.sim[count,]
  NsRand           <- sum(NsPerStratumUniv)
  
  progress("Generate cutoffs using a population of 10000")
  dat.tmp   <- GenerateX(N=25000, n=ni[2], prev.grp=prev.grp, c.parm=conf.param)
  dat.tmp$y <- GenerateY(X=cbind(1, dat.tmp$time, dat.tmp$grp, dat.tmp$conf, dat.tmp$time*dat.tmp$grp), Z=cbind(1, dat.tmp$time), id=dat.tmp$id,
                         beta=inits[1:5], sig.b0 = exp(inits[6]), sig.b1 =exp(inits[7]), rho = inits[8], sig.e = exp(inits[9]), RanefDist="Gaussian", ErrorDist="Gaussian")
  cutoffs  <- est.cutoffs(Y=dat.tmp$y, time=dat.tmp$time, id=dat.tmp$id, PropInCentralRegion=p.central)

  progress("Generate random data from truth")
  dat      <- GenerateX(N=N, n=ni[2], prev.grp=prev.grp, c.parm=conf.param)
  dat$y    <- GenerateY(X=cbind(1, dat$time, dat$grp, dat$conf, dat$time*dat$grp), Z=cbind(1, dat$time), id=dat$id,
                        beta=inits[1:5], sig.b0 = exp(inits[6]), sig.b1 =exp(inits[7]), rho = inits[8], sig.e = exp(inits[9]), RanefDist="Gaussian", ErrorDist="Gaussian")
  dat      <- dat[order(dat$id, dat$time),]
  droptime <- rep(sample(seq(ni[1], ni[2]), N, replace=TRUE), each=ni[2]) - 1
  dat      <- dat[dat$time<=droptime,]
  
  IntSlps <- CalcSSIntSlp( Y=dat$y, time=dat$time, id=dat$id)
  dat$Int <- IntSlps[[1]]
  dat$Slp <- IntSlps[[2]]

  ## identify stratum membership
  dat$StratInt <- identify.stratum(Y=dat$y, time=dat$time, id=dat$id, w.function="intercept", cutpoints=cutoffs$IntCutUniv, Int=dat$Int, Slp=dat$Slp)
  dat$StratSlp <- identify.stratum(Y=dat$y, time=dat$time, id=dat$id, w.function="slope",     cutpoints=cutoffs$SlpCutUniv, Int=dat$Int, Slp=dat$Slp)
 
  ## identify those sampled along with individual sampling probs and stratum sampling probs
  SampledRan <- random.sampling(id.long=dat$id, n=NsRand)
  SampledInt <- ods.sampling(id.long=dat$id, stratum.long=dat$StratInt, SamplingStrategy="IndepODS", NsPerStratum=NsPerStratumUniv)
  SampledSlp <- ods.sampling(id.long=dat$id, stratum.long=dat$StratSlp, SamplingStrategy="IndepODS", NsPerStratum=NsPerStratumUniv)
  
  ## add who will be sampled for each design to dat
  dat$SampledRan <- SampledRan
  dat$SampledInt <- SampledInt[[1]]
  dat$SampledSlp <- SampledSlp[[1]]
  dat$SampledMix1 <- dat$SampledInt*(dat$id <= N/3) + dat$SampledSlp*(dat$id > N/3)
  dat$SampledMix2 <- dat$SampledInt*(dat$id <= 2*N/3) + dat$SampledSlp*(dat$id > 2*N/3)
  
  ## Added subject specific sampling weights for each design if doing a weighted likelihood analysis
  dat$WeightsInt <- 1/SampledInt[[2]]
  dat$WeightsSlp <- 1/SampledSlp[[2]]
  dat$WeightsMix1 <- 1/(SampledInt[[2]]*(dat$id <= N/3) + SampledSlp[[2]]*(dat$id > N/3))
  dat$WeightsMix2 <- 1/(SampledInt[[2]]*(dat$id <= 2*N/3) + SampledSlp[[2]]*(dat$id > 2*N/3))
  
  dat$NoWeighting <-1
  
  dat$IntProbLow  <- SampledInt[[3]][1]
  dat$IntProbMid  <- SampledInt[[3]][2]
  dat$IntProbHigh <- SampledInt[[3]][3]
  dat$SlpProbLow  <- SampledSlp[[3]][1]
  dat$SlpProbMid  <- SampledSlp[[3]][2]
  dat$SlpProbHigh <- SampledSlp[[3]][3]
  dat$Mix1ProbLow  <- SampledInt[[3]][1]*(dat$id <= N/3) + SampledSlp[[3]][1]*(dat$id > N/3)
  dat$Mix1ProbMid  <- SampledInt[[3]][2]*(dat$id <= N/3) + SampledSlp[[3]][2]*(dat$id > N/3)
  dat$Mix1ProbHigh <- SampledInt[[3]][3]*(dat$id <= N/3) + SampledSlp[[3]][3]*(dat$id > N/3)
  dat$Mix2ProbLow  <- SampledInt[[3]][1]*(dat$id <= 2*N/3) + SampledSlp[[3]][1]*(dat$id > 2*N/3)
  dat$Mix2ProbMid  <- SampledInt[[3]][2]*(dat$id <= 2*N/3) + SampledSlp[[3]][2]*(dat$id > 2*N/3)
  dat$Mix2ProbHigh <- SampledInt[[3]][3]*(dat$id <= 2*N/3) + SampledSlp[[3]][3]*(dat$id > 2*N/3)
  
  dat$IntCutoff1 <- cutoffs$IntCutUniv[1]
  dat$IntCutoff2 <- cutoffs$IntCutUniv[2]
  dat$SlpCutoff1 <- cutoffs$SlpCutUniv[1]
  dat$SlpCutoff2 <- cutoffs$SlpCutUniv[2]
  dat$Mix1Cutoff1 <- dat$IntCutoff1*(dat$id <= N/3) + dat$SlpCutoff1*(dat$id > N/3)
  dat$Mix1Cutoff2 <- dat$IntCutoff2*(dat$id <= N/3) + dat$SlpCutoff2*(dat$id > N/3)
  dat$Mix2Cutoff1 <- dat$IntCutoff1*(dat$id <= 2*N/3) + dat$SlpCutoff1*(dat$id > 2*N/3)
  dat$Mix2Cutoff2 <- dat$IntCutoff2*(dat$id <= 2*N/3) + dat$SlpCutoff2*(dat$id > 2*N/3)
  
  dat$Int.w <- "intercept"
  dat$Slp.w <- "slope"
  dat$Mix1.w <- ifelse(dat$id <= N/3, dat$Int.w, dat$Slp.w)
  dat$Mix2.w <- ifelse(dat$id <= 2*N/3, dat$Int.w, dat$Slp.w)
  
  ## Stratum Sampling probabilities
  SampProbRan <- c(1,1,1)
  SampProbInt <- SampledInt[[3]]
  SampProbSlp <- SampledSlp[[3]]
  
  ## Datasets for sampled subjects
  datRan <- dat[dat$SampledRan==1,]
  datInt <- dat[dat$SampledInt==1,]
  datSlp <- dat[dat$SampledSlp==1,]
  datMix1 <- dat[dat$SampledMix1==1,]
  datMix2 <- dat[dat$SampledMix2==1,]
  
  cutpointsRan=cbind(datRan$IntCutoff1, datRan$IntCutoff2)  ## just need these numbers for the function, not used
  SampProbRan=matrix(1, ncol=3, nrow=length(datRan[,1])) 
  w.functionRan=datRan$Int.w                                ## just need these numbers for the function, not used
  
  cutpointsInt=cbind(datInt$IntCutoff1, datInt$IntCutoff2)
  SampProbInt=cbind(datInt$IntProbLow, datInt$IntProbMid, datInt$IntProbHigh)  
  w.functionInt=datInt$Int.w
  
  cutpointsSlp=cbind(datSlp$SlpCutoff1, datSlp$SlpCutoff2)
  SampProbSlp=cbind(datSlp$SlpProbLow, datSlp$SlpProbMid, datSlp$SlpProbHigh)  
  w.functionSlp=datSlp$Slp.w
  
  cutpointsMix1=cbind(datMix1$Mix1Cutoff1, datMix1$Mix1Cutoff2)
  SampProbMix1=cbind(datMix1$Mix1ProbLow, datMix1$Mix1ProbMid, datMix1$Mix1ProbHigh)  
  w.functionMix1=datMix1$Mix1.w
  
  cutpointsMix2=cbind(datMix2$Mix2Cutoff1, datMix2$Mix2Cutoff2)
  SampProbMix2=cbind(datMix2$Mix2ProbLow, datMix2$Mix2ProbMid, datMix2$Mix2ProbHigh)  
  w.functionMix2=datMix2$Mix2.w
  
  ## Datasets for unsampled subjects
  datNotRan  <- dat[dat$SampledRan==0,]
  datNotInt  <- dat[dat$SampledInt==0,]
  datNotSlp  <- dat[dat$SampledSlp==0,]
  datNotMix1 <- dat[dat$SampledMix1==0,]
  datNotMix2 <- dat[dat$SampledMix2==0,]
  
  cutpointsNotRan=cbind(datNotRan$IntCutoff1, datNotRan$IntCutoff2)  ## just need these numbers for the function, not used
  SampProbNotRan=matrix(1, ncol=3, nrow=length(datNotRan[,1])) 
  w.functionNotRan=datNotRan$Int.w                                ## just need these numbers for the function, not used
  
  cutpointsNotInt=cbind(datNotInt$IntCutoff1, datNotInt$IntCutoff2)
  SampProbNotInt=cbind(datNotInt$IntProbLow, datNotInt$IntProbMid, datNotInt$IntProbHigh)  
  w.functionNotInt=datNotInt$Int.w
  
  cutpointsNotSlp=cbind(datNotSlp$SlpCutoff1, datNotSlp$SlpCutoff2)
  SampProbNotSlp=cbind(datNotSlp$SlpProbLow, datNotSlp$SlpProbMid, datNotSlp$SlpProbHigh)  
  w.functionNotSlp=datNotSlp$Slp.w
  
  cutpointsNotMix1=cbind(datNotMix1$Mix1Cutoff1, datNotMix1$Mix1Cutoff2)
  SampProbNotMix1=cbind(datNotMix1$Mix1ProbLow, datNotMix1$Mix1ProbMid, datNotMix1$Mix1ProbHigh)  
  w.functionNotMix1=datNotMix1$Mix1.w
  
  cutpointsNotMix2=cbind(datNotMix2$Mix2Cutoff1, datNotMix2$Mix2Cutoff2)
  SampProbNotMix2=cbind(datNotMix2$Mix2ProbLow, datNotMix2$Mix2ProbMid, datNotMix2$Mix2ProbHigh)  
  w.functionNotMix2=datNotMix2$Mix2.w
    
    ###########################################################################
   ##
  progress("Model Fitting Begins")

  progress("Random Sampling")
  Fit.ran     <- acml.lmem2(formula.fixed=y~time*grp+conf, formula.random= ~time, id=id, data=datRan, InitVals=inits, ProfileCol=NA,
                            cutpoints=cutpointsRan, SampProb=SampProbRan, Weights=NoWeighting, w.function=w.functionRan)

    ###########################################################################
   ##
  progress("ACML")
  
  progress("... Intercept")
  Fit.int     <- acml.lmem2(formula.fixed=y~time*grp+conf, formula.random= ~time, id=id, data=datInt, InitVals=inits, ProfileCol=NA,
                            cutpoints=cutpointsInt, SampProb=SampProbInt, Weights=NoWeighting, w.function=w.functionInt)

  progress("... Slope")
  Fit.slp     <- acml.lmem2(formula.fixed=y~time*grp+conf, formula.random= ~time, id=id, data=datSlp, InitVals=inits, ProfileCol=NA,
                            cutpoints=cutpointsSlp, SampProb=SampProbSlp, Weights=NoWeighting, w.function=w.functionSlp)

  progress("... Mixture 1")  
  Fit.mix1    <- acml.lmem2(formula.fixed=y~time*grp+conf, formula.random= ~time, id=id, data=datMix1, InitVals=inits, ProfileCol=NA,
                            cutpoints=cutpointsMix1, SampProb=SampProbMix1, Weights=NoWeighting, w.function=w.functionMix1)
  
  progress("... Mixture 2")  
  Fit.mix2    <- acml.lmem2(formula.fixed=y~time*grp+conf, formula.random= ~time, id=id, data=datMix2, InitVals=inits, ProfileCol=NA,
                            cutpoints=cutpointsMix2, SampProb=SampProbMix2, Weights=NoWeighting, w.function=w.functionMix2)

      ###########################################################################
   ##
  progress("Weighted Likelihood")

  progress("... Intercept")
  Fit.int.wl  <- acml.lmem2(formula.fixed=y~time*grp+conf, formula.random= ~time, id=id, data=datInt, InitVals=inits, ProfileCol=NA,
                            cutpoints=cutpointsInt, SampProb=matrix(1, nrow=nrow(datInt), ncol=3), Weights=WeightsInt, w.function=w.functionInt)

  progress("... Slope")
  Fit.slp.wl  <- acml.lmem2(formula.fixed=y~time*grp+conf, formula.random= ~time, id=id, data=datSlp, InitVals=inits, ProfileCol=NA,
                            cutpoints=cutpointsSlp, SampProb=matrix(1, nrow=nrow(datSlp), ncol=3), Weights=WeightsSlp, w.function=w.functionSlp)
  
  progress("... Mixture 1")  
  Fit.mix1.wl <- acml.lmem2(formula.fixed=y~time*grp+conf, formula.random= ~time, id=id, data=datMix1, InitVals=inits, ProfileCol=NA,
                            cutpoints=cutpointsMix1, SampProb=matrix(1, nrow=nrow(datMix1), ncol=3), Weights=WeightsMix1, w.function=w.functionMix1)

  progress("... Mixture 2")  
  Fit.mix2.wl <- acml.lmem2(formula.fixed=y~time*grp+conf, formula.random= ~time, id=id, data=datMix2, InitVals=inits, ProfileCol=NA,
                            cutpoints=cutpointsMix2, SampProb=matrix(1, nrow=nrow(datMix2), ncol=3), Weights=WeightsMix2, w.function=w.functionMix2)
  
    ###########################################################################
   ##
  progress("Indirect Imputation Analyses")
  
  progress("... Random")
  Fit.ran.mi <- IndirectImputation(acml.fit=Fit.ran, datSampled=datRan, datNotSampled=datNotRan, 
                                   cutpointsNotSampled=cutpointsNotRan, w.functionNotSampled=w.functionNotRan, SampProbNotSampled=SampProbNotRan, n.imp=n.imp)
  progress("... Intercept")
  Fit.int.mi <- IndirectImputation(acml.fit=Fit.int, datSampled=datInt, datNotSampled=datNotInt, 
                                   cutpointsNotSampled=cutpointsNotInt, w.functionNotSampled=w.functionNotInt, SampProbNotSampled=SampProbNotInt, n.imp=n.imp)
  progress("... Slope")
  Fit.slp.mi <- IndirectImputation(acml.fit=Fit.slp, datSampled=datSlp, datNotSampled=datNotSlp, 
                                   cutpointsNotSampled=cutpointsNotSlp, w.functionNotSampled=w.functionNotSlp, SampProbNotSampled=SampProbNotSlp, n.imp=n.imp)
  progress("... Mixture 1")
  Fit.mix1.mi <- IndirectImputation(acml.fit=Fit.mix1, datSampled=datMix1, datNotSampled=datNotMix1, 
                                    cutpointsNotSampled=cutpointsNotMix1, w.functionNotSampled=w.functionNotMix1, SampProbNotSampled=SampProbNotMix1, n.imp=n.imp) 
  progress("... Mixture 2")
  Fit.mix2.mi <- IndirectImputation(acml.fit=Fit.mix2, datSampled=datMix2, datNotSampled=datNotMix2, 
                                   cutpointsNotSampled=cutpointsNotMix2, w.functionNotSampled=w.functionNotMix2, SampProbNotSampled=SampProbNotMix2, n.imp=n.imp) 

    ###########################################################################
   ##
  progress("Two Phase Reg")
  
  progress("... Random")
  Fit.ran.2p <- two_phase(datRan, datNotRan)

  progress("... Intercept")
  Fit.int.2p <- two_phase(datInt, datNotInt)

  progress("... Slope")
  Fit.slp.2p <- two_phase(datSlp, datNotSlp)

  progress("... Mixture 1")
  Fit.mix1.2p <- two_phase(datMix1, datNotMix1)
  
  progress("... Mixture 2")
  Fit.mix2.2p <- two_phase(datMix2, datNotMix2)

    ###########################################################################
   ##
  progress("Saving results.")
  
  fits <- list(
    Fit.ran,      Fit.int,      Fit.slp,      Fit.mix1,     Fit.mix2,
                  Fit.int.wl,   Fit.slp.wl,   Fit.mix1.wl,  Fit.mix2.wl,
    Fit.ran.mi,   Fit.int.mi,   Fit.slp.mi,   Fit.mix1.mi,  Fit.mix2.mi,
    Fit.ran.2p,   Fit.int.2p,   Fit.slp.2p,   Fit.mix1.2p,  Fit.mix2.2p
  )
  
  results <- data.frame(
    Job      = rep(run,   19),
    Scenario = rep(count, 19),
    Sampling = factor(c(
      "Random", "Intercept", "Slope", "Mix1", "Mix2",
                "Intercept", "Slope", "Mix1", "Mix2",
      "Random", "Intercept", "Slope", "Mix1", "Mix2",
      "Random", "Intercept", "Slope", "Mix1", "Mix2"
    )),
    Method   = factor(c(
      rep("ACML",     5),
      rep("WL",       4),
      rep("MI",       5),
      rep("TwoPhase", 5)
    ))
  )
  
  results <- cbind(results,
                   coverage(fits, inits),
                   diffs(fits, inits),
                   ests(fits),
                   vars(fits))
  
  save( results,
        file=paste0("output/run-", run, "-", count, ".RData")
      )
  
  progress(paste0("Run ", count, " complete."))
}