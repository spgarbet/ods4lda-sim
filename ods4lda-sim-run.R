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
source("ods4lda-library.R") # instead of library(ODS4LDA)
source("ods4lda-sim-fun.R")
source("ods4lda-sim-impute.R")
source("ods4lda-sim-setup.R")

options(width=200)

  #############################################################################
 ##
## Progress helper
progress   <- function(...)
{
  cat(date(), ' ')
  lapply(list(...), function(x) cat(x,'\n'))
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
  NsRand           <- NsRand.sim[count]
  
  progress("Generate cutoffs using a population of 10000")
  dat.tmp   <- GenerateX(N=20000, n=ni[2], prev.grp=prev.grp, c.parm=conf.param)
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
  #dat$StratBiv <- identify.stratum(Y=dat$y, time=dat$time, id=dat$id, w.function="bivar",     cutpoints=c(cutoffs$IntCutBiv,cutoffs$SlpCutBiv), Int=dat$Int, Slp=dat$Slp)

  ## identify those sampled along with individual sampling probs and stratum sampling probs
  SampledInt <- ods.sampling(id.long=dat$id, stratum.long=dat$StratInt, SamplingStrategy="IndepODS", NsPerStratum=NsPerStratumUniv)
  SampledSlp <- ods.sampling(id.long=dat$id, stratum.long=dat$StratSlp, SamplingStrategy="IndepODS", NsPerStratum=NsPerStratumUniv)
  #SampledBiv <- ods.sampling(id.long=dat$id, stratum.long=dat$StratBiv, SamplingStrategy="IndepODS", NsPerStratum=NsPerStratumBiv)
  SampledRan <- random.sampling(id.long=dat$id, n=NsRand)

  ## add who will be sampled for each design to dat
  dat$SampledInt <- SampledInt[[1]]
  dat$SampledSlp <- SampledSlp[[1]]
  #dat$SampledBiv <- SampledBiv[[1]]
  dat$SampledRan <- SampledRan
  dat$SampledMix <- dat$SampledInt*(dat$id <= N/2) + dat$SampledSlp*(dat$id > N/2)
  
  ## Added subject specific sampling probabilities for each design if doing a weighted likelihood analysis
  dat$SampProbiIntWL <- SampledInt[[2]]
  dat$SampProbiSlpWL <- SampledSlp[[2]]
  #dat$SampProbiBivWL <- SampledBiv[[2]]
  dat$SampProbiRanWL <- rep(1, length(dat[,1]))
  dat$SampProbiMixWL <- SampledInt[[2]]*(dat$id <= N/2) + SampledSlp[[2]]*(dat$id > N/2)
  
  dat$SampProbiRan <- dat$SampProbiBiv <- dat$SampProbiSlp <- dat$SampProbiInt <- dat$SampProbiMix <-rep(1, length(dat[,1]))
  
  dat$IntProbLow  <- SampledInt[[3]][1]
  dat$IntProbMid  <- SampledInt[[3]][2]
  dat$IntProbHigh <- SampledInt[[3]][3]
  dat$SlpProbLow  <- SampledSlp[[3]][1]
  dat$SlpProbMid  <- SampledSlp[[3]][2]
  dat$SlpProbHigh <- SampledSlp[[3]][3]
  dat$MixProbLow  <- SampledInt[[3]][1]*(dat$id <= N/2) + SampledSlp[[3]][1]*(dat$id > N/2)
  dat$MixProbMid  <- SampledInt[[3]][2]*(dat$id <= N/2) + SampledSlp[[3]][2]*(dat$id > N/2)
  dat$MixProbHigh <- SampledInt[[3]][3]*(dat$id <= N/2) + SampledSlp[[3]][3]*(dat$id > N/2)
  
  dat$IntCutoff1 <- cutoffs$IntCutUniv[1]
  dat$IntCutoff2 <- cutoffs$IntCutUniv[2]
  dat$SlpCutoff1 <- cutoffs$SlpCutUniv[1]
  dat$SlpCutoff2 <- cutoffs$SlpCutUniv[2]
  
  dat$MixCutoff1 <- dat$IntCutoff1*(dat$id <= N/2) + dat$SlpCutoff1*(dat$id > N/2)
  dat$MixCutoff2 <- dat$IntCutoff2*(dat$id <= N/2) + dat$SlpCutoff2*(dat$id > N/2)
  
  dat$Int.w <- "intercept"
  dat$Slp.w <- "slope"
  dat$Mix.w <- ifelse(dat$id <= N/2, dat$Int.w, dat$Slp.w)
  
  ## Stratum Sampling probabilities
  SampProbInt <- SampledInt[[3]]
  SampProbSlp <- SampledSlp[[3]]
  #SampProbBiv <- SampledBiv[[3]]
  SampProbRan <- c(1,1,1)

  ## Datasets for sampled subjects
  datInt <- dat[dat$SampledInt==1,]
  datSlp <- dat[dat$SampledSlp==1,]
  #datBiv <- dat[dat$SampledBiv==1,]
  datRan <- dat[dat$SampledRan==1,]
  datMix <- dat[dat$SampledMix==1,]

  ## Datasets for unsampled subjects
  datNotInt <- dat[dat$SampledInt==0,]
  datNotSlp <- dat[dat$SampledSlp==0,]
  #datNotBiv <- dat[dat$SampledBiv==0,]
  datNotRan <- dat[dat$SampledRan==0,]
  datNotMix <- dat[dat$SampledMix==0,]
  
  cutpointsRan.new=cbind(datRan$IntCutoff1, datRan$IntCutoff2)  ## just need these numbers for the function, not used
  SampProbRan.new=matrix(1, ncol=3, nrow=length(datRan[,1])) 
  w.functionRan.new=datRan$Int.w                                ## just need these numbers for the function, not used
  
  cutpointsInt.new=cbind(datInt$IntCutoff1, datInt$IntCutoff2)
  SampProbInt.new=cbind(datInt$IntProbLow, datInt$IntProbMid, datInt$IntProbHigh)  
  w.functionInt.new=datInt$Int.w
  
  cutpointsSlp.new=cbind(datSlp$SlpCutoff1, datSlp$SlpCutoff2)
  SampProbSlp.new=cbind(datSlp$SlpProbLow, datSlp$SlpProbMid, datSlp$SlpProbHigh)  
  w.functionSlp.new=datSlp$Slp.w
  
  cutpointsMix.new=cbind(datMix$MixCutoff1, datMix$MixCutoff2)
  SampProbMix.new=cbind(datMix$MixProbLow, datMix$MixProbMid, datMix$MixProbHigh)  
  w.functionMix.new=datMix$Mix.w
    
    ###########################################################################
   ##
  progress("Model Fitting Begins")
  
  progress("Random Sampling")
  Fit.ran     <- acml.lmem2(formula.fixed=y~time*grp+conf, formula.random= ~time, id=id, data=datRan, InitVals=inits, ProfileCol=NA,
                            cutpoints=cutpointsRan.new, SampProb=SampProbRan.new, SampProbiWL=SampProbiRan, w.function=w.functionRan.new)
  
  progress("Univariate intercept ACML")
  Fit.int     <- acml.lmem2(formula.fixed=y~time*grp+conf, formula.random= ~time, id=id, data=datInt, InitVals=inits, ProfileCol=NA,
                            cutpoints=cutpointsInt.new, SampProb=SampProbInt.new, SampProbiWL=SampProbiInt, w.function=w.functionInt.new)

  progress("Univariate slope ACML")
  Fit.slp     <- acml.lmem2(formula.fixed=y~time*grp+conf, formula.random= ~time, id=id, data=datSlp, InitVals=inits, ProfileCol=NA,
                            cutpoints=cutpointsSlp.new, SampProb=SampProbSlp.new, SampProbiWL=SampProbiSlp, w.function=w.functionSlp.new)

  progress("Mixture ACML")  
  Fit.mix     <- acml.lmem2(formula.fixed=y~time*grp+conf, formula.random= ~time, id=id, data=datMix, InitVals=inits, ProfileCol=NA,
                            cutpoints=cutpointsMix.new, SampProb=SampProbMix.new, SampProbiWL=SampProbiMix, w.function=w.functionMix.new)

    ###########################################################################
   ##
  progress("Indirect Imputation Analyses", "... ran.mi")
  Fit.ran.mi <- IndirectImputation(acml.fit=Fit.ran, datSampled=datRan, datNotSampled=datNotRan, n.imp=n.imp)
  progress("... int.mi")
  Fit.int.mi <- IndirectImputation(acml.fit=Fit.int, datSampled=datInt, datNotSampled=datNotInt, n.imp=n.imp)
  progress("... slp.mi")
  Fit.slp.mi <- IndirectImputation(acml.fit=Fit.slp, datSampled=datSlp, datNotSampled=datNotSlp, n.imp=n.imp)
  progress("... mix.mi")
  Fit.mix.mi <- IndirectImputation(acml.fit=Fit.mix, datSampled=datMix, datNotSampled=datNotMix, n.imp=n.imp) 

  progress(paste0("Run ", count, " complete. Saving data."))
  save( Fit.ran,      Fit.int,      Fit.slp,      Fit.mix,
        Fit.ran.mi,   Fit.int.mi,   Fit.slp.mi,   Fit.mix.mi,
        file=paste0("output/run-", run, "-", count, ".RData")
      ) 
}

  #############################################################################
 ##
## Use ods4lda-sim-accre.R to run on ACCRE
## Use ods4lda-sim-local.R to run locally
