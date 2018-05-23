  #############################################################################
 ##
## Perform multiple simulations from a set of known parameters of ods4lda
## Compares random, ACML, WL, and Imputation for int, slp and bivar
##  
.libPaths(c("~/R/rlib-3.4.0", .libPaths()))  # Append ACCRE local directory

library(lme4)
library(ODS4LDA)
library(reshape)
library(rms)
library(nlme)
library(mitools)
library(MASS)

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
  ni               <- ni.sim[count,]
  inits            <- inits.sim[count,]
  NsPerStratumUniv <- NsPerStratumUniv.sim[count,]
  NsPerStratumBiv  <- NsPerStratumBiv.sim[count,]
  NsRand           <- NsRand.sim[count]
  
  progress("Generate cutoffs using a population of 10000")
  dat.tmp   <- GenerateX(N=10000, n=ni[2], prev.grp=prev.grp, c.parm=conf.param)
  dat.tmp$y <- GenerateY(X=cbind(1, dat.tmp$time, dat.tmp$grp, dat.tmp$conf, dat.tmp$time*dat.tmp$grp), Z=cbind(1, dat.tmp$time), id=dat.tmp$id,
                         beta=inits[1:5], sig.b0 = exp(inits[6]), sig.b1 =exp(inits[7]), rho = inits[8], sig.e = exp(inits[9]), RanefDist="Gaussian", ErrorDist="Gaussian")
  cutoffs   <- est.cutoffs(Y=dat.tmp$y, time=dat.tmp$time, id=dat.tmp$id, PropInCentralRegion=p.central)

  progress("Generate random data from truth")
  dat      <- GenerateX(N=N, n=ni[2], prev.grp=prev.grp, c.parm=conf.param)
  dat$y    <- GenerateY(X=cbind(1, dat$time, dat$grp, dat$conf, dat$time*dat$grp), Z=cbind(1, dat$time), id=dat$id,
                        beta=inits[1:5], sig.b0 = exp(inits[6]), sig.b1 =exp(inits[7]), rho = inits[8], sig.e = exp(inits[9]), RanefDist="Gaussian", ErrorDist="Gaussian")
  dat      <- dat[order(dat$id, dat$time),]
  droptime <- rep(sample(seq(ni[1], ni[2]), N, replace=TRUE), each=ni[2]) - 1
  dat      <- dat[dat$time<=droptime,]
  
  ## for the imputation analysis
  dat$ymean  <- cluster.summary(dat$id, dat$y, mean)
  dat$ytmean <- cluster.summary(dat$id, dat$y*dat$time, mean)
  dat$tmean  <- cluster.summary(dat$id, dat$time, mean)
  dat$t2mean <- cluster.summary(dat$id, dat$time^2, mean)
  
  ## Calculate subject specific intercepts and slopes
  IntSlps <- CalcSSIntSlp( Y=dat$y, time=dat$time, id=dat$id)
  dat$Int <- IntSlps[[1]]
  dat$Slp <- IntSlps[[2]]
  
  ## identify stratum membership
  dat$StratInt <- identify.stratum(Y=dat$y, time=dat$time, id=dat$id, w.function="intercept", cutpoints=cutoffs$IntCutUniv, Int=dat$Int, Slp=dat$Slp)
  dat$StratSlp <- identify.stratum(Y=dat$y, time=dat$time, id=dat$id, w.function="slope",     cutpoints=cutoffs$SlpCutUniv, Int=dat$Int, Slp=dat$Slp)
  dat$StratBiv <- identify.stratum(Y=dat$y, time=dat$time, id=dat$id, w.function="bivar",     cutpoints=c(cutoffs$IntCutBiv,cutoffs$SlpCutBiv), Int=dat$Int, Slp=dat$Slp)
  
  ## identify those sampled along with individual sampling probs and stratum sampling probs
  SampledInt <- ods.sampling(id.long=dat$id, stratum.long=dat$StratInt, SamplingStrategy="IndepODS", NsPerStratum=NsPerStratumUniv)
  SampledSlp <- ods.sampling(id.long=dat$id, stratum.long=dat$StratSlp, SamplingStrategy="IndepODS", NsPerStratum=NsPerStratumUniv)
  SampledBiv <- ods.sampling(id.long=dat$id, stratum.long=dat$StratBiv, SamplingStrategy="IndepODS", NsPerStratum=NsPerStratumBiv)
  SampledRan <- random.sampling(id.long=dat$id, n=NsRand)
  
  ## add who will be sampled for each design to dat
  dat$SampledInt <- SampledInt[[1]]
  dat$SampledSlp <- SampledSlp[[1]]
  dat$SampledBiv <- SampledBiv[[1]]
  dat$SampledRan <- SampledRan
  
  ## Added subject specific sampling probabilities for each design if doing a weighted likelihood analysis
  dat$SampProbiIntWL <- SampledInt[[2]]
  dat$SampProbiSlpWL <- SampledSlp[[2]]
  dat$SampProbiBivWL <- SampledBiv[[2]]
  dat$SampProbiRanWL <- rep(1, length(dat[,1]))
  
  dat$SampProbiRan <- dat$SampProbiBiv <- dat$SampProbiSlp <- dat$SampProbiInt <- rep(1, length(dat[,1]))
  
  ## Stratum Sampling probabilities
  SampProbInt <- SampledInt[[3]]
  SampProbSlp <- SampledSlp[[3]]
  SampProbBiv <- SampledBiv[[3]]
  SampProbRan <- c(1,1,1)
  
  ## Datasets for sampled subjects
  datInt <- dat[dat$SampledInt==1,]
  datSlp <- dat[dat$SampledSlp==1,]
  datBiv <- dat[dat$SampledBiv==1,]
  datRan <- dat[dat$SampledRan==1,]
  
  ## Datasets for unsampled subjects
  datNotInt <- dat[dat$SampledInt==0,]
  datNotSlp <- dat[dat$SampledSlp==0,]
  datNotBiv <- dat[dat$SampledBiv==0,]
  datNotRan <- dat[dat$SampledRan==0,]
    
    ###########################################################################
   ##
  progress("Model Fitting Begins")
  
  progress("Random Sampling")
  Fit.ran     <- acml.lmem(formula.fixed=y~time*grp+conf, formula.random= ~time, id=id, data=datRan, InitVals=inits, ProfileCol=NA,
                          cutpoints=cutoffs$SlpCutUniv, SampProb=SampProbRan, SampProbiWL=SampProbiRan, w.function="slope")
  progress("Univariate intercept ACML")
  Fit.int     <- acml.lmem(formula.fixed=y~time*grp+conf, formula.random= ~time, id=id, data=datInt, InitVals=inits, ProfileCol=NA,
                           cutpoints=cutoffs$IntCutUniv, SampProb=SampProbInt, SampProbiWL=SampProbiInt, w.function="intercept")
  progress("Univariate intercept WL")
  Fit.int.wl  <- acml.lmem(formula.fixed=y~time*grp+conf, formula.random= ~time, id=id, data=datInt, InitVals=Fit.int$coefficients, ProfileCol=NA,
                           cutpoints=cutoffs$IntCutUniv, SampProb=c(1,1,1), SampProbiWL=SampProbiIntWL, w.function="intercept")

  progress("Univariate slope ACML")
  Fit.slp     <- acml.lmem(formula.fixed=y~time*grp+conf, formula.random= ~time, id=id, data=datSlp, InitVals=inits, ProfileCol=NA,
                           cutpoints=cutoffs$SlpCutUniv, SampProb=SampProbSlp, SampProbiWL=SampProbiSlp, w.function="slope")
  
  progress("Univariate slope WL")
  Fit.slp.wl  <- acml.lmem(formula.fixed=y~time*grp+conf, formula.random= ~time, id=id, data=datSlp, InitVals=Fit.slp$coefficients, ProfileCol=NA,
                           cutpoints=cutoffs$SlpCutUniv, SampProb=c(1,1,1), SampProbiWL=SampProbiSlpWL, w.function="slope")
    
  progress("Bivariate ACML")
  Fit.biv     <- acml.lmem(formula.fixed=y~time*grp+conf, formula.random= ~time, id=id, data=datBiv, InitVals=inits, ProfileCol=NA,
                           cutpoints=c(cutoffs$IntCutBiv, cutoffs$SlpCutBiv), SampProb=SampProbBiv, SampProbiWL=SampProbiBiv, w.function="bivar")
  
  progress("Bivariate WL")
  Fit.biv.wl  <- acml.lmem(formula.fixed=y~time*grp+conf, formula.random= ~time, id=id, data=datSlp, InitVals=Fit.biv$coefficients, ProfileCol=NA,
                           cutpoints=c(cutoffs$IntCutBiv, cutoffs$SlpCutBiv), SampProb=c(1,1), SampProbiWL=SampProbiBivWL, w.function="bivar")

    ###########################################################################
   ##
  progress("Indirect Imputation Analyses", "... ran.mi2")
  Fit.ran.mi2 <- IndirectImputation(acml.fit=Fit.ran, datSampled=datRan, datNotSampled=datNotRan, n.imp=n.imp)
  progress("... int.mi2")
  Fit.int.mi2 <- IndirectImputation(acml.fit=Fit.int, datSampled=datInt, datNotSampled=datNotInt, n.imp=n.imp)
  progress("... slp.mi2")
  Fit.slp.mi2 <- IndirectImputation(acml.fit=Fit.slp, datSampled=datSlp, datNotSampled=datNotSlp, n.imp=n.imp)
  progress("... biv.mi2")
  Fit.biv.mi2 <- IndirectImputation(acml.fit=Fit.biv, datSampled=datBiv, datNotSampled=datNotBiv, n.imp=n.imp)

    ###########################################################################
   ##
  progress("Direct Imputation Analyses")

  datIntMiss <- datSlpMiss <-datBivMiss <- datRanMiss <- dat

  datIntMiss$grp[datIntMiss$SampledInt==0] <- NA
  datSlpMiss$grp[datSlpMiss$SampledSlp==0] <- NA
  datBivMiss$grp[datBivMiss$SampledBiv==0] <- NA
  datRanMiss$grp[datRanMiss$SampledRan==0] <- NA

  datIntMiss <- datIntMiss[,c("id","y","time","grp","conf","ymean","ytmean","tmean","t2mean")]
  datSlpMiss <- datSlpMiss[,c("id","y","time","grp","conf","ymean","ytmean","tmean","t2mean")]
  datBivMiss <- datBivMiss[,c("id","y","time","grp","conf","ymean","ytmean","tmean","t2mean")]
  datRanMiss <- datRanMiss[,c("id","y","time","grp","conf","ymean","ytmean","tmean","t2mean")]

  progress("... ran.mi")
  Fit.ran.mi <- DirectImputation(new.dat=datRanMiss, n.imp=n.imp)
  progress("... int.mi")
  Fit.int.mi <- DirectImputation(new.dat=datIntMiss, n.imp=n.imp)
  progress("... slp.mi")
  Fit.slp.mi <- DirectImputation(new.dat=datSlpMiss, n.imp=n.imp)
  progress("... biv.mi")
  Fit.biv.mi <- DirectImputation(new.dat=datBivMiss, n.imp=n.imp)
  

  progress(paste0("Run ", count, " complete. Saving data."))
  save( Fit.ran,      Fit.int,      Fit.slp,      Fit.biv,
                      Fit.int.wl,   Fit.slp.wl,   Fit.biv.wl,
        Fit.ran.mi,   Fit.int.mi,   Fit.slp.mi,   Fit.biv.mi,
        Fit.ran.mi2,  Fit.int.mi2,  Fit.slp.mi2,  Fit.biv.mi2,
        file=paste0("output/run-", run, "-", count, ".RData")
      ) 
}

  #############################################################################
 ##
## For exection on local desktop
# library(parallel)
# 
# mclapply(1:4, mc.cores=8, function(x)
# {
#  set.seed(x)
#  sapply(1:4, function(y) simulation(x, y))
# })

  #############################################################################
 ##
## ACCRE batch run
args <- commandArgs(trailingOnly=TRUE)
x    <- as.numeric(args[1])
set.seed(x)
sapply(1:4, function(y) simulation(x,y))

