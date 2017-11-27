  #############################################################################
 ##
## Perform multiple simulations from a set of known parameters of ods4lda
## Compares random, ACML, WL, and Imputation for int, slp and bivar
##  
.libPaths("~/R/rlib-3.4.0")


library(lme4)
library(ODS4LDA)
library(reshape)
library(rms)
library(nlme)
library(mitools)
library(MASS)

source("ods4lda-sim-fun.R")

  #############################################################################
 ##
expit     <- function(x) exp(x)/(1+exp(x))
Odds2Prob <- function(odds) odds/(1+odds)

cluster.summary <- function( id, x, fun )
{
  xlist <- split( x, id )
  nj    <- unlist( lapply( xlist, length ) )
  xj    <- unlist( lapply( xlist, fun) )
  xsummary <- rep( xj, nj )
  xsummary
}

options(width=200)

  #############################################################################
 ##
DirectImputation <- function(new.dat, n.imp)
{
    new.dat$grp[is.na(new.dat$grp)] <- 9999
    new.dat.wide <- reshape(new.dat, idvar=c("id", "conf","grp","ymean","ytmean"), timevar=c("time"), direction="wide")
    new.dat.wide$grp[new.dat.wide$grp==9999] <- NA
    #imp          <- aregImpute(~grp+conf+y.0+y.1+y.2+y.3+y.4, data=new.dat.wide, n.impute=n.imp, type="regression")
    imp          <- aregImpute(~grp+conf+ymean+ytmean+y.0+y.1+y.2+y.3+y.4, data=new.dat.wide, n.impute=n.imp, nk=0, type="regression", pr=FALSE)

    Ests.Imp <- Covs.Imp <- list()
    cat("... ... Imputation # ")
    for (IMP in 1:n.imp)
    {
        cat(IMP, " ")

        ## get this complete imputation dataset
        imputed         <- as.data.frame(impute.transcan(imp, imputation=IMP, data=new.dat.wide, list.out=TRUE, pr=FALSE, check=FALSE))
        imputed$new.grp <- rbinom(length(imputed$grp),1, imputed$grp)
        imputed$grp     <- imputed$new.grp
        completed       <- imputed[,-grep("new.grp", names(imputed))]
        #tmp             <- completed[names(imputed)] <- cbind.data.frame(imputed)

        ## reshape back to long format
        long.dat <- reshape(completed, varying=c("y.0","y.1","y.2","y.3","y.4"), direction="long", idvar="id")
        long.dat <- long.dat[order(long.dat$id, long.dat$time),]

        ## do an lme on the fully imputed data
        lme.mod.imp <- lmer(y~time*grp+conf+(time | id), data=long.dat)

        ## get estimates and covariances
        Ests.Imp[[IMP]] <- c(fixef(lme.mod.imp))#, as.numeric(VarCorr(lme.mod.imp)[,"StdDev"][1:2]),
        Covs.Imp[[IMP]] <- as.matrix(vcov(lme.mod.imp))
    }
    cat("\n")
    out.tmp <- MIcombine(Ests.Imp, Covs.Imp)
    
    list(coefficients = out.tmp$coefficients, covariance = out.tmp$variance)
}

  #############################################################################
 ##
IndirectImputation <- function(acml.fit, datSampled, datNotSampled, n.imp)
{
  #datSampled <-datSampled
  #datNotSampled <- datNotSampled
  #n.imp <- n.imp
  
  Fit <- acml.fit

  ProfileCol   <- attr(Fit, "args")$ProfileCol
  SampProb     <- attr(Fit, "args")$SampProb
  SampProbiWL  <- attr(Fit, "args")$SampProbiWL
  w.function   <- attr(Fit, "args")$w.function
  cutpoints    <- attr(Fit, "args")$cutpoints
  fixed.form   <- attr(Fit, "args")$formula.fixed
  rand.form    <- attr(Fit, "args")$formula.random
  SampProbiVar <- attr(Fit, "args")$SampProbiVar

  ## Get estimates and variance-covariance matrix: Profiling is used right now only for the correlation parameter in the covariance matrix for the random effeects distribution
  if (is.na(ProfileCol))
  {
    n.par  <- length(Fit$coefficients)
    params <- Fit$coefficients
    vcovs  <- Fit$covariance
    for (m1 in 2:n.par)
    { 
      for (m2 in 1:(m1-1)) { vcovs[m1,m2] <- vcovs[m2,m1] }
    }
  } else {
    n.par  <- length(Fit$coefficients)
    params <- c(Fit$coefficients[1:(ProfileCol-1)], 0, Fit$coefficients[ProfileCol:n.par])
    vcovs  <- cbind(Fit$covariance[,1:(ProfileCol-1) ], 0, Fit$covariance[,ProfileCol:n.par])
    vcovs  <- rbind(vcovs[1:(ProfileCol-1), ], 0, vcovs[ProfileCol:n.par,])
    for (m1 in 2:(n.par+length(ProfileCol)))
    {
      for (m2 in 1:(m1-1)) {vcovs[m1,m2] <- vcovs[m2,m1] }
    }
  }
  Est.mi <- Cov.mi <- list()
  cat("... ... Imputation # ")
  for (j in 1:n.imp)
  {
    cat(j, " ")

      #######################################################################################
     ## Likelihood ratio piece of the imputation model
    ## Draw a sample from the parameter estimate distribution for the model Y | Xe, Xo, S=1
    Y.Xe.Xo.Seq1.MI.params <- rmvnorm(1, params, vcovs)

    ## Build exposure model piece of the imputation model based on the insample people.  Based on an offsetted logistic regression
    datSampled.0     <- datSampled.1 <- datSampled
    datSampled.0$grp <- 0
    datSampled.1$grp <- 1
    mod.y            <- model.response(model.frame(fixed.form, data=datSampled.0))
    fixed.matrix0    <- model.matrix(fixed.form, data=datSampled.0)
    fixed.matrix1    <- model.matrix(fixed.form, data=datSampled.1)
    rand.matrix      <- model.matrix(rand.form, data=datSampled.0)
    tmp.0            <- LogLikeCAndScore(params=Y.Xe.Xo.Seq1.MI.params, y=mod.y, x=fixed.matrix0, z=rand.matrix, id=datSampled.0$id,
                                         w.function=w.function, cutpoints=cutpoints, SampProb=SampProb, SampProbi=datSampled.0[,SampProbiVar], ProfileCol=ProfileCol, Keep.liC=TRUE)
    tmp.1            <- LogLikeCAndScore(params=Y.Xe.Xo.Seq1.MI.params, y=mod.y, x=fixed.matrix1, z=rand.matrix, id=datSampled.1$id,
                                         w.function=w.function, cutpoints=cutpoints, SampProb=SampProb, SampProbi=datSampled.1[,SampProbiVar], ProfileCol=ProfileCol, Keep.liC=TRUE)
    AC.Xe0.S1        <- exp(tmp.0$logACi)
    AC.Xe1.S1        <- exp(tmp.1$logACi)

    dup                        <- duplicated(datSampled$id)
    datSampled.firstobs        <- datSampled[!dup,]
    datSampled.firstobs$Offset <- log(AC.Xe1.S1/AC.Xe0.S1)
    Xe.Xo.Seq1.mod             <- glm(grp ~ conf+ offset(Offset), family=binomial, data=datSampled.firstobs)

    ## Draw a sample from the parameter estimate distribution for the model Xe | Xo, S=1
    Xe.Xo.Seq1.MI.params <- rmvnorm(1, Xe.Xo.Seq1.mod$coef, summary(Xe.Xo.Seq1.mod)$cov.unscaled)

    ## Apply to the unsampled people.

    ## Marginal exposure prevalence in unsampled.  Need to calculate an offset first
    datNotSampled.0 <- datNotSampled.1 <- datNotSampled
    datNotSampled.0$grp <- 0
    datNotSampled.1$grp <- 1
    mod.y            <- model.response(model.frame(fixed.form, data=datNotSampled.0))
    fixed.matrix0    <- model.matrix(fixed.form, data=datNotSampled.0)
    fixed.matrix1    <- model.matrix(fixed.form, data=datNotSampled.1)
    rand.matrix      <- model.matrix(rand.form, data=datNotSampled.0)

    tmp.0 <- LogLikeCAndScore(params=Y.Xe.Xo.Seq1.MI.params, y=mod.y, x=fixed.matrix0, z=rand.matrix, id=datNotSampled.0$id,
                              w.function=w.function, cutpoints=cutpoints, SampProb=SampProb, SampProbi=datNotSampled.0[,SampProbiVar], ProfileCol=ProfileCol, Keep.liC=TRUE)
    tmp.1 <- LogLikeCAndScore(params=Y.Xe.Xo.Seq1.MI.params, y=mod.y, x=fixed.matrix1, z=rand.matrix, id=datNotSampled.1$id,
                              w.function=w.function, cutpoints=cutpoints, SampProb=SampProb, SampProbi=datNotSampled.1[,SampProbiVar], ProfileCol=ProfileCol, Keep.liC=TRUE)
    prY.Xo.S0.Xe0 <- exp(tmp.0$liC)
    prY.Xo.S0.Xe1 <- exp(tmp.1$liC)

    AC.Xe0.S0  <- exp(tmp.0$logACi)
    AC.Xe1.S0  <- exp(tmp.1$logACi)

    dup                           <- duplicated(datNotSampled$id)
    datNotSampled.firstobs        <- datNotSampled[!dup,]
    datNotSampled.firstobs$Offset <- log(AC.Xe1.S0/AC.Xe0.S0)

    ## Marginal exposure odds for each subbject: pr(Xe=1 | Xo, S=0) / pr(Xe=0 | Xo, S=0)
    prXeEq1.Xo.S0 <- predict(Xe.Xo.Seq1.mod, newdata=datNotSampled.firstobs, type="response")
    ExposureOdds.S0 <- prXeEq1.Xo.S0/(1-prXeEq1.Xo.S0)

    ## Likelihood ratio for each subbject: pr(Y | Xe=1, Xo, S=0) / pr(Y | Xe=0, Xo, S=0)
    LR.S0           <- exp(prY.Xo.S0.Xe1)/exp(prY.Xo.S0.Xe0)

    ## Conditional exposure odds and then probability:  pr(Xe=1 | Y, Xo, S=0)
    odds            <- LR.S0*ExposureOdds.S0
    probs           <- Odds2Prob(odds)

    ## Do the imputation and add the imputed valued into the non-sampled subjects' data
    Impute.Xe.1       <- rbinom(length(probs), 1, probs)
    ni                <- c(unlist(tapply(datNotSampled$id, datNotSampled$id, length)))
    datNotSampled$grp <- rep(Impute.Xe.1, ni)

    ## Analysis of imputed data
    datAll.Imp <- rbind(datSampled, datNotSampled)
    datAll.Imp <- datAll.Imp[order(datAll.Imp$id, datAll.Imp$time),]
    ImpMod <- paste("lmer(", paste(fixed.form[c(2,3)], collapse="~"), paste("+(", rand.form[2], "| id) , data=datAll.Imp, REML=FALSE)", collapse=""), collapse="")
    lme.mi     <- eval(parse(text=ImpMod))
    #lmer(y~time*grp+conf+(time | id), data=datAll.Imp, REML=FALSE)

    ## Save results
    Est.mi[[j]] <- c(fixef(lme.mi))
    Cov.mi[[j]] <- as.matrix(vcov(lme.mi))
  }
  cat("\n")
  out.tmp <- MIcombine(Est.mi, Cov.mi)
  
  list(coefficients = out.tmp$coefficients, covariance   = out.tmp$variance)
}

  #############################################################################
 ##
## Known True Parameters
N                <- 1200
ni               <- 5
prev.grp         <- 0.4
conf.param       <- c(-0.5, 1)
n.imp            <- 200

tmp.inits <- rbind(c(75, -1, -1, -0, -.3, log(9),log(1),0,log(3.5)),
                   c(75, -1, -1, -1, -.3, log(9),log(1),0,log(3.5)),
                   c(75, -1, -1, -5, -.3, log(9),log(1),0,log(3.5)),
                   c(75, -1, -1, -1, -.3, log(4),log(1),0,log(3.5)))
tmpNsStratUniv <- rbind(c(100,50,100),c(100,50,100), c(100,50,100),c(100,50,100))
tmpNsStratBiv  <- rbind(c(50,200),c(50,200),c(50,200),c(50,200))

inits.sim            <- tmp.inits
NsPerStratumUniv.sim <- tmpNsStratUniv
NsPerStratumBiv.sim  <- tmpNsStratBiv

p.central        <- 0.8
NsRand           <- 250

  #############################################################################
 ##
## Progress helper
progress   <- function(...)
{
  lapply(list(...), function(x) cat(x,'\n'))
  cat(date(), '\n')
}

  #############################################################################
 ##
## Main Simulation Loop
simulation <- function(run, count)
{
  progress("Generating random data from known parameters")
  inits <- inits.sim[count,]
  NsPerStratumUniv <- NsPerStratumUniv.sim[count,]
  NsPerStratumBiv  <- NsPerStratumBiv.sim[count,]
    
    ##################################################
   ## Generate poplation with MAR missingness
  # dat          <- GenerateX(N=25000, n=10, prev.grp=.25, c.parm=c(-1, 2))
  # dat$y        <- GenerateY(X=cbind(1, dat$time, dat$grp, dat$conf, dat$time*dat$grp), Z=cbind(1, dat$time), id=dat$id,
  #                           beta=c(100, -1, -.5, -4, -.3), sig.b0 = 6, sig.b1 =1, rho = 0, sig.e = 2, RanefDist="Gaussian", ErrorDist="Gaussian")
  # dat$keep     <- GenerateMAR(Y=dat$y, X=dat$time, id=dat$id, param=c(2, -2, -.1), cutpoint.drop=65)
  # dat          <- dat[dat$keep==1,]
  # #dat$y[dat$keep==0] <- NA
  # cutoffs.popn <- est.cutoffs(Y=dat$y, time=dat$time, id=dat$id, PropInCentralRegion=0.8)
  ###################################################

  dat        <- GenerateX(N=N, n=ni, prev.grp=prev.grp, c.parm=conf.param)
  dat$y      <- GenerateY(X=cbind(1, dat$time, dat$grp, dat$conf, dat$time*dat$grp), Z=cbind(1, dat$time), id=dat$id,
                        beta=inits[1:5], sig.b0 = exp(inits[6]), sig.b1 =exp(inits[7]), rho = inits[8], sig.e = exp(inits[9]), RanefDist="Gaussian", ErrorDist="Gaussian")
  dat$ymean  <- cluster.summary(dat$id, dat$y, mean)
  dat$ytmean <- cluster.summary(dat$id, dat$y*dat$time, mean)
    
  ## deal with missingness soon.
  #dat$keep <- GenerateMAR(Y=dat$y, X=dat$time, id=dat$id, param=c(100, -2, -.1), cutpoint.drop=-100)
    
  dat$keep <- 1
  dat.tmp  <- dat[dat$keep==1,]
  cutoffs  <- est.cutoffs(Y=dat.tmp$y, time=dat.tmp$time, id=dat.tmp$id, PropInCentralRegion=p.central)
  
    #####################################################
   ## Generate cohort with MAR missingness
  #dat$y[dat$keep==0] <- NA
  #####################################################
  # PopnMod0    <- lm(y~time+conf, data=dat )
  # dat$ystar   <- PopnMod0$resid
  # IntSlps     <- CalcSSIntSlp( Y=dat$ystar, time=dat$time, id=dat$id)
  # dat$Intstar <- IntSlps[[1]]
  # dat$Slpstar <- IntSlps[[2]]
  
  IntSlps <- CalcSSIntSlp( Y=dat$y, time=dat$time, id=dat$id)
  dat$Int <- IntSlps[[1]]
  dat$Slp <- IntSlps[[2]]

  dat$StratInt <- identify.stratum(Y=dat$y, time=dat$time, id=dat$id, w.function="intercept", cutpoints=cutoffs$IntCutUniv, Int=dat$Int, Slp=dat$Slp)
  dat$StratSlp <- identify.stratum(Y=dat$y, time=dat$time, id=dat$id, w.function="slope",     cutpoints=cutoffs$SlpCutUniv, Int=dat$Int, Slp=dat$Slp)
  dat$StratBiv <- identify.stratum(Y=dat$y, time=dat$time, id=dat$id, w.function="bivar",     cutpoints=c(cutoffs$IntCutBiv,cutoffs$SlpCutBiv), Int=dat$Int, Slp=dat$Slp)

  SampledInt <- ods.sampling(id.long=dat$id, stratum.long=dat$StratInt, SamplingStrategy="IndepODS", NsPerStratum=NsPerStratumUniv)
  SampledSlp <- ods.sampling(id.long=dat$id, stratum.long=dat$StratSlp, SamplingStrategy="IndepODS", NsPerStratum=NsPerStratumUniv)
  SampledBiv <- ods.sampling(id.long=dat$id, stratum.long=dat$StratBiv, SamplingStrategy="IndepODS", NsPerStratum=NsPerStratumBiv)
  SampledRan <- random.sampling(id.long=dat$id, n=NsRand)

  dat$SampledInt <- SampledInt[[1]]
  dat$SampledSlp <- SampledSlp[[1]]
  dat$SampledBiv <- SampledBiv[[1]]
  dat$SampledRan <- SampledRan

  ## Use these if doing weighted likelihood
  dat$SampProbiIntWL <- SampledInt[[2]]
  dat$SampProbiSlpWL <- SampledSlp[[2]]
  dat$SampProbiBivWL <- SampledBiv[[2]]
  dat$SampProbiRanWL <- rep(1, length(dat[,1]))

  dat$SampProbiRan <- dat$SampProbiBiv <- dat$SampProbiSlp <- dat$SampProbiInt <- rep(1, length(dat[,1]))

  SampProbInt <- SampledInt[[3]]
  SampProbSlp <- SampledSlp[[3]]
  SampProbBiv <- SampledBiv[[3]]
  SampProbRan <- c(1,1,1)

  datInt <- dat[dat$SampledInt==1,]
  datSlp <- dat[dat$SampledSlp==1,]
  datBiv <- dat[dat$SampledBiv==1,]
  datRan <- dat[dat$SampledRan==1,]

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
  progress("Bivariate ACML Analysis")
  Fit.biv     <- acml.lmem(formula.fixed=y~time*grp+conf, formula.random= ~time, id=id, data=datBiv, InitVals=inits, ProfileCol=NA,
                           cutpoints=c(cutoffs$IntCutBiv, cutoffs$SlpCutBiv), SampProb=SampProbBiv, SampProbiWL=SampProbiBiv, w.function="bivar")
  progress("Bivariate WL Analysis")
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

  datIntMiss <- datIntMiss[,c("id","y","time","grp","conf","ymean","ytmean")]
  datSlpMiss <- datSlpMiss[,c("id","y","time","grp","conf","ymean","ytmean")]
  datBivMiss <- datBivMiss[,c("id","y","time","grp","conf","ymean","ytmean")]
  datRanMiss <- datRanMiss[,c("id","y","time","grp","conf","ymean","ytmean")]

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
# mclapply(1:200, mc.cores=8, function(x)
# {
#   set.seed(x)
#   sapply(1:4, function(y) simulation(x, y))
# })

  #############################################################################
 ##
## ACCRE batch run
args <- commandArgs(trailingOnly=TRUE)
x    <- as.numeric(args[1])
set.seed(x)
sapply(1:4, function(y) simulation(x,y))

