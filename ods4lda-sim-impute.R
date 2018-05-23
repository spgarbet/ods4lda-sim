  #############################################################################
 ##
expit     <- function(x) exp(x)/(1+exp(x))
Odds2Prob <- function(odds) odds/(1+odds)

cluster.summary <- function(id, x, fun)
{
  xlist    <- split( x, id )
  nj       <- unlist( lapply( xlist, length ) )
  xj       <- unlist( lapply( xlist, fun) )
  xsummary <- rep( xj, nj )
  xsummary
}

  #############################################################################
 ##

#####################
DirectImputation <- function(new.dat, n.imp)
{
    new.dat$grp[is.na(new.dat$grp)] <- 9999
    new.dat.wide <- reshape(new.dat, idvar=c("id", "conf","grp","ymean","ytmean","tmean","t2mean"), timevar=c("time"), direction="wide")
    new.dat.wide.merge <- new.dat.wide[,c("id",grep("y\\.", names(new.dat.wide), value=TRUE))]
    new.dat.wide$grp[new.dat.wide$grp==9999] <- NA
    #imp          <- aregImpute(~grp+conf+y.0+y.1+y.2+y.3+y.4, data=new.dat.wide, n.impute=n.imp, type="regression")
    #imp          <- aregImpute(~grp+conf+ymean+ytmean+y.0+y.1+y.2+y.3+y.4, data=new.dat.wide, n.impute=n.imp, nk=0, type="regression")
    imp          <- aregImpute(~grp+conf+ymean+ytmean+tmean+t2mean, data=new.dat.wide, n.impute=n.imp, nk=0, type="pmm")
    Ests.Imp <- Covs.Imp <- list()
    for (IMP in 1:n.imp)
    {
        ## get this complete imputation dataset
        imputed         <- as.data.frame(impute.transcan(imp, imputation=IMP, data=new.dat.wide, list.out=TRUE, pr=FALSE, check=FALSE))
        imputed$new.grp <- rbinom(length(imputed$grp),1, imputed$grp)
        imputed$grp     <- imputed$new.grp
        completed       <- imputed[,-grep("new.grp", names(imputed))]
        #tmp             <- completed[names(imputed)] <- cbind.data.frame(imputed)
        
        ## merge and reshape back to long format
        Completed <- cbind(new.dat.wide.merge, completed)
        #long.dat <- reshape(Completed, varying=c("y.0","y.1","y.2","y.3","y.4"), direction="long", idvar="id")
        long.dat <- reshape(Completed, varying=c(grep("y\\.", names(Completed), value=TRUE)), direction="long", idvar="id")
        long.dat <- long.dat[order(long.dat$id, long.dat$time),]
        
        ## do an lme on the fully imputed data
        #lme.mod.imp <- lme(Y~time*grp+conf,random=~time | id, data=long.dat)
        lme.mod.imp <- lmer(y~time*grp+conf+(time | id), data=long.dat)
        
        ## get estimates and covariances
        Ests.Imp[[IMP]] <- c(fixef(lme.mod.imp))#, as.numeric(VarCorr(lme.mod.imp)[,"StdDev"][1:2]),
        Covs.Imp[[IMP]] <- as.matrix(vcov(lme.mod.imp))
        
    }
    out.tmp <- MIcombine(Ests.Imp, Covs.Imp)
    list(coefficients = out.tmp$coefficients,
         covariance   = out.tmp$variance)
}


  #############################################################################
 ##
IndirectImputation <- function(acml.fit, datSampled, datNotSampled, n.imp)
{
    #datSampled <-datSlp
    #datNotSampled <- datNotSlp
    #n.imp <- n.imp
    #Fit <- Fit.slp
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
    if (is.na(ProfileCol)){ n.par  <- length(Fit$coefficients)
    params <- Fit$coefficients
    vcovs  <- Fit$covariance
    for (m1 in 2:n.par){ for (m2 in 1:(m1-1)){ vcovs[m1,m2] <- vcovs[m2,m1] }}
    }else{ n.par  <- length(Fit$coefficients)
    params <- c(Fit$coefficients[1:(ProfileCol-1)], 0, Fit$coefficients[ProfileCol:n.par])
    vcovs  <- cbind(Fit$covariance[,1:(ProfileCol-1) ], 0, Fit$covariance[,ProfileCol:n.par])
    vcovs  <- rbind(vcovs[1:(ProfileCol-1), ], 0, vcovs[ProfileCol:n.par,])
    for (m1 in 2:(n.par+length(ProfileCol))){ for (m2 in 1:(m1-1)){ #print(c(m1,m2))
        vcovs[m1,m2] <- vcovs[m2,m1] }}
    }
    Est.mi <- Cov.mi <- list()
    for (j in 1:n.imp){ #cat(j)
        
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
        #LR.S0           <- exp(prY.Xo.S0.Xe1)/exp(prY.Xo.S0.Xe0) #### is this right?
        LR.S0           <- prY.Xo.S0.Xe1/prY.Xo.S0.Xe0
        
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
    out.tmp <- MIcombine(Est.mi, Cov.mi)
    
    list(coefficients = out.tmp$coefficients,
         covariance   = out.tmp$variance)
}
