#library(MASS, lib.loc="/usr/local/biostat_r/lib/R")
#library(MASS)

expit <- function(x) exp(x)/(1+exp(x))
## Standardized Gamma Random effect
gen.std.gam <- function(n.obs, shape, rate){
    (rgamma(n.obs, shape, rate) - shape/rate)/(sqrt(shape/(rate^2)))
}

## Standardized T random effect
gen.std.t <- function(n.obs, df){
	rt(n.obs, df)*sqrt((df-2)/df)
	}

## Standardized binormal random effect
gen.std.binormal <- function(n.obs, mn0, mn1, sd0, sd1, p1){
	group <- rbinom(n.obs, 1, p1)
	val <- rnorm(n.obs, mn0, sd0)*(1-group)+rnorm(n.obs, mn1, sd1)*group
	grp.tmp <- rbinom(50000, 1, p1)
	val.tmp <- rnorm(50000, mn0, sd0)*(1-grp.tmp)+rnorm(50000, mn1, sd1)*grp.tmp

	out <- (val-mean(val.tmp))/(sqrt(var(val.tmp)))
    out
	}

## Standardized mixture of normals with difference variances
gen.std.mixnorm <- function(n.obs, mn0, mn1, sd0, sd1, p1, group){
	val     <- rnorm(n.obs, mn0, sd0)*(1-group)+rnorm(n.obs, mn1, sd1)*group
	grp.tmp <- rbinom(50000, 1, p1)
	val.tmp <- rnorm(50000, mn0, sd0)*(1-grp.tmp)+rnorm(50000, mn1, sd1)*grp.tmp
	out     <- (val-mean(val.tmp))/(sqrt(var(val.tmp)))
    out
}
GenerateX <- function(N, n, prev.grp, c.parm){
    id   <- rep(1:N, each=n)
    time <- rep(c(0:(n-1)), N)
    grp.tmp <- rbinom(N,1, prev.grp)
    conf.tmp <- rnorm(N, c.parm[1]+grp.tmp*c.parm[2], 1) 
    grp  <- rep(grp.tmp, each=n)
    conf <- rep(conf.tmp, each=n)
    out <- data.frame(id=id, time=time, grp=grp, conf=conf)
    out
}

GenerateY <- function(X, Z, id, beta, sig.b0 = 0.25, sig.b1 = 0.25, rho = 0, sig.e = 0.5, RanefDist, ErrorDist){
    lp <- X %*% beta
    cov.mat  <- matrix(c(sig.b0^2, rho*sig.b0*sig.b1, rho*sig.b0*sig.b1, sig.b1^2),2,2)
    ni       <- c(unlist(tapply(id, id, length)))
    N        <- length(unique(id))
    sum.ni   <- length(id)
    if (RanefDist=="Gaussian"){bi <- mvrnorm(N, mu=c(0,0), Sigma=cov.mat)
    }else if (RanefDist=="Gamma5"){b0i <- gen.std.gam(n.obs=N, shape=5, rate=sqrt(3))*sig.b0
                              b1i <- gen.std.gam(n.obs=N, shape=5, rate=sqrt(3))*sig.b1
                              bi  <- cbind(b0i, b1i)
    }else if (RanefDist=="Binormal1"){b0i <- gen.std.binormal(n.obs=N, mn0=0, mn1=1, sd0=1, sd1=1, p1=.25)*sig.b0
                                 b1i <- gen.std.binormal(n.obs=N, mn0=0, mn1=1, sd0=1, sd1=1, p1=.25)*sig.b1
                                 bi  <- cbind(b0i, b1i)
    }else if (RanefDist=="MixNorm2"){b0i <- gen.std.mixnorm(n.obs=N, mn0=0, mn1=0, sd0=1, sd1=sqrt(2), p1=pr.grp.1.overall, group=grp.i)*sig.b0
                                b1i <- gen.std.mixnorm(n.obs=N, mn0=0, mn1=0, sd0=1, sd1=sqrt(2), p1=pr.grp.1.overall, group=grp.i)*sig.b1
                                bi  <- cbind(b0i, b1i)
    }
    b        <- cbind(rep(bi[,1], ni),rep(bi[,2], ni))

    if (ErrorDist=="Gaussian"){error <- mvrnorm(sum.ni, mu=0, Sigma=sig.e)
    }else if (ErrorDist=="Gamma5"){error <- gen.std.gam(n.obs=sum.ni, shape=5, rate=sqrt(3))*sig.e
    }else if (RanefDist=="Binormal1"){error <- gen.std.binormal(n.obs=sum.ni, mn0=0, mn1=1, sd0=1, sd1=1, p1=.25)*sig.e
    }else if (RanefDist=="MixNorm2"){error <- gen.std.mixnorm(n.obs=sum.ni, mn0=0, mn1=0, sd0=1, sd1=sqrt(2), p1=pr.grp.1.overall, group=grp.i)*sig.e
    }

    Y <- X %*% beta + Z[,1]*b[,1] + Z[,2]*b[,2] + error
    return(Y)
}
## generate a missing at random mechanism
GenerateMAR <- function(Y, X, id, param, cutpoint.drop){
    ## set param[2] to 0 for MCAR
    keep <- rep(1, length(Y))
    L    <- length(Y)
    for (j in 2:L){
        if (id[j]==id[j-1] & keep[j-1]==1){ keep[j] <- rbinom(1,1,expit(param[1]+param[2]*I(Y[j-1]<cutpoint.drop) + param[3]*X[j]))
        }else{                              keep[j] <- 1}
    }
    keep
}

## do subject-specific linear regressions. Important to note that data must contain variables Y and time
LinRegFn <- function(data){  X  <- cbind(1, data$time)
                             Xt <- t(X)
                             Y  <- data$Y
                             solve(Xt %*% X) %*% Xt %*% Y}

## calc subject specific intercepts and slopes and output them with the same length as the longitudinal data
CalcSSIntSlp <- function( Y, time, id){
    data.tmp  <- data.frame(id=id, Y=Y, time=time)
    data.list <- split(data.tmp, id)
    L.id      <- c(unlist(tapply(id,id,length)))
    mtx       <- matrix(unlist(lapply(data.list, LinRegFn)), byrow=TRUE, ncol=2)
    out       <- list(Int = rep(mtx[,1], L.id), Slp = rep(mtx[,2], L.id))}

## Calculate the cutpoints for defining sampling strata for the int- slope- and biv-based sampling
## based on the proportion we want in the central region
## If you have a problem with this function it is likely due to the search for the central region
## under bivariate sampling.  You may have to adjust Del.  I think it is due to the discreteness of
## due to insufficient sample size.
est.cutoffs <- function(Y, time, id, PropInCentralRegion){
    p         <- PropInCentralRegion
    data.tmp  <- data.frame(id=id, Y=Y, time=time)
    data.list <- split(data.tmp, id)
    #print("Running individual regressions")
    out      <- matrix(unlist(lapply(data.list, LinRegFn)), byrow=TRUE, ncol=2)

    Ints <- quantile(out[,1], c((1-p)/2,(1+p)/2))
    Slps <- quantile(out[,2], c((1-p)/2,(1+p)/2))

    q1 <- .99
    Del <- 1
    while (Del>0.003){
        q1 <- q1-.001
        Del <- abs(mean(out[,1] > quantile(out[,1], probs=1-q1) & out[,1] < quantile(out[,1], probs=q1) &
                        out[,2] > quantile(out[,2], probs=1-q1) & out[,2] < quantile(out[,2], probs=q1)) - PropInCentralRegion)
        }

    out <- list(IntCutUniv = Ints,
                SlpCutUniv = Slps,
                IntCutBiv   = c(quantile(out[,1], probs=1-q1), quantile(out[,1], probs=q1)),
                SlpCutBiv   = c(quantile(out[,2], probs=1-q1), quantile(out[,2], probs=q1)))
    out
}

identify.stratum <- function(Y, time, id, w.function, cutpoints, Int, Slp){
    ## under bivar sampling cutpoints should be of the form (int.lower, int.upper, slp.lower, slp.upper)

    if (w.function == "intercept") stratum <- ifelse( Int <  cutpoints[1], 1,
                                              ifelse( Int >= cutpoints[2], 3, 2))
    if (w.function == "slope")     stratum <- ifelse( Slp <  cutpoints[1], 1,
                                              ifelse( Slp >= cutpoints[2], 3, 2))
    if (w.function=="bivar") stratum <- ifelse(Int >= cutpoints[1] & Int < cutpoints[2] & Slp >=cutpoints[3] & Slp <cutpoints[4], 1, 2)
    stratum
}


ods.sampling <- function(id.long,          # id in long format
                         stratum.long,     # stratum identifier in long format)
                         SamplingStrategy, # "IndepODS" or "DepODS"
                         NsPerStratum){    # target number of subjects sampled per stratum.  This is exact if using Dependent sampling
                                           # This should be of length 3 for univariate sampling and of length 2 for bivariate sampling

    strat.1       <- c(unlist(tapply(stratum.long, id.long, unique)))
    id.1          <- c(unlist(tapply(id.long, id.long, unique)))
    ni            <- c(unlist(tapply(id.long, id.long, length)))
    N             <- length(id.1)
    NPerStratum   <- c(unlist(tapply(id.1, strat.1, length)))
    SampleTooMany <- any(NsPerStratum>NPerStratum)
    if (SampleTooMany){ print("Warning: You want to sample more people than you have in one of the strata.  Sampling from that stratum with probability 1")
                        WhichStratum <- which(NsPerStratum>NPerStratum)
                        NsPerStratum[WhichStratum] <- NPerStratum[WhichStratum]}
    SampProbs <- NsPerStratum/NPerStratum
    SampProb.1 <- ifelse(strat.1==1, SampProbs[1],
                  ifelse(strat.1==2, SampProbs[2], SampProbs[3]))
    if (SamplingStrategy=="IndepODS"){Samp <- rbinom(N, 1, SampProb.1)}
    if (SamplingStrategy=="DepODS"){  Sampled.ids <- NULL
                                      for (mm in 1:length(NPerStratum)){
                                          Sampled.ids <- c( Sampled.ids, c(sample(id.1[strat.1==mm], NsPerStratum[mm], replace=FALSE)))}
                                      Samp <- ifelse(id.1 %in% Sampled.ids, 1, 0)}
    Sampled <- list(Sampled=rep(Samp, ni),SampProbi=rep(SampProb.1, ni), SampProbs=SampProbs)
    Sampled
}

## note that we really do not need quants, PopnQuants, and w.function but to make the fitting function work
random.sampling <- function(id.long, n=225){
    s <- sample(unique(id.long), n)
    Sampled <- as.integer(id.long %in% s)
    Sampled
    }


# 
# 
# 
# ODS.Sampling <- function(dat,                      ## a list generated from the GenPopnData() function
#                          PopnQuants,               ## a matrix from est.quants function
#                          w.function,               ## Response summary to sample on ("mean","intercept", or "slope")
#                          quants,                   ## Population quantiles to define the theoretical sampling strata
#                          TargetNSampledPerStratum, ## Theoretical (and possibly observed) number sampled per stratum
#                          SamplingStrategy,         ## Options are "IndepODS" and "DepODS"
#                          Univariate=FALSE){
# 
# 
#     dat        = dat
#     w.function = "intercept"
#     TargetNSampledPerStratum <- c(100,100,100)
# 
#     NCohort                    <- length(unique(dat$id))
#     NStratumThry               <- round(NCohort*c(quants[1], quants[2]-quants[1], 1-quants[2]))
#     SampProbThry               <- TargetNSampledPerStratum / NStratumThry
#     SampProbThry[1] <- ifelse(SampProbThry[1]>1, 1, SampProbThry[1])
#     SampProbThry[2] <- ifelse(SampProbThry[2]>1, 1, SampProbThry[2])
#     SampProbThry[3] <- ifelse(SampProbThry[3]>1, 1, SampProbThry[3])
# 
#     C1 <- ifelse(w.function=="mean",      PopnQuants[2,match(quants[1], PopnQuants[1,])],
#           ifelse(w.function=="intercept", PopnQuants[3,match(quants[1], PopnQuants[1,])],
#           ifelse(w.function=="slope",     PopnQuants[4,match(quants[1], PopnQuants[1,])])))
#     C2 <- ifelse(w.function=="mean",      PopnQuants[2,match(quants[2], PopnQuants[1,])],
#           ifelse(w.function=="intercept", PopnQuants[3,match(quants[2], PopnQuants[1,])],
#           ifelse(w.function=="slope",     PopnQuants[4,match(quants[2], PopnQuants[1,])])))
# 
#     uid <- unique(dat$id)
#     ni  <- c(unlist(tapply(dat$id, dat$id, length)))
#     SampVar <- NULL
#     for(i in uid){ yi      <- dat$Y[dat$id==i]
#     xi      <- dat$X[dat$id==i,1:2]
#     SampVar <- rbind(SampVar, c(mean(yi), solve(t(xi)%*%xi) %*% t(xi) %*% yi))
#     }
#     SampVar <- (w.function=="mean")*SampVar[,1] +
#                (w.function=="intercept")*SampVar[,2] +
#                (w.function=="slope")*SampVar[,3]
# 
#     SampStratum  <- ifelse(SampVar<C1, 1,
#                            ifelse(SampVar<C2, 2, 3))
#     NperStratum  <- unlist(tapply(uid, SampStratum, length))
# 
#     SampProbiThry <- ifelse(SampVar<C1, SampProbThry[1],
#                             ifelse(SampVar<C2, SampProbThry[2], SampProbThry[3]))
# 
#     Sampled     <- rbinom(length(SampProbiThry), 1, SampProbiThry)
#     SampProbObs <- c(tapply(Sampled, SampStratum, mean))
# 
#     SampProbiObs  <- ifelse(SampVar<C1, SampProbObs[1],
#                             ifelse(SampVar<C2, SampProbObs[2], SampProbObs[3]))
#     #print(rbind(SampProbThry, SampProbObs))
#     #print(cbind(SampProbiThry,SampProbiObs))
# 
# 
#     ## Independent Sampling
#     if (SamplingStrategy=="IndepODS") InODS <- uid[ Sampled==1]
# 
#     TheSample <- dat$id %in% InODS
#     X.ods     <- dat$X[TheSample,]
#     Y.ods     <- dat$Y[TheSample]
#     Z.ods     <- dat$Z[TheSample,]
#     id.ods    <- dat$id[TheSample]
#     if (SamplingStrategy=="IndepODS"){
#         SampProbi.ods <- rep(SampProbiThry, ni)[TheSample]
#         SampProb.ods  <- SampProbThry
#     }
#     SampStrat.ods <- SampStratum[InODS]
#     Qi <- SampVar[InODS]
#     dup.id <- duplicated(dat$id)
#     dat.univariate <- dat$X[!dup.id,]
#     dat.univ.ods   <- dat.univariate[InODS,]
# 
#     SampProb  <- SampProbThry
#     SampProbi <- rep(SampProbiThry, ni)
#     Qi        <- SampVar
# 
#     list(X=dat$X, Y=dat$Y, Z=dat$Z, id=dat$id,
#          SampProb=SampProb, SampProbi=SampProbi,
#          N=dat$N, n=dat$n, beta=dat$beta, sig.b0=dat$sig.b0,
#          sig.b1=dat$sig.b1, rho=dat$rho, sig.e=dat$sig.e,
#          prob.grp=dat$prob.grp, w.function=w.function,
#          cutpoint=c(C1,C2), SampStratum=SampStratum, Qi=Qi, InSample=TheSample)
# 
# }
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# GenPopnData <- function(N = 1000, n = 11, beta = c(1, 0.25, 0, 0.25),
#                      sig.b0 = 0.25, sig.b1 = 0.25, rho = 0, sig.e = 0.5,
#                      pr.conf=0.25, pr.grp.0 = 0.1, pr.grp.1 = 0.2, n.low=11, n.high=11,
#                      dist){
#     id       <- rep(1:N, each=n)
#     conf.i   <- rbinom(N, 1, pr.conf)
#     pr.grp.i <- ifelse(conf.i==1, pr.grp.1, pr.grp.0)
#     grp.i    <- rbinom(N,1,pr.grp.i)
#     conf     <- rep(conf.i, each=n)
#     grp      <- rep(grp.i, each=n)
#     #grp      <- rep(rbinom(N, 1, prob.grp), each=n)
#     ###### Add between subject variation in the time variable
#     ###### Note that we are creating ICC equal to 0.5
#     ###### var(sample(seq(-2,2,length.out=10), 100000, replace=TRUE))
# 
#     bl.age   <- rep(rep(0,N), each=n)
#     #time     <- rep(seq(-2,2, length.out=n), N)
#     time     <- rep(c(0:(n-1)), N)
#     age      <- bl.age + time
#     pr.grp.1.overall <- pr.conf*pr.grp.1 + (1-pr.conf)*pr.grp.0
# 
#     cov.mat  <- matrix(c(sig.b0^2, rho*sig.b0*sig.b1, rho*sig.b0*sig.b1, sig.b1^2),2,2)
#     if (dist %in% c("Gaussian","GaussianBal","Gamma5err","T5err","MixNorm3err")){bi <- mvrnorm(N, mu=c(0,0), Sigma=cov.mat)
#     }else if (dist=="Gamma"){b0i <- gen.std.gam(n.obs=N, shape=3, rate=sqrt(3))*sig.b0
#     	               b1i <- gen.std.gam(n.obs=N, shape=3, rate=sqrt(3))*sig.b1
#     	               bi  <- cbind(b0i, b1i)
#     }else if (dist=="Gamma5"){b0i <- gen.std.gam(n.obs=N, shape=5, rate=sqrt(3))*sig.b0
#     	               b1i <- gen.std.gam(n.obs=N, shape=5, rate=sqrt(3))*sig.b1
#     	               bi  <- cbind(b0i, b1i)
#     }else if (dist=="Gamma10"){b0i <- gen.std.gam(n.obs=N, shape=10, rate=sqrt(3))*sig.b0
#     	               b1i <- gen.std.gam(n.obs=N, shape=10, rate=sqrt(3))*sig.b1
#     	               bi  <- cbind(b0i, b1i)
#     }else if (dist=="Gamma15"){b0i <- gen.std.gam(n.obs=N, shape=15, rate=sqrt(3))*sig.b0
#     	               b1i <- gen.std.gam(n.obs=N, shape=15, rate=sqrt(3))*sig.b1
#     	               bi  <- cbind(b0i, b1i)
#     }else if (dist=="T5"){b0i <- gen.std.t(n.obs=N, df=5)*sig.b0
#     	           b1i <- gen.std.t(n.obs=N, df=5)*sig.b1
#     	           bi  <- cbind(b0i, b1i)
#     }else if (dist=="T10"){b0i <- gen.std.t(n.obs=N, df=10)*sig.b0
#     	           b1i <- gen.std.t(n.obs=N, df=10)*sig.b1
#     	           bi  <- cbind(b0i, b1i)
#     }else if (dist=="T15"){b0i <- gen.std.t(n.obs=N, df=15)*sig.b0
#     	           b1i <- gen.std.t(n.obs=N, df=15)*sig.b1
#     	           bi  <- cbind(b0i, b1i)
#     }else if (dist=="Binormal1"){b0i <- gen.std.binormal(n.obs=N, mn0=0, mn1=1, sd0=1, sd1=1, p1=.25)*sig.b0
#     	                  b1i <- gen.std.binormal(n.obs=N, mn0=0, mn1=1, sd0=1, sd1=1, p1=.25)*sig.b1
#     	                  bi  <- cbind(b0i, b1i)
#     }else if (dist=="Binormal2"){b0i <- gen.std.binormal(n.obs=N, mn0=0, mn1=2, sd0=1, sd1=1, p1=.25)*sig.b0
#     	                  b1i <- gen.std.binormal(n.obs=N, mn0=0, mn1=2, sd0=1, sd1=1, p1=.25)*sig.b1
#     	                  bi  <- cbind(b0i, b1i)
#     }else if (dist=="MixNorm1.25"){b0i <- gen.std.mixnorm(n.obs=N, mn0=0, mn1=0, sd0=1, sd1=sqrt(1.25), p1=pr.grp.1.overall, group=grp.i)*sig.b0
#     	                   b1i <- gen.std.mixnorm(n.obs=N, mn0=0, mn1=0, sd0=1, sd1=sqrt(1.25), p1=pr.grp.1.overall, group=grp.i)*sig.b1
#     	                   bi  <- cbind(b0i, b1i)
#     }else if (dist=="MixNorm1.5"){b0i <- gen.std.mixnorm(n.obs=N, mn0=0, mn1=0, sd0=1, sd1=sqrt(1.5), p1=pr.grp.1.overall, group=grp.i)*sig.b0
#     	                   b1i <- gen.std.mixnorm(n.obs=N, mn0=0, mn1=0, sd0=1, sd1=sqrt(1.5), p1=pr.grp.1.overall, group=grp.i)*sig.b1
#     	                   bi  <- cbind(b0i, b1i)
#     }else if (dist=="MixNorm2"){b0i <- gen.std.mixnorm(n.obs=N, mn0=0, mn1=0, sd0=1, sd1=sqrt(2), p1=pr.grp.1.overall, group=grp.i)*sig.b0
#     	                   b1i <- gen.std.mixnorm(n.obs=N, mn0=0, mn1=0, sd0=1, sd1=sqrt(2), p1=pr.grp.1.overall, group=grp.i)*sig.b1
#     	                   bi  <- cbind(b0i, b1i)
#     }else if (dist=="MixNorm2.25"){b0i <- gen.std.mixnorm(n.obs=N, mn0=0, mn1=0, sd0=1, sd1=sqrt(2.25), p1=pr.grp.1.overall, group=grp.i)*sig.b0
#     	                   b1i <- gen.std.mixnorm(n.obs=N, mn0=0, mn1=0, sd0=1, sd1=sqrt(2.25), p1=pr.grp.1.overall, group=grp.i)*sig.b1
#     	                   bi  <- cbind(b0i, b1i)
#     }else if (dist %in% c("MixNorm3","MixNorm3err3")){b0i <- gen.std.mixnorm(n.obs=N, mn0=0, mn1=0, sd0=1, sd1=sqrt(3), p1=pr.grp.1.overall, group=grp.i)*sig.b0
#     	                   b1i <- gen.std.mixnorm(n.obs=N, mn0=0, mn1=0, sd0=1, sd1=sqrt(3), p1=pr.grp.1.overall, group=grp.i)*sig.b1
#     	                   bi  <- cbind(b0i, b1i)
#     }
#     b        <- cbind(rep(bi[,1], each=n),rep(bi[,2], each=n))
# 
#     if (dist == "Gamma5err"){
#     	    error <- gen.std.gam(n.obs=N*n, shape=5, rate=sqrt(3))*sig.e
#     	}else if (dist == "T5err"){
#     	    error <- gen.std.t(n.obs=N*n, df=5)*sig.e
#     	}else if (dist %in% c("MixNorm3err","MixNorm3err3")){
#     		error <- gen.std.mixnorm(n.obs=N*n, mn0=0, mn1=0, sd0=1, sd1=sqrt(3), p1=pr.grp.1.overall, group=grp)*sig.e
#      }else{
#      	error    <- rnorm(N*n, 0, sig.e)
#      }
# 
# 
#     X <- as.matrix(cbind(1, age, grp, age*grp, conf))
#     Z <- X[,c(1:2)]
#     Y <- X %*% beta + Z[,1]*b[,1] + Z[,2]*b[,2] + error
# 
#     ## Induce MCAR dropout
#     n.obs <- rep(sample(rep(c(n.low:n.high),2), N, replace=TRUE), each=n)
#     obs.num <- rep(c(1:n), N)
#     X <- X[obs.num<= n.obs,]
#     Z <- Z[obs.num<= n.obs,]
#     Y <- Y[obs.num<= n.obs]
#     id <- id[obs.num<= n.obs]
#     list(id=id, X=X, Y=Y, Z=Z, ni=c(unlist(tapply(Y,id,length))),
#          N=N, n=n, beta=beta,sig.b0=sig.b0, sig.b1=sig.b1, rho=rho, sig.e=sig.e, pr.grp.1=pr.grp.1, pr.grp.0=pr.grp.0)
# }
# 
# 
# est.quants <- function(N, n, beta, sig.b0, sig.b1, rho, sig.e, pr.grp.0, pr.grp.1, pr.conf, quant, n.low=11, n.high=11, dist){
#     d        <- GenPopnData(N=N, n=n, beta=beta, sig.b0=sig.b0, sig.b1=sig.b1, rho=rho, sig.e=sig.e,
#                            pr.conf=pr.conf, pr.grp.0 = pr.grp.0, pr.grp.1 = pr.grp.1, dist=dist)
#     data.tmp <- data.frame(id=d$id, Y=d$Y, X.time=d$X[,2])
#     out      <- matrix(unlist(lapply(split(data.tmp, data.tmp$id), LinRegFn)), byrow=TRUE, ncol=2)
#     out      <- cbind(c(tapply(data.tmp$Y, data.tmp$id, mean)), out)
#     ## outputs 4 row with quantiles (row 1), and then means, intercepts, and slopes (rows)
#     r        <- rbind( quant,
#                        quantile(out[,1], probs=quant),
#                        quantile(out[,2], probs=quant),
#                        quantile(out[,3], probs=quant))
#     rownames(r) <- c("quant", "mean", "int", "slp")
#     r}
# 
# ## Do a search to find the quantiles that correspond to the central rectangle that contains 60 and 80 percent of
# ## the subject specific intercepts and slopes.  This is not necessary if slope and intercepts are independent
# ## but with unequal followup they were positiviely correlated.  Searches for the smallest 'rectangle' defined
# ## by quantiles that contains 60 and 80 percent of the data
# est.bivar.lims <- function(N, n, beta, sig.b0, sig.b1, rho, sig.e, pr.grp.0, pr.grp.1, pr.conf, quants, n.low=11, n.high=11, dist){
# 	   d        <- GenPopnData(N=N, n=n, beta=beta, sig.b0=sig.b0, sig.b1=sig.b1, rho=rho, sig.e=sig.e,
#                                 pr.conf=pr.conf, pr.grp.0 = pr.grp.0, pr.grp.1 = pr.grp.1, dist=dist)
#         #d        <- GenPopnData(N, n, beta, sig.b0, sig.b1, rho, sig.e, prob.grp, n.low=n.low, n.high=n.high)
#         data.tmp <- data.frame(id=d$id, Y=d$Y, X.time=d$X[,2])
#         out      <- matrix(unlist(lapply(split(data.tmp, data.tmp$id), LinRegFn)), byrow=TRUE, ncol=2)
#         out      <- cbind(c(tapply(data.tmp$Y, data.tmp$id, mean)), out)
#         print(cor(out[,2], out[,3]))
# 
#         q1 <- .99
#         Del <- 1
#         while (Del>0.001){ q1 <- q1-.00025
#                            Del <- abs(mean(out[,2] > quantile(out[,2], probs=1-q1) &
#                                            out[,2] < quantile(out[,2], probs=q1) &
#                                            out[,3] > quantile(out[,3], probs=1-q1) &
#                                            out[,3] < quantile(out[,3], probs=q1)) - quants[1])
#                            #print(c(q1,Del))
#                            q1}
#         q2 <- .99
#         Del <- 1
#         while (Del>0.001){ q2 <- q2-.00025
#                            Del <- abs(mean(out[,2] > quantile(out[,2], probs=1-q2) &
#                                            out[,2] < quantile(out[,2], probs=q2) &
#                                            out[,3] > quantile(out[,3], probs=1-q2) &
#                                            out[,3] < quantile(out[,3], probs=q2)) - quants[2])
#                            q2}
# 
#         rbind( c( quantile(out[,2], probs=1-q1), quantile(out[,2], probs=q1), quantile(out[,3], probs=1-q1), quantile(out[,3], probs=q1)),
#                c( quantile(out[,2], probs=1-q2), quantile(out[,2], probs=q2), quantile(out[,3], probs=1-q2), quantile(out[,3], probs=q2)))
# 
# }
# 
# est.cutoffs <- function(N, n, beta, sig.b0, sig.b1, rho, sig.e, pr.grp.0, pr.grp.1, pr.conf, quant, n.low=11, n.high=11, dist){
#     d        <- GenPopnData(N=N, n=n, beta=beta, sig.b0=sig.b0, sig.b1=sig.b1, rho=rho, sig.e=sig.e,
#                             pr.conf=pr.conf, pr.grp.0 = pr.grp.0, pr.grp.1 = pr.grp.1, dist=dist)
#     data.tmp <- data.frame(id=d$id, Y=d$Y, X.time=d$X[,2])
#     out      <- matrix(unlist(lapply(split(data.tmp, data.tmp$id), LinRegFn)), byrow=TRUE, ncol=2)
#     out      <- cbind(c(tapply(data.tmp$Y, data.tmp$id, mean)), out)
#     ## outputs 4 row with quantiles (row 1), and then means, intercepts, and slopes (rows)
#     r        <- rbind( quant,
#                        quantile(out[,1], probs=quant),
#                        quantile(out[,2], probs=quant),
#                        quantile(out[,3], probs=quant))
#     rownames(r) <- c("quant", "mean", "int", "slp")
#     r}
# 
# 
# ODS.Sampling.Bivar <- function(dat,                      ## a list generated from the GenPopnData() function
#                                PopnQuantsBivar,          ## a matrix from est.bivar.lims function
#                                PopnPropInRectangle,      ## Proportion of subjects in the central rectangle
#                                TargetNSampledPerStratum, ## Theoretical (and possibly observed) number sampled per stratum
#                                SamplingStrategy){        ## Options are "IndepODS" and "DepODS"
# 
#     NCohort                    <- length(unique(dat$id))
#     NStratumThry               <- round(NCohort*c(PopnPropInRectangle, (1-PopnPropInRectangle)))
#     SampProbThry               <- TargetNSampledPerStratum / NStratumThry
#     SampProbThry[1] <- ifelse(SampProbThry[1]>1, 1, SampProbThry[1])
#     SampProbThry[2] <- ifelse(SampProbThry[2]>1, 1, SampProbThry[2])
# 
#     Lims <- (PopnPropInRectangle==.6)*PopnQuantsBivar[1,] + (PopnPropInRectangle==.8)*PopnQuantsBivar[2,]
# 
#     uid <- unique(dat$id)
#     ni  <- c(unlist(tapply(dat$id, dat$id, length)))
#     SampVar <- NULL
#     for(i in uid){ yi      <- dat$Y[dat$id==i]
#                    xi      <- dat$X[dat$id==i,1:2]
#                    SampVar <- rbind(SampVar, c(mean(yi), solve(t(xi)%*%xi) %*% t(xi) %*% yi))
#     }
#    print(sum(SampVar[,2]>Lims[1] & SampVar[,2]<Lims[2]))
#    print(sum(SampVar[,3]>Lims[3] & SampVar[,3]<Lims[4]))
# 
#     SampStratum  <- ifelse(SampVar[,2]>Lims[1] & SampVar[,2]<Lims[2] &
#                            SampVar[,3]>Lims[3] & SampVar[,3]<Lims[4], 1,2)
#     print(table(SampStratum))
# 
#     NperStratum  <- unlist(tapply(uid, SampStratum, length))
# 
# 
#     SampProbiThry <- ifelse(SampStratum==1, SampProbThry[1],
#                      ifelse(SampStratum==2, SampProbThry[2], NA))
# 
# 
#     Sampled     <- rbinom(length(SampProbiThry), 1, SampProbiThry)
#     SampProbObs <- c(tapply(Sampled, SampStratum, mean))
# 
#     SampProbiObs  <- ifelse(SampStratum==1, SampProbObs[1],
#                      ifelse(SampStratum==2, SampProbObs[2], NA))
# 
#     ## Independent Sampling
#     if (SamplingStrategy=="IndepODS") InODS <- uid[ Sampled==1]
# 
#     TheSample <- dat$id %in% InODS
#     X.ods     <- dat$X[TheSample,]
#     Y.ods     <- dat$Y[TheSample]
#     Z.ods     <- dat$Z[TheSample,]
#     id.ods    <- dat$id[TheSample]
#     if (SamplingStrategy=="IndepODS"){
#         SampProbi.ods <- rep(SampProbiThry, ni)[TheSample]
#         SampProb.ods  <- SampProbThry
#     }
#     #if (SamplingStrategy=="DepODS"){
#     #	   SampProbi.ods <- rep(SampProbiObs, ni)[TheSample]
#     #	   SampProb.ods  <- SampProbObs
#     #}
#     SampStrat.ods <- SampStratum[InODS]
#     Qi <- SampVar[InODS,2:3]
#     dup.id <- duplicated(dat$id)
#     dat.univariate <- dat$X[!dup.id,]
#     dat.univ.ods   <- dat.univariate[InODS,]
# 
#     list(#X=X.ods, Y=Y.ods, Z=Z.ods, id=id.ods,
#     X=dat$X, Y=dat$Y, Z=dat$Z, id=dat$id,
#       SampProb=SampProb.ods, SampProbi=SampProbi.ods,
#       N=dat$N, n=dat$n, beta=dat$beta, sig.b0=dat$sig.b0,
#       sig.b1=dat$sig.b1, rho=dat$rho, sig.e=dat$sig.e,
#       prob.grp=dat$prob.grp,w.function="bivar",
#       #cutpoint=c(C1,C2),
#       SampStratum=SampStrat.ods, Qi=Qi, dat.univ.ods=dat.univ.ods,
#       cutpoint=Lims, InSample=TheSample)
# }
# 
# ODS.Sampling <- function(dat,                      ## a list generated from the GenPopnData() function
#                          PopnQuants,               ## a matrix from est.quants function
#                          w.function,               ## Response summary to sample on ("mean","intercept", or "slope")
#                          quants,                   ## Population quantiles to define the theoretical sampling strata
#                          TargetNSampledPerStratum, ## Theoretical (and possibly observed) number sampled per stratum
#                          SamplingStrategy,         ## Options are "IndepODS" and "DepODS"
#                          Univariate=FALSE){
# 
#     NCohort                    <- length(unique(dat$id))
#     NStratumThry               <- round(NCohort*c(quants[1], quants[2]-quants[1], 1-quants[2]))
#     SampProbThry               <- TargetNSampledPerStratum / NStratumThry
#     SampProbThry[1] <- ifelse(SampProbThry[1]>1, 1, SampProbThry[1])
#     SampProbThry[2] <- ifelse(SampProbThry[2]>1, 1, SampProbThry[2])
#     SampProbThry[3] <- ifelse(SampProbThry[3]>1, 1, SampProbThry[3])
# 
#     C1 <- ifelse(w.function=="mean",      PopnQuants[2,match(quants[1], PopnQuants[1,])],
#           ifelse(w.function=="intercept", PopnQuants[3,match(quants[1], PopnQuants[1,])],
#           ifelse(w.function=="slope",     PopnQuants[4,match(quants[1], PopnQuants[1,])])))
#     C2 <- ifelse(w.function=="mean",      PopnQuants[2,match(quants[2], PopnQuants[1,])],
#           ifelse(w.function=="intercept", PopnQuants[3,match(quants[2], PopnQuants[1,])],
#           ifelse(w.function=="slope",     PopnQuants[4,match(quants[2], PopnQuants[1,])])))
# 
#     uid <- unique(dat$id)
#     ni  <- c(unlist(tapply(dat$id, dat$id, length)))
#     SampVar <- NULL
#     for(i in uid){ yi      <- dat$Y[dat$id==i]
#                    xi      <- dat$X[dat$id==i,1:2]
#                    SampVar <- rbind(SampVar, c(mean(yi), solve(t(xi)%*%xi) %*% t(xi) %*% yi))
#     }
#     SampVar <- (w.function=="mean")*SampVar[,1] +
#                (w.function=="intercept")*SampVar[,2] +
#                (w.function=="slope")*SampVar[,3]
# 
#     SampStratum  <- ifelse(SampVar<C1, 1,
#                     ifelse(SampVar<C2, 2, 3))
#     NperStratum  <- unlist(tapply(uid, SampStratum, length))
# 
#     SampProbiThry <- ifelse(SampVar<C1, SampProbThry[1],
#                      ifelse(SampVar<C2, SampProbThry[2], SampProbThry[3]))
# 
#     Sampled     <- rbinom(length(SampProbiThry), 1, SampProbiThry)
#     SampProbObs <- c(tapply(Sampled, SampStratum, mean))
# 
#     SampProbiObs  <- ifelse(SampVar<C1, SampProbObs[1],
#                      ifelse(SampVar<C2, SampProbObs[2], SampProbObs[3]))
#     #print(rbind(SampProbThry, SampProbObs))
#     #print(cbind(SampProbiThry,SampProbiObs))
# 
# 
#     ## Independent Sampling
#     if (SamplingStrategy=="IndepODS") InODS <- uid[ Sampled==1]
# 
#     TheSample <- dat$id %in% InODS
#     X.ods     <- dat$X[TheSample,]
#     Y.ods     <- dat$Y[TheSample]
#     Z.ods     <- dat$Z[TheSample,]
#     id.ods    <- dat$id[TheSample]
#     if (SamplingStrategy=="IndepODS"){
#         SampProbi.ods <- rep(SampProbiThry, ni)[TheSample]
#         SampProb.ods  <- SampProbThry
#     }
#     SampStrat.ods <- SampStratum[InODS]
#     Qi <- SampVar[InODS]
#     dup.id <- duplicated(dat$id)
#     dat.univariate <- dat$X[!dup.id,]
#     dat.univ.ods   <- dat.univariate[InODS,]
# 
#     SampProb  <- SampProbThry
#     SampProbi <- rep(SampProbiThry, ni)
#     Qi        <- SampVar
# 
#     list(X=dat$X, Y=dat$Y, Z=dat$Z, id=dat$id,
#       SampProb=SampProb, SampProbi=SampProbi,
#       N=dat$N, n=dat$n, beta=dat$beta, sig.b0=dat$sig.b0,
#       sig.b1=dat$sig.b1, rho=dat$rho, sig.e=dat$sig.e,
#       prob.grp=dat$prob.grp, w.function=w.function,
#       cutpoint=c(C1,C2), SampStratum=SampStratum, Qi=Qi, InSample=TheSample)
# 
# }
# #########################################################
# #########################################################
# #########################################################
# ## note that we really do not need quants, PopnQuants, and w.function but to make the fitting function work
# Random.Sampling <- function(d, quants=c(.1, .9), PopnQuants, w.function="mean", n=225){
#     s <- sample(unique(d$id), n)
#     TheSample <- d$id %in% s
#     X.rand     <- d$X
#     Y.rand     <- d$Y
#     Z.rand     <- d$Z
#     id.rand    <- d$id
# 
#     C1 <- PopnQuants[2,match(quants[1], PopnQuants[1,])]
#     C2 <- PopnQuants[2,match(quants[2], PopnQuants[1,])]
# 
#     list(X=X.rand, Y=Y.rand, Z=Z.rand, id=id.rand, InSample=TheSample,
#          N=d$N, n=d$n, beta=d$beta, sig.b0=d$sig.b0,
#          sig.b1=d$sig.b1, rho=d$rho, sig.e=d$sig.e,
#          prob.grp=d$prob.grp,
#          ## output below is only needed for the fitting function to run.  They are not used.
#          SampProb=c(1,1,1), cutpoint=c(C1,C2), SampProbi=rep(1, length(Y.rand)), w.function=w.function)
# }
