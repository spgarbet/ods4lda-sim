

expit <- function(x) exp(x)/(1+exp(x))

## Standardized Gamma Random effect
gen.std.gam <- function(n.obs, shape, rate) (rgamma(n.obs, shape, rate) - shape/rate)/(sqrt(shape/(rate^2)))

## Standardized T random effect
gen.std.t <- function(n.obs, df) rt(n.obs, df)*sqrt((df-2)/df)

## Standardized binormal random effect
gen.std.binormal <- function(n.obs, mn0, mn1, sd0, sd1, p1)
{
	group <- rbinom(n.obs, 1, p1)
	val <- rnorm(n.obs, mn0, sd0)*(1-group)+rnorm(n.obs, mn1, sd1)*group
	grp.tmp <- rbinom(50000, 1, p1)
	val.tmp <- rnorm(50000, mn0, sd0)*(1-grp.tmp)+rnorm(50000, mn1, sd1)*grp.tmp

	(val-mean(val.tmp))/(sqrt(var(val.tmp)))
}

## Standardized mixture of normals with difference variances
gen.std.mixnorm <- function(n.obs, mn0, mn1, sd0, sd1, p1, group)
{
	val     <- rnorm(n.obs, mn0, sd0)*(1-group)+rnorm(n.obs, mn1, sd1)*group
	grp.tmp <- rbinom(50000, 1, p1)
	val.tmp <- rnorm(50000, mn0, sd0)*(1-grp.tmp)+rnorm(50000, mn1, sd1)*grp.tmp
	
	(val-mean(val.tmp))/(sqrt(var(val.tmp)))
}


#################################

GenerateX <- function(N, n, prev.grp, c.parm)
{
    id   <- rep(1:N, each=n)
    time <- rep(c(0:(n-1)), N)
    
    # This is for properly specified data
    # grp.tmp  <- rbinom(N,1, prev.grp)
    # conf.tmp <- rnorm(N, c.parm[1]+grp.tmp*c.parm[2], 1) 
    
    # This is for model misspecification
    conf.tmp <- rnorm(N, 0, 1)
    alpha    <- c.parm[1]
    beta     <- c.parm[2]
    grp.tmp  <- rbinom(N, 1, inv.logit(alpha+beta*(conf.tmp^2)))
      
    # Start data frame
    grp  <- rep(grp.tmp, each=n)
    conf <- rep(conf.tmp, each=n)
    
    data.frame(id=id, time=time, grp=grp, conf=conf)
}

# Order of beta is 
# B_0, B_t, B_g, B_c, B_gt
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
random.sampling <- function(id.long, n=225)
{
    s <- sample(unique(id.long), n)
    Sampled <- as.integer(id.long %in% s)
    Sampled
}

