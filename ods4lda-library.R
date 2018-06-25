#' Calculate V_i = Z_i D t(Z_i) + sig_e^2 I_{n_i}
#'
#' Calculate V_i = Z_i D t(Z_i) + sig_e^2 I_{n_i}
#' @param zi n_i by 2 design matrix for the random effects
#' @param sigma0 std dev of the random intercept distribution
#' @param sigma1 std dev of the random slope distribution
#' @param rho correlation between the random intercept and slope
#' @param sigmae std dev of the measurement error distribution
#' @return V_i
#' @export
#'
vi.calc <- function(zi, sigma0, sigma1, rho, sigmae){
    rho.sig0.sig1 <- rho*sigma0*sigma1
    zi %*% matrix(c(sigma0^2, rho.sig0.sig1,rho.sig0.sig1,sigma1^2), nrow=2) %*% t(zi) +
        sigmae*sigmae*diag(length(zi[,1]))
}
#' Ascertainment correction piece for univariate sampling
#'
#' Calculate the (not yet log transformed) ascertainment correction under a univariate Q_i
#'
#' @param cutpoints cutpoints defining the sampling regions. (a vector of length 2)
#' @param SampProb Sampling probabilities from within each region (vector of length 3).
#' @param mu_q a scalar for the mean value of the Q_i distribution
#' @param sigma_q a scalar for the standard deviation of the Q_i distribution
#' @return Not yet log transformed ascertainment correction
#' @export
ACi1q <- function(cutpoints, SampProb, mu_q, sigma_q){
    CDFs <- pnorm(c(-Inf, cutpoints, Inf), mu_q, sigma_q)
    sum( SampProb*(CDFs[2:length(CDFs)] - CDFs[1:(length(CDFs)-1)]) )
}

#' Ascertainment correction piece for bivariate sampling
#'
#' Calculate the (not yet log transformed) ascertainment correction under a bivariate Q_i
#'
#' @param cutpoints cutpoints defining the sampling regions. (a vector of length 4: c(xlow, xhigh, ylow, yhigh))
#' @param SampProb Sampling probabilities from within each of two sampling regions; central region and outlying region (vector of length 2).
#' @param mu_q a 2-vector for the mean value of the bivariate Q_i distribution
#' @param sigma_q a 2 by 2 covariance matrix for the bivariate Q_i distribution
#' @return Not yet log transformed ascertainment correction
#' @export
ACi2q <- function(cutpoints, SampProb, mu_q, sigma_q){
    (SampProb[1]-SampProb[2])*pmvnorm(lower=c(cutpoints[c(1,3)]), upper=c(cutpoints[c(2,4)]), mean=mu_q, sigma=sigma_q)[[1]] + SampProb[2]
}

#' Log of the Ascertainment correction for univariate sampling
#'
#' Calculate the log transformed ascertainment correction under a univariate Q_i.  Also return vi
#'
#' @param yi n_i-response vector
#' @param xi n_i by p design matrix for fixed effects
#' @param zi n_i by 2 design matric for random effects (intercept and slope)
#' @param wi the pre-multiplier of yi to generate the sampling variable q_i
#' @param beta mean model parameter vector
#' @param sigma0 std dev of the random intercept distribution
#' @param sigma1 std dev of the random slope distribution
#' @param rho correlation between the random intercept and slope
#' @param sigmae std dev of the measurement error distribution
#' @param cutpoints cutpoints defining the sampling regions. (a vector of length 2)
#' @param SampProb Sampling probabilities from within each region (vector of length 3).
#' @return log transformed ascertainment correction
#' @export
logACi1q <- function(yi, xi, zi, wi, beta, sigma0, sigma1, rho, sigmae, cutpoints, SampProb){
    vi      <- vi.calc(zi, sigma0, sigma1, rho, sigmae)
    mu      <- xi %*% beta
    mu_q    <- (wi %*% mu)[,1]
    sigma_q <- sqrt((wi %*% vi %*% t(wi))[1,1])
    return(list(vi=vi, logACi=log(ACi1q(cutpoints, SampProb, mu_q, sigma_q))))
}

#' Log of the Ascertainment correction piece for bivariate sampling
#'
#' Calculate the log transformed ascertainment correction under a bivariate Q_i.  Also return vi
#'
#' @param yi n_i-response vector
#' @param xi n_i by p design matrix for fixed effects
#' @param zi n_i by 2 design matric for random effects (intercept and slope)
#' @param wi the pre-multiplier of yi to generate the sampling variable q_i
#' @param beta mean model parameter p-vector
#' @param sigma0 std dev of the random intercept distribution
#' @param sigma1 std dev of the random slope distribution
#' @param rho correlation between the random intercept and slope
#' @param sigmae std dev of the measurement error distribution
#' @param cutpoints cutpoints defining the sampling regions. (a vector of length 4 c(xlow, xhigh, ylow, yhigh))
#' @param SampProb Sampling probabilities from within each region (vector of length 2 c(central region, outlying region)).
#' @return log transformed ascertainment correction
#' @export
logACi2q <- function(yi, xi, zi, wi, beta, sigma0, sigma1, rho, sigmae, cutpoints, SampProb){
    vi      <- vi.calc(zi, sigma0, sigma1, rho, sigmae)
    mu      <- xi %*% beta
    mu_q    <- as.vector(wi %*% mu)
    sigma_q <- wi %*% vi %*% t(wi)
    sigma_q[2,1] <- sigma_q[1,2]
    return(list(vi=vi, logACi= log( ACi2q(cutpoints=cutpoints, SampProb=SampProb, mu_q=mu_q, sigma_q=sigma_q))))
}


#' Calculate a subject-specific contribution to a log-likelihood for longitudinal normal data
#'
#' Calculate a subject-specific contribution to a log-likelihood for longitudinal normal data
#' @param yi n_i-response vector
#' @param xi n_i by p design matrix for fixed effects
#' @param beta mean model parameter vector
#' @param vi the variance covariance matrix (ZDZ+Sige2*I)
#' @return subject specific contribution to the log-likelihood
#' @export
#'
li.lme <- function(yi, xi, beta, vi){
    resid <- yi - xi %*% beta
    -(1/2) * (length(xi[,1])*log(2*pi) + log(det(vi)) + t(resid) %*% solve(vi) %*% resid )[1,1]
}

#' Calculate the conditional likelihood for the univariate and bivariate sampling cases across all subjects (Keep.liC=FALSE) or the subject specific contributions to the conditional likelihood along with the log-transformed ascertainment correction for multiple imputation (Keep.liC=TRUE).
#'
#' Calculate the conditional likelihood for the univariate and bivariate sampling cases across all subjects (Keep.liC=FALSE) or the subject specific contributions to the conditional likelihood along with the log-transformed ascertainment correction for multiple imputation (Keep.liC=TRUE).
#'
#' @param y response vector
#' @param x sum(n_i) by p design matrix for fixed effects
#' @param z sum(n_i) by 2 design matric for random effects (intercept and slope)
#' @param w.function sum(n_i) vector with possible values that include "mean" "intercept" "slope" and "bivar."  There should be one unique value per subject
#' @param id sum(n_i) vector of subject ids
#' @param beta mean model parameter p-vector
#' @param sigma0 std dev of the random intercept distribution
#' @param sigma1 std dev of the random slope distribution
#' @param rho correlation between the random intercept and slope
#' @param sigmae std dev of the measurement error distribution
#' @param cutpoints A matrix with the first dimension equal to sum(n_i).  These cutpoints define the sampling regions [bivariate Q_i: each row is a vector of length 4 c(xlow, xhigh, ylow, yhigh); univariate Q_i: each row is a vector of length 2 c(k1,k2) to define the sampling regions, i.e., low, middle, high].  Each subject should have n_i rows of the same values.
#' @param SampProb A matrix with the first dimension equal to sum(n_i).   Sampling probabilities from within each region [bivariate Q_i: each row is a vector of length 2 c(central region, outlying region); univariate Q_i: each row is a vector of length 3 with sampling probabilities for each region]. Each subject should have n_i rows of the same values.
#' @param SampProbi Subject specific sampling probabilities.  A vector of length sum(n_i).  Not used unless using weighted Likelihood
#' @param Keep.liC If FALSE, the function returns the conditional log likelihood across all subjects.  If TRUE, subject specific contributions and exponentiated subject specific ascertainment corrections are returned in a list.
#' @return If Keep.liC=FALSE, conditional log likelihood.  If Keep.liC=TRUE, a two-element list that contains subject specific likelihood contributions and exponentiated ascertainment corrections.
#' @export
#'
LogLikeC2 <- function(y, x, z, w.function, id, beta, sigma0, sigma1, rho, sigmae, cutpoints, SampProb, SampProbi, Keep.liC=FALSE){

    subjectData <- CreateSubjectData(id=id,y=y,x=x,z=z,SampProbi=SampProbi,SampProb=SampProb,cutpoints=cutpoints,w.function=w.function)
    
    # id.tmp        <- split(id,id)
    # y.tmp         <- split(y,id)
    # x.tmp         <- split(x,id)
    # z.tmp         <- split(z,id)
    # SampProbi.tmp <- split(SampProbi,id)
    # 
    # subjectData <- vector('list', length=length(unique(id)))
    # subjectData <- list()
    # uid <- as.character(unique(id))
    # for(j in seq(along=uid)){
    #     i <- uid[j]
    #     subjectData[[j]] <- list(idi=as.character(unique(id.tmp[[i]])),
    #                              xi=matrix(x.tmp[[i]], ncol=ncol(x)),
    #                              zi=matrix(z.tmp[[i]], ncol=ncol(z)),
    #                              yi=y.tmp[[i]],
    #                              SampProbi.i=SampProbi.tmp[[i]])
    # }
    # names(subjectData) <- uid
    liC.and.logACi <- lapply(subjectData, LogLikeiC2, beta=beta, sigma0=sigma0, sigma1 = sigma1, rho = rho, sigmae = sigmae)

    if (Keep.liC == FALSE){out <- -1*Reduce('+', liC.and.logACi)[1]  ## sum ss contributions to liC
    }else{ out <- list(liC    = c(unlist(sapply(liC.and.logACi, function(x) x[1]))), ## ss contributions to liC
                       logACi = c(unlist(sapply(liC.and.logACi, function(x) x[2]))))} ## ss contributions to ACi
    out
}



#' Calculate the ss contributions to the conditional likelihood for the univariate and bivariate sampling cases.
#'
#' Calculate the ss contributions to the conditional likelihood for the univariate and bivariate sampling cases.
#'
#' @param subjectData a list containing: yi, xi, zi, SampProbi.i
#' @param w.function options include "mean" "intercept" "slope" and "bivar"
#' @param beta mean model parameter p-vector
#' @param sigma0 std dev of the random intercept distribution
#' @param sigma1 std dev of the random slope distribution
#' @param rho correlation between the random intercept and slope
#' @param sigmae std dev of the measurement error distribution
#' @param cutpoints cutpoints defining the sampling regions. (a vector of length 4 c(xlow, xhigh, ylow, yhigh))
#' @param SampProb Sampling probabilities from within each region (vector of length 2 c(central region, outlying region)).
#' @return ss contributions to the conditional log likelihood.  This is an internal function used by LogLikeC2
#' @export
#'
#'
LogLikeiC2 <- function(subjectData, beta, sigma0, sigma1, rho, sigmae){
        yi          <- subjectData[["yi"]]
        xi          <- subjectData[["xi"]]
        zi          <- subjectData[["zi"]]
        SampProbi.i <- subjectData[["SampProbi.i"]]
        w.function  <- subjectData[["w.function.i"]]
        SampProb    <- subjectData[["SampProb.i"]]
        cutpoints   <- subjectData[["cutpoints.i"]]
        ni          <- length(yi)
        t.zi        <- t(zi)
        if (w.function != "bivar"){
            if (w.function=="mean")      wi <- t(rep(1/ni, ni))
            if (w.function=="intercept") wi<- (solve(t.zi %*% zi) %*% t.zi)[1,]
            if (w.function=="slope")     wi<- (solve(t.zi %*% zi) %*% t.zi)[2,]
            wi         <- matrix(wi, 1, ni)
            IPWi       <- 1/ unique(SampProbi.i)
            tmp        <- logACi1q(yi, xi, zi, wi, beta, sigma0, sigma1, rho, sigmae, cutpoints, SampProb)
            logACi     <- tmp[["logACi"]]
            liC        <- li.lme(yi, xi, beta, tmp[["vi"]])*IPWi - logACi

        }else{
            wi         <- solve(t.zi %*% zi) %*% t.zi
            IPWi       <- 1/ unique(SampProbi.i)
            tmp        <- logACi2q(yi, xi, zi, wi, beta, sigma0, sigma1, rho, sigmae, cutpoints, SampProb)
            logACi     <- tmp[["logACi"]]
            liC        <- li.lme(yi, xi, beta, tmp[["vi"]])*IPWi - logACi

        }
        return(c(liC, logACi))
}

#' Gradient of the log of the ascertainment correction piece for sampling based on bivariate Q_i
#'
#' Calculate the gradient of the log transformed ascertainment correction under designs that sample based on a bivariate Q_i (numerically)
#'
#' @param subjectData a list containing: yi, xi, zi
#' @param w.function Sampling variable q_i function "mean", "intercept", "slope", "bivar".  This only gets called if w.function="bivar"
#' @param beta mean model parameter p-vector
#' @param sigma0 std dev of the random intercept distribution
#' @param sigma1 std dev of the random slope distribution
#' @param rho correlation between the random intercept and slope
#' @param sigmae std dev of the measurement error distribution
#' @param cutpoints cutpoints defining the sampling regions. (a vector of length 4 c(xlow, xhigh, ylow, yhigh))
#' @param SampProb Sampling probabilities from within each region (vector of length 2 c(central region, outlying region)).
#' @return gradient of the log transformed ascertainment correction under the bivariate sampling design
#' @importFrom numDeriv grad
#' @importFrom mvtnorm pmvnorm
#' @export
logACi2q.score2 <- function(subjectData, beta, sigma0, sigma1, rho, sigmae){

    yi          <- subjectData[["yi"]]
    xi          <- subjectData[["xi"]]
    zi          <- subjectData[["zi"]]
    SampProbi.i <- subjectData[["SampProbi.i"]]  ## not used
    w.function  <- subjectData[["w.function.i"]]
    SampProb    <- subjectData[["SampProb.i"]]
    cutpoints   <- subjectData[["cutpoints.i"]]
    
    # yi          <- subjectData[["yi"]]
    # xi          <- subjectData[["xi"]]
    # zi          <- subjectData[["zi"]]
    t.zi        <- t(zi)
    wi          <- solve(t.zi %*% zi) %*% t.zi
    t.wi        <- t(wi)

    param   <- c(beta, sigma0, sigma1, rho, sigmae)
    npar    <- length(param)
    Deriv <- sapply(1:npar,  function(rr)
                             {
                                grad(function(x) { new.param <- param
                                                   new.param[rr] <- x
                                                   vi      <- vi.calc(zi, new.param[(npar-3)], new.param[(npar-2)], new.param[(npar-1)], new.param[npar])
                                                   mu_q    <- as.vector(wi %*% (xi %*% new.param[1:(npar-4)]))
                                                   sigma_q <- wi %*% vi %*% t.wi
                                                   ## sometimes the off diagonals different in the 7th or 8th decimal place.  This seeks to fix that.  Not sure why it happened.
                                                   sigma_q[2,1] <- (sigma_q[2,1] + sigma_q[1,2]) / 2
                                                   sigma_q[1,2] <- sigma_q[2,1]
                                                   pmvnorm(lower=c(cutpoints[c(1,3)]), upper=c(cutpoints[c(2,4)]), mean=mu_q, sigma=sigma_q)[[1]]
                                                 },
                                     param[rr],
                                     method="simple")
                             }
    )

    vi      <- vi.calc(zi, sigma0, sigma1, rho, sigmae)
    mu_q    <- as.vector(wi %*% (xi %*% beta))
    sigma_q <- wi %*% vi %*% t.wi
    sigma_q[2,1] <- (sigma_q[2,1] + sigma_q[1,2]) / 2
    sigma_q[1,2] <- sigma_q[2,1]

    (SampProb[1]-SampProb[2])*Deriv / ACi2q(cutpoints, SampProb, mu_q, sigma_q)
}


#' Gradient of the log of the ascertainment correction piece for sampling based on univariate Q_i
#'
#' Calculate the gradient of the log transformed ascertainment correction for sampling based on univariate Q_i
#'
#' @param subjectData a list containing: yi, xi, zi
#' @param w.function Sampling variable q_i function "mean", "intercept", "slope", "bivar".  This only gets called if w.function in mean, slope, intercept
#' @param beta mean model parameter p-vector
#' @param sigma0 std dev of the random intercept distribution
#' @param sigma1 std dev of the random slope distribution
#' @param rho correlation between the random intercept and slope
#' @param sigmae std dev of the measurement error distribution
#' @param cutpoints cutpoints defining the sampling regions. (a vector of length 2 to define low, medium and high values of $Q_i$).
#' @param SampProb Sampling probabilities from within each region (vector of length 3 to define sampling probabilities within sampling regions
#' @return gradient of the log transformed ascertainment correction under univariate $Q_i$
#' @export
#' @importFrom stats dnorm
logACi1q.score2 <- function(subjectData, beta, sigma0, sigma1, rho, sigmae){

    yi          <- subjectData[["yi"]]
    xi          <- subjectData[["xi"]]
    zi          <- subjectData[["zi"]]
    SampProbi.i <- subjectData[["SampProbi.i"]]  ## not used
    w.function  <- subjectData[["w.function.i"]]
    SampProb    <- subjectData[["SampProb.i"]]
    cutpoints   <- subjectData[["cutpoints.i"]]
    
    #yi          <- subjectData[["yi"]]
    #xi          <- subjectData[["xi"]]
    #zi          <- subjectData[["zi"]]
    t.zi        <- t(zi)
    ni          <- length(yi)
    if (w.function=="mean")      wi <- t(rep(1/ni, ni))
    if (w.function=="intercept") wi<- (solve(t.zi %*% zi) %*% t.zi)[1,]
    if (w.function=="slope")     wi<- (solve(t.zi %*% zi) %*% t.zi)[2,]
    wi      <- matrix(wi, 1, ni)

    vi        <- vi.calc(zi, sigma0, sigma1, rho, sigmae)
    t.wi      <- t(wi)
    wi.zi     <- wi %*% zi
    t.wi.zi   <- t(wi.zi)
    rho.sig0  <- rho*sigma0
    rho.sig1  <- rho*sigma1
    sig0.sig1 <- sigma0*sigma1

    mu      <- xi %*% beta
    mu_q    <- (wi %*% mu)[,1]
    sigma_q <- sqrt((wi %*% vi %*% t.wi)[1,1])

    l <- ACi1q(cutpoints, SampProb, mu_q, sigma_q)
    p <- SampProb[1:(length(SampProb)-1)] - SampProb[2:(length(SampProb))]
    f <- dnorm(cutpoints, mu_q, sigma_q)

    d_li_beta <- (wi %*% xi) * sum(p*f) / l

    #f_alpha_k <- sum(p*f*(cutpoints - mu_q)) / (l * sigma_q *2 * sigma_q)
    f_alpha_k <- sum(p*f*(cutpoints - mu_q)) / (l * 2* sigma_q^2 )

    a5        <- (wi.zi %*% matrix(c(2*sigma0, rho.sig1,  rho.sig1,  0),        nrow=2) %*% t.wi.zi)[1,1]
    a6        <- (wi.zi %*% matrix(c(0,        rho.sig0,  rho.sig0,  2*sigma1), nrow=2) %*% t.wi.zi)[1,1]
    a7        <- (wi.zi %*% matrix(c(0,        sig0.sig1, sig0.sig1, 0),        nrow=2) %*% t.wi.zi)[1,1]
    a8        <- (wi %*% (2 * sigmae * diag(length(yi))) %*% t.wi)[1,1]
    c(d_li_beta, c(f_alpha_k * c(a5, a6, a7, a8)))
}

#' Subject specific contribution to the lme model score (also returns marginal Vi=Cov(Y|X))
#'
#' Subject specific contribution to the lme model score (also returns marginal Vi=Cov(Y|X))
#'
#' @param subjectData a list that contains yi, xi, zi, SampProbi.i.  Note that SampProbi.i is used for IPW only.
#' @param beta mean model parameter p-vector
#' @param sigma0 std dev of the random intercept distribution
#' @param sigma1 std dev of the random slope distribution
#' @param rho correlation between the random intercept and slope
#' @param sigmae std dev of the measurement error distribution
#' @return Subject specific contribution to the log-likelihood score (also returns marginal Vi=Cov(Y|X))
#' @export
li.lme.score <- function(subjectData, beta, sigma0, sigma1, rho, sigmae){
    yi          <- subjectData[["yi"]]
    xi          <- subjectData[["xi"]]
    zi          <- subjectData[["zi"]]
    t.zi        <- t(zi)
    IPW.i       <- 1/ unique(subjectData[["SampProbi.i"]])

    resid         <- yi - xi %*% beta
    vi            <- vi.calc(zi, sigma0, sigma1, rho, sigmae)
    inv.v         <- solve(vi)
    t.resid       <- t(resid)
    t.resid.inv.v <- t.resid %*% inv.v
    inv.v.resid   <- inv.v %*% resid

    rho.sig0  <- rho*sigma0
    rho.sig1  <- rho*sigma1
    sig0.sig1 <- sigma0*sigma1
    t.zi      <- t(zi)

    dvk5 <- zi %*% matrix(c(2*sigma0, rho.sig1,  rho.sig1,  0),nrow=2) %*% t.zi
    dvk6 <- zi %*% matrix(c(0,        rho.sig0,  rho.sig0,  2*sigma1),nrow=2) %*% t.zi
    dvk7 <- zi %*% matrix(c(0,        sig0.sig1, sig0.sig1, 0),nrow=2) %*% t.zi
    dvk8 <- 2 * sigmae  * diag(length(yi))

    l14 <- t(xi) %*% inv.v.resid # for Beta

    l5  <- -0.5*(sum(diag(inv.v %*% dvk5)) - t.resid.inv.v %*% dvk5 %*% inv.v.resid)[1,1]
    l6  <- -0.5*(sum(diag(inv.v %*% dvk6)) - t.resid.inv.v %*% dvk6 %*% inv.v.resid)[1,1]
    l7  <- -0.5*(sum(diag(inv.v %*% dvk7)) - t.resid.inv.v %*% dvk7 %*% inv.v.resid)[1,1]
    l8  <- -0.5*(sum(diag(inv.v %*% dvk8)) - t.resid.inv.v %*% dvk8 %*% inv.v.resid)[1,1]
    list(gr=append(l14, c(l5,l6,l7,l8))*IPW.i,
         vi=vi)
}

#' Calculate the gradient of the conditional likelihood for the univariate and bivariate sampling cases across all subjects (CheeseCalc=FALSE) or the cheese part of the sandwich estimator if CheeseCalc=TRUE.
#'
#' Calculate the gradient of the conditional likelihood for the univariate and bivariate sampling cases across all subjects (CheeseCalc=FALSE) or the cheese part of the sandwich estimator if CheeseCalc=TRUE.
#'
#' @param y response vector
#' @param x sum(n_i) by p design matrix for fixed effects
#' @param z sum(n_i) by 2 design matric for random effects (intercept and slope)
#' @param w.function sum(n_i) vector with possible values that include "mean" "intercept" "slope" and "bivar."  There should be one unique value per subject
#' @param id sum(n_i) vector of subject ids
#' @param beta mean model parameter p-vector
#' @param sigma0 std dev of the random intercept distribution
#' @param sigma1 std dev of the random slope distribution
#' @param rho correlation between the random intercept and slope
#' @param sigmae std dev of the measurement error distribution
#' @param cutpoints A matrix with the first dimension equal to sum(n_i).  These cutpoints define the sampling regions [bivariate Q_i: each row is a vector of length 4 c(xlow, xhigh, ylow, yhigh); univariate Q_i: each row is a vector of length 2 c(k1,k2) to define the sampling regions, i.e., low, middle, high].  Each subject should have n_i rows of the same values.
#' @param SampProb A matrix with the first dimension equal to sum(n_i).   Sampling probabilities from within each region [bivariate Q_i: each row is a vector of length 2 c(central region, outlying region); univariate Q_i: each row is a vector of length 3 with sampling probabilities for each region]. Each subject should have n_i rows of the same values.
#' @param SampProbi Subject specific sampling probabilities.  A vector of length sum(n_i).  Not used unless using weighted Likelihood
#' @param CheeseCalc If FALSE, the function returns the gradient of the conditional log likelihood across all subjects.  If TRUE, the cheese part of the sandwich esitmator is calculated.
#' @return If CheeseCalc=FALSE, gradient of conditional log likelihood.  If CheeseCalc=TRUE, the cheese part of the sandwich estimator is calculated.
#' @export
LogLikeC.Score2 <- function(y, x, z, w.function, id, beta, sigma0, sigma1, rho, sigmae, cutpoints, SampProb, SampProbi, CheeseCalc=FALSE){
    param.vec <- c(beta, log(sigma0),log(sigma1),log((1+rho)/(1-rho)),log(sigmae))
    n.par     <- length(param.vec)

    subjectData <- CreateSubjectData(id=id,y=y,x=x,z=z,SampProbi=SampProbi,SampProb=SampProb,cutpoints=cutpoints,w.function=w.function)
    # id.tmp        <- split(id,id)
    # y.tmp         <- split(y,id)
    # x.tmp         <- split(x,id)
    # z.tmp         <- split(z,id)
    # SampProbi.tmp <- split(SampProbi,id)
    # subjectData   <- vector('list', length=length(unique(id)))
    # subjectData   <- list()
    # uid           <- as.character(unique(id))
    # for(j in seq(along=uid)){
    #     i <- uid[j]
    #     subjectData[[j]] <- list(idi=as.character(unique(id.tmp[[i]])),
    #                              xi=matrix(x.tmp[[i]], ncol=ncol(x)),
    #                              zi=matrix(z.tmp[[i]], ncol=ncol(z)),
    #                              yi=y.tmp[[i]],
    #                              SampProbi.i=SampProbi.tmp[[i]])
    # }
    # names(subjectData) <- uid
    UncorrectedScorei <- lapply(subjectData, li.lme.score, beta=beta, sigma0=sigma0, sigma1=sigma1, rho=rho, sigmae=sigmae)
    Gradienti         <- lapply(UncorrectedScorei, function(x) x[['gr']]) ## create a list of ss contributions to gradient
    UncorrectedScore  <- Reduce('+', Gradienti)  ## Note if using IPW this is actually a corrected score (corrected by the IPW)

    ## NOTE HERE: I used the first element of w.function in this call.  This means, for now, we cannot mix bivar with other
    ## sampling schemes.  This also applies to the cheese calculation
    if (w.function[[1]] != "bivar"){
        logACi.Score <- lapply(subjectData, logACi1q.score2, beta=beta, sigma0=sigma0, sigma1=sigma1, rho=rho, sigmae=sigmae)
        logAC.Score  <- Reduce('+', logACi.Score)
        CorrectedScore <- UncorrectedScore + logAC.Score
    }else{
        logACi.Score <- lapply(subjectData, logACi2q.score2, beta=beta, sigma0=sigma0, sigma1=sigma1, rho=rho, sigmae=sigmae)
        logAC.Score  <- Reduce('+', logACi.Score)
        CorrectedScore <- UncorrectedScore - logAC.Score  ## notice this has the opposite sign compared to above.  Remember to check
    }

    if (CheeseCalc==TRUE){
        if (w.function[[1]] != "bivar"){ GradiMat <- mapply("+", Gradienti, logACi.Score)
        }else{                      GradiMat <- mapply("-", Gradienti, logACi.Score)} ## notice this has the opposite sign compared to above.  Remember to check
        ## Need to use the chain rule: note that param,vec is on the unconstrained scale but Gradi was calculated on the constrained parameters
        GradiMat[c((n.par-3):n.par),] <- GradiMat[c((n.par-3):n.par),]*c(exp(param.vec[(n.par-3)]),
                                                                         exp(param.vec[(n.par-2)]),
                                                                         2*exp(param.vec[(n.par-1)])/((exp(param.vec[(n.par-1)])+1)^2),
                                                                         exp(param.vec[(n.par)]))
        cheese <- matrix(0,  n.par, n.par)
        for (mm in 1:ncol(GradiMat)) cheese <- cheese + outer(GradiMat[,mm], GradiMat[,mm])
    }
    out <- -CorrectedScore
    if (CheeseCalc==TRUE) out <- cheese
    out
}

# # ###################################
# y=odsInt$Y
# x=as.matrix(cbind(1, odsInt[,c("time","snp","snptime","confounder")]))
# z=as.matrix(cbind(1, odsInt$time))
# id=odsInt$id
# beta=c(5,  1, -2.5,  0.75,  0)
# sigma0= 1.6094379
# sigma1 = 0.2231436
# rho = -0.5108256
# sigmae = 1.6094379
# cutpoints=c(-2.569621, 9.992718)
# SampProb=c(1, 0.1228, 1)
# w.function="intercept"
# SampProbi = rep(1, length(y))
#
# tmp1 <- gradient.nll.lme(y, x, z, w.function, id, beta, sigma0, sigma1, rho, sigmae, cutpoints, SampProb, SampProbi, CheeseCalc=FALSE)
# tmp2 <- LogLikeC.Score2(y, x, z, w.function, id, beta, sigma0, sigma1, rho, sigmae, cutpoints, SampProb, SampProbi, CheeseCalc=FALSE)
# tmp1 <- gradient.nll.lme(y, x, z, w.function, id, beta, sigma0, sigma1, rho, sigmae, cutpoints, SampProb, SampProbi, CheeseCalc=TRUE)
# tmp2 <- LogLikeC.Score2(y, x, z, w.function, id, beta, sigma0, sigma1, rho, sigmae, cutpoints, SampProb, SampProbi, CheeseCalc=TRUE)
# tmp1
# tmp2
# #
# y=odsSlp$Y
# x=as.matrix(cbind(1, odsSlp[,c("time","snp","snptime","confounder")]))
# z=as.matrix(cbind(1, odsSlp$time))
# id=odsSlp$id
# beta=c(5,  1, -2.5,  0.75,  0)
# sigma0= 1.6094379
# sigma1 = 0.2231436
# rho = -0.5108256
# sigmae = 1.6094379
# cutpoints=c(-0.7488912,  3.4557775)
# SampProb=c(1, 0.1228, 1)
# w.function="slope"
# SampProbi = rep(1, length(y))
#
# tmp1 <- gradient.nll.lme(y, x, z, w.function, id, beta, sigma0, sigma1, rho, sigmae, cutpoints, SampProb, SampProbi, CheeseCalc=TRUE)
# tmp2 <- LogLikeC.Score2(y, x, z, w.function, id, beta, sigma0, sigma1, rho, sigmae, cutpoints, SampProb, SampProbi, CheeseCalc=TRUE)
# tmp1
# tmp2
# #
# y=odsBiv$Y
# x=as.matrix(cbind(1, odsBiv[,c("time","snp","snptime","confounder")]))
# z=as.matrix(cbind(1, odsBiv$time))
# id=odsBiv$id
# beta=c(5,  1, -2.5,  0.75,  0)
# sigma0= 1.6094379
# sigma1 = 0.2231436
# rho = -0.5108256
# sigmae = 1.6094379
# cutpoints=c(-4.413225, 11.935188, -1.390172, 4.084768)
# SampProb=c(0.122807, 1)
# w.function="bivar"
# SampProbi = rep(1, length(y))
#
# date()
# tmp1 <- gradient.nll.lme(y, x, z, w.function, id, beta, sigma0, sigma1, rho, sigmae, cutpoints, SampProb, SampProbi, CheeseCalc=TRUE)
# date()
# tmp2 <- LogLikeC.Score2(y, x, z, w.function, id, beta, sigma0, sigma1, rho, sigmae, cutpoints, SampProb, SampProbi, CheeseCalc=TRUE)
# date()
# tmp1
# tmp2


#' Calculate the ascertainment corrected log likelihood and score
#'
#' Calculate the ascertainment corrected log likelihood and score
#'
#' @param params parameter vector c(beta, log(sigma0), log(sigma1), rho, sigmae)
#' @param y response vector
#' @param x sum(n_i) by p design matrix for fixed effects
#' @param z sum(n_i) by 2 design matric for random effects (intercept and slope)
#' @param id sum(n_i) vector of subject ids
#' @param w.function sum(n_i) vector with possible values that include "mean" "intercept" "slope" and "bivar."  There should be one unique value per subject
#' @param cutpoints A matrix with the first dimension equal to sum(n_i).  These cutpoints define the sampling regions [bivariate Q_i: each row is a vector of length 4 c(xlow, xhigh, ylow, yhigh); univariate Q_i: each row is a vector of length 2 c(k1,k2) to define the sampling regions, i.e., low, middle, high].  Each subject should have n_i rows of the same values.
#' @param SampProb A matrix with the first dimension equal to sum(n_i).   Sampling probabilities from within each region [bivariate Q_i: each row is a vector of length 2 c(central region, outlying region); univariate Q_i: each row is a vector of length 3 with sampling probabilities for each region]. Each subject should have n_i rows of the same values.
#' @param SampProbi Subject specific sampling probabilities.  A vector of length sum(n_i).  Not used unless using weighted Likelihood
#' @param ProfileCol the column number(s) for which we want fixed at the value of param.  Maimizing the log likelihood for all other parameters
#'                   while fixing these columns at the values of params[ProfileCol]
#' @param Keep.liC If TRUE outputs subject specific conditional log lileihoods to be used for the imputation procedure described in the AOAS paper keep z sum(n_i) by 2 design matric for random effects (intercept and slope)
#' @return The conditional log likelihood with a "gradient" attribute (if Keep.liC=FALSE) and subject specific contributions to the conditional likelihood if Keep.liC=TRUE).
#' @export
LogLikeCAndScore2 <- function(params, y, x, z, id, w.function, cutpoints, SampProb, SampProbi, ProfileCol=NA, Keep.liC=FALSE){
    npar   <- length(params)
    beta   <- params[1:(npar-4)]
    sigma0 <- exp(params[(npar-3)])
    sigma1 <- exp(params[(npar-2)])
    rho    <- (exp(params[(npar-1)])-1) / (exp(params[(npar-1)])+1)
    sigmae <- exp(params[npar])

    out     <- LogLikeC2( y=y, x=x, z=z, w.function=w.function, id=id, beta=beta, sigma0=sigma0,
                                sigma1=sigma1, rho=rho, sigmae=sigmae, cutpoints=cutpoints, SampProb=SampProb, SampProbi=SampProbi,
                                Keep.liC=Keep.liC)
    GRAD    <- LogLikeC.Score2(y=y, x=x, z=z, w.function=w.function, id=id, beta=beta, sigma0=sigma0,
                                sigma1=sigma1, rho=rho, sigmae=sigmae, cutpoints=cutpoints, SampProb=SampProb, SampProbi=SampProbi)

    ## Need to use the chain rule: note that params is on the unconstrained scale but GRAD was calculated on the constrained parameters
    GRAD[(npar-3)] <- GRAD[(npar-3)]*exp(params[(npar-3)])
    GRAD[(npar-2)] <- GRAD[(npar-2)]*exp(params[(npar-2)])
    GRAD[(npar-1)] <- GRAD[(npar-1)]*2*exp(params[(npar-1)])/((exp(params[(npar-1)])+1)^2)
    GRAD[npar]     <- GRAD[npar]*exp(params[npar])
    ## Force the gradient of the fixed parameter to be zero, so that it does not move
    if (!is.na(ProfileCol)) GRAD[ProfileCol] <- 0
    attr(out,"gradient") <- GRAD

    out
}


## If you do not want to use the ascertainment correction term in the conditional likelihood
## set all SampProb values equal to each other.  This would be the case if you were doing
## straightforward maximum likelihood (albeit computationally inefficient) or weighted likelihood.

#' Fitting function: ACML for a linear mixed effects model (random intercept and slope)
#'
#' Fitting function: ACML or WL for a linear mixed effects model (random intercept and slope)
#'
#' @param formula.fixed formula for the fixed effects (of the form y~x)
#' @param formula.random formula for the random effects (of the form ~z)
#' @param data data frame (should contain everything in formula.fixed, formula.random, id, and SampProbiWL)
#' @param id sum(n_i) vector of subject ids
#' @param w.function sum(n_i) vector with possible values that include "mean" "intercept" "slope" and "bivar."  There should be one unique value per subject
#' @param InitVals starting values for c(beta, log(sigma0), log(sigma1), rho, log(sigmae))
#' @param cutpoints A matrix with the first dimension equal to sum(n_i).  These cutpoints define the sampling regions [bivariate Q_i: each row is a vector of length 4 c(xlow, xhigh, ylow, yhigh); univariate Q_i: each row is a vector of length 2 c(k1,k2) to define the sampling regions, i.e., low, middle, high].  Each subject should have n_i rows of the same values.
#' @param SampProb A matrix with the first dimension equal to sum(n_i).   Sampling probabilities from within each region [bivariate Q_i: each row is a vector of length 2 c(central region, outlying region); univariate Q_i: each row is a vector of length 3 with sampling probabilities for each region]. Each subject should have n_i rows of the same values.
#' @param SampProbiWL Subject specific sampling probabilities.  A vector of length sum(n_i).  Not used unless using weighted Likelihood.  To only be used if doing IPWL.  Note if doing IPWL, only use robcov (robust variances) and not covar.  This should be contained in data.
#' @param ProfileCol the column number(s) for which we want fixed at the value of param.  Maximizing the log likelihood for all other parameters
#'                   while fixing these columns at the values of params[ProfileCol]
#' @return Ascertainment corrected Maximum likelihood: Ests, covar, LogL, code, robcov
#' @importFrom stats model.frame
#' @importFrom stats model.matrix
#' @importFrom stats model.response
#' @importFrom stats na.omit
#' @importFrom stats nlm
#' @importFrom stats pnorm
#' @export
acml.lmem2 <- function(formula.fixed, 
                          formula.random, 
                          data, 
                          id,   
                          w.function="mean",           
                          InitVals,                    
                          cutpoints = c(0,5),         
                          SampProb = NA,             
                          SampProbiWL=NA,             
                          ProfileCol=NA){              ## Columns to be held fixed while doing profile likelihood.  It is fixed at its initial value.
    
    if(is.null(formula.random)) {stop('Specify the random effects portion of the model.  It is currently NULL.')}
    if(!is.data.frame(data)) {
        data <- as.data.frame(data)
        warning('data converted to data.frame.')
    }
    terms = unique( c(all.vars(formula.fixed), all.vars(formula.random),as.character(substitute(id)), as.character(substitute(SampProbiWL))) )
    data  = data[,terms]
    
    if(any(is.na(data))) data = na.omit(data)

    id0   =  as.character(substitute(id))
    id    = data$id = data[ , id0 ]

    fixed.f = model.frame(formula.fixed, data)
    fixed.t = attr(fixed.f, "terms")
    y      = model.response(fixed.f,'numeric')
    uy     = unique(y)
    x      = model.matrix(formula.fixed, fixed.f)

    rand.f = model.frame(formula.random, data)
    z      = model.matrix(formula.random, fixed.f)

    #if (is.na(SampProb[1])) SampProb = c(1,1,1)
    SampProbi0   =  as.character(substitute(SampProbiWL))
    SampProbi    = data$SampProbi = data[ , SampProbi0 ]

    acml.fit <- nlm(LogLikeCAndScore2, InitVals, y=y, x=x, z=z, id=id, w.function=w.function, cutpoints=cutpoints, SampProb=SampProb,
               SampProbi=SampProbi, ProfileCol=ProfileCol,
               stepmax=4, iterlim=250, check.analyticals = TRUE, print.level=0)

    ## Calculate the observed information and then invert to get the covariance matrix
    npar <- length(acml.fit$estimate)
    Hessian.eps <- 1e-7
    eps.mtx     <- diag(rep(Hessian.eps, npar))
    grad.at.max <- acml.fit$gradient
    ObsInfo.tmp <- ObsInfo <- matrix(NA, npar, npar)

    ## Observed Information
    for (j in 1:npar){
        temp        <- LogLikeCAndScore2(acml.fit$estimate+eps.mtx[j,], y=y, x=x, z=z, id=id,
                                        w.function=w.function, cutpoints=cutpoints,
                                        SampProb=SampProb,SampProbi=SampProbi, ProfileCol=ProfileCol)
        ObsInfo.tmp[j,] <- (attr(temp,"gradient")-grad.at.max)/(Hessian.eps)
    }
    for (m in 1:npar){
        for (n in 1:npar){ ObsInfo[m,n] <-  (ObsInfo.tmp[m,n]+ObsInfo.tmp[n,m])/2}}
    ## Cheese part of the sandwich estimator
    Cheese <- LogLikeC.Score2(y=y, x=x, z=z, w.function=w.function,
                             id=id, beta=acml.fit$estimate[c(1:(npar-4))],
                             sigma0=exp(acml.fit$estimate[(npar-3)]),
                             sigma1=exp(acml.fit$estimate[(npar-2)]),
                             rho=   (exp(acml.fit$estimate[(npar-1)])-1) / (exp(acml.fit$estimate[(npar-1)])+1),
                             sigmae=exp(acml.fit$estimate[npar]),
                             cutpoints=cutpoints,
                             SampProb=SampProb,
                             SampProbi=SampProbi,
                             CheeseCalc=TRUE)

    if (!is.na(ProfileCol)){
        acml.fit$estimate <- acml.fit$estimate[-ProfileCol]
        ObsInfo <- ObsInfo[-ProfileCol, -ProfileCol]
        Cheese  <- Cheese[-ProfileCol, -ProfileCol]
    }


    out              <- NULL
    out$call         <- match.call()
    out$coefficients <- acml.fit$estimate
    out$covariance   <- solve(ObsInfo)
    out$robcov       <- solve(ObsInfo)%*%Cheese%*%solve(ObsInfo)
    out$logLik       <- -acml.fit$minimum
    attr(out,'args') <- list(formula.fixed=formula.fixed,
                             formula.random=formula.random,
                             id=id,
                             w.function=w.function,
                             cutpoints = cutpoints,
                             SampProb = SampProb,
                             SampProbiWL=SampProbi,
                             SampProbiVar = as.character(substitute(SampProbiWL)),
                             ProfileCol=ProfileCol)
    if(kappa(out$covar) > 1e5) warning("Poorly Conditioned Model")
    out
}

#' Create a list of subject-specific data
#'
#' @param id sum(n_i) vector of subject ids
#' @param y response vector
#' @param x sum(n_i) by p design matrix for fixed effects
#' @param z sum(n_i) by 2 design matric for random effects (intercept and slope)
#' @param SampProbi Subject specific sampling probabilities.  A vector of length sum(n_i).  Not used unless using weighted Likelihood
#' @param SampProb A matrix with the first dimension equal to sum(n_i).   Sampling probabilities from within each region [bivariate Q_i: each row is a vector of length 2 c(central region, outlying region); univariate Q_i: each row is a vector of length 3 with sampling probabilities for each region]. Each subject should have n_i rows of the same values.
#' @param w.function sum(n_i) vector with possible values that include "mean" "intercept" "slope" and "bivar."  There should be one unique value per subject
#' @param cutpoints A matrix with the first dimension equal to sum(n_i).  These cutpoints define the sampling regions [bivariate Q_i: each row is a vector of length 4 c(xlow, xhigh, ylow, yhigh); univariate Q_i: each row is a vector of length 2 c(k1,k2) to define the sampling regions, i.e., low, middle, high].  Each subject should have n_i rows of the same values.
#' @export
CreateSubjectData <- function(id,y,x,z,SampProbi,SampProb,cutpoints,w.function){
    id.tmp        <- split(id,id)
    y.tmp         <- split(y,id)
    x.tmp         <- split(x,id)
    z.tmp         <- split(z,id)
    SampProbi.tmp <- split(SampProbi,id)
    SampProb.tmp  <- split(SampProb,id)
    cutpoints.tmp  <- split(cutpoints,id)
    w.function.tmp  <- split(w.function,id)
    
    subjectData <- vector('list', length=length(unique(id)))
    subjectData <- list()
    uid <- as.character(unique(id))
    for(j in seq(along=uid)){
        i <- uid[j]
        subjectData[[j]] <- list(idi          = as.character(unique(id.tmp[[i]])),
                                 xi           = matrix(x.tmp[[i]], ncol=ncol(x)),
                                 zi           = matrix(z.tmp[[i]], ncol=ncol(z)),
                                 yi           = y.tmp[[i]],
                                 SampProbi.i  = SampProbi.tmp[[i]],
                                 SampProb.i   = matrix(SampProb.tmp[[i]], ncol=ncol(SampProb))[1,],
                                 w.function.i = as.character(unique(w.function.tmp[[i]])),
                                 cutpoints.i  = matrix(cutpoints.tmp[[i]], ncol=ncol(cutpoints))[1,])
    }
    names(subjectData) <- uid
    subjectData
}