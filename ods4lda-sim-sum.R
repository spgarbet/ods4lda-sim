# Summarizing

Cases <- rbind(c(0,log(9)), c(-1, log(9)),c(-5,log(9)), c(-1, log(4)) )
sqrt.diag <- function(x) sqrt(diag(x))

CalcSums <- function(est.parms, est.cov, params){

    avg.est<-var.est<-avg.var <- coverage <- NULL
    for (i in 1:1){

        scenario.num <- which(simplify2array(params)[4,]==Cases[i,1] & simplify2array(params)[4,]==Cases[i,2])
        avg.est <- rbind(avg.est, c(Cases[i,], apply(simplify2array(est.parms[scenario.num]), 1, mean)))
        var.est <- rbind(var.est, c(Cases[i,], apply(simplify2array(est.parms[scenario.num]), 1, var)))
        avg.var <- rbind(avg.var, c(Cases[i,], diag(apply(simplify2array(est.cov[scenario.num]), 1:2, mean))))
        #print(est.parms[scenario.num])
        pars <- simplify2array(params[scenario.num])
        ests <- simplify2array(est.parms[scenario.num])
        ses <- apply(simplify2array(est.cov[scenario.num]), 3, sqrt.diag)
        #print(ests)
        #print(ses)
        lci <- ests - qnorm(.975)*ses
        uci <- ests + qnorm(.975)*ses
        Covered <- (lci <= pars & uci>= pars)
        coverage <- rbind(coverage, c(Cases[i,], apply(Covered, 1, mean)))
    }
    out <- list(avg.est=avg.est, var.est=var.est, avg.var=avg.var, coverage=coverage)
    out
}


params4mi <- lapply(params, function(x) x[1:5])
rs.sum <- CalcSums(rs.est, rs.cov, params)
rs.mi.sum <- CalcSums(rs.mi.est, rs.mi.cov, params4mi)
rs.mi2.sum <- CalcSums(rs.mi2.est, rs.mi2.cov, params4mi)
int.sum <- CalcSums(int.est, int.cov, params)
int.wl.sum <- CalcSums(int.wl.est, int.wl.cov, params)
int.mi.sum <- CalcSums(int.mi.est, int.mi.cov, params4mi)
int.mi2.sum <- CalcSums(int.mi2.est, int.mi2.cov, params4mi)

slp.sum <- CalcSums(slp.est, slp.cov, params)
slp.wl.sum <- CalcSums(slp.wl.est, slp.wl.cov, params)
slp.mi.sum <- CalcSums(slp.mi.est, slp.mi.cov, params4mi)
slp.mi2.sum <- CalcSums(slp.mi2.est, slp.mi2.cov, params4mi)

biv.sum <- CalcSums(biv.est, biv.cov, params)
biv.wl.sum <- CalcSums(biv.wl.est, biv.wl.cov, params)
biv.mi.sum <- CalcSums(biv.mi.est, biv.mi.cov, params4mi)
biv.mi2.sum <- CalcSums(biv.mi2.est, biv.mi2.cov, params4mi)


round(rs.sum$var.est/int.sum$var.est, 4)[1:7]
round(rs.sum$var.est[1:7]/int.mi.sum$var.est, 4)
round(int.sum$var.est[1:7]/int.mi.sum$var.est, 4)

round(rbind(rs.sum$avg.est[1:7], rs.mi.sum$avg.est, rs.mi2.sum$avg.est,
      int.sum$avg.est[1:7],  int.wl.sum$avg.est[1:7], int.mi.sum$avg.est, int.mi2.sum$avg.est,
      slp.sum$avg.est[1:7],  slp.wl.sum$avg.est[1:7], slp.mi.sum$avg.est, slp.mi2.sum$avg.est,
      biv.sum$avg.est[1:7],  biv.wl.sum$avg.est[1:7], biv.mi.sum$avg.est, biv.mi2.sum$avg.est),3)

round(sqrt(rbind(rs.sum$avg.var[1:7], rs.mi.sum$avg.var, rs.mi2.sum$avg.var,
            int.sum$avg.var[1:7],  int.wl.sum$avg.var[1:7], int.mi.sum$avg.var, int.mi2.sum$avg.var,
            slp.sum$avg.var[1:7],  slp.wl.sum$avg.var[1:7], slp.mi.sum$avg.var, slp.mi2.sum$avg.var,
            biv.sum$avg.var[1:7],  biv.wl.sum$avg.var[1:7], biv.mi.sum$avg.var, biv.mi2.sum$avg.var)),3)

round(rbind(rs.sum$var.est[1:7], rs.mi.sum$var.est, rs.mi2.sum$var.est,
            int.sum$var.est[1:7],  int.wl.sum$var.est[1:7], int.mi.sum$var.est, int.mi2.sum$var.est,
            slp.sum$var.est[1:7],  slp.wl.sum$var.est[1:7], slp.mi.sum$var.est, slp.mi2.sum$var.est,
            biv.sum$var.est[1:7],  biv.wl.sum$var.est[1:7], biv.mi.sum$var.est, biv.mi2.sum$var.est),5)

rs.var.mat <- matrix(rep(rs.sum$var.est[1:7],15), byrow=TRUE, ncol=7)
round(rs.var.mat /
          rbind(rs.sum$var.est[1:7], rs.mi.sum$var.est, rs.mi2.sum$var.est,
            int.sum$var.est[1:7],  int.wl.sum$var.est[1:7], int.mi.sum$var.est, int.mi2.sum$var.est,
            slp.sum$var.est[1:7],  slp.wl.sum$var.est[1:7], slp.mi.sum$var.est, slp.mi2.sum$var.est,
            biv.sum$var.est[1:7],  biv.wl.sum$var.est[1:7], biv.mi.sum$var.est, biv.mi2.sum$var.est),3)

rs.var.mat <- matrix(rep(rs.sum$avg.var[1:7],15), byrow=TRUE, ncol=7)
round(rs.var.mat /
          rbind(rs.sum$avg.var[1:7], rs.mi.sum$avg.var, rs.mi2.sum$avg.var,
                int.sum$avg.var[1:7],  int.wl.sum$avg.var[1:7], int.mi.sum$avg.var, int.mi2.sum$avg.var,
                slp.sum$avg.var[1:7],  slp.wl.sum$avg.var[1:7], slp.mi.sum$avg.var, slp.mi2.sum$avg.var,
                biv.sum$avg.var[1:7],  biv.wl.sum$avg.var[1:7], biv.mi.sum$avg.var, biv.mi2.sum$avg.var),3)


