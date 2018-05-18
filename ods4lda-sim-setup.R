   #############################################################################
 ##
## Known True Parameters
   
N           <- 1600
prev.grp    <- 0.25
conf.param  <- c(-0.5, 1)
n.imp       <- 75
p.central   <- 0.8

ni.sim               <- rbind(c(4, 6), 
                              c(4, 6), 
                              c(7,13), 
                              c(7,13))
inits.sim            <- rbind(c(75, -1, -1, -1, 1, log(9),log(1),0,log(3.5)),
                              c(75, -1, -1, -5, 1, log(9),log(1),0,log(3.5)),
                              c(75, -1, -1, -1, 1, log(9),log(1),0,log(3.5)),
                              c(75, -1, -1, -5, 1, log(9),log(1),0,log(3.5)))
NsPerStratumUniv.sim <- matrix( rep(c(150,100,150),4), ncol=3, byrow=TRUE)
NsPerStratumBiv.sim  <- matrix( rep(c(100,300),4), ncol=2, byrow=TRUE)
NsRand.sim           <- rep(400,4)

# params <- NsCentral <- list()
# rs.est <- int.est <- slp.est <- biv.est <- list()
# rs.mi.est <- int.mi.est <- slp.mi.est <- biv.mi.est <- list()
# rs.mi2.est <- int.mi2.est <- slp.mi2.est <- biv.mi2.est <- list()
# int.wl.est <- slp.wl.est <- biv.wl.est <- list()
# rs.cov <- int.cov <- slp.cov <- biv.cov <- list()
# rs.mi.cov <- int.mi.cov <- slp.mi.cov <- biv.mi.cov <- list()
# rs.mi2.cov <- int.mi2.cov <- slp.mi2.cov <- biv.mi2.cov <- list()
# int.wl.cov <- slp.wl.cov <- biv.wl.cov <- list()