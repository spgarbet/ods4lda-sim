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
inits.sim            <- rbind(c(75, -1, -.5, -1, -.5, log(9), log(1.25), 0, log(3.5)),
                              c(75, -1, -.5, -4, -.5, log(9), log(1.25), 0, log(3.5)),
                              c(75, -1, -.5, -1, -.5, log(9), log(1.25), 0, log(3.5)),
                              c(75, -1, -.5, -4, -.5, log(9), log(1.25), 0, log(3.5)))
NsPerStratumUniv.sim <- matrix( rep(c(100,100,100),4), ncol=3, byrow=TRUE)
NsPerStratumBiv.sim  <- matrix( rep(c(100,200),4), ncol=2, byrow=TRUE)
NsRand.sim           <- rep(300,4)
