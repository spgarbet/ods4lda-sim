   #############################################################################
 ##
## Known True Parameters
N           <- 2000
prev.grp    <- c(0.3,
                 0.3,
                 0.3,
                 0.3)
   
n.imp       <- 50
p.central   <- c(0.8,
                 0.8,
                 0.8,
                 0.8)
conf.param.sim       <- rbind( c(-0.25, 0.5),
                               c(-0.25, 0.5),
                               c(-0.25, 0.5),
                               c(-0.25, 0.5))
ni.sim               <- rbind(c(4, 6), 
                              c(4, 6), 
                              c(4, 6), 
                              c(4, 6))
inits.sim            <- rbind(c(75, -1, -.5, -6, -.5, log(9),log(1.25),0,log(3.5)),
                              c(75, -1, -.5, -6, -.5, log(9),log(1.25),0,log(3.5)),
                              c(75, -1, -.5, -6, -.5, log(9),log(1.25),0,log(3.5)),
                              c(75, -1, -.5, -6, -.5, log(9),log(1.25),0,log(3.5)))
NsPerStratumUniv.sim <- matrix( c(150,100,150,
                                  150,100,150,
                                  150,100,150,
                                  150,100,150), ncol=3, byrow=TRUE)
# Control TwoPhase RanTao
hn_scale=c(1.5, 2, 0.5, 0.1)
