   #############################################################################
 ##
## Known True Parameters
N           <- 2000
prev.grp    <- 0.25
conf.param  <- c(-0.25, 0.5)
n.imp       <- 75
p.central   <- 0.8

ni.sim               <- rbind(c(4, 6), 
                              c(4, 6), 
                              c(4, 6), 
                              c(4, 6))
inits.sim            <- rbind(c(75, -1, -.5, -1, -.5, log(9),log(1.25),0,log(3.5)),
                              c(75, -1, -.5, -5, -.5, log(9),log(1.25),0,log(3.5)),
                              c(75, -1, -.5, -1, -.5, log(9),log(1.25),0,log(3.5)),
                              c(75, -1, -.5, -1, -.5, log(1.5),log(0.5),0,log(3.5)))
NsPerStratumUniv.sim <- matrix( c(125,125,125,
                                  125,125,125,
                                  175,175,175,
                                  125,125,125), ncol=3, byrow=TRUE)
NsPerStratumBiv.sim  <- matrix( c(125,250,
                                  125,250,
                                  175,350,
                                  125,250), ncol=2, byrow=TRUE)
NsRand.sim           <- c(375,375,525,375)
