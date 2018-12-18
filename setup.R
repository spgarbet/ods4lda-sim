##################################
# For misspecification
dmisspecified <- function(x, alpha, beta) dchisq( (x-alpha) / beta, 1)/beta
ex.misspecified <- Vectorize(function(alpha, beta)
  integrate(function(x) inv.logit(x) * dmisspecified(x, alpha, beta),
            alpha,
            Inf)[1][[1]]
)
find.alpha <- function(beta, zeta)
{
  optimize(function(y) (ex.misspecified(y, beta) - zeta)^2, c(-100, 100))$minimum
}

   
  #############################################################################
 ##
## Known True Parameters
N                    <- 2000
prev.grp             <- c(0.3,
                          0.3,
                          0.3)
n.imp                <- 50
p.central            <- c(0.8, # 10th to 90th percentile
                          0.8,
                          0.8)
conf.param.sim       <- rbind( c(find.alpha(0.5, prev.grp[1]), 0.5),
                               c(find.alpha(1, prev.grp[1]), 1),
                               c(find.alpha(2, prev.grp[1]), 2))

ni.sim               <- rbind(c(4, 6), 
                              c(4, 6), 
                              c(4, 6))

# Ordering
# B_0, B_t, B_g, B_c, B_gt, log(s_b0), log(s_b1), F_z(rho), log(s_e)
inits.sim            <- rbind(c(75, -1, -0.5, -2, -0.5, log(9),log(1.25),0,log(3.5)),
                              c(75, -1, -0.5, -2, -0.5, log(9),log(1.25),0,log(3.5)),
                              c(75, -1, -0.5, -2, -0.5, log(9),log(1.25),0,log(3.5)))
                             # c(75, -0.5, -0.5, -1, -6, log(9),log(1.25),0,log(3.5)),
                             # c(75, -0.5, -0.5, -1, -2, log(9),log(1.25),0,log(3.5)),
                             # c(75, -0.5, -0.5, -1, -2, log(9),log(1.25),0,log(3.5)))
NsPerStratumUniv.sim <- matrix( c(150,100,150,
                                  150,100,150,
                                  150,100,150), ncol=3, byrow=TRUE)
# Control TwoPhase RanTao
hn_scale             <- rep(2, 3)

