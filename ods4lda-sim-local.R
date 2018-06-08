source("ods4lda-sim-run.R")

  #############################################################################
 ##
## For exection on local desktop
library(parallel)
 
mclapply(1:4, mc.cores=8, function(x)
{
  set.seed(x)
  sapply(1:4, function(y) simulation(x, y))
})