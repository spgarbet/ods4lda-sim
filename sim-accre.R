source("run.R")

  #############################################################################
 ##
## ACCRE batch run
args <- commandArgs(trailingOnly=TRUE)
x    <- as.numeric(args[1])
set.seed(x)
sapply(1:8, function(y) simulation(x,y))
