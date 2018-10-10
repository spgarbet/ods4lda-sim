two_phase <- function(sampled, notsampled, hn_scale)
{
  interest     <- c("id", "time", "grp", "conf", "y")
  x            <- sampled[, interest]
  y            <- notsampled[, interest]
  y$grp        <- NA
  x            <- rbind(x, y)
  simZ         <- unique(x$conf)
  n            <- length(simZ) # number of subjects in the dataset
  N_SIEVE      <- 10 # number of breaks in the histogram, set for 10 now but may need to fine tune it later
  simBspline_Z <- matrix(NA, nrow=n, ncol=N_SIEVE)
  cut_z        <- cut(simZ, breaks=quantile(simZ, probs=seq(0, 1, 1/N_SIEVE)), include.lowest = TRUE) # split the range of Z into evenly N_SIEVE spaced quantiles
  
  # create the indicators for each quantile
  for (i in 1:N_SIEVE) 
    simBspline_Z[,i] <- as.numeric(cut_z == names(table(cut_z))[i])
  colnames(simBspline_Z) <- paste("bs", 1:N_SIEVE, sep="")
  
  expand       <- match(x[,'id'], unique(x[,'id']))
  x            <- cbind(x, as.data.frame(simBspline_Z[expand,]))
  
  fit <- smle_lmm(ID="id", Y="y", Time="time", X="grp", Z="conf",
                  Bspline_Z = colnames(simBspline_Z), data = x, hn_scale=hn_scale)
  
  # Reformat results for use in summarization code
  reorder              <- c(1,4,2,3,5) # Reorder to ordering of other fits
  coef                 <- fit$coefficients[reorder, 1]
  serr                 <- fit$coefficients[reorder, 2]
  
  fit$coefficients     <- coef
  fit$covariance       <- matrix(rep(0, 25), nrow=5)
  diag(fit$covariance) <- serr^2
  
  fit
}