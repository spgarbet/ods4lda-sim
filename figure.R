n <- 2000

load("results-2018-06-11.RData")

x <- data.frame(
  factor     = c("I", "T", "G", "C", "GT"),
  baseline   = rep(375/n, 5),
  efficiency = c(3.473, 3.766, 1.781, 4.891, 1.928)
)

barplot(x$efficiency*x$baseline, names=x$factor, ylim=c(0,1),
        col=c('grey', 'grey', 'yellow', 'grey', 'yellow'), 
        space=0.1,
        oma=c(0,0,0,0),
        labels=FALSE)
abline(h=375/n, col='red', lwd=2)

library(tidyr)
library(dplyr)
library(reshape2)

var.sums <- 
  subset(results, method %in% c("acml", "mi2")) %>%
  group_by(scenario, sampling, method) %>%
  summarise(I=sum(var.int),
            T=sum(var.time),
            G=sum(var.grp),
            C=sum(var.conf),
            TG=sum(var.tg)) %>%
  gather(factor, variance, -scenario, -sampling, -method)

var.sums$factor <- factor(var.sums$factor, levels=c("I", "T", "G", "C", "TG"))

ref <- var.sums %>% filter(sampling == "random" & method=="acml")

N                    <- 2000
NsRand.sim           <- c(375,375,525,375)

reference <- Vectorize(function(sc, fc)
{
  subset(ref, scenario==sc & factor==fc)$variance
})

normalize <- Vectorize(function(sc) NsRand.sim[sc] / N)

efficiency <-
  filter(var.sums, !(sampling == "random" & method=="acml")) %>%
  group_by(scenario, sampling, method, factor) %>%
  summarise(efficiency = reference(scenario, factor)/variance,
            percent    = efficiency*normalize(scenario) )

par(mfrow=c(5, 6), mar=c(1, 0, 1, 0))
for(smp in c("Intercept", "Slope", "Bivariate"))
{
  for(mth in c("ACML", "MI"))
  {
    plot(0:1, 0:1, typ="n", axes=FALSE, ann=FALSE)
    text(0.5, 0, mth, cex=2)
    text(0.5, 0.5, smp, cex=2.2)
  }
}

for(scn in 1:4)
{
  for(smp in c("intercept", "slope", "bivariate"))
  {
    for(mth in c("acml", "mi2"))
    {
      x <- subset(efficiency, sampling == smp &
                              method   == mth &
                              scenario == scn)
      xx <- barplot(x$percent, names=x$factor, ylim=c(0,1),
        col=c('grey', 'grey', 'yellow', 'grey', 'yellow'), 
        space=0.1,
        oma=c(0,0,0,0),
        axes=FALSE,
        axisnames = FALSE)
      #if(scn == 4)
        mtext(c("I", "T", "G", "C", "TG"), side=1, at=xx,line=0.1, cex=0.67)
      text(x= xx, y=x$percent, label=sprintf("%04.2f", round(x$efficiency,2)), pos=3, cex=0.8, col='red')
      abline(h=NsRand.sim[scn]/N, col='red', lwd=2)
    }
  }
}
