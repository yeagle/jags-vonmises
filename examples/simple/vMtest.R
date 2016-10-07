#Hierarchical 2-parameter Mixture model for JAGS using the rjags package

rm(list=ls(all=TRUE))
graphics.off()

library(rjags)
library(circular)
library(coda)
load.module("vonmises")

N <- 100
Kappa <- 10

memset <- rvonmises(N,0,0)
x <- NULL
for (i in 1:N) {
  x[i] <- rvonmises(1, memset[i], Kappa)
}

data = list (
  "m" = memset,
  "x" = x,
  "N" = N  )

try (
  if (1 == 1) {
    # fitting with rjags
    jagsmodel<-jags.model("vMtest.bug", data=data, n.chains=4, n.adapt=5000)
    codasamples <- coda.samples(jagsmodel, variable.names=c("precision", "x"),n.iter=10000, thin=2)  # turn object of class "jags" into "mcmc.list"
    samples <- as.matrix(codasamples)
  }
)

x.post <- samples[,1:N]
x11()
hist(x.post,100)
