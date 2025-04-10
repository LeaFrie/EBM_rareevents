library(R.matlab)
library(tidyverse)
library(splines)
library(ggplot2)
library(patchwork)
library(KSD)
library(matlib)
library(CompQuadForm)
library(numDeriv)
library(metRology)
library(ordinal)
library(SpatialExtremes)

R_fct <- function(x1,x2){
  a <- 3 + 0.1*(x1-x2)^2 - (x1+x2)/sqrt(2)
  b <- 3 + 0.1*(x1-x2)^2 + (x1+x2)/sqrt(2)
  c <- (x1-x2) + 6/sqrt(2)
  d <- (x2-x1) + 6/sqrt(2)
  return(-min(a,b,c,d))
}
#################################################################################

nr_it = 4000
p <- NaN
for(l in 1:50){
  rr <- NaN
  for(i in 1:nr_it){
    rr[i] = R_fct(rnorm(1,0,1),rnorm(1,0,1))
  }

  p[l] <- length(which(rr>=2))/nr_it
}

mean(p)
sd(p)/mean(p)




