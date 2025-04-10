library(R.matlab)
library(tidyverse)
library(dplyr)
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
library(philentropy)

R_fct <- function(x1,x2){
  a <- 3 + 0.1*(x1-x2)^2 - (x1+x2)/sqrt(2)
  b <- 3 + 0.1*(x1-x2)^2 + (x1+x2)/sqrt(2)
  c <- (x1-x2) + 6/sqrt(2)
  d <- (x2-x1) + 6/sqrt(2)
  return(-min(a,b,c,d))
}

rr <- NaN
for(i in 1:1000000){
  rr[i] = R_fct(rnorm(1,0,1),rnorm(1,0,1))
}
hist(rr, breaks=seq(-10,100,length.out=50))
abline(v=-6/sqrt(2))

hist(rgumbel(1000000,2,3), breaks=seq(-10,60,length.out=50))
# hist(rr, add=TRUE, col='blue', breaks=seq(-10,60,length.out=50))
abline(v=-6/sqrt(2))

c1 <- rgb(255,192,203, max = 255, alpha = 80)
c2 <- rgb(173,250,230, max = 255, alpha = 100)

hist(rgumbel(1000000,0,3.5), breaks=seq(-10,400,length.out=400), xlim=c(-10,30))
hist(rgev(1000000, loc = 0, scale = 2, shape = 1/(6/sqrt(2))), breaks=seq(-10,400,length.out=400), add=TRUE, col=c1)
hist(rr, add=TRUE, col=c2, breaks=seq(-10,400,length.out=400))


################################################################################
## functions
U_fct <- function(x1,x2) {
  logprior <- sum(log(dnorm(x1))) + log(dnorm(x2))
  return(-logprior)
}

gk <- function(x){
  return(exp(-(eps*x)^2)) # gaussian 
}

rbf <- function(x){
  y <- NaN
  for(ii in 1:length(x)) {
    y[ii] <- sum(weights * gk(abs(x[ii] - centers)))
  }
  return(y)
}

# probability
ff <- function(x) {
  
  a <- NaN
  for(l in 1:length(x)){
    
    x_xxx <- xxx[abs(x[l]-xxx)==min(abs(x[l]-xxx))][1]
    if(p_dist == 'gev'){
      a[l] <- - V[abs(x[l]-xxx)==min(abs(x[l]-xxx))][1] - log(dgev(x_xxx,p_m,p_sd,shape=sh))
    }
    if(p_dist == 'normal'){
      a[l] <- - V[abs(x[l]-xxx)==min(abs(x[l]-xxx))][1] - log(dnorm(x_xxx,p_m,p_sd))
    }
  }
  
  return(exp(-a))
}

score_f1 <- function(x){
  return(-(p_df+1)/2*1/(1+x^2/p_df)*2*x/p_df)
}

score_f2 <- function(x){
  return(-1/p_sd^2*(x-p_m))
}
  
find_median_distance <- function(Z){
  
  if(is.data.frame(Z)){
    Z = data.matrix(Z)
  }else{
    Z = as.array(Z)
  }
  size1 <- dim(Z)[1]
  size2 <- dim(Z)[2]
  
  # if size of Z is greater than 100, randomly sample 100 points
  if(size1 > 100){
    if(is.na(size2)){
      Zmed <- Z[sample(size1,100)]
    }else{
      Zmed <- Z[sample(size1,100),]
    }
    size1 = 100
  }else{
    Zmed <- Z
  }
  
  Zmedsq <- Zmed * Zmed;
  if(is.na(dim(Z)[2]))
    G <- Zmedsq
  else
    G <- rowSums(Zmedsq)
  
  # Create row/col repeated matrices
  Q <- rep.col(G,size1)
  R <- rep.row(t(G),size1)
  
  dists <- Q + R - 2 * Zmed %*% t(Zmed)
  dists[lower.tri(dists, diag = TRUE)] = 0
  dists <- array(dists,dim=c(size1^2,1))
  median_dist <- median(dists[dists > 0 ])
  
  return(median_dist)
}


ratio_median_heuristic <- function(Z, score_function){
  
  Z = as.array(Z)
  size1 <- dim(Z)[1]
  size2 <- dim(Z)[2]
  
  # if size of Z is greater than 100, randomly sample 100 points
  if(size1 > 100){
    if(is.na(size2)){
      Zmed <- Z[sample(size1,100)]
    }else{
      Zmed <- Z[sample(size1,100),]
    }
    size1 = 100
    print('Sampled (Heuristic)')
  }else{
    Zmed <- Z
    print('Original 100 dataset used (Heuristic)')
  }
  
  Zmedsq <- Zmed * Zmed;
  if(is.na(dim(Z)[2])){
    G <- Zmedsq
  } else{
    G <- rowSums(Zmedsq)
  }
  
  # Create row/col repeated matrices
  Q <- rep.col(G,size1)
  R <- rep.row(t(G),size1)
  
  dists <- Q + R - 2 * Zmed %*% t(Zmed)
  dists[lower.tri(dists, diag = TRUE)] = 0
  dists <- array(dists,dim=c(size1^2,1))
  median_dist <- median(dists[dists > 0 ])
  
  if(is.na(size2)){
    Zmed = as.double(Zmed)
  }
  sqx <- score_function(Zmed)
  sqxx <- sqx %*% t(sqx)
  sqxx[lower.tri(sqxx, diag = TRUE)] = 0
  sqxx <- array(sqxx,dim=c(size1^2,1))
  median_sqxx <- median(sqxx[sqxx > 0])
  
  h <- (median_dist / median_sqxx)^(1/4)
  medInfo <- list("h"=h, "median_dist" = median_dist, "median_sqxx" = median_sqxx)
}

rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}


rep.col<-function(x,n){
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}


repmat <- function(X,a=1,b=1){
  rows <- dim(X)[1]
  cols <- dim(X)[2]
  if(is.null(cols))
    cols <- 1
  rowRep <- matrix(rep(t(X),a),ncol = cols, byrow = TRUE)
  newX <- matrix(rep(rowRep, b),ncol=cols*b)
  return(newX)
}

# Cleans up workspace
# Removes all current variables

cleanup <- function(){
  rm(list=ls())
}


# Returns size of x
#
# Returns a list, which stores the dimension of the matrix.
#
# @param x   A matrix of data
# @return    List where
# n = number of rows of x;
# dim = number of columns of x
#
# @examples
# x <- matrix(c(1,2,3,4,5,6,7,8),ncol=2)
# dim <- getDim(x)
getDim <- function(x){
  if(is.array(x)){
    n <- dim(x)[1]; dimen <- dim(x)[2]
  }else{
    x <- array(x)
    n <- dim(x)[1]; dimen <- 1
  }
  
  result <- list("n" = n, "dim" = dimen)
  return(result)
}

# Score function for gamma distribution
gamma_score <- function(x, shape, rate=1,scale=1/rate){
  return ((shape-1)/x - 1/scale)
}

# Function that can be used to retain legend of a plot
# http://www.sthda.com/
get_legend<-function(myggplot){
  tmp <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# Function that rounds up given number to three significant digits
custround <- function(x){
  round(x,3)
}

ksd_corr <- function(x, score_function, nboot, a, kernel = "rbf", width = -1){
  #cleanup()
  # set.seed(1)
  
  ######## NEED TO IMPLEMENT#######
  #list[kernel, h, nboot] = process_varargin(varargin, 'kernel', 'rbf', 'width', -1, 'nboot', 0)
  
  h <- width
  
  # If score_function = 'gaussian', calculate the score function for standard norm
  if(typeof(score_function) == 'character'){
    if(tolower(score_function) == 'gaussian'){
      score_function <- function(val){
        return (-1 * val)
      }
    }
  }else{
    
  }
  
  # Set up the bandwidth of RBF Kernel
  if((typeof(h)== 'character' & tolower(h) == 'median') | h == -1){
    h = sqrt(0.5 * find_median_distance(x))
  }else if(h == -2){
    #h = 1
    medInfo = ratio_median_heuristic(x, score_function)
    h = medInfo$h
  }
  
  # Get parameters
  dimensions <-getDim(x)
  n <- dimensions$n; dimen <- dimensions$dim
  if(is.numeric(score_function) | is.vector(score_function)){
    Sqx = score_function
  }else{
    Sqx = score_function(x)
  }
  
  
  if(tolower(kernel) == 'rbf'){
    XY <- x %*% t(x)
    if(is.null(dim(x)[2])){
      x2 <- x^2
      sumx <- x * Sqx
    } else{
      x2 <- array(rowSums(x^2),dim=c(n,1))
      sumx <- rowSums(x * Sqx)
    }
    X2e <- repmat(x2,1,n)
    
    H <- (X2e + t(X2e) - 2*XY)
    Kxy <- exp(-H/(2*h^2))
    
    sqxdy <- -(Sqx %*% t(x) - repmat(sumx,1,n))/h^2
    dxsqy <- t(sqxdy)
    dxdy <- (-H/h^4 + dimen/h^2)
    
    M = (Sqx %*% t(Sqx) + sqxdy + dxsqy + dxdy) * Kxy
    
  }else{
    warning('Wrong Kernel')
  }
  
  
  M2 <- M - diag(diag(M))
  ksd <- sum(M2)/(n^2)
  ksdV <- sum(M) / (n^2)
  
  ##---------------------Bootstrap methods---------------------------##
  
  bootstrapSamples <-  rep(NA,nboot)
  for(i in 1:nboot){
    # mimic markov process
    w = matrix(0,nr_samples1,1)
    w[1,1] = 1
    for(ii in 2:nr_samples1){
      u = runif(1)
      if(u > a){
        w[ii,1] = w[ii-1,1]
      }
      if(u <= a){
        w[ii,1] = -w[ii-1,1]
      }
    }
    
    bootstrapSamples[i] <- 1/nr_samples1^2*t(w)%*%M2%*%w
  }
  p <- mean(bootstrapSamples >= ksd)
  
  info <- list("bandwidth"=h, "M" = M, "nboot" = nboot, "ksd_V" = ksdV)
  result <- list("ksd" = ksd, "p"=p, "bootStrapSamples"=bootstrapSamples, "info"=info)
  
  return(result)
}

tf <- function(x){
  if(sh==0){
    return(exp(-x))
  }
  if(sh>0){
    return((1+sh*x)^(-1/sh))
  }
}

score_gev <- function(x){
  if(sh>0){
    delt_t <- -1/sh*(1+sh*x)^(-1/sh-1)
  }
  if(sh==0){
    delt_t <- -exp(-x)
  }
  y <- (sh+1)*1/tf(x)*delt_t - delt_t
  return(y)
}



################################################################################
pV_dens <- function(x, samples, adj=1){
  y <- density(samples, bw = "nrd", adjust = adj, 
               kernel = "gaussian", n=length(xxx), from=x[1], to=x[length(x)])
  return(y$y)
}

adj <- 1

plott='yes'

# p(r)
p_dist = 'normal'
if(p_dist=='student'){
  p_m <- 20
  p_sd <- 7
  p_df <- 1
  score_f <- score_f1
}
if(p_dist=='normal'){
  p_m <- 0
  p_sd <- 2
  score_f <- score_f2
}
if(p_dist=='gumbel'){
  p_m <- 2
  p_sd <- 3.0
  score_f <- score_gumbel
}
if(p_dist=='gev'){
  p_m <- 2
  p_sd <- 2
  sh <- 0.33
  score_f <- score_gev
}

hist(rgev(1000,p_m,p_sd,0))

a_ksd <- 0.5
samp_burn=5
nr_samples1 <- round(max((0.1/a_ksd)*500,100))
th <- 4
nr_samples <- th*nr_samples1+samp_burn

# stop
stop_crit <- 'NO'
stop_con <- 1.0
nr_tt <- 200

xxx <- seq(-10,100, 0.1)
V <- seq(0,0, length.out=length(xxx))

conf <- matrix(0,6,2)
forward_mean <- matrix(0,6,1)
nr_optim <- matrix(0,6,2)

ppp <- 1
prob_vec <- NaN
prob_vec1 <- NaN
ksd_vec <- NaN
nr_vec <- NaN


while(ppp <= 50) {

a=as.numeric(Sys.time())
set.seed(a) 

x1_old = rnorm(1)
x2_old = rnorm(1)

loss = NaN
prob = NaN
prob1 = NaN
ksd <- NaN
ksd_p <- NaN
ar <- NaN

jr = 2.5
learn_rate <- 6.5
learn_rate_dc <- 1.0
learn_rate_end <- 0
momentum <- 0.5
velocity = seq(0,0,length.out=length(xxx))

V <- seq(0,0, length.out=length(xxx))

it <- 1
stopp <- 0

while(stopp < stop_con & it < nr_tt){
  
  set.seed(it*ppp*runif(1))

  if(p_dist=='student'){
    samples1 <- rt.scaled(nr_samples1, p_df, mean = p_m, sd = p_sd)
  }
  if(p_dist=='normal'){
    samples1 <- rnorm(nr_samples1,p_m,p_sd)
  }
  if(p_dist=='gumbel'){
    samples1 <- rgumbel(nr_samples1,p_m,p_sd)
  }
  if(p_dist=='gev'){
    samples1 <- rgev(nr_samples1,p_m,p_sd, shape=sh)
  }

arr = 0
the = matrix(NA,nr_samples,2)
outt <- matrix(NA,nr_samples,1)

u <- R_fct(x1_old, x2_old)
out_old = U_fct(x1_old, x2_old) + V[abs(u-xxx)==min(abs(u-xxx))]

for(t in 1:nr_samples){
  
  x1_new = x1_old + jr*rnorm(1)
  x2_new = x2_old + jr*rnorm(1)
  
  u <- R_fct(x1_new, x2_new)
  out_new = U_fct(x1_new, x2_new) + V[abs(u-xxx)==min(abs(u-xxx))]
  
  acce = exp(-out_new + out_old)

  if (acce > runif(1)) {
    x1_old = x1_new
    x2_old = x2_new
    arr = arr + 1
    out_old = out_new
  }
  
  the[t,1] = x1_old
  the[t,2] = x2_old
  outt[t,1] = out_old
}

ar[it] = arr/nr_samples

pV_samples = the[round(seq(samp_burn, nr_samples-1, length.out=nr_samples1)),]
samples2 = matrix(0,nr_samples1,1)
for (ii in 1:nr_samples1){
  samples2[ii] = R_fct(pV_samples[ii,1],pV_samples[ii,2])
}

# loss[it] <- KL(rbind(samples1,samples2[,1]), est.prob="empirical")

# ## proba
i <- integrate(ff,-10,100, subdivisions=100, stop.on.error = FALSE)
ii <- integrate(ff,0,100, subdivisions=100, stop.on.error = FALSE)
ii1 <- integrate(ff,2,100, subdivisions=100, stop.on.error = FALSE)
prob[it] <- ii$value/i$value
prob1[it] <- ii1$value/i$value

if(p_dist=='student'){
  xxxx <- c((samples2-p_m)/p_sd %>% as.vector())
}
if(p_dist=='normal'){
  xxxx <- c(samples2%>% as.vector())
}
if(p_dist=='gumbel'){
  xxxx <- c((samples2-p_m)/p_sd %>% as.vector())
}
if(p_dist=='gev'){
  xxxx <- c((samples2-p_m)/p_sd %>% as.vector())
}

kksd <- ksd_corr(xxxx, score_f, nboot=1000, a=a_ksd)

ksd[it] <- kksd$ksd
ksd_p[it] <- kksd$p

stopp <- ksd_p[it]

if(plott=='yes'){
par(mfrow=c(3,2))
x <- seq(-20,50,0.05)
plot(xxx,V, type='l', xlab='r', ylab='', ylim=c(-20,20),xlim=c(-10,30), lwd=4, cex.main=1.5, main='V(r)', cex.lab=1.5)
ax <- seq(-10,1000,length.out=800)
hist(samples1, breaks=ax, xlab='r', ylab='', main='Samples p(r) and pV(r)', ylim=c(0,60), cex.main=1.5, cex.lab=1.5, xlim=c(-10,30))
hist(samples2, breaks=ax,add=TRUE, col='cornflowerblue')
# plot(loss, xlim=c(1,nr_tt), ylim=c(-10,10), xlab='Iteration', ylab='', cex.main=1.5, cex.lab=1.5, main='Loss')
# abline(h=0)
plot(prob, log="y", xlim=c(1,nr_tt), ylim=c(10^(-5),10^(-1)), xlab='Iteration', ylab='', cex.main=1.5, cex.lab=1.5, main='Risk probability')
abline(h=4e-03)
plot(prob1, log="y", xlim=c(1,nr_tt), ylim=c(10^(-7),10^(-3)), xlab='Iteration', ylab='', cex.main=1.5, cex.lab=1.5, main='Risk probability')
abline(h=1e-05)
plot(ksd, xlim=c(1,nr_tt), ylim=c(-0.05,0.5), xlab='Iteration', ylab='', cex.main=1.5, cex.lab=1.5, main='KSD')
abline(h=0)
plot(ksd_p[1:it], xlim=c(1,nr_tt), ylim=c(0,0.1), xlab='Iteration', ylab='', cex.main=1.5, cex.lab=1.5, main='P-value')
abline(h=stop_con) }
#plot(xxx,pV_dens(xxx, samples2), xlim=c(-10,30))  }
# plot(the[,1], the[,2], xlim=c(-6,6), ylim=c(-6,6))  }
# hist(the[samp_burn:nr_samples,1], xlim= c(-6,6))
# hist(the[samp_burn:nr_samples,2], xlim= c(-6,6))}

if(ar[it] < 0.3){
  jr = jr/1.1
}
if (ar[it] >= 0.3){
  jr = jr*1.1
}

if(learn_rate > learn_rate_end){
  learn_rate <- learn_rate*learn_rate_dc
}

it <- it + 1

if(p_dist == 'gev'){
  velocity <- momentum*velocity + (1-momentum)*(dgev(xxx,p_m,p_sd, shape=sh) - pV_dens(xxx, samples2, adj))
}
if(p_dist == 'normal'){
  velocity <- momentum*velocity + (1-momentum)*(dnorm(xxx,p_m,p_sd) - pV_dens(xxx, samples2, adj))
}

V <- V - learn_rate*velocity

}

if(stop_crit=='YES'){
if(stopp >= stop_con){  
  prob_vec[ppp] = prob[it-1]
  prob_vec1[ppp] = prob1[it-1]
  ksd_vec[ppp] = ksd_p[it-1] 
  nr_vec[ppp] = it
  ppp <- ppp + 1
  }
}

if(stop_crit=='NO'){
  prob_vec[ppp] = prob[it-1]
  prob_vec1[ppp] = prob1[it-1]
  ksd_vec[ppp] = ksd_p[it-1] 
  nr_vec[ppp] = it
  ppp <- ppp + 1
  }

}

