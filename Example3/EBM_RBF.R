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

nr_R <- 100

gumbel_scale <- sqrt(1*6)/pi
gumbel_loc <- 2-gumbel_scale*0.5772156

mean_param = 12
std_dev_param <- 2
lognorm_mean <- log(mean_param^2 / sqrt(std_dev_param^2 + mean_param^2))
lognorm_sd <- sqrt(log(1 + (std_dev_param^2 / mean_param^2)))

ZS <- rnorm(1)
ZR <- rnorm(nr_R)

S = qgumbel(pnorm(ZS), gumbel_loc,gumbel_scale)
RR = qlnorm(pnorm(ZR), lognorm_mean/nr_R, lognorm_sd/sqrt(nr_R))
R = prod(RR)

# likelihood
data <- 8^{1/nr_R}
sigma_data <- 0.05

################################################################################
## functions
U_fct <- function(ZS, ZR, data) {
  logprior <- 0 # pCN sum(log(dnorm(ZR))) + log(dnorm(ZS))
  R = qlnorm(pnorm(ZR), lognorm_mean/nr_R, lognorm_sd/sqrt(nr_R))
  loglh <- -0.5*sum(((log(R)-log(data))/sigma_data)^2)
  return(-loglh-logprior)
}

R_fct <- function(ZS, ZR) {
  S = qgumbel(pnorm(ZS), gumbel_loc,gumbel_scale)
  RR = qlnorm(pnorm(ZR), lognorm_mean/nr_R, lognorm_sd/sqrt(nr_R))
  R = prod(RR)
  return(S-R)
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

ff <- function(x) {
  if(p_dist == 'normal'){
    a <- - rbf(x) - log(dnorm(x,p_m,p_sd))
  }
  if(p_dist == 'student'){
    a <- - rbf(x) - log(dt.scaled(x, p_df, mean = p_m, sd = p_sd))
  }
  if(p_dist=='gev'){
    a <- - rbf(x) - log(dgev(x, p_m, p_sd, shape = sh))
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

score_gev <- function(x){
  s <- sh
  y <- -(s+1)/(s*x+1) + (1 + s*x)^(-1/s-1)
  return(y)
}



################################################################################
plott='yes'

# p(r)
p_dist = 'gev'
if(p_dist=='student'){
  p_m <- 20
  p_sd <- 7
  p_df <- 1
  score_f <- score_f1
}
if(p_dist=='normal'){
  p_m <- 0
  p_sd <- 5
  score_f <- score_f2
}
if(p_dist=='gev'){
  p_m <- 0
  p_sd <- 7
  sh <- 0
  score_f <- score_gev
}

hist(rnorm(1000,0,8), breaks=seq(-100,3000,0.5), xlim=c(-20,40), col='blue')
hist(rgev(1000,p_m,p_sd,sh), breaks=seq(-100,3000,0.5), xlim=c(-20,40), add=TRUE)

# stop
stop_crit <- 'NO'
stop_con <- 1.0
nr_tt <- 27 # 7700/nr_samples
a_ksd <- 0.5

samp_burn=10
nr_samples1 <- round(max((0.1/a_ksd)*500,100))/2
th <- 6
nr_samples <- th*nr_samples1+samp_burn

conf <- matrix(0,6,2)
forward_mean <- matrix(0,6,1)
nr_optim <- matrix(0,6,2)

ppp <- 1
prob_vec <- NaN
prob_vec1 <- NaN
ksd_vec <- NaN
nr_vec <- NaN

Vsafe <- matrix(0,1601, 50)

while(ppp < 50) {

a=as.numeric(Sys.time())
set.seed(a) 

## initialize RBF basis
weights <- seq(0,0,length.out=500) 
centers <- seq(-100,100,length.out=500)
eps <- 0.5

ZS_old = rnorm(1)
ZR_old = rnorm(nr_R)

loss = NaN
prob = NaN
ksd <- NaN
ksd_p <- NaN
ar <- NaN

jr = 0.5
learn_rate <- 0.5 # 50
learn_rate_end <- 0
momentum <- 0.5
velocity = seq(0,0,length.out=500)

it <- 1
stopp <- 0

while(it < nr_tt){   # stopp < stop_con 
  
  set.seed(it*ppp*runif(1))

  if(p_dist=='student'){
    samples1 <- rt.scaled(nr_samples1, p_df, mean = p_m, sd = p_sd)
  }
  if(p_dist=='normal'){
    samples1 <- rnorm(nr_samples1,p_m,p_sd)
  }
  if(p_dist=='gev'){
    samples1 <- rgev(nr_samples1,p_m,p_sd, shape=sh)
  }

arr = 0
the = matrix(NA,nr_samples,nr_R+1)
outt <- matrix(NA,nr_samples,1)

u <- R_fct(ZS_old, ZR_old)
out_old = U_fct(ZS_old, ZR_old, data) + rbf(u)

for(t in 1:nr_samples){
  
  ZS_new = sqrt(1-jr^2)*ZS_old + jr*rnorm(1)
  ZR_new = sqrt(1-(0.7*jr)^2)*ZR_old + 0.7^2*jr*rnorm(nr_R)
  
  u <- R_fct(ZS_new, ZR_new)
  out_new = U_fct(ZS_new, ZR_new, data) + rbf(u)
  
  acce = exp(-out_new + out_old)

  if (acce > runif(1)) {
    ZS_old = ZS_new
    ZR_old = ZR_new
    arr = arr + 1
    out_old = out_new
  }
  
  the[t,1] = ZS_old
  the[t,2:(1+nr_R)] = ZR_old
  outt[t,1] = out_old
}

ar[it] = arr/nr_samples

pV_samples = the[round(seq(samp_burn, nr_samples-1, length.out=nr_samples1)),]
samples2 = matrix(0,nr_samples1,1)
for (ii in 1:nr_samples1){
  samples2[ii] = R_fct(pV_samples[ii,1],pV_samples[ii,2:(1+nr_R)])
}

## loss
loss[it] <- mean(rbf(samples1)-rbf(samples2))

# probability
i <- integrate(ff,-100,100)
ii <- integrate(ff,0,100)
prob[it] <- ii$value/i$value

if(p_dist=='student'){
  xxx <- c((samples2-p_m)/p_sd %>% as.vector())
}
if(p_dist=='normal'){
  xxx <- c(samples2%>% as.vector())
}
if(p_dist=='gev'){
  xxx <- c((samples2-p_m)/p_sd %>% as.vector())
}

kksd <- ksd_corr(xxx, score_f, nboot=1000, a=a_ksd)

ksd[it] <- kksd$ksd
ksd_p[it] <- kksd$p

stopp <- ksd_p[it]

if(plott=='yes'){
par(mfrow=c(3,2))
x <- seq(-40,40,0.05)
plot(x,rbf(x), type='l', xlab='r', ylab='', ylim=c(-20,20), lwd=4, cex.main=1.5, main='V(r)', cex.lab=1.5)
ax <- seq(-40,2000,0.5)
hist(samples1, breaks=ax, xlab='r', ylab='', main='Samples p(r) and pV(r)', ylim=c(0,20), cex.main=1.5, cex.lab=1.5, xlim=c(-20,20))
hist(samples2, breaks=ax,add=TRUE, col='cornflowerblue')
plot(loss, xlim=c(1,nr_tt), ylim=c(-10,10), xlab='Iteration', ylab='', cex.main=1.5, cex.lab=1.5, main='Loss')
abline(h=0)
plot(prob, log="y", xlim=c(1,nr_tt), ylim=c(10^(-7),10^(-2)), xlab='Iteration', ylab='', cex.main=1.5, cex.lab=1.5, main='Risk probability')
abline(h=2.1e-05)
plot(ksd, xlim=c(1,nr_tt), ylim=c(-0.05,0.5), xlab='Iteration', ylab='', cex.main=1.5, cex.lab=1.5, main='KSD')
abline(h=0)
plot(ksd_p[1:it], xlim=c(1,nr_tt), ylim=c(0,0.5), xlab='Iteration', ylab='', cex.main=1.5, cex.lab=1.5, main='P-value')
abline(h=0.05) }

for(ii in 1:length(centers)) {
  velocity[ii] <- momentum*velocity[ii] + (1-momentum)*(mean(gk(samples1-centers[ii]))-mean(gk(samples2-centers[ii])))
  weights[ii] <- weights[ii] - learn_rate * velocity[ii]
}

if(ar[it] < 0.3){
  jr = max(min(jr/1.1,1),0)
}
if (ar[it] >= 0.3){
  jr = max(min(jr*1.1,1),0)
}

if(learn_rate > learn_rate_end){# before 0.05
  learn_rate <- learn_rate*exp(-0.005*it)  # learn_rate/learn_rate_dc
}

it <- it + 1

}

if(stop_crit=='YES'){
if(stopp >= stop_con){  
  prob_vec[ppp] = prob[it-1]
  prob_vec1[ppp] = prob[7700/nr_samples]
  ksd_vec[ppp] = ksd_p[it-1] 
  nr_vec[ppp] = it
  ppp <- ppp + 1
  }
}

if(stop_crit=='NO'){
  prob_vec[ppp] = prob[it-1]
  ksd_vec[ppp] = ksd_p[it-1] 
  nr_vec[ppp] = it
  Vsafe[,ppp] <- rbf(x)
  ppp <- ppp + 1
}

print(prob_vec[ppp-1])

}
