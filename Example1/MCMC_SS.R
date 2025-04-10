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
library(combinat)
library(Metrics)

data <- readMat('Data_example1.mat')

y_true = data$y
theta_true = data$theta.true
lh_true = data$lh.true
R_true = data$R.true

nr_cell = 9
obs_noise = 0.05
obs_cov = diag(3)*obs_noise^2
J = matrix(0,3,nr_cell); J[1,3] = 1; J[2,7] = 1; J[3,9] = 1
prior_mean = 1
prior_sd = 0.3

## analytical solution
posterior_cov = inv(inv(prior_sd^2*diag(nr_cell)) + t(J)%*%inv(obs_cov)%*%J)
posterior_mean = posterior_cov%*%(t(J)%*%inv(obs_cov)%*%t(y_true) + 1/prior_sd^2 * prior_mean*matrix(1,nr_cell,1))
prob_true <- davies(20, lambda = diag(posterior_cov), delta = posterior_mean^2/diag(posterior_cov), lim=50000, acc=10^(-10))$Qq

#################################################################################
a=as.integer(Sys.time())
set.seed(a)

statistics <- matrix(0,20,10)

indd <- 1

for(i in 2:2){

# number of samples per subset
if(i<=2){
  nr_it = 60
}
if(i<=4 & i > 2){
    nr_it = 60
}

# number of MH steps per representative of subset   
  if(i%%2 == 1){
    mh_steps = 10
  }
  if(i%%2 == 0){
    mh_steps = 20
  }

budget <- 72000
steps_burn = 100 # for first sample posterior
thin = 500 # for first sample posterior
nr_it_0 = steps_burn+thin*nr_it

nr_risk = floor((budget - nr_it_0)/(nr_it*mh_steps))
x = seq(0,nr_risk,1)
alpha2_fix = log10(x)/log10(nr_risk)*15 + 5
alpha2_fix[1] = 0

#####
pesti = matrix(0,20,1)
forwards = matrix(0,20,1)

for(ppp in 1:50){

# posterior

jr = 0.3
posterior <- matrix(0,nr_it_0,11)
R_samp <- matrix(0,nr_it_0,1)
count = matrix(0,nr_it_0,1)

theta_old = prior_mean + prior_sd*rnorm(9,0,1)
y_old = theta_old[c(3,7,9)]
lh_old = sum(log(dnorm(y_old, y_true, obs_noise)))
prior_old = sum(log(dnorm(theta_old,prior_mean,prior_sd)))

for(it in 1:nr_it_0){

  theta_new = NaN
  for (c in c(3,7,9)){
    theta_new[c] = theta_old[c] + 1/10*jr*rnorm(1,0,1)
  }
  for(c in c(1,2,4,5,6,8)) {
    theta_new[c] = theta_old[c] + jr*rnorm(1,0,1)
  }
  
  y_new = theta_new[c(3,7,9)]
  lh_new = sum(log(dnorm(y_new, y_true, obs_noise)))
  prior_new = sum(log(dnorm(theta_new,prior_mean,prior_sd)))
  
  arr = exp(lh_new + prior_new - lh_old - prior_old)
  
  if(arr > runif(1)){
    theta_old = theta_new
    lh_old = lh_new
    prior_old = prior_new
    count[it,1] = 1
  }
  
  posterior[it,1:9] = theta_old
  posterior[it,10] = lh_old
  posterior[it,11] = prior_old
  
  R_samp[it,1] = sum(theta_old^2)

}

#######
# subset sampling
T = matrix(0,1000,1)

states = array(NA,c(1000,nr_it,9))
risks = array(NA,c(1000,nr_it,1))
prob = matrix(1,1000,1)
states[1,,] = posterior[seq(steps_burn,nr_it_0,length.out=nr_it),1:9]
risks[1,,] = rowSums(states[1,,]^2)
T[2] = alpha2_fix[2]
prob[2] = sum(risks[1,,]>=T[2])/nr_it

th <- 2

while(T[th] < 20) {

  ind = which(risks[th-1,,]>=T[th])
  count1 = matrix(0,mh_steps,nr_it)
  
  for(it in 1:nr_it){
    
    theta_old = states[th-1,sample(ind,1),]
    y_old = theta_old[c(3,7,9)]
    lh_old = sum(log(dnorm(y_old, y_true, obs_noise)))
    prior_old = sum(log(dnorm(theta_old,prior_mean,prior_sd)))

    for(st in 1:mh_steps){
    
      theta_new = NaN
      for (c in c(3,7,9)){
        theta_new[c] = theta_old[c] + 1/10*jr*rnorm(1,0,1)
      }
      for(c in c(1,2,4,5,6,8)) {
        theta_new[c] = theta_old[c] + jr*rnorm(1,0,1)
      }
      
      y_new = theta_new[c(3,7,9)]
      lh_new = sum(log(dnorm(y_new, y_true, obs_noise)))
      prior_new = sum(log(dnorm(theta_new,prior_mean,prior_sd)))
      
      arr = exp(lh_new + prior_new - lh_old - prior_old)
      
      if(arr > runif(1) & sum(theta_new^2) >= T[th]){
        theta_old = theta_new
        lh_old = lh_new
        prior_old = prior_new
        count1[st,it] = 1
      }
      
    states[th,it,] = theta_old
    risks[th,it,] = sum(theta_old^2)
      
  
    }
  }
  
  arr = sum(count1)/(mh_steps*nr_it)
    if(arr < 0.3){
      jr = jr/1.01
    }
    if (arr >= 0.3){
      jr = jr*1.01
    }
  
  th = th + 1
  T[th] = alpha2_fix[th]
  prob[th] = sum(risks[th-1,,]>=T[th])/nr_it
  
  if(T[th] >= 20){
    T[th] = 20
    prob[th] = sum(risks[th-1,,]>=T[th])/nr_it
    rounds <- th - 2
  }

}

pesti[ppp] = prod(prob)
forwards[ppp] = nr_it_0 + nr_risk*nr_it*mh_steps

}

statistics[indd,1] <- mean(forwards)
statistics[indd,2] <-mean(pesti)
statistics[indd,3] <-sd(pesti)/mean(pesti)
statistics[indd,4] <-rmse(pesti,prob_true)

indd <- indd + 1

}

##################################################################################
