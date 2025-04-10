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
a=as.integer(Sys.time())
set.seed(a)

statistics <- matrix(0,20,10)

for(indd in 5:5){

# number of samples per subset
if(indd<=2){
  nr_it = 40
}
if(indd<=4 & indd > 2){
    nr_it = 60
}
if(indd<=6 & indd > 4){
    nr_it = 80
}

# number of MH steps per representative of subset   
  if(indd%%2 == 1){
    mh_steps = 5
  }
  if(indd%%2 == 0){
    mh_steps = 10
  }

budget <- 10000
nr_risk = floor((budget - nr_it)/(nr_it*mh_steps))
x = seq(0,nr_risk,1)
alpha2_fix = log10(x)/log10(nr_risk)*5 - 3
alpha2_fix[abs(alpha2_fix - 0)==min(abs(alpha2_fix - 0))] = 0

#####
pesti = matrix(0,20,1)
pesti1 = matrix(0,20,1)
forwards = matrix(0,20,1)

for(ppp in 1:500){

# posterior
  
for(p in 1:nr_it){
  rr <- NaN
  posterior <- matrix(0,nr_it,2)
  for(i in 1:nr_it){
    theta1 = rnorm(1,0,1)
    theta2 = rnorm(1,0,1)
    rr[i] = R_fct(theta1, theta2)
    posterior[i,1]= theta1
    posterior[i,2]= theta2
  }
}

#######
# subset sampling
TTT = matrix(0,1000,1)

states = array(NA,c(1000,nr_it,2))
risks = array(NA,c(1000,nr_it,1))
prob = matrix(1,1000,1)
states[1,,] = posterior[,]
risks[1,,] = rr
TTT[2] = alpha2_fix[2]
prob[2] = sum(risks[1,,]>=TTT[2])/nr_it

th <- 2
jr <- 0.5

while(TTT[th] < 2) {

  ind = which(risks[th-1,,]>=TTT[th])
  count1 = matrix(0,mh_steps,nr_it)
  
  for(it in 1:nr_it){
    
    theta_old = states[th-1,sample(ind,1),]
    prior_old = sum(log(dnorm(theta_old)))
   
    for(st in 1:mh_steps){
      
      theta_new = sqrt(1-jr^2)*theta_old + jr*rnorm(2)
      
      if(R_fct(theta_new[1], theta_new[2]) >= TTT[th]){
        theta_old = theta_new
        count1[st,it] = 1
      }
      
    states[th,it,] = theta_old
    risks[th,it,] = R_fct(theta_old[1], theta_old[2])
      
  
    }
  }

  arr = sum(count1)/(mh_steps*nr_it)
    if(arr < 0.3){
      jr = jr/1.01
    }
    if (arr >= 0.3){
      jr = jr*1.01
    }
  jr = max(min(jr,1),0)
  
  th = th + 1
  TTT[th] = alpha2_fix[th]
  prob[th] = sum(risks[th-1,,]>=TTT[th])/nr_it
  
  if(TTT[th] >= 2){
    TTT[th] = 2
    prob[th] = sum(risks[th-1,,]>=TTT[th])/nr_it
    rounds <- th - 2
  }

}

pesti[ppp] = prod(prob)
pesti1[ppp] = prod(prob[1:which(abs(alpha2_fix - 0)==min(abs(alpha2_fix - 0)))])
forwards[ppp] = (nr_risk+1)*nr_it*mh_steps

}

statistics[indd,1] <- mean(forwards)
statistics[indd,2] <-mean(pesti)
statistics[indd,3] <-sd(pesti)/mean(pesti)
statistics[indd,4] <-mean(pesti1)
statistics[indd,5] <-sd(pesti1)/mean(pesti1)

}


