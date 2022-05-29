

rm(list=ls())

library(doParallel)

# function

###############################################################################

newton<-function (ftn, x0, tol = 1e-09, max.iter = 50) {
  x <- x0
  fx <- ftn(x)
  iter <- 0
  
  while ((abs(fx[1]) > tol) && (iter < max.iter)) {
    x <- x - fx[1]/fx[2]
    fx <- ftn(x)
    iter <- iter + 1
  }
  if (abs(fx[1]) > tol) {
    return(NULL)
  }
  else {
    return(x)
  }
}

###############################################################################

nr_ftn<-function(x){
  
  fx <- sum((Z-mean(Z))*Y*exp(-x*A))
  dfx <- sum((Z-mean(Z))*Y*exp(-x*A)*-A)
  return(c(fx, dfx))
  
}

###############################################################################

gestimation_ftn <- function(Z,A,U,Y) {
  # returns function value and its derivative at x
  
  newton(nr_ftn,x0=2)
  
}

###############################################################################

n=10^3;iter=10^3;boot=10^3

ratio_est<-g_coef<-c()
ratio_coverage<-g_coverage<-c()

for (i in 1:iter) {
  

  # generate data
  Z1<-rbinom(n=n,size = 1,prob = 0.5)
  U1<-rbinom(n=n,size = 1,prob = 0.5)
  A1<-rbinom(n=n,size=1,prob=exp(-1+Z1+U1)/(1+exp(-1+Z1+U1)))
  psi<-1
  Y1<-rpois(n,exp(A1+U1))
  

  # Bootstrap : 2sls
  
  ratio_est2 <- list() 
  
  for(i in 1:n){
    
    idx<-sample(1:n,replace = T)
    Z<-Z1[idx]
    U<-U1[idx]
    A<-A1[idx]
    Y<-Y1[idx]
    
    
    denom<-lm(A~Z)
    numer<-glm(Y~Z,family = 'poisson')
    

    ratio_est2[i] <- numer$coefficients[2]/denom$coefficients[2]
  }

  
  ###############################################################################    
  

  # 2) bootstrap : g-estimation
  
  g_coef2 <- list() 
  
  for(i in 1:n){
    
    idx<-sample(1:n,replace = T)
    Z<-Z1[idx]
    U<-U1[idx]
    A<-A1[idx]
    Y<-Y1[idx]
    
    g_coef2[i] <- gestimation_ftn(Z=Z,A=A,U=U,Y=Y)
    
  }
  
  
  # list -> vector
  ratio_est2<-unlist(ratio_est2)
  g_coef2<-unlist(g_coef2)
  
  # 2. Estimator 
  # 1) 2sls
  Z<-Z1;U<-U1;A<-A1;Y<-Y1
  
  denom<-lm(A~Z)
  numer<-glm(Y~Z,family = 'poisson')
  ratio_est[i]<-numer$coefficients[2]/denom$coefficients[2]
  
  # 2) G-estimation
  g_coef[i]<-gestimation_ftn(Z=Z,A=A,U=U,Y=Y)
  
  # 3.Coverage
  ratio_coverage[i]<-ifelse((psi>=quantile(ratio_est2,.025))&(psi<=quantile(ratio_est2,.975)),1,0)
  g_coverage[i]<-ifelse((psi>=quantile(g_coef2,.025))&(psi<=quantile(g_coef2,.975)),1,0)
  
  
  
  
  # 4.Check results
  if(i%%10==0){
    
    print(paste0(i,'-th iteration'))
    print(paste0('mean(2sls) : ',round(mean(ratio_est),4)))
    print(paste0('mean(gestiation) : ', round(mean(g_coef),4)))
    print(paste0('coverage(2sls) : ', sum(ratio_coverage)))
    print(paste0('coverage(gestiation) : ',sum(g_coverage)))
    
  }
  
}
