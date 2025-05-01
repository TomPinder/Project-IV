#PMMH for SIR using MDB
library(MASS)
library(mvtnorm)
alpha <- function(theta, Xt){
  c(-theta[1]*Xt[1]*Xt[2],
    theta[1]*Xt[1]*Xt[2]-theta[2]*Xt[2])
}

betaXt <- function(theta, Xt){
  matrix(c(-sqrt(theta[1]*Xt[1]*Xt[2]), 0,
           sqrt(theta[1]*Xt[1]*Xt[2]), -sqrt(theta[2]*Xt[2])),
         ncol=2, nrow=2, byrow=TRUE)
}

#SIR simulation
simSIR <- function(x10, x20, deltat, T, beta, gamma){
  len <- T/deltat
  xmat <- matrix(0, nrow=len+1, ncol=2)
  Xt <- c(x10, x20)
  xmat[1,] <- Xt
  for(i in 2:(len+1)){
    muXt <- c(-beta*Xt[1]*Xt[2],
              beta*Xt[1]*Xt[2]-gamma*Xt[2])
    sigmaXt <- matrix(c(-sqrt(beta*Xt[1]*Xt[2]), 0,
                        sqrt(beta*Xt[1]*Xt[2]), -sqrt(gamma*Xt[2])),
                      ncol=2, nrow=2, byrow=TRUE)
    Wt <- c(rnorm(1, 0, sqrt(deltat)), rnorm(1, 0, sqrt(deltat)))
    deltaXt <- muXt*deltat + sigmaXt%*%Wt
    Xt <- Xt + deltaXt
    Xt <- abs(Xt)
    xmat[i,] <- Xt
  }
  return(xmat)
}

#Function to generate MV bridge for SIR model
mvdg <- function(x0, xT, theta, deltatau, T){
  beta <- theta[1]
  gamma <- theta[2]
  m <- T/deltatau
  tau <- 0
  xmat <- matrix(0, nrow=m+1, ncol=2)
  Xtau <- x0
  xmat[1,] <- Xtau
  for(i in 1:(m-1)){
    muXtau <- (xT-Xtau)/(T-tau)
    psiXtau <- sqrt((T-tau-deltatau)/(T-tau)*deltatau)*betaXt(theta, Xtau)
    Wtau <- mvrnorm(1, c(0,0), diag(c(1, 1)))
    Xtau <- Xtau + muXtau*deltatau + psiXtau%*%Wtau
    Xtau <- abs(Xtau)
    xmat[(i+1),] <- Xtau
    tau <- tau + deltatau
  }
  xmat[(m+1),] <- xT
  return(xmat)
}

#Evaluate E-M likelihood for SIR model (on log scale)
SIRemlik <- function(x, deltat, theta){
  beta <- theta[1]
  gamma <- theta[2]
  len <- dim(x)[1]
  xnow <- x[1:(len-1),]
  xnext <- x[2:len,]
  lik <- 0
  for(i in 1:(len-1)){
    lik <- lik + dmvnorm(xnext[i,], 
                         mean = xnow[i,]+c(-beta*xnow[i,1]*xnow[i,2],
                                           beta*xnow[i,1]*xnow[i,2]-gamma*xnow[i,2])*deltat,
                         sigma = matrix(c(beta*xnow[i,1]*xnow[i,2], -beta*xnow[i,1]*xnow[i,2],
                                          -beta*xnow[i,1]*xnow[i,2], beta*xnow[i,1]*xnow[i,2]+gamma*xnow[i,2]),
                                        nrow=2, ncol=2, byrow=TRUE)*deltat, log = TRUE)
  }
  return(lik)
}

#Evaluate proposal density (on log scale)
SIRpropdens <- function(x, deltatau, theta, T){
  beta <- theta[1]
  gamma <- theta[2]
  len <- dim(x)[1]
  xnow <- x[1:(len-2),]
  xnext <- x[2:(len-1),]
  xT <- x[len,]
  lik <- 0
  tau <- 0
  for(i in 1:(len-2)){
    lik <- lik + dmvnorm(xnext[i,],
                         mean = xnow[i,]+(xT-xnow[i,])/(T-tau)*deltatau,
                         sigma = (T-tau-deltatau)/(T-tau)*deltatau*
                           matrix(c(beta*xnow[i,1]*xnow[i,2], -beta*xnow[i,1]*xnow[i,2],
                                    -beta*xnow[i,1]*xnow[i,2], beta*xnow[i,1]*xnow[i,2]+gamma*xnow[i,2]),
                                  nrow=2, ncol=2, byrow=TRUE), log=TRUE)
    tau <- tau + deltatau
  }
  return(lik)
}

SIRlogwt <- function(x, theta, deltatau, T){
  return(SIRemlik(x=x, deltat=deltatau, theta=theta) -
           SIRpropdens(x=x, deltatau=deltatau, theta=theta, T=T))
}

#Evaluate p-hat (again on log scale)
SIRwr <- function(x0, xT, theta, deltatau, T, N){
  logwtvec <- vector("numeric", N)
  for(i in 1:N){
    bridge <- mvdg(x0, xT, theta, deltatau, T)
    logwtvec[i] <- SIRlogwt(bridge, theta, deltatau, T)
  }
  logp <- log(mean(exp(logwtvec)))
  return(logp)
}

#Phatu function
SIRphatu <- function(x, theta, deltatau, N, inter_obs){
  len <- dim(x)[1] - 1 
  logpvec <- vector("numeric", len)
  for(i in 1:len){
    logp <- SIRwr(x0=x[i,], xT=x[i+1,], theta = theta, deltatau=deltatau, T=inter_obs, N=N)
    logpvec[i] <- logp
  }
  return(sum(logpvec))
}

lprior <- function(psi){
  return(sum(dmvnorm(psi, c(0,0), diag(c(1,1)), log = TRUE)))
}

mh = function(x,iter,sigma_tune,deltatau,N,inter_obs)
{
  mat <- matrix(0,ncol=2,nrow=iter) #store theta samples in each row
  psi <- log(c(exp(-6), 0.5)) #initialise
  mat[1,] <- psi
  a <- 0
  current_lik <- -9999
  for (i in 2:iter)
  {
    #Propose a candidate value psi
    can <- mvrnorm(n=1, mu = c(psi[1], psi[2]), Sigma = sigma_tune)
    #loglikeihood at candidate
    llikecan <- SIRphatu(theta = exp(can), x = x, deltatau = deltatau, N = N, inter_obs=inter_obs)
    #loglikelihood at current psi
    llikepsi <- current_lik
    #log of prior at candidate
    lpriorcan <- lprior(can)
    #log of prior at current psi
    lpriorpsi <- lprior(psi)
    #log acceptance probability
    laprob <- llikecan + lpriorcan - llikepsi - lpriorpsi
    u <- runif(1)
    if (log(u) < laprob)
    {
      psi <- can #Accept candidate
      current_lik <- llikecan
      a <- a + 1
    }
    mat[i,] <- psi #Update chain
  }
  print(a/iter)
  return(mat)
}

Nvar <- function(x, theta, deltatau, N, inter_obs, samplesize){
  probs <- vector("numeric", samplesize)
  for(i in 1:samplesize){
    probs[i] <- SIRphatu(x, theta, deltatau, N, inter_obs)
  }
  var(probs)
}

Nvar(x=obs, theta=c(exp(-6),0.5), deltatau=0.2, N=25, inter_obs=1, samplesize=1000)

set.seed(1)
obs <- simSIR(x10=762, x20=1, deltat=0.01, T=15, beta=exp(-6), gamma=0.5)
plot(ts(obs[,1],start=0,deltat=0.01))
plot(ts(obs[,2],start=0,deltat=0.01))
obs<-obs[1+(0:15)*100,]
plot(ts(obs[,1],start=0,deltat=1))
plot(ts(obs[,2],start=0,deltat=1))

out_dg <- mh(x=obs, iter=1000, sigma_tune=diag(c(0.01,0.01)), deltatau=0.2, N=25, inter_obs=1)
V <- var(out_dg)
system.time(out_dg_2 <- mh(x=obs, iter=10000, sigma_tune=V, deltatau=0.2, N=25, inter_obs=1))
plot(ts(out_dg_2[,1]), ylab="")
abline(h=-6, col="red", lwd=2)
plot(ts(out_dg_2[,2]), ylab="")
abline(h=log(0.5), col="red", lwd=2)
library(coda)
effectiveSize(out_dg_2)