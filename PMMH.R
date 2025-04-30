#Applying the PMMH algorithm to SIR model with residual bridge
library(mvtnorm)
library(MASS)

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

detSIR <- function(deltat, x0, T, theta){
  m <- T/deltat
  xmat <- matrix(0, nrow=m+1, ncol=2)
  Xt <- x0
  xmat[1,] <- Xt
  for(i in 2:(T/deltat+1)){
    Xt <- Xt + alpha(theta, Xt)*deltat
    xmat[i,] <- Xt
  }
  return(xmat)
}

rb <- function(x0, xT, theta, deltatau, T){
  beta <- theta[1]
  gamma <- theta[2]
  m <- T/deltatau
  tau <- 0
  xmat <- matrix(0, nrow=m+1, ncol=2)
  Xtau <- x0
  eta <- detSIR(deltat=deltatau, x0=x0, T=T, theta=theta)
  xmat[1,] <- Xtau
  for(i in 1:(m-1)){
    muXtau <- (eta[(i+1),]-eta[i,])/deltatau + ((xT-Xtau)-(eta[(m+1),]-eta[i,]))/(T-tau)
    psiXtau <- sqrt((T-tau-deltatau)/(T-tau)*deltatau)*betaXt(theta, Xtau)
    Wtau <- mvrnorm(1, c(0,0), diag(c(1, 1)))
    Xtau <- Xtau+muXtau*deltatau + psiXtau%*%Wtau
    Xtau <- abs(Xtau)
    xmat[(i+1),] <- Xtau
    tau <- tau + deltatau
  }
  xmat[(m+1),] <- xT
  return(xmat)
}

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

SIRpropdens_rb <- function(x, deltatau, theta, T){
  beta <- theta[1]
  gamma <- theta[2]
  len <- dim(x)[1]
  xnow <- x[1:(len-2),]
  xnext <- x[2:(len-1),]
  xT <- x[len,]
  x0 <- x[1,]
  lik <- 0
  tau <- 0
  eta <- detSIR(deltat=deltatau, x0=x0, T=T, theta=theta)
  for(i in 1:(len-2)){
    muXtau <- (eta[(i+1),]-eta[i,])/deltatau + ((xT-xnow[i,])-(eta[(len),]-eta[i,]))/(T-tau)
    lik <- lik + dmvnorm(xnext[i,],
                         mean = xnow[i,] + muXtau*deltatau,
                         sigma = (T-tau-deltatau)/(T-tau)*deltatau*
                           matrix(c(beta*xnow[i,1]*xnow[i,2], -beta*xnow[i,1]*xnow[i,2],
                                    -beta*xnow[i,1]*xnow[i,2], beta*xnow[i,1]*xnow[i,2]+gamma*xnow[i,2]),
                                  nrow=2, ncol=2, byrow=TRUE), log=TRUE)
    tau <- tau + deltatau
  }
  return(lik)
}

SIRlogwt_rb <- function(x, theta, deltatau, T){
  return(SIRemlik(x=x, deltat=deltatau, theta=theta) -
           SIRpropdens_rb(x=x, deltatau=deltatau, theta=theta, T=T))
}

#Evaluate p-hat (again on log scale)
SIRwr_rb <- function(x0, xT, theta, deltatau, T, N){
  logwtvec <- vector("numeric", N)
  for(i in 1:N){
    bridge <- rb(x0, xT, theta, deltatau, T)
    logwt <- SIRlogwt_rb(bridge, theta, deltatau, T)
    logwtvec[i] <- logwt
  }
  logp <- log(mean(exp(logwtvec)))
  return(logp)
}

SIRphatu_rb <- function(x, theta, deltatau, N, inter_obs){
  len <- dim(x)[1] - 1 
  logpvec <- vector("numeric", len)
  for(i in 1:len){
    logp <- SIRwr_rb(x0=x[i,], xT=x[i+1,], theta = theta, deltatau=deltatau, T=inter_obs, N=N)
    logpvec[i] <- logp
  }
  return(sum(logpvec))
}

lprior <- function(psi){
  return(sum(dmvnorm(psi, c(0,0), diag(c(1,1)), log = TRUE)))
}


mh_rb = function(x,iter,sigma_tune,deltatau,N,inter_obs)
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
    llikecan <- SIRphatu_rb(theta = exp(can), x = x, deltatau = deltatau, N = N, inter_obs=inter_obs)
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

Nvar_rb <- function(x, theta, deltatau, N, inter_obs, samplesize){
  probs <- vector("numeric", samplesize)
  for(i in 1:samplesize){
    probs[i] <- SIRphatu_rb(x, theta, deltatau, N, inter_obs)
  }
  var(na.omit(probs))
}

set.seed(1)
obs <- simSIR(x10=762, x20=1, deltat=0.01, T=15, beta=exp(-6), gamma=0.5)
plot(ts(obs[,1],start=0,deltat=0.01))
plot(ts(obs[,2],start=0,deltat=0.01))
obs<-obs[1+(0:15)*100,]
plot(ts(obs[,1],start=0,deltat=1))
plot(ts(obs[,2],start=0,deltat=1))

set.seed(1)
out_part_synth <- mh_rb(x=obs, iter=1000, sigma_tune=diag(c(0.001,0.001)), deltatau=0.2, N=5, inter_obs=1)
plot(ts(out_part_synth[,1]), ylab="")
abline(h=-6, col="red", lwd=2)
plot(ts(out_part_synth[,2]), ylab="")
abline(h=log(0.5), col="red", lwd=2)
betahat_part_synth <- exp(mean(out_part_synth[,1]))
gammahat_part_synth <- exp(mean(out_part_synth[,2]))
quantile(exp(out_part_synth[,1]), probs=c(0.025,0.975))
quantile(exp(out_part_synth[,2]), probs=c(0.025,0.975))
hist(exp(out_part_synth[,1]), main="", xlab=expression(beta), col = "darkslategray3")
abline(v = exp(-6), col='red', lwd=2)
hist(exp(out_part_synth[,2]), main="", xlab=expression(gamma), col="darkseagreen3")
abline(v=0.5, col='red', lwd=2)
#For R0
betahat_part_synth/gammahat_part_synth*763
hist(exp(out_part_synth[,1])/exp(out_part_synth[,2])*763, main="", xlab="R0", col="darkgoldenrod2")
abline(v=exp(-6)/0.5*763, lwd=2, col="red")
quantile(exp(out_part_synth[,1])/exp(out_part_synth[,2])*763, probs=c(0.025,0.975))

V <- var(out_part_synth)
out_part_synth <- mh_rb(x=obs, iter=10000, sigma_tune=V, deltatau=0.2, N=5, inter_obs=1)
hist(exp(out[,1]), col="purple3")

Nvec <- c(1, 2, 3, 4, 5)
varvec <- vector("numeric", 5)
for(i in 1:length(Nvec)){
  varvec[i] <- Nvar_rb(x=obs, theta=c(exp(-6),0.5), deltatau=0.2, N=Nvec[i], inter_obs=1, samplesize=5000)
  print(varvec[i])
}