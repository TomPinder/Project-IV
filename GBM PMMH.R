#PMMH for GBM
library(MASS)
library(mvtnorm)

EM <- function(x0, T, deltat, mu, sigma){
  N <- T/deltat
  vec <- vector("numeric", N)
  x <- x0
  vec[1] <- x
  for(n in 1:N){
    x <- x + mu*x*deltat + sigma*sqrt(abs(x)*deltat)*rnorm(1)
    vec[n+1] <- x
  }
  return(vec)
}

exact <- function(x0=1, T=10, deltat=0.1, mu, sigma){
  N <- T/deltat
  vec <- vector("numeric", N)
  x <- x0
  t <- 0
  vec[1] <- x
  for(n in 1:N){
    t <- t + deltat
    x <- rlnorm(1, log(x) + (mu - sigma^2/2)*deltat, sigma*sqrt(deltat))
    vec[n+1] <- x
  }
  return(vec)
}

dg <- function(T, deltatau, x0, xT, sigma){
  m <- T/deltatau
  tau <- 0
  xvec <- vector("numeric", m+1)
  xvec[1] <- x0
  for(i in 1:(m-1)){
    x <- rnorm(1, xvec[i] + (xT-xvec[i])/(T-tau)*deltatau, 
               sqrt((T-tau-deltatau)/(T-tau)*deltatau*xvec[i]^2*sigma^2))
    tau <- tau + deltatau
    xvec[i+1] <- x
  }
  xvec[m+1] <- xT
  return(xvec)
}

em_lik <- function(psi, x, deltat){
  xnow <- x[1:(length(x)-1)]
  xnext <- x[2:length(x)]
  mu <- exp(psi[1])
  sigma <- exp(psi[2])
  return(sum(dnorm(xnext, xnow + mu*xnow*deltat, 
                   sigma*sqrt(xnow^2*deltat),
                   log = TRUE)))
}

prop_dens <- function(x, sigma, deltatau, T){
  xT <- x[length(x)]
  xnow <- x[1:(length(x)-2)] 
  xnext <- x[2:(length(x)-1)]
  tau <- seq(0, T-2*deltatau, deltatau)
  
  return(sum(dnorm(xnext, xnow + (xT-xnow)/(T-tau)*deltatau, 
                   sqrt((T-tau-deltatau)/(T-tau)*deltatau*xnow^2*sigma^2),
                   log = TRUE)))
}

logwt <- function(theta, x, deltatau, T){
  psi <- log(theta)
  sigma <- theta[2]
  return(em_lik(psi, x, deltatau) - prop_dens(x, sigma, deltatau, T))
}

wr <- function(N=100, T=1, deltatau=0.1, x0, xT, theta){
  xmat <- matrix(nrow = T/deltatau+1, ncol = N)
  logwtvec <- vector("numeric", N)
  for(i in 1:N){
    x <- dg(T, deltatau, x0, xT, theta[2])
    xmat[,i] <- x
    logwt <- logwt(theta, x, deltatau, T)
    logwtvec[i] <- logwt
  }
  xsampleindex <- sample(1:N, N, TRUE, exp(logwtvec))
  xsample <- xmat[,xsampleindex]
  logp <- log(mean(exp(logwtvec)))
  return(list(xmat = xmat, logwtvec = logwtvec, 
              xsample = xsample, xsampleindex = xsampleindex,
              logp = logp))
}

phatu <- function(x, theta, deltatau, N){
  len <- length(x) - 1 #-1 since x includes x0
  logpvec <- vector("numeric", len)
  for(i in 1:len){
    logp <- wr(N=N, T=1, deltatau, x0=x[i], xT=x[i+1], theta = theta)$logp
    logpvec[i] <- logp
  }
  return(sum(logpvec))
}

lprior <- function(psi){
  return(sum(dmvnorm(psi, c(0,0), diag(c(1,1)), log = TRUE)))
}

u_lik <- function(x, psi, deltatau = 0.2, N){
  theta <- exp(psi)
  return(phatu(x, theta, deltatau, N))
}

mh = function(x,iter,sigma_tune,deltatau,N)
{
  mat <- matrix(0,ncol=2,nrow=iter) #store theta samples in each row
  psi <- log(c(0.5,0.5)) #initialise
  mat[1,] <- psi
  a <- 0
  current_lik <- -9999
  for (i in 2:iter)
  {
    #Propose a candidate value psi
    can <- mvrnorm(n=1, mu = c(psi[1], psi[2]), Sigma = sigma_tune)
    #loglikeihood at candidate
    llikecan <- u_lik(psi = can, x = x, deltatau = deltatau, N = N)
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

set.seed(1)
obs <- exact(1, 15, 0.01, 0.3, 0.5)
plot(ts(obs,start=0,deltat=0.01))
obs<-obs[1+(0:15)*100]
plot(ts(obs,start=0,deltat=1), lwd=2, ylab="Xt")

Nvar_GBM <- function(x, theta, deltatau, N, samplesize){
  probs <- vector("numeric", samplesize)
  for(i in 1:samplesize){
    probs[i] <- phatu(x, theta, deltatau, N)
  }
  var(probs)
}
Nvar_GBM(obs, theta=c(0.3,0.5), deltatau=0.2, N=3, samplesize=10000)

out <- mh(obs, 10000, sigma_tune = diag(c(0.1, 0.05)), deltatau = 0.2, N=3)
V <- var(out)
out <- mh(obs, 10000, sigma_tune = 0.5*V, deltatau = 0.2, N=3)
mean(exp(out[,1]))
mean(exp(out[,2]))
quantile(exp(out[,1]), probs = c(0.025, 0.5, 0.975))
quantile(exp(out[,2]), probs = c(0.025, 0.5, 0.975))
plot(ts(exp(out[,2])), ylab=expression(sigma))
plot(ts(exp(out[,1])), ylab=expression(mu))
hist(exp(out[,1]), main='', xlab=expression(mu), col="darkslategray3")
abline(v=0.3, col="red", lwd=2)
hist(exp(out[,2]), main='', xlab=expression(sigma), col="darkseagreen3")
abline(v=0.5, col = "red", lwd=2)