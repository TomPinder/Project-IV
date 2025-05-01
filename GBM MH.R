#GBM Metropolis-Hastings
library(MASS)
#Simulation functions
EM <- function(x0, T, deltat, mu, sigma){
  N <- T/deltat
  vec <- vector("numeric", N)
  x <- x0
  vec[1] <- x
  for(n in 1:N){
    x <- x + mu*x*deltat + sigma*sqrt(x*deltat)*rnorm(1)
    x <- abs(x)
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


#Log likelihoods
exact_lik <- function(psi, x = out, deltat){
  xnow <- x[1:(length(x)-1)]
  xnext <- x[2:length(x)]
  mu <- exp(psi[1])
  sigma <- exp(psi[2])
  return(sum(dlnorm(xnext, 
                    log(xnow) + (mu-sigma^2/2)*deltat,
                    sigma*sqrt(deltat), 
                    log = TRUE)))
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

lprior <- function(psi){
  return(sum(dnorm(psi, 0, 1, log = TRUE)))
}


mh = function(x,N,sigma,deltat,lik)
{
  mat <- matrix(0,ncol=2,nrow=N) #store theta samples in each row
  psi <- c(0,0) #initialise
  mat[1,] <- psi
  a <- 0
  for (i in 2:N)
  {
    #Propose a candidate value psi
    can <- mvrnorm(n=1, mu = c(psi[1], psi[2]), Sigma = sigma)
    #loglikeihood at candidate
    llikecan <- lik(psi = can, x = x, deltat = deltat)
    #loglikelihood at current psi
    llikepsi <- lik(psi = psi, x = x, deltat = deltat)
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
      a <- a + 1
    }
    mat[i,] <- psi #Update chain
  }
  print(a/N)
  return(mat)
}

#Using exact likelihood
set.seed(2)
out <- exact(x0 = 1, T = 20, deltat = 1, mu = 0.3, sigma = 0.5)
plot(ts(out, start=0), ylab = "Xt", lwd=2)
out_exact<- mh(x = out, N = 10000, sigma = diag(c(1,1)), 
               deltat = 1, lik = exact_lik)
V <- var(out_exact)
sqrt(V)
out_exact <- mh(x = out, N = 10000, sigma = 3*V, 
                deltat = 1, lik = exact_lik)
exp(mean(out_exact[,1]))
exp(mean(out_exact[,2]))
par(mfrow = c(1,2))
plot(ts(out_exact[,1]))
plot(ts(out_exact[,2]))
par(mfrow = c(1,1))
quantile(exp(out_exact[,1]), probs = c(0.025, 0.975))
quantile(exp(out_exact[,2]), probs = c(0.025, 0.975))
hist(exp(out_exact[,2]), col='darkseagreen3', main = "Density of sigma", xlab= "(Exact Likelihood)")
abline(v = 0.5, col = "red")
hist(exp(out_exact[,1]), col='darkslategray3', main = "Density of mu", xlab= "(Exact Likelihood)", breaks=18)
abline(v = 0.3, col = "red")

#Using EM likelihood
set.seed(2)
out <- exact(x0 = 1, T = 20, deltat = 1, mu = 0.3, sigma = 0.5)
plot(ts(out))
out2 <- mh(x = out, N = 10000, sigma = diag(c(1,1)), 
           deltat = 1, lik = em_lik)
V <- var(out2)
sqrt(V)
out2 <- mh(x = out, N = 10000, sigma = 3*V, 
           deltat = 1, lik = em_lik)
exp(mean(out2[,1]))
exp(mean(out2[,2]))
par(mfrow = c(1,2))
plot(ts(out2[,1]), ylab = "log(mu)", main = "Mixing for E-M Likelihood")
plot(ts(out2[,2]), ylab = "log(sigma)", main = "Mixing for E-M Likelihood")
par(mfrow = c(1,1))
quantile(exp(out2[,1]), probs = c(0.025, 0.975))
quantile(exp(out2[,2]), probs = c(0.025, 0.975))
hist(exp(out2[,2]), col='darkseagreen3', main = "Density of sigma", xlab= "(Approximate Likelihood)")
abline(v = 0.5, col = "red")
hist(exp(out2[,1]), col='darkslategray3', main = "Density of mu", xlab= "(Approximate Likelihood)", breaks=18)
abline(v = 0.3, col = "red")

#Smaller timestep
set.seed(1)
out <- exact(x0 = 1, T = 20, deltat = 0.1, mu = 0.3, sigma = 0.5)
plot(ts(out, start=0, deltat=0.1), xlab="Time", ylab='Xt', lwd=2)
outsmall <- mh(x = out, N = 10000, sigma = diag(c(1,1)), 
               deltat = 0.1, lik = em_lik)
V <- var(outsmall)
outsmall2 <- mh(x = out, N = 10000, sigma = 3*V, 
                deltat = 0.1, lik = em_lik)
mean(exp(outsmall2[,1]))
mean(exp(outsmall2[,2]))
hist(exp(outsmall2[,1]), col="darkslategray3", main = "Density of mu", xlab="")
abline(v=0.3,col=2)
hist(exp(outsmall2[,2]), col="darkseagreen3", main="Density of sigma", xlab="", breaks=45, xlim=c(0.37,0.6))
abline(v=0.5, col=2)
quantile(exp(outsmall2[,1]), probs=c(0.025,0.975))
quantile(exp(outsmall2[,2]), probs=c(0.025,0.975))