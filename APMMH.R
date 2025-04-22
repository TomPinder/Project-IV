library(mvtnorm)
library(MASS)
#rb into the inference
alpha <- function(theta, Xt){
  c(-theta[1]*Xt[1]*Xt[2],
    theta[1]*Xt[1]*Xt[2]-theta[2]*Xt[2])
}


betaXt <- function(theta, Xt){
  matrix(c(-sqrt(theta[1]*Xt[1]*Xt[2]), 0,
           sqrt(theta[1]*Xt[1]*Xt[2]), -sqrt(theta[2]*Xt[2])),
         ncol=2, nrow=2, byrow=TRUE)
}


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




#Updating the mh function to handle the unobserved case
mh_rb_unobs = function(x_obs, iter, sigma_theta, sigma_x_u, deltatau, N, inter_obs,
                       x_u_init)
{
  M <- length(x_obs)
  mat <- matrix(0,ncol=2,nrow=iter) #store theta samples in each row
  x_u_mat <- matrix(0, ncol=M, nrow=iter) #store x_u in each row
  psi <- log(c(exp(-6), 0.5)) #initialise
  mat[1,] <- psi
  a <- 0
  a_x_u <- 0
  current_lik_theta <- -9999
  current_phat_x_u <- c(0, rep(-9999, M-1))
  phat_x_u_can <- vector("numeric", M)
  phat_x_u_can[1] <- 0
  x_u <- x_u_init
  x_u_mat[1,] <- x_u
  x <- cbind(x_u, x_obs)
  for (i in 2:iter)
  {
    #Propose a candidate value psi
    can <- mvrnorm(n=1, mu = c(psi[1], psi[2]), Sigma = sigma_theta)
    #loglikeihood at candidate
    for(m in 2:M){
      phat_x_u_can[m] <- SIRphatu_rb(x=x[(m-1):m,], theta=exp(can), deltatau=deltatau, N=N, inter_obs=inter_obs)
    }
    llikecan <- sum(phat_x_u_can)
    #loglikelihood at current psi
    llikepsi <- current_lik_theta
    #log of prior at candidate
    lpriorcan <- lprior(can)
    #log of prior at current psi
    lpriorpsi <- lprior(psi)
    #log acceptance probability
    laprob <- llikecan + lpriorcan - llikepsi - lpriorpsi
    accept <- runif(1)
    if (log(accept) < laprob)
    {
      psi <- can #Accept candidate
      current_lik_theta <- llikecan
      current_phat_x_u <- phat_x_u_can
      a <- a + 1
    }
    mat[i,] <- psi #Update chain
    
    #Step for unobserved x for times 1 to T-1
    for(m in 2:(M-1)){
      #Propose x_u at observation point m
      x_u_can <- abs(rnorm(1, mean = x_u[m], sd=sigma_x_u))
      #Make a matrix of required x's with x_u_can inserted at point m 
      x_can <- cbind(c(x_u[m-1], x_u_can, x_u[m+1]), x_obs[(m-1):(m+1)])
      #Likelihood for x_u_can
      phat_x_u_can[m] <- SIRphatu_rb(theta=exp(psi), x=x_can[1:2,], deltatau=deltatau, N=N, inter_obs=inter_obs)
      phat_x_u_can[m+1] <- SIRphatu_rb(theta=exp(psi), x=x_can[2:3,], deltatau=deltatau, N=N, inter_obs=inter_obs)
      llike_x_u_can <- sum(phat_x_u_can[m:(m+1)])
      #Likelihood for x_u
      llike_x_u <- sum(current_phat_x_u[m:(m+1)])
      laprob_x_u <- llike_x_u_can - llike_x_u
      accept <- runif(1)
      if(log(accept)<laprob_x_u){
        x_u[m] <- x_u_can
        current_phat_x_u[m] <- phat_x_u_can[m]
        current_phat_x_u[m+1] <- phat_x_u_can[m+1]
        current_lik_theta <- sum(current_phat_x_u)
        a_x_u <- a_x_u + 1
      }
    }
    #For unobserved x at time T
    x_u_can <- abs(rnorm(1, mean=x_u[M], sd=sigma_x_u))
    x_can <- cbind(c(x_u[M-1], x_u_can), x_obs[(M-1):(M)])
    llike_x_u_can <- SIRphatu_rb(theta=exp(psi), x=x_can, deltatau=deltatau, N=N, inter_obs=inter_obs)
    laprob_x_u <- llike_x_u_can - current_phat_x_u[M]
    accept <- runif(1)
    if(log(accept)<laprob_x_u){
      x_u[M] <- x_u_can
      current_phat_x_u[M] <- llike_x_u_can
      current_lik_theta <- sum(current_phat_x_u)
      a_x_u <- a_x_u + 1
    }
    x_u_mat[i,] <- x_u
  }
  print(a/iter)
  print(a_x_u/(iter*(M-1)))
  return(list(psi=mat, x_u=x_u_mat))
}

set.seed(1)
sim <- simSIR(762, 1, 0.01, 15, exp(-6), 0.5)
plot(ts(sim))
sim <- sim[1+(0:15)*100,]
plot(ts(sim))
x_u_prop <- simSIR(762, 1, 0.01, 15, exp(-6), 0.5)[,1]
x_u_prop <- x_u_prop[1+(0:15)*100]
plot(ts(x_u_prop, start=0), lwd=2, col=2, ylab="Number of Susceptibles")
lines(ts(x_u_prop2, start=0), lwd=2, col="orange")
lines(ts(sim[,1], start=0), lwd=2, col=3)

out_part1 <- mh_rb_unobs(x_obs=sim[,2], iter=500, sigma_theta=diag(c(0.001, 0.001)),
                   sigma_x_u=4, deltatau=0.2, N=5, inter_obs=1,
                   x_u_init=x_u_prop)
x_u_prop2 <- colMeans(out_part1$x_u[200:500,])
out_part1 <- mh_rb_unobs(x_obs=sim[,2], iter=1000, sigma_theta=diag(c(0.001, 0.001)),
                         sigma_x_u=4, deltatau=0.2, N=5, inter_obs=1,
                         x_u_init=x_u_prop2)
V <- var(out_part1$psi)
var(out_part1$x_u)
out_part2 <- mh_rb_unobs(x_obs=sim[,2], iter=10000, sigma_theta=3.8*V,
                         sigma_x_u=8, deltatau=0.2, N=5, inter_obs=1,
                         x_u_init=x_u_prop2)
out_part3 <- mh_rb_unobs(x_obs=sim[,2], iter=1000, sigma_theta=3.8*V,
                         sigma_x_u=8, deltatau=0.2, N=5, inter_obs=1,
                         x_u_init=x_u_prop2)
plot(ts(out_part3$psi[,1]), ylab="log(beta)")
abline(h=-6, col=2, lwd=2)

plot(ts(exp(out_part2$psi)))
plot(ts(out_part2$psi[,1]), ylab="log(beta)")
abline(h=-6, col=2, lwd=2)
plot(ts(out_part2$psi[,2]), ylab="log(gamma)")
abline(h=log(0.5), col=2, lwd=2)

hist(exp(out_part2$psi[,1]), main="", xlab=expression(beta), col="darkslategray3")
abline(v=exp(-6), col=2, lwd=2)
hist(exp(out_part2$psi[,2]), main="", xlab=expression(gamma), col="darkseagreen3")
abline(v=0.5, col=2, lwd=2)
hist(exp(out_part2$psi[,1])/exp(out_part2$psi[,2])*763, main="", xlab="R0", col="darkgoldenrod2")
abline(v=exp(-6)/0.5*763, lwd=2, col="red")

x_u_synth_cred <- x_credint(out_part2$x_u)
plot(ts(x_u_synth_cred[,1], start=0), lwd=2, ylab="Number of susceptibles")
lines(ts(x_u_synth_cred[,2], start=0), lwd=2)
lines(ts(x_u_synth_cred[,3], start=0), lwd=2)
lines(ts(sim[,1], start=0), lwd=2, col=2)

quantile(exp(out_part2$psi[,1]), probs=c(0.025,0.975))
quantile(exp(out_part2$psi[,2]), probs=c(0.025,0.975))
quantile(exp(out_part2$psi[,1])/exp(out_part2$psi[,2])*763, probs=c(0.025,0.975))
mean(exp(out_part2$psi[,1]))
mean(exp(out_part2$psi[,2]))
mean(exp(out_part2$psi[,1])/exp(out_part2$psi[,2])*763)

plot(ts(out_part2$x_u[,2:11]))
plot(ts(out_part2$x_u[,12:16]))
mean(out_part2$psi[,1])
mean(exp(out_part2$psi[,2]))
plot(hist(exp(out_part2$psi[,1])))
abline(v=exp(-6), lwd=2, col=2)
plot(hist(exp(out$psi[,2])))
abline(v=0.5, lwd=2, col=2)

V <- var(out$psi)
out <- mh_rb_unobs(sim[,2], iter=10000, sigma_theta=2*V,
                   sigma_x_u=5, deltatau=0.2, N=10, inter_obs=1,
                   x_u_init=sim[,1])
plot(ts(exp(out$psi)))
plot(ts(out$x_u[,2:11]))
plot(ts(out$x_u[,12:16]))
mean(out$psi[,1])
mean(exp(out$psi[,2]))
plot(hist(exp(out$psi[,1])))
abline(v=exp(-6), lwd=2, col=2)
plot(hist(exp(out$psi[,2])))
abline(v=0.5, lwd=2, col=2)


x_credint <- function(x){
  ci <- matrix(0, ncol=3, nrow=dim(x)[2])
  for(i in 1:dim(x)[2]){
    ci[i,] <- quantile(x[,i], probs = c(0.025, 0.5, 0.975))
  }
  return(ci)
}
x_u_credint <- x_credint(out$x_u)
x_u_credint
plot(ts(x_u_credint[,3]))
lines(ts(x_u_credint[,2]))
lines(ts(x_u_credint[,1]))

lines(ts(sim[,1]))


#Now using Fuchs influenza dataset
flu <- c(1, 3, 6, 25, 73, 221, 294, 257, 236, 189, 125,
         67, 26, 10, 3)
plot(ts(flu), lwd=2, ylab="Number confined to bed", col=2)
lines(ts(sim[,2]), lwd=2, col=3)


set.seed(1)
sim <- simSIR(762, 1, 0.01, 14, exp(-6), 0.5)
sim <- sim[1+(0:14)*100,]
plot(ts(sim))

plot(ts(sim[,2], start=0))
lines(ts(flu, start=0))
x_propflu <- sim[,1]
outflu <- mh_rb_unobs(x_obs=flu, iter=2000, sigma_theta=diag(c(0.01, 0.01)),
                      sigma_x_u=5, deltatau=0.2, N=5, inter_obs=1,
                      x_u_init=x_propflu)

x_propflu2 <- colMeans(outflu$x_u[800:2000,])
outflu2 <- mh_rb_unobs(x_obs=flu, iter=1000, sigma_theta=diag(c(0.01, 0.01)),
                       sigma_x_u=5, deltatau=0.2, N=5, inter_obs=1,
                       x_u_init=x_propflu2)
V <- var(outflu2$psi)
outflu3 <- mh_rb_unobs(x_obs=flu, iter=10000, sigma_theta=3.8*V,
                       sigma_x_u=7, deltatau=0.2, N=5, inter_obs=1,
                       x_u_init=x_propflu2)
var(outflu3$x_u)

Nvar_rb(x=cbind(x_propflu2, flu), theta=c(exp(-6),0.5), deltatau=0.2, N=5, inter_obs=1, samplesize=5000)

plot(ts(outflu3$psi[,1]), ylab="log(beta)")
plot(ts(outflu3$psi[,2]), ylab="log(gamma)")
plot(ts(outflu3$x_u[,2:11]))
plot(ts(outflu3$x_u[,12:15]))

hist(exp(outflu3$psi[,1]), xlab=expression(beta), col="darkslategray3", main="")
hist(exp(outflu3$psi[,2]), xlab=expression(gamma), col="darkseagreen3", main="")
hist(exp(outflu3$psi[,1])/exp(outflu3$psi[,2])*763, col="darkgoldenrod2", xlab="R0", main="")
plot(density(outflu3$psi[,1]))

quantile(exp(outflu3$psi[,1]), probs=c(0.025,0.975))
quantile(exp(outflu3$psi[,2]), probs=c(0.025,0.975))
quantile(exp(outflu3$psi[,1])/exp(outflu3$psi[,2])*763, probs=c(0.025,0.975))

mean(exp(outflu3$psi[,1]))
mean(exp(outflu3$psi[,2]))
mean(exp(outflu3$psi[,1])/exp(outflu3$psi[,2])*763)

x_u_flu_credint <- x_credint(outflu3$x_u)
x_u_flu_credint
plot(ts(x_u_flu_credint[,3]), ylim = c(0,765), ylab="Number of susceptibles", lwd=2)
lines(ts(x_u_flu_credint[,2]), lwd=2)
lines(ts(x_u_flu_credint[,1]), lwd=2)

mean(exp(outflu3$psi[,1]))
mean(exp(outflu3$psi[,2]))
hist(exp(outflu3$psi[,1]))
abline(v=1.85/763)
hist(exp(outflu3$psi[,2]))
abline(v=0.49)

V <- var(outflu$psi)
outflu2 <- mh_rb_unobs(x_obs=flu, iter=10000, sigma_theta=2*V,
                      sigma_x_u=5, deltatau=0.2, N=10, inter_obs=1,
                      x_u_init=x_prop)
plot(ts(exp(outflu2$psi)))
plot(ts(outflu2$x_u[,2:11]))
plot(ts(outflu2$x_u[,12:15]))
betahat <- mean(exp(outflu2$psi[,1]))
gammahat <- mean(exp(outflu2$psi[,2]))
x_u_flu_credint <- x_credint(outflu3$x_u)
x_u_flu_credint
plot(ts(x_u_flu_credint[,3]), ylim = c(0,765))
lines(ts(x_u_flu_credint[,2]))
lines(ts(x_u_flu_credint[,1]))
lines(flu)
hist(exp(outflu2$psi[,1]), main = "Density of beta")
hist(exp(outflu2$psi[,2]), main = "Density of gamma")
betahat/gammahat*763
betahat*763

quantile(exp(outflu2$psi[,1]), probs=c(0.025, 0.975))
quantile(exp(outflu3$psi[,2]), probs=c(0.025, 0.975))

#Multiplying my estimate by 763 so comparable with the Fuchs estimate:
quantile(763*exp(outflu2$psi[,1]), probs=c(0.025, 0.975))
#Next would like to initialise theta somewhere else and see if I converge on the same values


