#PMMH and cPMMH toy examples
#PMMH
mh<-function(iters=10000,v=1) 
{
  vec=vector("numeric", iters)
  x=0
  targ_old=target(x)
  vec[1]=x
  a <- 0
  for (i in 2:iters) {
    innov=rnorm(1,0,v)
    can=x+innov
    targ_new=target(can)
    aprob=targ_new/targ_old
    u=runif(1)
    if (u < aprob) { 
      x=can
      targ_old=targ_new
      a <- a+1
    }
    vec[i]=x
  }
  print(a/iters)
  return(vec)
}

#Lognormal noise
target <- function(x){
  u <- rlnorm(1, -5/2, sqrt(5))
  return(dnorm(x)*u)
}

out=mh()
plot(ts(out), ylab=expression(theta))
hist(out,freq=FALSE, main="", xlab="", ylim = c(0,0.4), xlim=c(-4,4))
lines(seq(-4,4,0.001),dnorm(seq(-4,4,0.001)),type="l")

#Correlated pseudo-marginal algorithm which targets N(0,1)
mhnew<-function(iters=10000,v=1,rho=0.95,sig=1) 
{
  #Initialise u
  u <- rnorm(1,0,1)
  mat <- matrix(nrow=iters, ncol=2, dimnames=list(c(),c('x','u')))
  x=0
  targ_old=target(x,u,sig)
  mat[1,]=c(x,u)
  a <- 0
  for (i in 2:iters) {
    innov=rnorm(1,0,v)
    can=x+innov
    ucan=rnorm(1,rho*u,sqrt(1-rho^2))
    targ_new=target(can,ucan,sig)
    aprob=targ_new/targ_old
    r=runif(1)
    if (r < aprob) { 
      x=can
      u <- ucan
      targ_old=targ_new
      a <- a+1
    }
    mat[i,]=c(x,u)
  }
  print(a/iters)
  return(mat)
}


target<-function(x,u,sig)
{
  return(dnorm(x)*dnorm(u)*exp(-sig^2/2+sig*u))
}

out=mhnew(rho=0.95,sig=5)

plot(ts(out[,1]), ylab=expression(theta))
hist(out[,1], freq=FALSE, main="", xlab="", ylim = c(0,0.4), xlim=c(-4,4))
lines(seq(-4,4,0.001),dnorm(seq(-4,4,0.001)),type="l")