#Runs for different numbers of bridges
set.seed(1)
obs <- simSIR(x10=762, x20=1, deltat=0.01, T=15, beta=exp(-6), gamma=0.5)
plot(ts(obs[,1],start=0,deltat=0.01))
plot(ts(obs[,2],start=0,deltat=0.01))
obs<-obs[1+(0:15)*100,]
plot(ts(obs[,1],start=0,deltat=1))
plot(ts(obs[,2],start=0,deltat=1))

#N=1
set.seed(1)
out_part_synth_1 <- mh_rb(x=obs, iter=1000, sigma_tune=diag(c(0.001,0.001)), deltatau=0.2, N=1, inter_obs=1)

V <- var(out_part_synth_1)
out_part_synth_1 <- mh_rb(x=obs, iter=1000, sigma_tune=3*V, deltatau=0.2, N=1, inter_obs=1)
system.time(out_part_synth_1 <- mh_rb(x=obs, iter=10000, sigma_tune=3*V, deltatau=0.2, N=1, inter_obs=1))

plot(ts(out_part_synth_1[,1]), ylab="")
abline(h=-6, col="red", lwd=2)
plot(ts(out_part_synth_1[,2]), ylab="")
abline(h=log(0.5), col="red", lwd=2)
betahat_part_synth_1 <- exp(mean(out_part_synth_1[,1]))
gammahat_part_synth_1 <- exp(mean(out_part_synth_1[,2]))
quantile(exp(out_part_synth_1[,1]), probs=c(0.025,0.975))
quantile(exp(out_part_synth_1[,2]), probs=c(0.025,0.975))
hist(exp(out_part_synth_1[,1]), main="", xlab=expression(beta), col = "darkslategray3")
abline(v = exp(-6), col='red', lwd=2)
hist(exp(out_part_synth_1[,2]), main="", xlab=expression(gamma), col="darkseagreen3")
abline(v=0.5, col='red', lwd=2)
#For R0
betahat_part_synth_1/gammahat_part_synth_1*763
hist(exp(out_part_synth_1[,1])/exp(out_part_synth_1[,2])*763, main="", xlab="R0", col="darkgoldenrod2")
abline(v=exp(-6)/0.5*763, lwd=2, col="red")
quantile(exp(out_part_synth_1[,1])/exp(out_part_synth_1[,2])*763, probs=c(0.025,0.975))

library(coda)
effectiveSize(out_part_synth_1)

#N=2
set.seed(1)
out_full_synth_2 <- mh_rb(x=obs, iter=1000, sigma_tune=diag(c(0.001,0.001)), deltatau=0.2, N=2, inter_obs=1)

V <- var(out_full_synth_2)
out_full_synth_2 <- mh_rb(x=obs, iter=1000, sigma_tune=3*V, deltatau=0.2, N=2, inter_obs=1)
system.time(out_full_synth_2 <- mh_rb(x=obs, iter=10000, sigma_tune=3*V, deltatau=0.2, N=2, inter_obs=1))

plot(ts(out_full_synth_2[,1]), ylab="")
abline(h=-6, col="red", lwd=2)
plot(ts(out_full_synth_2[,2]), ylab="")
abline(h=log(0.5), col="red", lwd=2)
betahat_full_synth_2 <- exp(mean(out_full_synth_2[,1]))
gammahat_full_synth_2 <- exp(mean(out_full_synth_2[,2]))
quantile(exp(out_full_synth_2[,1]), probs=c(0.025,0.975))
quantile(exp(out_full_synth_2[,2]), probs=c(0.025,0.975))
hist(exp(out_full_synth_2[,1]), main="", xlab=expression(beta), col = "darkslategray3")
abline(v = exp(-6), col='red', lwd=2)
hist(exp(out_full_synth_2[,2]), main="", xlab=expression(gamma), col="darkseagreen3")
abline(v=0.5, col='red', lwd=2)
#For R0
betahat_full_synth_2/gammahat_full_synth_2*763
hist(exp(out_full_synth_2[,1])/exp(out_full_synth_2[,2])*763, main="", xlab="R0", col="darkgoldenrod2")
abline(v=exp(-6)/0.5*763, lwd=2, col="red")
quantile(exp(out_full_synth_2[,1])/exp(out_full_synth_2[,2])*763, probs=c(0.025,0.975))

effectiveSize(out_full_synth_2)


#N=3
set.seed(1)
out_full_synth_3 <- mh_rb(x=obs, iter=1000, sigma_tune=diag(c(0.001,0.001)), deltatau=0.2, N=3, inter_obs=1)

V <- var(out_full_synth_3)
out_full_synth_3 <- mh_rb(x=obs, iter=1000, sigma_tune=2*V, deltatau=0.2, N=3, inter_obs=1)
system.time(out_full_synth_3 <- mh_rb(x=obs, iter=10000, sigma_tune=2*V, deltatau=0.2, N=3, inter_obs=1))

plot(ts(out_full_synth_3[,1]), ylab="")
abline(h=-6, col="red", lwd=2)
plot(ts(out_full_synth_3[,2]), ylab="")
abline(h=log(0.5), col="red", lwd=2)
betahat_full_synth_3 <- exp(mean(out_full_synth_3[,1]))
gammahat_full_synth_3 <- exp(mean(out_full_synth_3[,2]))
quantile(exp(out_full_synth_3[,1]), probs=c(0.025,0.975))
quantile(exp(out_full_synth_3[,2]), probs=c(0.025,0.975))
hist(exp(out_full_synth_3[,1]), main="", xlab=expression(beta), col = "darkslategray3")
abline(v = exp(-6), col='red', lwd=2)
hist(exp(out_full_synth_3[,2]), main="", xlab=expression(gamma), col="darkseagreen3")
abline(v=0.5, col='red', lwd=2)
#For R0
betahat_full_synth_3/gammahat_full_synth_3*763
hist(exp(out_full_synth_3[,1])/exp(out_full_synth_3[,2])*763, main="", xlab="R0", col="darkgoldenrod2")
abline(v=exp(-6)/0.5*763, lwd=2, col="red")
quantile(exp(out_full_synth_3[,1])/exp(out_full_synth_3[,2])*763, probs=c(0.025,0.975))

effectiveSize(out_full_synth_3)

#N=4
set.seed(1)
out_full_synth_4 <- mh_rb(x=obs, iter=1000, sigma_tune=diag(c(0.001,0.001)), deltatau=0.2, N=4, inter_obs=1)

V <- var(out_full_synth_4)
out_full_synth_4 <- mh_rb(x=obs, iter=1000, sigma_tune=V, deltatau=0.2, N=4, inter_obs=1)
system.time(out_full_synth_4 <- mh_rb(x=obs, iter=10000, sigma_tune=V, deltatau=0.2, N=4, inter_obs=1))

plot(ts(out_full_synth_4[,1]), ylab="")
abline(h=-6, col="red", lwd=2)
plot(ts(out_full_synth_4[,2]), ylab="")
abline(h=log(0.5), col="red", lwd=2)
betahat_full_synth_4 <- exp(mean(out_full_synth_4[,1]))
gammahat_full_synth_4 <- exp(mean(out_full_synth_4[,2]))
quantile(exp(out_full_synth_4[,1]), probs=c(0.025,0.975))
quantile(exp(out_full_synth_4[,2]), probs=c(0.025,0.975))
hist(exp(out_full_synth_4[,1]), main="", xlab=expression(beta), col = "darkslategray3")
abline(v = exp(-6), col='red', lwd=2)
hist(exp(out_full_synth_4[,2]), main="", xlab=expression(gamma), col="darkseagreen3")
abline(v=0.5, col='red', lwd=2)
#For R0
betahat_full_synth_4/gammahat_full_synth_4*763
hist(exp(out_full_synth_4[,1])/exp(out_full_synth_4[,2])*763, main="", xlab="R0", col="darkgoldenrod2")
abline(v=exp(-6)/0.5*763, lwd=2, col="red")
quantile(exp(out_full_synth_4[,1])/exp(out_full_synth_4[,2])*763, probs=c(0.025,0.975))

effectiveSize(out_full_synth_4)


#N=5
set.seed(1)
out_full_synth_5 <- mh_rb(x=obs, iter=1000, sigma_tune=diag(c(0.001,0.001)), deltatau=0.2, N=5, inter_obs=1)
plot(ts(out_full_synth_5[,1]), ylab="")
abline(h=-6, col="red", lwd=2)
plot(ts(out_full_synth_5[,2]), ylab="")
abline(h=log(0.5), col="red", lwd=2)
V <- var(out_full_synth_5)
out_full_synth_5 <- mh_rb(x=obs, iter=1000, sigma_tune=V, deltatau=0.2, N=5, inter_obs=1)
system.time(out_full_synth_5 <- mh_rb(x=obs, iter=10000, sigma_tune=V, deltatau=0.2, N=5, inter_obs=1))
effectiveSize(out_full_synth_5)


#N=10
set.seed(1)
out_full_synth_10 <- mh_rb(x=obs, iter=1000, sigma_tune=diag(c(0.001,0.001)), deltatau=0.2, N=10, inter_obs=1)

V <- var(out_full_synth_10)
out_full_synth_10 <- mh_rb(x=obs, iter=1000, sigma_tune=3*V, deltatau=0.2, N=10, inter_obs=1)

system.time(out_full_synth_10 <- mh_rb(x=obs, iter=10000, sigma_tune=3*V, deltatau=0.2, N=10, inter_obs=1))

plot(ts(out_full_synth_10[,1]), ylab="")
abline(h=-6, col="red", lwd=2)
plot(ts(out_full_synth_10[,2]), ylab="")
abline(h=log(0.5), col="red", lwd=2)
effectiveSize(out_full_synth_10)
