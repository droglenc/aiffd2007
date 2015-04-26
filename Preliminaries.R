library(FSA)

options(show.signif.stars=FALSE,contrasts=c("contr.sum","contr.poly"))

set.seed(983452)

setwd("c:/aaaWork/web/fishR/BookVignettes/AIFFD/")

d6 <- read.table("data/box5_6.txt",header=TRUE)

library(car)

X <- rnorm(50,5)
Y <- 3*X+25+rnorm(5,2)

lm1 <- lm(Y~X)
coef(lm1)
anova(lm1)

options(contrasts=c("contr.sum","contr.poly"))

SSE <- function(theta,obs) {
  sum((obs-theta)^2)
}

normLH <- function(pars,obs) {
  mu <- pars[1]
  sig <- pars[2]
  sum(log(dnorm(obs,mean=mu,sd=sig)))
}

dat <- rnorm(1000)
mean(dat)
opt1 <- optimize(SSE,c(-2,2),obs=dat)
opt1$minimum
opt1$objective

opt2 <- optim(c(0,1),normLH,control=list(fnscale=-1),obs=dat)
opt2$par
opt2$value

opt3 <- optim(c(0,1),normLH,method="L-BFGS-B",lower=c(NA,0),upper=c(NA,2),
              control=list(fnscale=-1),obs=dat)
opt3$par
opt3$value

ages <- 0:9
Nt <- round(10000*exp(-0.2*ages)+rnorm(length(ages),800,300),0)
plot(Nt~ages,pch=19,xlab="Ages",ylab="Abundance")

edSSE <- function(pars,t,y) {  # t=age, y=abundance
  N0 <- pars[1]                # isolate N0 param from first position
  Z <- pars[2]                 # isolate Z param from second position
  sum((y-N0*exp(-Z*t))^2)      # compute and return SSE
}

par.start <- c(15000,0.3)
par.scale <- c(1e4,1e-1)
opt4 <- optim(par.start,edSSE,control=list(parscale=par.scale),t=ages,y=Nt)
opt4$par
opt4$value

plot(Nt~ages,pch=19,xlab="Ages",ylab="Abundance")
curve(opt4$par[1]*exp(-opt4$par[2]*x),c(0,9),lwd=2,lty=2,add=TRUE)


# Script created at 2015-04-26 11:23:31
