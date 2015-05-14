library(FSA)          # Subset, capHistSum, mrOpen, removal
library(Rcapture)     # desc, closedp, profileCI, openp

setwd("C:/aaaWork/Web/fishR/BookVignettes/aiffd2007")

options(width=90,continue=" ",show.signif.stars=FALSE,
        contrasts=c("contr.sum","contr.poly"))
par(mar=c(3.5,3.5,1,1),mgp=c(2.1,0.4,0),tcl=-0.2)

y <- c(0.7,0.1,0.6,0.3,0.4,0.1,3.2,0.4,0.6,1.4,0.2,0.1,2.5,0.4,4.6,2.2,0.5,1.6,
       0.4,0.4,1.5,0.8,0.0,0.2,2.1,0.4,0.4,0.1,1.1,0.6)
( n <- length(y) )
L <- 100
A <- 500

( w <- sum(y)/(n-1) )
( D <- n/(2*L*w) )
( N <- D*A )
( V <- n/((n/N)^2)*(1-(n/N)+(n/(n-2))) )
( CI <- N + qnorm(c(0.025,0.975))*sqrt(V) )

( w <- 1/sqrt(2/(pi*sum((y^2)/n))) )
( D <- n/(2*L*w) )
( N <- D*A )
( V <- n/((n/N)^2)*(1-(n/N)+(n/(n-2))) )
( CI <- N + qnorm(c(0.025,0.975))*sqrt(V) )

dstnc <- function(y,L,A,type=c("Exponential","Half-Normal"),conf.level=0.95) {
  type <- match.arg(type)
  if (type=="Exponential") { w <- sum(y)/(n-1) }
    else { w <- 1/sqrt(2/(pi*sum((y^2)/n))) }
  D <- n/(2*L*w)
  N <- D*A
  V <- n/((n/N)^2)*(1-(n/N)+(n/(n-2)))
  CI <- N + qnorm(c((1-conf.level)/2,1-(1-conf.level)/2))*sqrt(V)
  list(type=type,w=w,D=D,N=N,V=V,CI=CI)
}

dstnc(y,L,A)                      # Exponential decline in sightability
d1 <- dstnc(y,L,A,"Half-Normal")  # Half-normal decline in sightability
round(d1$N,0)
round(d1$CI,0)

N <- 60
p <- 0.4
dbinom(20,size=N,prob=p)

catch <- seq(0,60)                         # create all values between 0 and 60
pr <- dbinom(catch,size=N,prob=p)          # find probabilities for all possible catches
( max.catch <- catch[which(pr==max(pr))] ) # find catch with maximum probability
plot(pr~catch,type="l",xlab="Number Caught",ylab="Probablity")
abline(v=max.catch,lty=2,col="red")        # mark the maximum catch
axis(1,max.catch,max.catch,col="red")

possN <- seq(0,100)                        # create possible values of N between 0 and 100
lh <- dbinom(20,size=possN,prob=p)         # find probabilities for possible Ns
( max.lh <- possN[which(lh==max(lh))] )    # find possN with maximum lh
plot(lh~possN,type="l",xlab="Population Size",ylab="Likelihood")
abline(v=max.lh,lty=2,col="red")           # mark the maximum possN

d3 <- read.table("data/Box8_3.txt",header=TRUE)
str(d3)

( d3ch <- d3[,-1] )

( desc <- descriptive(d3ch) )
plot(desc)

capHistSum(d3ch)

# for Sexpr below
res1 <- closedp(d3ch)

( res1 <- closedp(d3ch) )
res1$results

( res1a <- res1$results[c("M0","Mt","Mb"),"AIC"] )
res1a-min(res1a)+26.598

( CI1 <- profileCI(d3ch,m="M0") )

d4 <- read.table("data/Box8_4.txt",header=TRUE)
str(d4)

chsum <- capHistSum(d4,cols2use=-1)
chsum$methodB.top
chsum$methodB.bot

res1 <- mrOpen(chsum)
summary(res1,parm="N")
confint(res1,parm="N")

p <- 0.2
q <- 1-p
catch <- c(24,17,8)
n <- sum(catch)
denom <- p+p*q+p*q*q
( lh <- catch[1]*log(p/denom)+catch[2]*log(p*q/denom)+catch[3]*log(p*q*q/denom) )

crmvlLH <- function(p,catch) {
  t <- length(catch)
  geoms <- dgeom(0:(t-1),p)
  denom <- sum(geoms) 
  sum(catch*log(geoms/denom))
}

crmvlLH(0.2,catch)

ps <- seq(0.01,0.99,0.01)
plhs <- sapply(ps,crmvlLH,catch=catch)
round(plhs,1)                           # for display only

( maxplh <- max(plhs) )
( phat <- ps[which(plhs==maxplh)] )

plot(plhs~ps,type="l",xlab="Probability of Capture (p)",ylab="Log-Likelihood")
lines(rep(phat,2),c(min(plhs),maxplh),lty=2,col="red")
plot(plhs~ps,type="l",xlab="Probablity of Capture (p)",ylab="Log-Likelihood",ylim=c(-55,-49))
lines(rep(phat,2),c(min(plhs),maxplh),lty=2,col="red")  

logLikCI <- function(logLiks,vals,conf.level=0.95) {
  rel.lhs <- max(logLiks)-logLiks
  CIrng <- vals[which(rel.lhs<qchisq(conf.level,1))]
  range(CIrng)
}

( phat.CI <- logLikCI(plhs,ps) )

( Nhat <- n/(phat+phat*(1-phat)+phat*((1-phat)^2)) )
Ns <- n/(ps+ps*(1-ps)+ps*((1-ps)^2))
( Nhat.CI <- logLikCI(plhs,Ns) )
plot(plhs~Ns,type="l",xlab="Population Size",ylab="Log-Likelihood")
lines(rep(Nhat,2),c(min(plhs),maxplh),lty=2,col="red")
plot(plhs~Ns,type="l",xlab="Population Size",ylab="Log-Likelihood",
     ylim=c(-55,-49),xlim=c(50,400))
lines(rep(Nhat,2),c(min(plhs),maxplh),lty=2,col="red")

optim1 <- optimize(crmvlLH,interval=c(0,1),maximum=TRUE,catch=catch)
optim1$maximum      # p corresponding to maximum likelihood value
optim1$objective    # maximum likelihood value

optim2 <- optim(0.2,crmvlLH,catch=catch,method="L-BFGS-B",lower=0.01,upper=0.99,
                control=list(fnscale=-1))
optim2$par          # p corresponding to maximum likelihood value
optim2$value        # maximum likelihood value

rem <- removal(catch,method="Seber3")
summary(rem)
confint(rem)           # normal approximation CI

d6 <- read.table("data/Box8_6.txt",header=TRUE)
str(d6)

M <- 0.2                              # constant M
q <- 0.0001                           # initial value for q
etas <- rep(1,length(d6$catch))       # use 1s as initial values for etas
deltas <- rep(1,length(d6$catch))     # use 1s as initial values for deltas
par <- c(q,etas,deltas)               # put all parameters into one vector

ABsse <- function(par,catch,A,R,M,sse.only=TRUE)  {
  n <- length(catch)                  # get number of years of data
  q <- par[1]                         # isolate q parameter
  eta <- par[2:(n+1)]                 # isolate eta parameters (begin 2nd, n long)
  delta <- par[(n+2):(2*n+1)]         # isolate delta params (begin after eta, goto end)
  nhat <- A*eta                       # "True" adult index
  rhat <- R*delta                     # "True" recruitment index
  ntilde <- c(NA,(nhat[-n]-q*catch[-n]+rhat[-n])*exp(-M)) # compute N-tilde
  sse.1 <- sum(log10(eta)^2)          # Sum of log etas
  sse.2 <- sum(log10(delta)^2)        # Sum of log deltas
  sse.3 <- sum((A-ntilde)^2,na.rm=TRUE) # Sum of squared deviations Adult CPE
  sse <- sse.1+sse.2+sse.3            # compute overall SSE
  if (sse.only) sse
    else list(sse=sse,parts=c(sse.1,sse.2,sse.3),q,
              vals=data.frame(nhat=nhat,ntilde=ntilde,rhat=rhat,eta,delta))
}

ABsse(par,catch=d6$catch,A=d6$aIndex,R=d6$rIndex,M=M)

ABoptim1 <- optim(par,ABsse,catch=d6$catch,A=d6$aIndex,R=d6$rIndex,M=M)
ABoptim1$value                        # "minimum" SSE
ABoptim1$counts[1]                    # number of iterations
ABoptim1$convergence                  # a "1" means that maximum iterations limit was met

ABoptim2 <- optim(par,ABsse,catch=d6$catch,A=d6$aIndex,R=d6$rIndex,M=M,
                  control=list(maxit=25000))
ABoptim2$value                        # "minimum" SSE
ABoptim2$counts[1]                    # number of iterations
ABoptim2$convergence                  # a "0" means completed successfully

par.scale <- c(1e-4,rep(1,length(d6$catch)),rep(1,length(d6$catch)))  # Set parameter scales
ABoptim3 <- optim(par,ABsse,catch=d6$catch,A=d6$aIndex,R=d6$rIndex,M=M,
                  control=list(parscale=par.scale,maxit=25000))
ABoptim3$value                        # "minimum" SSE
ABoptim3$counts[1]                    # number of iterations
ABoptim3$convergence                  # a "0" means completed successfully

ABoptim4 <- optim(par,ABsse,method="BFGS",catch=d6$catch,A=d6$aIndex,R=d6$rIndex,M=M,
                  control=list(parscale=par.scale,maxit=25000))
ABoptim4$value
ABoptim4$counts
ABoptim4$convergence

ABoptim5 <- optim(par,ABsse,method="CG",catch=d6$catch,A=d6$aIndex,R=d6$rIndex,M=M,
                  control=list(parscale=par.scale,maxit=25000))
ABoptim5$value
ABoptim5$counts
ABoptim5$convergence

( qhat <- ABoptim4$par[1] )
( etahat <- ABoptim4$par[2:(length(d6$catch)+1)] )
( deltahat <- ABoptim4$par[(length(d6$catch)+2):(2*length(d6$catch)+1)] )

( Nhat <- d6$aIndex*etahat/qhat )
( Rhat <- d6$rIndex*deltahat/qhat )

res <- read.table("data/Box8_6res.txt",header=TRUE)
str(res)

parbox <- c(0.0001147,res$etahat,res$deltahat)
ABsse(parbox,catch=d6$catch,A=d6$aIndex,R=d6$rIndex,M=M)

# for Sexpr below
compTbl <- data.frame(year=d6$year,NhatBFGS=round(Nhat,0),NhatBook=res$Nhat,
                      Npdiff=(round(Nhat,0)-res$Nhat)/res$Nhat*100,RhatBFGS=round(Rhat,0),
                      RhatBook=res$Rhat,Rpdiff=(round(Rhat,0)-res$Rhat)/res$Rhat*100)

compTbl <- data.frame(year=d6$year,NhatBFGS=round(Nhat,0),NhatBook=res$Nhat,
                      Npdiff=(round(Nhat,0)-res$Nhat)/res$Nhat*100,RhatBFGS=round(Rhat,0),
                      RhatBook=res$Rhat, Rpdiff=(round(Rhat,0)-res$Rhat)/res$Rhat*100)
round(compTbl,1)

d7 <- read.table("data/Box8_7.txt",header=TRUE)
str(d7)
d7$cpe <- d7$catch/d7$effort

B0 <- 800000                                  # initial value of B0
K <- 1000000                                  # initial value of K
q <- 0.0001                                   # initial value of q
r <- 0.17                                     # initial value of r
pars <- c(B0,K,q,r)                           # put all parameters into one vector
names(pars) <- c("B0","K","q","r")            # name the parameters

  n <- length(B)                              # get number of years of data
  B0 <- par["B0"]                             # isolate B0 parameter
  K <- par["K"]                               # isolate K parameter
  q <- par["q"]                               # isolate q parameter
  r <- par["r"]                               # isolate r parameter
  predB <- numeric(n)
  predB[1] <- B0
  for (i in 2:n) predB[i] <- predB[i-1]+r*predB[i-1]*(1-predB[i-1]/K)-B[i-1]
  predCPE <- q*predB
  sse <- sum((CPE-predCPE)^2)
  if (SSE.only) sse
    else list(sse=sse,predB=predB,predCPE=predCPE)
}

( res1 <- SPsse(pars,d7$catch,d7$cpe) )

pars.scale <- c(1e5,1e6,1e-4,1e-1)            # set rescale values for parameters
SPoptim1 <- optim(pars,SPsse,control=list(maxit=10000,parscale=pars.scale),B=d7$catch,CPE=d7$cpe)
SPoptim1$value                                # "minimum" SSE
SPoptim1$counts[1]                            # number of iterations
SPoptim1$convergence                          # a "0" means completed successfully

SPoptim2 <- optim(pars,SPsse,method="BFGS",control=list(maxit=10000,parscale=pars.scale),
                  B=d7$catch,CPE=d7$cpe)
SPoptim2$value                                # "minimum" SSE
SPoptim2$counts                               # number of iterations
SPoptim2$convergence                          # a "0" means completed successfully

SPoptim3 <- optim(pars,SPsse,method="CG",control=list(maxit=10000,parscale=pars.scale),
                  B=d7$catch,CPE=d7$cpe)

SPoptim1$par                                  # the optimal parameter estimates

res3 <- SPsse(SPoptim1$par,d7$catch,d7$cpe,SSE.only=FALSE)
str(res3)
plot(cpe~year,data=d7,pch=19,xlab="Year",ylab="CPE")
lines(d7$year,res3$predCPE,lwd=2,col="red")

parsbox <- c(732506,1160771,0.0001484,0.4049) # put all parameters into one vector
names(parsbox) <- c("B0","K","q","r")         # name the parameters
res2 <- SPsse(parsbox,d7$catch,d7$cpe,SSE.only=FALSE)
res2$sse
cbind(SPoptim1$par,parsbox)

compTbl <- data.frame(year=d7$year,BNM=round(res3$predB,0),BBox=round(res2$predB,0),
                      CPENM=round(res3$predCPE,1),CPEBox=round(res2$predCPE,1))
compTbl
pdiffTbl <- data.frame(year=d7$year,Bpdiff=(compTbl$BNM-compTbl$BBox)/compTbl$BBox*100,
                       CPEpdiff=(compTbl$CPENM-compTbl$CPEBox)/compTbl$CPEBox*100)
round(pdiffTbl,1)

d8 <- read.table("data/Box8_8.txt",header=TRUE)
str(d8)
d8$logmeanw <- log(d8$meanw)
d8$var.logmeanw <- d8$var.meanw/(d8$meanw^2)

( d8a <- Subset(d8,date=="8Mar74") )
( d8b <- Subset(d8,date=="29Jul74") )

meanB <- (d8a$B+d8b$B)/2
var.meanB <- (d8a$var.B+d8b$var.B)/4             # equation 8.49
G <- d8b$logmeanw-d8a$logmeanw                   # given below equation 8.47
var.G <- d8a$var.logmeanw + d8b$var.logmeanw     # equation 8.50
P1 <- meanB*G                                    # equation 8.47
var.P1 <- var.meanB*(G^2)+var.G*(meanB^2)        # equation 8.48
cbind(d8a$age,meanB,var.meanB,G,var.G,P1,var.P1) # make a table for display purposes only

( P <- P1/0.181/1000 )
var.P <- var.P1/(0.181^2)/(1000^2)

conf.level <- 0.95
( z <- qnorm(1-(1-conf.level)/2) )
me <- z*sqrt(var.P)
P.lci <- P-me
p.uci <- P+me
round(cbind(d1$age,P,var.P,me,P.lci,p.uci),3)

Ptotal <- sum(P)
var.Ptotal <- sum(var.P)
me.Ptotal <- z*sqrt(var.Ptotal)
Ptotal.lci <- Ptotal-me.Ptotal
Ptotal.uci <- Ptotal+me.Ptotal
round(cbind(Ptotal,var.Ptotal,me.Ptotal,Ptotal.lci,Ptotal.uci),3)

plot(dens~meanw,data=d8a,type="b",pch=19,col="red",xlab="Mean Individual Weight",
     ylab="Density",xlim=c(0,max(d8a$meanw,d8b$meanw)),ylim=c(0,max(d8a$dens,d8b$dens)))
points(dens~meanw,data=d8b,type="b",pch=19,col="blue")
legend("topright",legend=c("8Mar74","29Jul74"),pch=19,col=c("red","blue"),lwd=1)

d9 <- read.table("data/Box8_9.txt",header=TRUE)
str(d9)

c <- length(d9$meanN)
CPI <- 3                                                              

temp1 <- d9$meanN[1:(c-2)]
temp2 <- d9$meanN[3:c]
Ndiff2 <- temp1-temp2
Ndiff1 <- d9$meanN[1]-d9$meanN[2]
Ndiff3 <- d9$meanN[c-1]-d9$meanN[c]
( Ndiff <- c(Ndiff1,Ndiff2,Ndiff3) )

( Phat1 <- 0.5*c*sum(d9$meanw*Ndiff)/CPI )

 # note the "plus" in the first term below
wdiff <- c(d9$meanw[1]+d9$meanw[2],d9$meanw[1:(c-2)]-d9$meanw[3:c],d9$meanw[c-1]-d9$meanw[c])
( var.Phat1 <- ((0.5*c)^2)*sum((wdiff^2)*d9$varN+d9$varw*(Ndiff^2))/(CPI^2) )

Phat <- Phat1/1000
var.Phat <- var.Phat1/(1000^2)

conf.level <- 0.95
( z <- qnorm(1-(1-conf.level)/2) )
me <- z*sqrt(var.Phat)
P.lci <- Phat-me
p.uci <- Phat+me
round(cbind(Phat,var.Phat,me,P.lci,p.uci),3)


# Script created at 2015-05-13 21:31:15
