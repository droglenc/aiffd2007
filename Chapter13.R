library(FSA)          # fitPlot, hoCoef,lagratio, popSizesPlot, rcumsum, srFuns, srStarts
library(TTR)          # SMA
library(popbio)       # eigenanalysis

setwd("C:/aaaWork/Web/fishR/BookVignettes/aiffd2007")

options(width=90,continue=" ",show.signif.stars=FALSE, contrasts=c("contr.sum","contr.poly"))
par(mar=c(3.5,3.5,1,1),mgp=c(2.1,0.4,0),tcl=-0.2)

( x <- 0:4 )
( k <- length(x) )
( nx <- c(52000,2118,779,159,13) )
( mx <- c(0,0,43,122.6,346.2) )

( dx <- -1*diff(c(nx,0)) )
( qx <- dx/nx )
( px <- 1-qx )
( lx <- nx/nx[1] )

( Lx <- SMA(c(nx,0),n=2)[-1] )

( Tx <- rcumsum(Lx) )
( ex <- Tx/nx )
# rcumsum includes the full last value, subtract half to effectively have added half
( Sx <- rcumsum(lx)-0.5*lx[k] )

ax <- rep(0,k)
# last entry below should be NA
( Wx <- c((qx[1:(k-1)]*Sx[2:k]^2)/(px[1:(k-1)]*(nx[1:(k-1)]-0.5*ax[1:(k-1)])),NA) )

cbind(x,nx,dx,lx,qx,px,mx,Lx,Tx,Sx,Wx,ex)  # for display only

( GRR <- sum(mx) )
( R0 <- sum(lx*mx) )
( G <- sum(x*lx*mx)/R0 )
( var.e0 <- sum(Wx,na.rm=TRUE) )
( CI.e0 <- ex[1]+c(-1,1)*qt(0.975,nx[1]-1)*sqrt(var.e0) )

EL <- function(r,x,lx,mx) sum(exp(-r*x)*lx*mx)-1
( r.init <- log(R0)/G )
( r.final <- uniroot(EL,c(0.03,0.05),x,lx,mx)$root )

rs <- seq(0.03,0.05,0.00001)
ELs <- sapply(rs,EL,x,lx,mx)
plot(ELs~rs,type="l",xlab="Possible r",ylab="Euler-Lotka Value")
abline(h=0,lty=3)
abline(v=r.final,lwd=2,lty=3,col="red")
plot(ELs~rs,type="l",xlab="Possible r",ylab="Euler-Lotka Value",xlim=c(0.039,0.041))
abline(h=0,lty=3)
abline(v=r.final,lwd=2,lty=3,col="red")

( lambda <- exp(r.final) )
( DT <- log(2)/r.final )
Cx.int <- lambda^(-x)*lx
( Cx <- Cx.int/sum(Cx.int) )

lt2lm <- function(x,lx,mx) {
  k <- length(x)  
  px <- c(lagratio(lx),0)
  fx <- px[1:(k-1)]*mx[2:k]
  M <- cbind(rbind(fx,diag(px[-k],nrow=(k-1))),rep(0,k))
  rownames(M) <- colnames(M) <- x
  M
}

( M <- lt2lm(x,lx,mx) )
M.eig <- eigen.analysis(M)
M.eig$lambda1
M.eig$stable.stage
log(M.eig$lambda1)  #shows this lambda is similar to r from above

M <- matrix(c(0,5,10,0.5,0,0,0,0.2,0),nrow=3,byrow=TRUE)
colnames(M) <- rownames(M) <- 1:3       # 1,2,3 for ages 1-3
M

n0 <- c(0,0,10)
names(n0) <- 1:3     # not really needed, but may help
n0

( n1 <- M %*% n0 )
( n2 <- M %*% n1 )

exp0 <- cbind(n0,n1,n2)                            # first three values
for (i in 3:7) exp0 <- cbind(exp0,M %*% exp0[,i])  # loop through next five
colnames(exp0) <- 1:8                              # name for time steps
exp0                                               # take a look

( N_exp0 <- apply(exp0,2,sum) )

M.eig <- eigen(M)
M.eig$values[1]
M.eig$vectors[,1]/sum(M.eig$vectors[,1])

exp1 <- pop.projection(M,n0,8)
exp1$stage.vectors
exp1$lambda
exp1$stable.stage

M.eig1 <- eigen.analysis(M)
M.eig1$lambda1
M.eig1$stable.stage

stable.stage(M)
generation.time(M)
net.reproductive.rate(M)

( M.elast <- elasticity(M) )
image2(M.elast)

exp1 <- pop.projection(M,n0,30)
exp1$lambda
exp1$stable.stage

popSizesPlot(exp1)

stage.vector.plot(exp1$stage.vector)
abline(h=stable.stage(M),lwd=2,lty=3,col="red")

# for Sexpr below
( M2 <- matrix(c(0,5,10,0.4,0,0,0,0.2,0),nrow=3,byrow=TRUE) )
M2.eig <- eigen.analysis(M2)
M2.eig$lambda1

( M2 <- matrix(c(0,5,10,0.4,0,0,0,0.2,0),nrow=3,byrow=TRUE) )
M2.eig <- eigen.analysis(M2)
M2.eig$lambda1
M2.eig$stable.stage

( M3 <- matrix(c(0,15,12,0.5,0,0,0,0.2,0),nrow=3,byrow=TRUE) )
M3.eig <- eigen.analysis(M3)
M3.eig$lambda1

( M3 <- matrix(c(0,15,12,0.5,0,0,0,0.2,0),nrow=3,byrow=TRUE) )
M3.eig <- eigen.analysis(M3)
M3.eig$lambda1
M3.eig$stable.stage

popSizesPlot(exp1)                                            # base case
expM2 <- pop.projection(M2,n0,30)
popSizesPlot(expM2,add=TRUE,col="red")                        # pollution scenario
expM3 <- pop.projection(M3,n0,30)
popSizesPlot(expM3,add=TRUE,col="blue")                       # stocking scenario
legend("topleft",legend=c("base","pollution","stocking"),col=c("black","red","blue"),lwd=2)
popSizesPlot(exp1,use.log=TRUE)
popSizesPlot(expM2,add=TRUE,use.log=TRUE,col="red")
popSizesPlot(expM3,add=TRUE,use.log=TRUE,col="blue")

( survs <- seq(0.2,0.5,0.02) )
iters <- length(survs)                         # get the number of iterations
rs <- numeric(iters)                           # initiate vector to hold r values
for (i in 1:iters) {                           # loop i from 1 to number of iterations
  Mtemp <- M                                   # put M matrix into Mtemp
  Mtemp[3,2] <- survs[i]                       # replace with current survival rate
  rs[i] <- log(eigen.analysis(Mtemp)$lambda1)  # find r and put in vector of rs
}
rs                                             # see if vector is filled

plot(rs~survs,type="b",pch=19,xlab="Age-1 survival probability",
     ylab="Rate of population increase (r)")
abline(h=log(eigen.analysis(M)$lambda1),lty=2,col="red")

year <- c(NA,1986:1990)                        # year labels
ct <- c(NA,90000,113300,155860,181128,198584)  # catches per year
iters <- length(year)-1                        # number of iterations required

K <- 1160771                                   # Carrying capacity
q <- 0.0001484                                 # Catchability coef (given but not used)
( B <- c(732506,rep(NA,length(year)-1)) )      # An initialized biomass vector
r <- 0.4049                                    # Base intrinsic growth rate

B1 <- B
for (i in 1:iters) B1[i+1] <- B1[i]+r*B1[i]*(1-B1[i]/K)-ct[i+1]
B1

( r <- 0.4049*0.98 )
B2 <- B
for (i in 1:iters) B2[i+1] <- B2[i]+r*B2[i]*(1-B2[i]/K)-ct[i+1]

ct1 <- ct*0.984
B3 <- B
for (i in 1:iters) B3[i+1] <- B3[i]+r*B3[i]*(1-B3[i]/K)-ct1[i+1]

ct2 <- ct*0.989
B4 <- B
for (i in 1:iters) B4[i+1] <- B4[i]+r*B4[i]*(1-B4[i]/K)-ct2[i+1]

res <- round(cbind(year,ct,B1,B2,B3,B4),0)  # round for display only
print(res,na.print="")                      # don't shown NAs

d4 <- read.table("data/Box13_4.txt",header=TRUE)
str(d4)
d4$logR <- log(d4$recruits)

# for Sexpr below
sr1 <- srFuns("Ricker",param=1)
( sr1s <- srStarts(recruits~stock,data=d4,type="Ricker",param=1) )
sr1r <- nls(logR~log(sr1(stock,a,b)),data=d4,start=sr1s)

sr1 <- srFuns("Ricker",param=1)
( sr1s <- srStarts(recruits~stock,data=d4,type="Ricker",param=1) )
sr1r <- nls(logR~log(sr1(stock,a,b)),data=d4,start=sr1s)
coef(sr1r)

plot(recruits~stock,data=d4,pch=19,xlim=c(0,2000),ylim=c(0,2000))
curve(sr1(x,coef(sr1r)[1],coef(sr1r)[2]),from=0,to=2000,col="red",lwd=2,add=TRUE)
abline(coef=c(0,1),col="blue",lwd=2,lty=2)

( alpha <- coef(sr1r)[1] )
( beta <- coef(sr1r)[2] )

E <- 1000
( Zi <- log(E) - log(alpha) )

( Zi.new <- Zi*1.10 )
( alpha.new <- E*exp(-Zi.new) )

plot(recruits~stock,data=d4,pch=19,xlim=c(0,2000),ylim=c(0,2000))
curve(sr1(x,alpha,beta),from=0,to=2000,col="red",lwd=2,add=TRUE)
abline(coef=c(0,1),col="blue",lty=2,lwd=2)
curve(sr1(x,alpha.new,beta),from=0,to=2000,col="blue",lwd=2,add=TRUE)

stocksize <- 0:max(d4$stock)                 # create stock sizes to model over
recruits1 <- sr1(stocksize,alpha,beta)       # predict recruits using base model
recruits2 <- sr1(stocksize,alpha.new,beta)   # predict recruits with increased mortality
diffrecruits <- recruits2-recruits1          # find difference in predicted recruits
pdiffrecruits <- diffrecruits/recruits1*100  # express difference as a percentage of base
plot(diffrecruits~stocksize,type="l",lwd=2,xlab="Stock Size",ylab="Loss of Recruits")
plot(pdiffrecruits~stocksize,type="l",lwd=2,xlab="Stock Size",ylab="Percentage Loss of Recruits")

d5 <- read.table("data/Box13_5.txt",header=TRUE)
str(d5)
head(d5,n=10)   # first 10 rows

d5a <- reshape(d5,v.names="abund",timevar="cohort",idvar="age",direction="wide")
str(d5a)

d5a[,3]             # ask for 3rd column, or
d5a[,"abund.1953"]  # by the name of the varible
head(d5a)           # first 10 rows

( m <- dim(d5a)[1] )

( w.1953 <- lagratio(d5a[,"abund.1953"]) )
( k.1953 <- abs(log10(w.1953)) )

( K.1953 <- sum(k.1953) )

( n.coh <- ncol(d5a)-1 )             # number of columns in d5a, -1 (for age column)
( k <- matrix(rep(NA,n.coh*(m-1)),nrow=(m-1)) )             # initialize matrix
for (i in 1:n.coh) k[,i] <- abs(log10(lagratio(d5a[,i+1]))) # need +1 to ignore age col
rownames(k) <- c("k0","k1","k2","k3")                       # rename rows and col labels
colnames(k) <- 1949:1959
round(k,4)                                                  # rounded for display only

( K <- apply(k,2,sum) )
( kt <- t(k) )
yrs <- as.numeric(rownames(kt))
matplot(yrs,kt,type="l",lwd=2,ylim=c(0,2),xlab="Cohort",ylab="k(i) values")
legend("topleft",c("k0","k1","k2","k3"),lty=1:4,lwd=2,col=1:4,cex=0.8)
plot(K~yrs,type="l",lwd=2,ylim=c(0,max(K)),xlab="Cohort",ylab="K")

lm1 <- lm(kt[,1]~K)
summary(lm1)

slps <- slpp <- rsqs <- numeric(4)             # initialize vector to hold results
for (i in 1:4) {                               # look through k values
  lm1 <- lm(kt[,i]~K)                          # fit regression
  slps[i] <- coef(lm1)[2]                      # extract slope
  slpp[i] <- anova(lm1)[1,"Pr(>F)"]            # extract slope p-value
  rsqs[i] <- cor(kt[,i],K)^2                   # extract r-squared
}
round(cbind(slps,slpp,rsqs),4)                 # summary table

( d5b <- d5a[,-1]  )                             # delete age column for simplicity
( Ntp1 <- log10(as.numeric(d5b[2,])) )
( Ntp <- log10(as.numeric(d5b[1,])) )

# for Sexpr below
res1 <- lm(Ntp~Ntp1)

res1 <- lm(Ntp~Ntp1)
coef(res1)
hoCoef(res1,2,1)

res2 <- lm(Ntp1~Ntp)
coef(res2)
hoCoef(res2,2,1)

slps1 <- slpp1 <- t1 <- numeric(4)      # initial vectors for 1st regression results
slps2 <- slpp2 <- t2 <- numeric(4)      # initial vectors for 2nd regression results
for (i in 1:4) {                        # loop through ages
  tmp2 <- log10(as.numeric(d5b[i+1,]))
  tmp1 <- log10(as.numeric(d5b[i,]))
  lm1 <- lm(tmp2~tmp1)                  # first regression
  slps1[i] <- coef(lm1)[2]              # extract results
  temp <- hoCoef(lm1,2,1)
  slpp1[i] <- temp[1,"p value"]
  t1[i] <- temp[1,"T"]
  lm2 <- lm(tmp1~tmp2)                  # second regression
  slps2[i] <- coef(lm2)[2]              # extract results
  temp <- hoCoef(lm2,2,1)
  slpp2[i] <- temp[1,"p value"]
  t2[i] <- temp[1,"T"] 
}
round(cbind(slps1,t1,slpp1,slps2,t2,slpp2),4)  # summary table

tmpN <- as.numeric(d5b[1,])                    # isolate abundance data for age-0
tmpK <- as.numeric(k[1,])                      # isolate k data for age-0
lm1 <- lm(tmpK~tmpN)                           # regression fitting
summary(lm1)                                   # extract results

fitPlot(lm1,xlab="Initial Age-0 Density",ylab="k-value from Age-0 to Age-1",main="")

slpsE <- slppE <- tE <- numeric(4)             # initialize results vectors
for (i in 1:4) {                               # loop through ages
  tmpN <- as.numeric(d5b[i,])
  tmpK <- as.numeric(k[i,])
  lm1 <- lm(tmpK~tmpN)                         # fit regression
  slpsE[i] <- coef(lm1)[2]                     # extract results
  temp <- hoCoef(lm1,2)
  slppE[i] <- temp[1,"p value"]
  tE[i] <- temp[1,"T"]
}
round(cbind(slpsE,tE,slppE),5)                 # summary table


# Script created at 2015-05-15 17:34:26
