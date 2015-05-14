library(FSA)       # Summarize, Subset, binCI, rcumsum, catchCurve, chapmanRobson, mrOpen
library(car)       # Anova

setwd("C:/aaaWork/Web/fishR/BookVignettes/aiffd2007")

options(width=90,continue=" ",show.signif.stars=FALSE,
        contrasts=c("contr.sum","contr.poly"))
par(mar=c(3.5,3.5,1,1),mgp=c(2.1,0.4,0),tcl=-0.2)

N0 <- 1000                       # Initial population
N12 <- 700                       # Final population
( d <- N0-N12 )                  # Interval absolute mortality (i.e., number of deaths)
A12 <- d/N0                      # Interval (12 month) mortality rate
100*A12
( Z12 <- -log(1-A12) )           # Instantaneous mortality rate

Z1 <- Z12/12                     # Monthly instantaneous mortality rate
( Z4 <- 4*Z1 )                   # Four-month Z estimate
( Z8 <- 8*Z1 )                   # Eight-month Z estimate

( A4 <- 1-exp(-Z4) )             # Interval (4-month) mortality rate
( A8 <- 1-exp(-Z8) )             # Interval (8-month) mortality rate

d4 <- N0*A4                      # Deaths in first four-months
round(d4,0)
N4 <- N0-d4                      # Estimated survivors after 4-months
round(N4,0)
d8 <- N0*A8                      # Deaths in first eight-months
round(d8,0)
N8 <- N0-d8                      # Estimated survivors after 4-months
round(N8,0)
d8.4 <- d8-d4                    # Deaths in second four-months
round(d8.4,0)
d12.8 <- d-d8                    # Deaths in third four-months
round(d12.8,0)

d4+d8.4+d12.8                    # Total deaths in the three four-month periods
(d4+d8.4+d12.8) == d             # Does this equal the first deaths calculation?

age <- 1:7
catch <- c(257,407,147,32,17,5,4)
rbind(age,catch)                        # for display purposes only.

catch[2]                         # 2nd position
catch[c(2,4)]                    # 2nd & 4th positions
catch[2:4]                       # 2nd through 4th positions
catch[-2]                        # all but the 2nd position

( N <- sum(catch[2:7]) )         # sum in positions 2 thru 7
( A <- 100*catch[2]/N )          # annual mortality estimate
( se.A <- sqrt(A*(100-A)/N) )    # SE of A

100*binCI(catch[2],N)            # default 95% CI for A

( N <- rcumsum(catch) )
100*catch/N
100*binCI(catch,N)
data.frame(age,catch,N,A=100*catch/N,100*binCI(catch,N))  # combine for display only

( d3 <- read.table("data/Box6_3.txt",header=TRUE) )

cr1 <- chapmanRobson(d3$age,d3$catch,3:10)
summary(cr1)
confint(cr1)
plot(cr1)

d3$logcatch <- log(d3$catch)
str(d3)

cc1 <- lm(logcatch~age,data=d3,subset=age>=3)
anova(cc1)
summary(cc1)
confint(cc1)

( Z <- -coef(cc1)[2] )        # 2nd coefficient of model results
( ZCI <- -confint(cc1)[2,] )  # 2nd row of confint results

( A <- 100*(1-exp(-Z)) )
( ACI <- 100*(1-exp(-ZCI)) )

d3$W <- predict(cc1,d3)
d3

cc2 <- lm(logcatch~age,data=d3,subset=age>=3,weights=W)
anova(cc1)
summary(cc2)
( Z <- -coef(cc2)[2] )
( ZCI <- -confint(cc2)[2,] )

cc3 <- catchCurve(catch~age,d3,3:10)                    # unweighted catch-curve
anova(cc3)
summary(cc3)
confint(cc3)
plot(cc3)

cc4 <- catchCurve(catch~age,d3,3:10,use.weights=TRUE)   # weighted catch-curve
summary(cc4)
confint(cc4)

d5 <- read.table("data/Box6_5.txt",header=TRUE)
str(d5)
d5$age <- 2001-d5$year.class         # find age
d5$ycs <- d5$age0/min(d5$age0)       # create year-class strength
d5$adjcatch <- d5$catch/d5$ycs       # adjust catch for year-class strength
d5                                   # now, what does the data frame look like

cc1 <- catchCurve(adjcatch~age,data=d5)
summary(cc1)
confint(cc1)
plot(cc1)

d5$yccf <- max(d5$age0)/d5$age0
d5$adjcatch2 <- d5$catch*d5$yccf
d5                                   # now, what does the data frame look like

cc2 <- catchCurve(adjcatch2~age,data=d5)
summary(cc2)

d5$yccf2 <- mean(d5$age0)/d5$age0
d5$adjcatch3 <- d5$catch*d5$yccf2
cc3 <- catchCurve(adjcatch3~age,data=d5)
summary(cc3)

d6 <- read.table("data/Box6_6.txt",header=TRUE)
str(d6)
d6$logcatch <- log(d6$catch)

lm1 <- lm(logcatch~age*lake,data=d6,subset=age>=2)
anova(lm1)
Anova(lm1,type="III")
Anova(lm1,type="II")
fitPlot(lm1,main="")

# for Sexpr below
lm2 <- lm(logcatch~age,data=d6,subset=age>=2)

lm2 <- lm(logcatch~age,data=d6,subset=age>=2)
anova(lm2)
coef(lm2)
confint(lm2)

( d7 <- data.frame(age=1:6,catch=c(65,66,27,6,4,1)) )

cr1 <- chapmanRobson(catch~age,d7,2:6)
summary(cr1)
confint(cr1)
plot(cr1)

cc1 <- catchCurve(catch~age,d7,2:5)
summary(cc1)
confint(cc1)
plot(cc1)

bins <- 9:52
freqs <- c(5,39,43,55,64,86,106,99,82,81,56,45,27,36,43,52,65,73,63,42,44,40,37,31,36,
           15,19,18,13,8,20,19,9,12,9,11,13,13,4,10,6,3,3,4)
length(bins)
length(freqs)

len.cm <- rep(bins,times=freqs)
table(len.cm)
len.mm <- len.cm*10

set.seed(2101)

rands <- sample(0:9,length(len.mm),replace=TRUE)
len.mm <- len.mm + rands

hbins <- seq(90,530,by=10)
h <- hist(len.mm,breaks=hbins,right=FALSE,col="gray90",main="")
h$mids
h$counts

write.table(data.frame(len=len.mm),"Box6_8.txt",quote=FALSE,row.names=FALSE,col.names=TRUE)

d8 <- read.table("data/Box6_8.txt",header=TRUE)
str(d8)
d8a <- Subset(d8,len>=150)
str(d8a)

Linf <- 636
K <- 0.226
( Lmean <- mean(d8a$len) )
( Lmin <- min(d8a$len) )
( Z1 <- K*(Linf-Lmean)/(Lmean-Lmin) )

( Lmedian <- median(d8a$len) )
Ymedian <- -log(1-Lmedian/Linf)
Ymin <- -log(1-Lmin/Linf)
( Z2 <- 0.693*K / (Ymedian-Ymin) )

bins <- seq(150,530,by=10)
h <- hist(d8a$len,breaks=bins,right=FALSE,xlab="Total Length (mm)",main="",col="gray90")
h$mids
h$counts

tprime <- -log(1-h$mids/Linf)
round(tprime,3)                # rounded only to compare to values in Box 6.8
logN <- log(h$counts)

lm1 <- lm(logN~tprime)
coef(lm1)
( b <- coef(lm1)[2] )
( Z3 <- K*(1-b) )

ord.len <- order(d8a$len)
ord.top3 <- ord.len[(length(d8a$len)-2):length(d8a$len)]
( len.top3 <- d8a$len[ord.top3] )
( Lmean.top3 <- mean(len.top3) )

( Linf1 <- Lmean.top3/0.95 )

( Lmax <- max(d8a$len) )
( Linf2 <- exp(0.044+0.984*log(Lmax)) )

( Lx <- h$breaks )                             # lower limit of each TL interval
Lmeans <- NULL                                 # initiate the Lmeans vector
for (i in 1:length(Lx)) {                      # loop through length intervals
  Lmeans <- c(Lmeans,mean(d8a$len[d8a$len>=Lx[i]]))
}
Lmeans                                         # mean for each interval & higher
Ldiff <- Lmeans-Lx

lm2 <- lm(Ldiff~Lx)
a <- coef(lm2)[1]      # isolate intercept
b <- coef(lm2)[2]      # isolate slope
( Linf3 <- -a/b )

s1 <- rep(NA,5)
s2 <- c(43,rep(NA,4))
s3 <- c(28,31,NA,NA,NA)
s4 <- c(12,16,48,NA,NA)
s5 <- c(3,9,19,37,NA)
( aiffd.top <- rbind(s1,s2,s3,s4,s5) )
( mb.top <- t(aiffd.top) )

n <- c(643,489,712,630,68)    # m_i_ from book
m <- c(0,43,59,76,68)         # r_.i_ from book
u <- n-m                      
R <- n                        # assumes no accidental deaths, thus m_i_ from book
( mb.bot <- rbind(m,u,n,R) )

mr1 <- mrOpen(mb.top, mb.bot)
summary(mr1,verbose=TRUE)

summary(mr1,parm="phi")
confint(mr1,parm="phi")


# Script created at 2015-05-13 19:49:45
