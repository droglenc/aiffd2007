library(FSA)          # Summarize, fitPlot, residPlot
library(car)          # Anova, recode, leveneTest
library(Hmisc)        # rcorr
library(multcomp)     # glht, mcp

setwd("C:/aaaWork/Web/fishR/BookVignettes/aiffd2007")

options(width=90,continue=" ",show.signif.stars=FALSE,
        contrasts=c("contr.sum","contr.poly"))
par(mar=c(3.5,3.5,1,1),mgp=c(2.1,0.4,0),tcl=-0.2)

year <- 0:14                   # years for simulation
N0 <- 10000                    # initial population size
A <- 0.95                      # annual survival rate
abund <- round(N0*A^year,0)    # abundance from exponential decay model
cpe <- round(abund/1000,1)     # simulated cpe data
cbind(year,abund)              # for display only

me1 <- c(-0.10,-0.05,0.05,0.10)
me1a <- sample(me1,length(cpe),replace=TRUE)
cpe.me1 <- round(cpe+cpe*me1a,1)

me2 <- c(-0.4,-0.2,0.1,0.4)
me2a <- sample(me2,length(cpe),replace=TRUE)
cpe.me2 <- round(cpe+cpe*me2a,1)
( d1 <- data.frame(year,abund,cpe,cpe.me1,cpe.me2) )     # for storage

plot(cpe~year,data=d1,ylim=c(0,15),xlab="Year",ylab="CPE",type="l",lwd=2,col="black")
points(d1$year,d1$cpe.me1,type="b",pch=19,lwd=1,col="blue")
points(d1$year,d1$cpe.me2,type="b",pch=19,lwd=1,col="red")
legend("topright",legend=c("Truth","5 or 10%","20 or 40%"),col=c("black","blue","red"),
       lwd=c(2,1,1),pch=c(NA,19,19))

cpe.sim <- function(year,cpe,me1,me2=NULL) {
  cpe1 <- cpe + cpe*sample(me1,length(cpe),replace=TRUE)
  if (!is.null(me2)) cpe2 <- cpe + cpe*sample(me2,length(cpe),replace=TRUE) 
  max.y <- max(c(cpe,cpe1,cpe2))
  plot(year,cpe,ylim=c(0,max.y),xlab="Year",ylab="CPE",type="l",lwd=2,col="black")
  points(year,cpe1,type="b",pch=19,lwd=1,col="blue")
  if (!is.null(me2)) points(year,cpe2,type="b",pch=19,lwd=1,col="red")
  if (is.null(me2)) legend("topright",legend=c("Truth","CPE1"),col=c("black","blue"),lwd=c(2,1),pch=19)
    else legend("topright",legend=c("Truth","CPE1","CPE2"),col=c("black","blue","red"),lwd=c(2,1,1),pch=19)
}

cpe.sim(d1$year,d1$cpe,me1,me2)

d2 <- read.table("data/Box7_2.txt",header=TRUE)
str(d2)
d2$logcpe.A <- log10(d2$cpe.A+1)
d2$logcpe.B <- log10(d2$cpe.B+1)
str(d2)

lapply(as.list(d2[,-1]),Summarize,digits=3)
hist(d2$cpe.A,main="")
hist(d2$cpe.B,main="")
hist(d2$logcpe.A,main="")
hist(d2$logcpe.B,main="")

( mn.B <- mean(d2$logcpe.B) )
( sd.B <- sd(d2$logcpe.B) )

power.t.test(n=NULL,delta=mn.B*0.10,power=0.9,sig.level=0.05,sd=sd.B,type="two.sample",
             alt="one.sided")
power.t.test(n=NULL,delta=mn.B*0.20,power=0.8,sig.level=0.10,sd=sd.B,type="two.sample",
             alt="one.sided")

perc.delta <- seq(0.05,0.25,0.01)
deltas <- mn.B*perc.delta
p1 <- sapply(deltas,power.t.test,n=NULL,power=0.9,sig.level=0.05,sd=sd.B,type="two.sample",
             alt="one.sided")
p1.n <- as.numeric(p1["n",])
plot(p1.n~perc.delta,xlab="Proportion Population Decline",ylab="Sample Size",type="l",
     ylim=c(0,max(p1.n)))

d3 <- read.table("data/Box7_3.txt",header=TRUE)
str(d3)

lm1 <- lm(cpue~lake+veg,data=d3)
anova(lm1)
fitPlot(lm1,xlab="Lake",main="",interval=FALSE)

lm2 <- lm(cpue~veg,data=d3)
anova(lm2)
fitPlot(lm2,xlab="Vegetation Treatment",main="")

d5 <- read.table("data/Box7_5.txt",header=TRUE)
str(d5)
d5$fDepth <- factor(d5$Depth)
d5$fMonth <- car::recode(d5$Month,"1='June'; 2='August'",as.factor.result=TRUE)
d5$fMonth <- factor(d5$fMonth,levels=c("June","August"))
str(d5)

lm1 <- lm(CPE~fMonth*fDepth,data=d5)
anova(lm1)

Summarize(CPE~fMonth*fDepth,data=d5,digits=2)
fitPlot(lm1,change.order=TRUE,xlab="Depth (m)",main="")

mc1 <- glht(lm1,mcp(fDepth="Tukey"))
summary(mc1)

residPlot(lm1)
leveneTest(lm1)

hist(lm1$residuals,main="")

d6 <- read.table("data/Box7_6.txt",header=TRUE)
str(d6)
d6$logcpe <- log10(d6$cpe+1)

cor(d6[,c("gravel","veg","slope")])
rcorr(as.matrix(d6[,c("gravel","veg","slope")]))

pairs(~gravel+veg+slope,data=d6,pch=19)

lm1 <- lm(cpe~gravel,data=d6)
lm2 <- lm(cpe~veg,data=d6)
lm3 <- lm(cpe~slope,data=d6)
lm1a <- lm(logcpe~gravel,data=d6)
lm2a <- lm(logcpe~veg,data=d6)
lm3a <- lm(logcpe~slope,data=d6)
anova(lm1a)
summary(lm1a)
confint(lm1a)
fitPlot(lm1a,xlab="Gravel (%)",ylab="log10(CPE+1)",main="")

# for Sexpr below
lm4 <- lm(logcpe~gravel*veg,data=d6)

lm4 <- lm(logcpe~gravel*veg,data=d6)
anova(lm4)


# Script created at 2015-05-13 19:49:59
