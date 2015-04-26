library(FSA)          # fitPlot, Subset
library(NCStats)      # addSigLetters
library(car)          # Anova, durbinWatsonTest, vif
library(Hmisc)        # rcorr
library(Kendall)      # kendall
library(lsmeans)      # lsmeans
library(multcomp)     # glht, mcp, cld 
library(nlstools)     # overview, nlsBoot
library(plotrix)      # thigmophobe.labels, rescale

setwd("c:/aaaWork/web/fishR/BookVignettes/AIFFD/")

options(show.signif.stars=FALSE,contrasts=c("contr.sum","contr.poly"))

set.seed(983452)

d1 <- read.table("data/box4_1.txt",header=TRUE)
str(d1)
d1$lcatch <- log(d1$catch)

d1$age <- factor(d1$age)
d1$yearcl <- factor(d1$yearcl)

xtabs(~age+yearcl,data=d1)

lm1 <- lm(lcatch~age+yearcl,data=d1)
Anova(lm1,type="III")

lsmeans(lm1,~age)
lsmeans(lm1,~yearcl)

mc1a <- glht(lm1,mcp(age="Tukey"))
summary(mc1a)
mc1yc <- glht(lm1,mcp(yearcl="Tukey"))
summary(mc1yc)
cld(mc1yc)

fitPlot(lm1,which="yearcl",xlab="Year Class",ylab="Loge(Catch)",main="")
addSigLetters(lm1,which="yearcl",lets=c("c","a","ab","bc","bc","ab","c","ab"), 
              pos=c(4,2,2,2,4,2,2,2),cex=0.75)

d2 <- read.table("data/box4_2.txt",header=TRUE)
str(d2)
d2a <- Subset(d2,year>=1972)
str(d2a)

plot(age0cpe~year,data=d2a,ylab="Age-0 CPE",xlab="Year",main="",pch=19)
plot(age0cpe~year,data=d2a,type="b",ylab="Age-0 CPE",xlab="Year",main="",pch=19)

Summarize(d2a$age0cpe,digits=4)
Summarize(d2a$year,digits=4)

rcorr(as.matrix(d2a[,c("age0cpe","year")]))

Kendall(d2a$age0cpe,d2a$year)

lm1 <- lm(age0cpe~year,data=d2a)
anova(lm1)
summary(lm1)

durbinWatsonTest(lm1,max.lag=5)

plot(lm1$residuals~d2a$year,ylab="Residuals",xlab="Year",main="",pch=19)
abline(h=0,lty=2)                          # adds horizontal reference line at 0
plot(lm1$residuals~d2a$year,type="b",ylab="Residuals",xlab="Year",main="",pch=19)
abline(h=0,lty=2)

fitPlot(lm1,ylab="Age-0 CPE",xlab="Year",main="")

d2a$logage0cpe <- log(d2a$age0cpe)
lm2 <- lm(logage0cpe~year,data=d2a)
plot(lm2$residuals~d2a$year,ylab="Residuals",xlab="Year",main="",pch=19)
abline(h=0,lty=2)                          # adds horizontal reference line at 0
plot(lm2$residuals~d2a$year,type="b",ylab="Residuals",xlab="Year",main="",pch=19)
abline(h=0,lty=2)

anova(lm2)
summary(lm2)
durbinWatsonTest(lm2,max.lag=5)
ccf(lm2$residuals,d2a$year)

plot(age0cpe~year,data=d2a,ylab="Age-0 CPE",xlab="Year",pch=19)
pCPE <- exp(predict(lm2))
lines(pCPE~d2a$year,lwd=2,col="red")

d3 <- read.table("data/box4_3.txt",header=TRUE)
str(d3)
d3$rep <- factor(d3$rep)
d3$time <- factor(d3$time)

rmsp2 <- function(object,type=c("III","II","I")) {
  type <- match.arg(type)
  # extract df and SS of appropriate rows of ANOVA table
  if (type=="I") { res <- anova(object)[1:2,1:3] }
    else if (type=="III") { res <- Anova(object,type=type)[2:4,2:1] }
      else { res <- Anova(object,type=type)[1:3,2:1] }
  # compute MSs
  res[,"Mean Sq"] <- res[,2]/res[,1]
  # MS in third position is the error MS
  errorMS <- res[3,"Mean Sq"]
  # compute F for first two positions (put NA in last position)
  res[,"F"] <- c(res[1:2,"Mean Sq"]/errorMS,NA)
  # convert Fs to p-values
  res[,"PR(>F)"] <- c(pf(res[1:2,"F"],res[1:2,"Df"],res[3,"Df"],lower.tail=FALSE),NA)  
  res
}

lm1 <- lm(terms(count~area+rep+area*rep+time+time*area,keep.order=TRUE),data=d3)
Anova(lm1,type="II")

rmsp2(lm1)

lm1a <- lm(terms(count~area+rep+area*rep+time,keep.order=TRUE),data=d3)
mc1 <- glht(lm1a,mcp(time="Tukey"))
summary(mc1)
cld(mc1)

d4 <- read.table("data/box4_4.txt",header=TRUE)
str(d4)
d4$lnum <- log(d4$num+1)
d4$lmeanret <- log10(d4$meanret)

lm1 <- lm(lnum~age,data=d4)
d4$plnum1 <- predict(lm1)

lm2 <- lm(lnum~age,weights=plnum1,data=d4)
anova(lm2)
summary(lm2)

( preds <- predict(lm2,se.fit=TRUE) )

ycs <- data.frame(age=d4$age, yearcl=d4$yearcl, lnum=d4$lnum, meanret=d4$meanret,
                  lmeanret=d4$lmeanret, plnum=preds$fit, pse=preds$se.fit,
                  resid=lm2$residuals, sresid=rstandard(lm2), cooksD=cooks.distance(lm2))
round(ycs,3)    # rounded for display purposes only

plot(lnum~age,data=d4,xlab="Year Class (19__)",xaxt="n", ylab="loge Number-at-Age",pch=19)
axis(1,at=ycs$age,labels=ycs$yearcl)  # Add 'new' x-axis
abline(lm2,lwd=2)                     # Add regression line
for (i in 1:nrow(ycs)) {              # Add residual lines
  with(ycs,lines(c(age[i],age[i]),c(lnum[i],plnum[i]),lty=2,lwd=2))
}

plot(sresid~yearcl,data=ycs,type="b",xlab="Year Class (19__)", ylab="Studentized Residual",
     pch=19,ylim=c(-2,2))
abline(h=0,lty=2)                     # Add horizontal line at 0
plot(sresid~yearcl,data=ycs,type="h",xlab="Year Class (19__)", ylab="Studentized Residual",
     pch=19,ylim=c(-2,2),lwd=3)
abline(h=0,lty=2)

d5 <- read.table("data/box4_5.txt",header=TRUE)
str(d5)
d5

lapply(as.list(d5[,-c(1,3)]),Summarize,digits=4)

pairs(~cpe0+winstage+winret+sprstage,data=d5,pch=19)

rcorr(as.matrix(d5[,c("cpe0","winstage","winret","sprstage")]))

lm1 <- lm(cpe0~winstage+winret+sprstage,data=d5)
summary(lm1)
vif(lm1)

lm2 <- lm(cpe0~winstage+sprstage,data=d5)
summary(lm2)

lm3 <- lm(cpe0~winstage,data=d5)
summary(lm3)

plot(cpe0~winstage,data=d5,pch=19,xlim=c(170.5,172),ylim=c(0,12),
     xlab="Mean Winter (Jan-Mar) Stage",ylab="Age-0 CPE")
thigmophobe.labels(d5$winstage,d5$cpe0,d5$yearcl-1900,cex=0.8)
abline(v=171.95,lty=3)
abline(lm3,lty=2)

plot(resid~meanret,data=ycs,pch=19,xlab="Mean Retention Time", ylab="Residual")
abline(h=0,lty=2)
plot(resid~lmeanret,data=ycs,pch=19,xlab="log10 Mean Retention Time", ylab="Residual")
abline(h=0,lty=2)

lm3 <- lm(lnum~age+lmeanret,weights=plnum1,data=d4)
summary(lm3)

Summarize(ycs$resid,digits=4)
Summarize(ycs$meanret,digits=4)
rcorr(cbind(ycs$resid,d4$meanret))

d7 <- read.table("data/box4_7.txt",header=TRUE)
str(d7)
d7$lrecruit <- log(d7$recruit)

plot(recruit~spawner,data=d7,ylab="Catch Rate of Age-0 Crappie",
     xlab="Biomass of Age-2+ Crappie",pch=19)

bhst <- list(a=0.03,b=0.002)                     # list of starting values
bhsr <- recruit~(a*spawner)/(1+b*spawner)        # B-H model as an R formula
nls1 <- nls(bhsr,data=d7,start=bhst,trace=TRUE)
overview(nls1)

bhst2 <- list(a=0.01,b=0.002)
bhsr2 <- lrecruit~log((a*spawner)/(1+b*spawner))
nls2 <- nls(bhsr2,data=d7,start=bhst2)
overview(nls2)

bhbc <- nlsBoot(nls2,niter=2000)
summary(bhbc)

d7$pred.lrec <- fitted(nls2)
d7$res.lrec <- residuals(nls2)
d7
sum(d7$res.lrec)

lm1 <- lm(pred.lrec~lrecruit,data=d7)
anova(lm1)
summary(lm1)

plot(recruit~spawner,data=d7,pch=19)
cnls1 <- coef(nls1)
curve((cnls1[1]*x)/(1+cnls1[2]*x),from=min(d7$spawner),to=max(d7$spawner),col="red",
      lwd=2,lty=2,add=TRUE)
cnls2 <- coef(nls2)
curve((cnls2[1]*x)/(1+cnls2[2]*x),from=min(d7$spawner),to=max(d7$spawner),col="blue",
      lwd=2,lty=2,add=TRUE)
legend("topleft",legend=c("Additive Error","Multiplicative Error"),col=c("red","blue"),
       lwd=2,lty=2,cex=0.5)

d8 <- read.table("data/box4_8.txt",header=TRUE)
str(d8)
d8$logR <- log(d8$recruit)
d8$logS <- log(d8$spawner)
d8$ratio <- d8$recruit/d8$spawner
d8$lratio <- log(d8$ratio)
view(d8)

par(mar=c(3.3,3.5,0.4,0.2),mgp=c(2.35,0.3,0))
plot(recruit~spawner,data=d8,ylab="Number of Age-0 Walleye",
     xlab="Number of Age-5 Walleye",pch=19)
with(d8,thigmophobe.labels(spawner,recruit,labels=year-1900,cex=0.8))
plot(recruit~spawner,data=d8,ylab="Number of Age-0 Walleye",
     xlab="Number of Age-5 Walleye",pch=19,cex=rescale(mtempcv,c(0.5,2)))

rst <- list(a=4,b=0)                           # Starting values
rsr <- logR~log(spawner*exp(a+b*spawner))      # Ricker model as an R formula
nls1 <- nls(rsr,data=d8,start=rst,trace=TRUE)
overview(nls1)

d8$pred1.lrec <- fitted(nls1)
d8$res1.lrec <- residuals(nls1)
sum(d8$res1.lrec)

lm1 <- lm(pred1.lrec~logR,data=d8)
summary(lm1)

lm2 <- lm(lratio~spawner,data=d8)
summary(lm2)
d8$pred2.lrec <- log(exp(fitted(lm2))*d8$spawner)
lm2a <- lm(pred2.lrec~logR,data=d8)
summary(lm2a)

rst2 <- list(a=4,b=0,c=-7)
rsr2 <- logR~log(spawner*exp(a+b*spawner+c*mtempcv))
nls2 <- nls(rsr2,data=d8,start=rst2,trace=TRUE)
overview(nls2)
d8$pred3.lrec <- fitted(nls2)
d8$res3.lrec <- residuals(nls2)
sum(d8$res2.lrec)
lm4 <- lm(pred3.lrec~logR,data=d8)
summary(lm4)

anova(nls1,nls2)

AIC(nls1,nls2)

par(mar=c(3.3,3.5,0.4,0.2),mgp=c(2.35,0.3,0))
rsr3 <- recruit~spawner*exp(a+b*spawner)
nls3 <- nls(rsr3,data=d8,start=rst)
plot(recruit~spawner,data=d8,pch=19,xlim=c(0,max(d8$spawner)))
cnls1 <- coef(nls1)
curve(x*exp(cnls1[1]+cnls1[2]*x),from=0,to=max(d8$spawner),col="red",lwd=2,lty=2,add=TRUE)
cnls3 <- coef(nls3)
curve(x*exp(cnls3[1]+cnls3[2]*x),from=0,to=max(d8$spawner),col="blue",lwd=2,lty=2,add=TRUE)
legend("topright",legend=c("Multiplicative Error","Additive Error"),col=c("red","blue"),
       lwd=2,lty=2,cex=0.5)

d8$logR <- log(d8$recruit)
rst <- list(a=4,b=0)
rsr <- logR~log(spawner*exp(a-b*spawner))
nls1 <- nls(rsr,data=d8,start=rst)
overview(nls1)

rbc <- nlsBoot(nls1,niter=2000)
summary(rbc)
confint(rbc)                  # default 95% CI
confint(rbc,conf.level=0.9)   # illustrative 90% CI

str(rbc)
view(rbc$coefboot)

rbc.d <- as.data.frame(rbc$coefboot)
Summarize(rbc.d$a)
Summarize(rbc.d$b)
t.test(rbc.d$a)
t.test(rbc.d$b)

quantile(rbc.d$a,probs=c(0,1,5,10,25,50,75,90,95,99,100)/100)
quantile(rbc.d$b,probs=c(0,1,5,10,25,50,75,90,95,99,100)/100)

hist(rbc.d$a,xlab="Bootstrap Estimates of a",main="")
abline(v=coef(nls1)["a"],lty=2,lwd=2,col="red")       # put nls estimate on hist
abline(v=confint(rbc)["a",],lty=3,col="blue")         #   and bootstrap CIs
hist(rbc.d$b,xlab="Bootstrap Estimates of b",main="")
abline(v=coef(nls1)["b"],lty=2,lwd=2,col="red")       # put nls estimate on hist
abline(v=confint(rbc)["b",],lty=3,col="blue")         #   and bootstrap CIs

par(mar=c(3.3,3.5,0.4,0.2),mgp=c(2.35,0.3,0))
plot(rbc)


# Script created at 2015-04-26 15:56:58
