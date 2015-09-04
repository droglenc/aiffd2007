library(FSA)       # Summarize, fitPlot, binCI, view
library(NCStats)   # sdCalc
library(arm)       # se.ranef
library(car)       # Anova, leveneTest, outlierTest
library(languageR) # pvals.fnc
library(lattice)   # bwplot, densityplot
library(lme4)      # lmer and related extractors
library(lsmeans)   # lsmeans
library(nortest)   # ad.test, cvm.test
library(sciplot)   # se
library(survey)    # svydesign, svymean

setwd("c:/aaaWork/web/fishR/BookVignettes/AIFFD/")

options(show.signif.stars=FALSE,contrasts=c("contr.sum","contr.poly"))

d1 <- read.table("data/box3_1.txt",header=TRUE)
str(d1)

Summarize(d1$catch)

summary(d1$catch)
sd(d1$catch)

se(d1$catch)

df <- length(d1$catch)-1             # compute df as n-1
C <- 0.95                            # set confidence level
a2 <- (1-C)/2                        # find upper tail probability
( ct <- qt(a2,df,lower.tail=FALSE) ) # find critical t-value

( mn <- mean(d1$catch) )              # find (and save) the mean
( semn <- se(d1$catch) )              # find (and save) the SE of the mean
( LCI <- mn-ct*semn )
( UCI <- mn+ct*semn )

t.test(d1$catch)

sdCalc(d1$catch)

( d2 <- read.table("data/box3_2.txt",header=TRUE) )

ttl <- sum(d2$n)                                      # compute (and save) total catch
( p <- d2$n/ttl )                                     # find proportion of "successes"

( seps <- sqrt(p*(1-p)/(ttl-1)) )

alpha <- 0.05
n <- ttl
a <- 55                                              # demonstrate for first group only
( f.lwr <- qf(alpha,2*(n-a+1),2*a,lower.tail=FALSE) )
( f.upr <- qf(alpha,2*a+2,2*(n-a+1)-2,lower.tail=FALSE) )

( LCI <- a / (a+(n-a+1)*f.lwr) )
( UCI <- ((a+1)*f.upr)/(n-a+(a+1)*f.upr) )

LCIs <- d2$n/(d2$n+(ttl-d2$n+1)*qf(alpha,2*(ttl-d2$n+1),2*d2$n,lower.tail=FALSE))
UCIs <- ((d2$n+1)*qf(alpha,2*d2$n+2,2*(ttl-d2$n+1)-2,lower.tail=FALSE))/
  (ttl-d2$n+(d2$n+1)*qf(alpha,2*d2$n+2,2*(ttl-d2$n+1)-2,lower.tail=FALSE))
data.frame(age=d2$age,p,seps,LCIs,UCIs)               # put in data frame just to show

binCI(d2$n,ttl)

d3 <- read.table("data/box3_3.txt",header=TRUE)
str(d3)
d3$ttlW <- d3$Zooplankton + d3$Benthos + d3$Fish
view(d3)

srs <- svydesign(id=~1,data=d3)
( res1 <- svyratio(~Zooplankton,~ttlW,srs) )
confint(res1)

( res2 <- svyratio(~Zooplankton+Benthos+Fish,~ttlW,srs) )
confint(res2)

d4 <- read.table("data/box3_4.txt",header=TRUE)
str(d4)
view(d4)

dstrat <- svydesign(id=~1,strata=~stratum,fpc=~cells,data=d4)
summary(dstrat)

( d.mns <- svymean(~catch,dstrat) )
confint(d.mns)

cv(d.mns)
svymean(~catch,dstrat,deff=TRUE)

d5 <- read.table("data/box3_5.txt",header=TRUE)
str(d5)
d5$net <- factor(d5$net)
view(d5)

( catch <- tapply(d5$catch1,d5$net,FUN=mean) )
( ttl.wt <- tapply(d5$weight,d5$net,FUN=sum) )

dclust <- svydesign(id=~net,data=d5)
summary(dclust)
( mn.wt <- svymean(~weight,dclust,na.rm=TRUE) )
confint(mn.wt)

dclust2 <- svydesign(id=~net+fish,data=d5)

d6 <- read.table("data/box3_6.txt",header=TRUE)
str(d6)

dsyst <- svydesign(id=~start,data=d6)
summary(dsyst)
( mn.wid <- svymean(~width,dsyst) )
confint(mn.wt)

d7 <- read.table("data/box3_7.txt",header=TRUE)
str(d7)
d7a <- d7[!is.na(d7$measured),]
str(d7a)

mn.est1 <- mean(d7a$estimated)
d7a$c.estimated <- d7a$estimated - mn.est1
lm1 <- lm(measured~c.estimated,data=d7a)

mn.est1 <- mean(d7a$estimated)
d7a$c.estimated <- d7a$estimated - mn.est1
lm1 <- lm(measured~c.estimated,data=d7a)
coef(lm1)

predict(lm1,data.frame(c.estimated=mean(d7$estimated)-mn.est1),interval="c",se.fit=TRUE)

d8 <- read.table("data/box3_8.txt",header=TRUE)
str(d8)
d8$fbag_limit <- factor(d8$bag_limit)

lm1 <- lm(catch~fbag_limit,data=d8)
anova(lm1)
Anova(lm1,type="III")

lsmeans(lm1,~fbag_limit)

fitPlot(lm1,xlab="Walleye Bag Limit",ylab="Walleye Catch Rate",main="")

Anova(lm1,type="II")

lm1$residuals               # for demonstration purposes only
ad.test(lm1$residuals)
cvm.test(lm1$residuals)
shapiro.test(lm1$residuals)

ks.test(lm1$residuals,"pnorm")

qqnorm(lm1$residuals,main="")
hist(lm1$residuals,main="")

bartlett.test(catch~fbag_limit,data=d8)
leveneTest(catch~fbag_limit,data=d8)
leveneTest(lm1)

residPlot(lm1)

studres(lm1)

outlierTest(lm1)

lm1 <- lm(catch~region*fbag_limit,data=d8)

Anova(lm1,type="III")
lsmeans(lm1,~region)
lsmeans(lm1,~fbag_limit)

fitPlot(lm1,change.order=TRUE,xlab="Walleye Bag Limit",ylab="Walleye Catch Rate",main="")

fitPlot(lm1,which="fbag_limit",xlab="Walleye Bag Limit",ylab="Walleye Catch Rate",main="")
fitPlot(lm1,which="region",xlab="Region",ylab="Walleye Catch Rate",main="")

lm2 <- lm(catch~region+fbag_limit,data=d8)
Anova(lm2,type="III")
lsmeans(lm2,~fbag_limit)
lsmeans(lm2,~region)

Anova(lm1,type="II")
Anova(lm2,type="II")

d11 <- read.table("data/box3_11.txt",header=TRUE)
str(d11)

lm1 <- lm(growth~egg_diameter*substrate,data=d11)
Anova(lm1,type="III")

fitPlot(lm1,xlab="Egg Diameter (mm)",ylab="Final length (mm)",legend="topleft",main="")

Anova(lm1,type="II")

d12 <- read.table("data/box3_12.txt",header=TRUE)
str(d12)
d12$year <- factor(d12$year)
d12$lakeID <- factor(d12$lakeID)
str(d12)

round(tapply(d12$bluegill,list(d12$herbicide,d12$year), FUN=mean),1) # round for display only
round(tapply(d12$bluegill,list(d12$herbicide,d12$year), FUN=sd),1)

bwplot(bluegill~herbicide|year,data=d12)

lme1 <- lmer(bluegill~1+herbicide+(1|year),data=d12,REML=TRUE)

summary(lme1)

ranef(lme1)
se.ranef(lme1)
anova(lme1)

lme1a <- lmer(bluegill~0+herbicide+(1|year),data=d12,REML=TRUE)
smry1a <- summary(lme1a)
coef(smry1a)

# needed for Sexpr below
lme2 <- lmer(bluegill~1+(1|year),data=d12,REML=TRUE)  # fit without fixed term

lme2 <- lmer(bluegill~1+(1|year),data=d12,REML=TRUE)  # fit without fixed term
lrt(lme2,com=lme1)                                    # perform LRT comparison

( fe1a <- fixef(lme1a) )
( re1a <- ranef(lme1a)$year )
pred.C <- fe1a[1]+re1a
pred.H <- fe1a[2]+re1a
pred <- t(cbind(pred.C,pred.H))
colnames(pred) <- rownames(re1a)
rownames(pred) <- c("Control","Herbicide")
round(pred,1)          # rounded for display purposes only.

round(tapply(d12$bluegill,list(d12$herbicide,d12$year),FUN=mean)-pred,1)

p1a <- pvals.fnc(lme1a,nsim=1000,withMCMC=TRUE)
p1a$fixed
p1a$random

mcmcM <- as.matrix(p1a$mcmc)
m <- data.frame(Value=mcmcM[,1],Predictor=rep(colnames(mcmcM)[1],nrow(mcmcM)))
for (i in 2:ncol(mcmcM)) {
  mtmp <- data.frame(Value=mcmcM[,i],Predictor=rep(colnames(mcmcM)[i],nrow(mcmcM)))
  m <- rbind(m, mtmp)
}
densityplot(~Value|Predictor,data=m,scales=list(relation="free"),
    par.strip.text=list(cex=0.75),xlab="Posterior Values",ylab="Density",pch=".")

d13 <- read.table("data/box3_13.txt",header=TRUE)
str(d13)
d13$arcsurv <- asin(d13$survival/100)
d13$arcsrsurv <- asin(sqrt(d13$survival/100))
str(d13)

lm1 <- lm(arcsurv~size*release,data=d13)
Anova(lm1,type="III")

lsmeans(lm1,~size)
lsmeans(lm1,~release)
lsmeans(lm1,~size*release)

# needed for Sexpr below
lm2 <- lm(arcsrsurv~size*release,data=d13)
p <- Anova(lm2,type="III")
p.int <- p["size:release","Pr(>F)"]
p.size <- p["size","Pr(>F)"] 
p.release <- p["release","Pr(>F)"]

lm2 <- lm(arcsrsurv~size*release,data=d13)
Anova(lm2,type="III")

fitPlot(lm2,xlab="Release Size",ylab="Arcsine Square Root of Survival",legend="topleft",main="")

fitPlot(lm2,which="size",xlab="Release Size",ylab="Arcsine Square Root of Survival",legend="topleft",main="",ylim=c(0.01,0.12))
fitPlot(lm2,which="release",xlab="Release Site",ylab="Arcsine Square Root of Survival",legend="topleft",main="",ylim=c(0.01,0.12))

Anova(lm1,type="II")
Anova(lm2,type="II")

d14 <- read.table("data/box3_14.txt",header=TRUE)
d14$lake <- factor(d14$lake)
str(d14)

bwplot(length~(herbicide:lake),data=d14)
Summarize(length~(lake:herbicide),data=d14,digits=1)

lme1 <- lmer(length~1+herbicide+(1|lake:herbicide),data=d14,REML=TRUE)

summary(lme1)
anova(lme1)
lme2 <- lmer(length~1+(1|lake:herbicide),data=d14,REML=TRUE) # fit without fixed term

lme2 <- lmer(length~1+(1|lake:herbicide),data=d14,REML=TRUE) # fit without fixed term
lrt(lme2,com=lme1)                                           # perform LRT comparison

lme1a <- lmer(length~0+herbicide+(1|lake:herbicide),data=d14,REML=TRUE)
smry1a <- summary(lme1a)
coef(smry1a)

ranef(lme1)
se.ranef(lme1)

p1a <- pvals.fnc(lme1a,nsim=1000,withMCMC=TRUE)
p1a$fixed
p1a$random

mcmcM <- as.matrix(p1a$mcmc)
m <- data.frame(Value=mcmcM[,1],Predictor=rep(colnames(mcmcM)[1],nrow(mcmcM)))
for (i in 2:ncol(mcmcM)) {
  mtmp <- data.frame(Value=mcmcM[,i],Predictor=rep(colnames(mcmcM)[i],nrow(mcmcM)))
  m <- rbind(m, mtmp)
}
densityplot(~Value|Predictor,data=m,scales=list(relation="free"),
    par.strip.text=list(cex=0.75),xlab="Posterior Values",ylab="Density",pch=".")

d15 <- read.table("data/box3_15.txt",header=TRUE)
str(d15)
d15$year <- factor(d15$year)
d15$time <- factor(d15$time)
d15$cove <- factor(d15$cove)
d15$logbio <- log10(d15$biomass)
str(d15)

rmsp2 <- function(object,type=c("III","II","I")) {
  type <- match.arg(type)
  if (type=="I") { res <- anova(object)[1:2,1:3] }       # extract df and SS of first three rows
    else if (type=="III") { res <- Anova(object,type=type)[2:4,2:1] }
      else { res <- Anova(object,type=type)[1:3,2:1] }
  res[,"Mean Sq"] <- res[,2]/res[,1]                     # compute MS
  errorMS <- res[3,"Mean Sq"]                            # MS in 3rd position is error MS
  res[,"F"] <- c(res[1:2,"Mean Sq"]/errorMS,NA)          # compute F for first 2 positions (put NA in las position)
  res[,"PR(>F)"] <- c(pf(res[1:2,"F"],res[1:2,"Df"],res[3,"Df"],lower.tail=FALSE),NA)  # convert to p-values
  res
}

lm1 <- lm(terms(logbio~cove+treat+cove:treat+time+treat:time,keep.order=TRUE),data=d15)
Anova(lm1,type="III")

rmsp2(lm1)

Anova(lm1,type="II")
rmsp2(lm1,type="II")


# Script created at 2015-09-04 15:08:19
