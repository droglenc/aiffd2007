library(FSA)       # Summarize, Subset, headtail, fitPlot, lencat, alkIndivAge
library(car)       # Anova, recode
library(lsmeans)   # lsmeans
library(nlstools)  # overview

setwd("C:/aaaWork/Web/fishR/BookVignettes/aiffd2007")

options(width=90,continue=" ",show.signif.stars=FALSE, contrasts=c("contr.sum","contr.poly"))
par(mar=c(3.5,3.5,1,1),mgp=c(2.1,0.4,0),tcl=-0.2)

d1.age <- read.table("data/box5_1_Aged.txt",header=TRUE)
str(d1.age)
d1.len <- read.table("data/box5_1_Length.txt",header=TRUE)
str(d1.len)

Summarize(d1.age$tl,digits=1)

d1.age1 <- lencat(~tl,data=d1.age,startcat=80,w=20)
headtail(d1.age1)      # to verify the correct addition of the length categories

d.raw <- xtabs(~LCat+age,data=d1.age1)
d.key <- prop.table(d.raw,margin=1)
round(d.key,3)                           # rounded for display only

d.key[which(is.na(d.key))] <- 0
round(d.key,3)                           # rounded for display only

Summarize(d1.len$tl,digits=2)

d1.len1 <- Subset(d1.len,tl<500)

d1.len2 <- lencat(~tl,data=d1.len1,startcat=80,w=20)
( len.freq <- xtabs(~LCat,data=d1.len2) )

key3003 <- d.key["300","3"]*100
lf300 <- len.freq["300"]
key3803 <- d.key["380","3"]*100
key3804 <- d.key["380","4"]*100
lf380 <- len.freq["380"]

length(len.freq)  # number of columns in length frequency vector
dim(d.key)        # dimensions (rows, columns) of the age-length key

0*0+0*0+0*0+0*0+0*0+0*0+0*0+0*0+0*0+0*0+0*0+3*1+6*1+12*1+12*1+30*0.333+28*0+48*0+
  51*0+61*0+83*0+0*0

( age.freq <- len.freq %*% d.key )

d1.len3 <- alkIndivAge(d.key,~tl,data=d1.len1)
headtail(d1.len3)

d1.comb <- rbind(d1.len3,d1.age[,c("tl","age")])
headtail(d1.comb)

xtabs(~age,data=d1.comb)
Summarize(tl~age,data=d1.comb,digits=2)

d2 <- read.table("data/box5_2.txt",header=TRUE)
str(d2)
head(d2)

d2$Si2 <- d2$Si*10/24   # convert radial measurements
d2$Sc2 <- d2$Sc*10/24   # convert total scale radius
head(d2)

d2$LiDL <- d2$Lc*d2$Si2/d2$Sc2
head(d2)
Summarize(LiDL~inc,data=d2)

lm1 <- lm(LiDL~Si2,data=d2)
anova(lm1)
summary(lm1)

a <- coef(lm1)[1]

a <- coef(lm1)[1]               # get the correction factor (first coefficient in lm1)
d2$LiFL <- a + (d2$Lc-a)*(d2$Si2/d2$Sc2)
head(d2)
Summarize(LiFL~inc,data=d2)

d2$incSqrd <- d2$inc^2
lm2 <- lm(LiFL~sex+inc+incSqrd+inc*sex+incSqrd*sex,data=d2)
anova(lm2)             # type-I SS
Anova(lm2,type="III")  # type-III SS

# base schematic for adding points to
plot(LiFL~inc,data=d2,col="white",xlab="Age",ylab="Back-Calculated Length")
# add points for males and females, with slight offsets so they can be seen
points(LiFL~I(inc-0.1),col="red",pch=16,data=Subset(d2,sex=="M"))
points(LiFL~I(inc+0.1),col="black",pch=16,data=Subset(d2,sex=="F"))

# create a sequence of increment numbers
incs <- seq(1,17,0.1)
# predict ages for males and females & then plot the lines
predM <- predict(lm2,data.frame(inc=incs,incSqrd=incs^2,sex=rep("M",length(incs))))
predF <- predict(lm2,data.frame(inc=incs,incSqrd=incs^2,sex=rep("F",length(incs))))
lines(predM~incs,col="red",lwd=2)
lines(predF~incs,col="black",lwd=2)

d4 <- read.table("data/box5_4.txt",header=TRUE)
str(d4)

vbStarts(tl~age,data=d4,type="typical",dynamicPlot=TRUE)

sv <- list(Linf=481,K=0.47,t0=0)

vbmdl <- tl~Linf*(1-exp(-K*(age-t0)))

ssvb <- nls(vbmdl,start=sv,data=d4)
overview(ssvb)

plot(tl~age,data=d4,pch=16,xlab="Age",ylab="Total Length (mm)",xlim=c(0,11),ylim=c(0,550))
ages <- seq(0,11,0.1)
preds <- predict(ssvb,data.frame(age=ages))
lines(preds~ages,lwd=2,col="red")

d5 <- read.table("data/box5_5.txt",header=TRUE)
str(d5)
d5$bcage <- factor(d5$bcage)

d5$group <- recode(d5$bcyear,"1995:1999='normal';else='dry'",as.factor.result=TRUE)
str(d5)

head(d5,n=10)

d5fish <- unique(d5[,1:6])
str(d5fish)

xtabs(~year+age,data=d5fish)

tbl2 <- xtabs(~group+bcage,data=d5)
addmargins(tbl2)

lm1 <- lm(growth~bcage*group,data=d5)

lm1 <- lm(growth~bcage*group,data=d5)
anova(lm1)

Anova(lm1,type="III")

lm2 <- lm(growth~bcage+group,data=d5)
anova(lm2)
Anova(lm2,type="III")

lsmeans(lm2,~group)
lsmeans(lm2,~bcage)

fitPlot(lm2,xlab="Age",ylab="Growth (mm)",legend="topright",main="")

# an alternative look
library(effects)
plot(effect("bcage",lm2),xlab="Age",ylab="Growth (mm)",main="",ylim=c(0,140))
plot(effect("group",lm2),xlab="Rainfall Group",ylab="Growth (mm)",main="",ylim=c(0,140))
plot(effect("bcage*group",lm2))

d6 <- read.table("data/box5_6.txt",header=TRUE)
str(d6)
d6$timeouty <- d6$timeoutd/365

svLinf <- max(d6$clrecap)
svK <- with(d6, mean((log(clrecap)-log(clmark))/timeouty))
( Fsv <- list(Linf=svLinf,K=svK) )

Fvbmdl <- clrecap ~ clmark+(Linf-clmark)*(1-exp(-K*timeouty))
tvb <- nls(Fvbmdl,start=Fsv,data=d6)
overview(tvb)


# Script created at 2015-05-13 19:49:12
